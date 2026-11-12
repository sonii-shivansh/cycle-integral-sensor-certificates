"""Exact-rational audit of pairwise multi-window CIFC constraints.

All supplied binary64 sample values, times, and parameter bounds are treated
as EXACTLY their rational binary64 values. Infeasibility yields a checkable
Farkas certificate for the resulting *linear relaxation*, not proof that the
physical measurement assumptions hold. Does not use SciPy or floating-point
arithmetic internally. Intended for small, independently reviewed examples.
"""
from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
import math


def _q(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral) or isinstance(value, str):
        return Fraction(value)
    v = float(value)
    if not math.isfinite(v):
        raise ValueError("samples and bounds must be finite")
    return Fraction.from_float(v)


@dataclass(frozen=True)
class ExactWindow:
    integrals: tuple[Fraction, ...]
    increments: tuple[Fraction, ...]
    length: Fraction
    max_gap: Fraction


def exact_window(times, readings):
    """Exact trapezoid integrals on increasing, possibly nonuniform time grids.

    Uses a COMMON reference shift before integration; this changes no
    differential or linear-feasibility identity in exact arithmetic.
    """
    t = tuple(map(_q, times))
    if len(t) < 2:
        raise ValueError("need at least two time samples")
    gaps = tuple(t[k + 1] - t[k] for k in range(len(t) - 1))
    if any(g <= 0 for g in gaps):
        raise ValueError("timestamps must be strictly increasing")
    z = tuple(tuple(map(_q, row)) for row in readings)
    if not z or any(len(row) != len(t) for row in z):
        raise ValueError("readings must have shape (n_sensors, len(times))")
    reference = z[0][0]
    integrals = []
    for row in z:
        values = [value - reference for value in row]
        integral = sum(
            (gap * (values[k] + values[k + 1]) / 2
             for k, gap in enumerate(gaps)), Fraction(0)
        )
        integrals.append(integral)
    return ExactWindow(tuple(integrals),
                       tuple(row[-1] - row[0] for row in z),
                       t[-1] - t[0], max(gaps))


@dataclass(frozen=True)
class Inequality:
    """Exact rational linear constraint a*tau_i + b*tau_j <= rhs."""
    a: Fraction
    b: Fraction
    rhs: Fraction
    label: str


@dataclass(frozen=True)
class FarkasCertificate:
    """Nonnegative multipliers lambda with lambda*A=0, lambda*rhs<0."""
    weights: tuple[Fraction, ...]
    combined_rhs: Fraction

    def check(self, rows):
        if len(self.weights) != len(rows) or any(x < 0 for x in self.weights):
            return False
        return (
            sum((w*r.a for w, r in zip(self.weights, rows)), Fraction(0)) == 0
            and sum((w*r.b for w, r in zip(self.weights, rows)), Fraction(0)) == 0
            and sum((w*r.rhs for w, r in zip(self.weights, rows)), Fraction(0)) == self.combined_rhs
            and self.combined_rhs < 0
        )


@dataclass(frozen=True)
class ExactVerdict:
    status: str  # exact_infeasible or feasible_relaxation
    inequalities: tuple[Inequality, ...]
    certificate: FarkasCertificate | None


def _sensor_bounds(value, count, name):
    if isinstance(value, (list, tuple)) or (
            hasattr(value, "ndim") and getattr(value, "ndim") != 0):
        values = tuple(map(_q, value))
        if len(values) != count:
            raise ValueError(name + " requires one value per sensor")
    else:
        values = (_q(value),) * count
    if any(v < 0 for v in values):
        raise ValueError(name + " must be nonnegative")
    return values


def pairwise_constraints(windows, indices, max_tau, eps, curvature):
    """Exactly eliminate shared J_w from the relaxed two-sensor LP."""
    if not windows:
        raise ValueError("need at least one window")
    if len(windows) > 24:
        raise ValueError("exact two-dimensional verifier limited to 24 windows")
    n = len(windows[0].integrals)
    if not n or any(
        len(w.integrals) != n or len(w.increments) != n
        or w.length <= 0 or w.max_gap <= 0 or w.max_gap > w.length
        for w in windows
    ):
        raise ValueError("invalid window summaries")
    if (len(indices) != 2 or indices[0] == indices[1] or any(
            not isinstance(i, Integral) or i < 0 or i >= n for i in indices)):
        raise ValueError("need two distinct valid sensor indices")
    i, j = indices
    M = _q(max_tau)
    if M < 0:
        raise ValueError("max_tau must be nonnegative")
    e = _sensor_bounds(eps, n, "eps")
    k = _sensor_bounds(curvature, n, "curvature")
    rows = []
    for wi, w in enumerate(windows):
        Ii, Ij = w.integrals[i], w.integrals[j]
        di, dj = w.increments[i], w.increments[j]
        Di = w.length*e[i] + w.length*w.max_gap**2*k[i]/12
        Dj = w.length*e[j] + w.length*w.max_gap**2*k[j]/12
        radius = Di + Dj
        # lower_i <= upper_j
        rows.append(Inequality(
            di - 2*e[i], -(dj + 2*e[j]), Ij - Ii + radius,
            f"window{wi}:lower_i<=upper_j"
        ))
        # lower_j <= upper_i
        rows.append(Inequality(
            -(di + 2*e[i]), dj - 2*e[j], Ii - Ij + radius,
            f"window{wi}:lower_j<=upper_i"
        ))
    rows.extend([
        Inequality(Fraction(-1), Fraction(0), Fraction(0), "tau_i>=0"),
        Inequality(Fraction(1), Fraction(0), M, "tau_i<=M"),
        Inequality(Fraction(0), Fraction(-1), Fraction(0), "tau_j>=0"),
        Inequality(Fraction(0), Fraction(1), M, "tau_j<=M"),
    ])
    return tuple(rows)


def _eliminate(rows, coordinate):
    """Fourier-Motzkin variable elimination, preserving positive multipliers."""
    pos, neg, zero = [], [], []
    for coeffs, rhs, weights in rows:
        a = coeffs[coordinate]
        if a > 0:
            pos.append((coeffs, rhs, weights))
        elif a < 0:
            neg.append((coeffs, rhs, weights))
        else:
            zero.append((coeffs, rhs, weights))
    output = list(zero)
    for pc, pr, pw in pos:
        for nc, nr, nw in neg:
            a, b = -nc[coordinate], pc[coordinate]
            output.append((
                tuple(a*p + b*n for p, n in zip(pc, nc)),
                a*pr + b*nr,
                tuple(a*p + b*n for p, n in zip(pw, nw)),
            ))
    return output


def exact_pairwise_verdict(windows, indices, max_tau, eps, curvature):
    """Exact feasibility and independently checkable rational Farkas proof.

    The claim is restricted to the *relaxed* two-variable LP on the exactly
    represented input samples, times and bounds. Feasibility is not a proof of
    sensor health. Worst-case runtime can grow quickly with window count.
    """
    constraints = pairwise_constraints(windows, indices, max_tau, eps, curvature)
    n = len(constraints)
    rows = [
        ((r.a, r.b), r.rhs,
         tuple(Fraction(int(i == j)) for i in range(n)))
        for j, r in enumerate(constraints)
    ]
    for coordinate in (0, 1):
        rows = _eliminate(rows, coordinate)
        for coefficients, rhs, weights in rows:
            if coefficients == (0, 0) and rhs < 0:
                cert = FarkasCertificate(weights, rhs)
                if not cert.check(constraints):
                    raise ArithmeticError("internal exact proof verification failed")
                return ExactVerdict("exact_infeasible", constraints, cert)
    if any(coefficients != (0, 0) or rhs < 0 for coefficients, rhs, _ in rows):
        raise ArithmeticError("exact variable elimination did not terminate correctly")
    return ExactVerdict("feasible_relaxation", constraints, None)


def certificate_summary(verdict):
    """Human-readable contradiction, suitable for inclusion in a review."""
    if verdict.certificate is None:
        return {"status": verdict.status, "proof": "no contradiction in the LP relaxation"}
    if not verdict.certificate.check(verdict.inequalities):
        raise ArithmeticError("Farkas certificate is invalid")
    return {
        "status": verdict.status,
        "weighted_constraints": [
            (row.label, str(weight))
            for row, weight in zip(verdict.inequalities, verdict.certificate.weights)
            if weight
        ],
        "combined_rhs": str(verdict.certificate.combined_rhs),
        "contradiction": f"0 <= {verdict.certificate.combined_rhs}",
    }
