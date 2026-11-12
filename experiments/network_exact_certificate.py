"""Exact-rational Farkas verifier for multi-sensor multi-window CIFC LPs.

This certifies infeasibility of the *finite linear outer relaxation* formed from
recorded numerical samples interpreted as exact rationals. It does not prove
correct physical assumptions, hardware failure, scientific novelty, or
feasibility of the full continuous-time model.

A floating-point SciPy dual LP proposes a sparse witness. The witness is
reconstructed by exact rational Gaussian elimination and accepted ONLY when
an independent rational equality/inequality check proves it valid.

If reconstruction fails, the result is NOT_CERTIFIED, never a false proof.
This research verifier is not a complete or guaranteed-fast exact LP solver.
"""
from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
import math

import numpy as np
from scipy.optimize import linprog


def _q(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, (Integral, str)):
        return Fraction(value)
    f = float(value)
    if not math.isfinite(f):
        raise ValueError("nonfinite input")
    return Fraction.from_float(f)


def _vector(value, count):
    if isinstance(value, (list, tuple, np.ndarray)) and np.ndim(value) > 0:
        values = tuple(_q(v) for v in value)
        if len(values) != count:
            raise ValueError("per-sensor bound length mismatch")
    else:
        values = (_q(value),) * count
    if any(v < 0 for v in values):
        raise ValueError("negative bound")
    return values


@dataclass(frozen=True)
class NetworkRow:
    coefficients: tuple[Fraction, ...]
    rhs: Fraction
    label: str


@dataclass(frozen=True)
class NetworkFarkasCertificate:
    weights: tuple[Fraction, ...]
    weighted_rhs: Fraction

    def check(self, rows):
        """Recheck the original exact system independently of HiGHS."""
        if not rows or len(self.weights) != len(rows):
            return False
        if any(weight < 0 for weight in self.weights):
            return False
        dim = len(rows[0].coefficients)
        if any(len(r.coefficients) != dim for r in rows):
            return False
        return (
            all(sum(
                (self.weights[i] * rows[i].coefficients[j]
                 for i in range(len(rows))), Fraction(0)
            ) == 0 for j in range(dim))
            and sum(
                (weight * row.rhs for weight, row in zip(self.weights, rows)),
                Fraction(0)
            ) == self.weighted_rhs
            and self.weighted_rhs < 0
        )


@dataclass(frozen=True)
class NetworkAudit:
    status: str  # certified_infeasible / not_certified
    certificate: NetworkFarkasCertificate | None
    message: str
    rows: tuple[NetworkRow, ...]


def exact_network_rows(windows, max_tau, eps, curvature, sensor_indices=None):
    """Construct exact Ax<=b for tau_i (bounded) and J_w (unrestricted).

    Windows should be produced by exact_pairwise_certificate.exact_window(),
    which computes trapezoid integrals from actual represented samples using
    exact Fraction arithmetic and one reference shift per window.
    """
    if not windows:
        raise ValueError("at least one window is necessary")
    n = len(windows[0].integrals)
    selected = tuple(range(n)) if sensor_indices is None else tuple(sensor_indices)
    if (not selected or len(selected) != len(set(selected))
            or any(not isinstance(i, Integral) or i < 0 or i >= n
                   for i in selected)):
        raise ValueError("invalid sensor subset")
    M = _q(max_tau)
    if M < 0:
        raise ValueError("negative response time bound")
    eb, kb = _vector(eps, n), _vector(curvature, n)
    W = len(windows)
    dim = len(selected) + W
    rows = []

    for w_index, win in enumerate(windows):
        if (len(win.integrals) != n or len(win.increments) != n
                or _q(win.length) <= 0 or _q(win.max_gap) <= 0
                or _q(win.max_gap) > _q(win.length)):
            raise ValueError("invalid window summary")
        L, h = _q(win.length), _q(win.max_gap)
        for local_index, i in enumerate(selected):
            I, delta = _q(win.integrals[i]), _q(win.increments[i])
            D = L * eb[i] + L * h**2 * kb[i] / 12

            # I + tau_i*(delta - 2*eps_i) - D <= J_w
            a = [Fraction(0)] * dim
            a[local_index] = delta - 2 * eb[i]
            a[len(selected) + w_index] = Fraction(-1)
            rows.append(NetworkRow(tuple(a), D - I,
                                   f"window{w_index}/sensor{i}:lower"))

            # J_w <= I + tau_i*(delta + 2*eps_i) + D
            a = [Fraction(0)] * dim
            a[local_index] = -(delta + 2 * eb[i])
            a[len(selected) + w_index] = Fraction(1)
            rows.append(NetworkRow(tuple(a), I + D,
                                   f"window{w_index}/sensor{i}:upper"))

    for local_index, i in enumerate(selected):
        a = [Fraction(0)] * dim
        a[local_index] = Fraction(-1)
        rows.append(NetworkRow(tuple(a), Fraction(0), f"sensor{i}:tau>=0"))
        a = [Fraction(0)] * dim
        a[local_index] = Fraction(1)
        rows.append(NetworkRow(tuple(a), M, f"sensor{i}:tau<=M"))
    return tuple(rows)


def _rational_linear_solve(matrix, rhs):
    """Exact Gauss-Jordan for over/under-determined systems; free vars=0."""
    m, n = len(matrix), len(matrix[0])
    table = [[_q(x) for x in matrix[i]] + [_q(rhs[i])]
             for i in range(m)]
    rank = 0
    pivots = []
    for col in range(n):
        pivot = next((j for j in range(rank, m) if table[j][col]), None)
        if pivot is None:
            continue
        table[rank], table[pivot] = table[pivot], table[rank]
        divisor = table[rank][col]
        table[rank] = [v / divisor for v in table[rank]]
        for j in range(m):
            if j == rank:
                continue
            multiplier = table[j][col]
            if multiplier:
                table[j] = [
                    a - multiplier * b
                    for a, b in zip(table[j], table[rank])
                ]
        pivots.append((col, rank))
        rank += 1
        if rank == m:
            break
    if any(
        all(table[i][j] == 0 for j in range(n)) and table[i][-1] != 0
        for i in range(m)
    ):
        return None
    result = [Fraction(0)] * n
    for col, row in pivots:
        result[col] = table[row][-1]
    return result


def certify_network_infeasibility(
    windows, max_tau, eps, curvature, sensor_indices=None
):
    """Return a rigorous exact contradiction, or explicitly NO CERTIFICATE.

    The candidate HiGHS dual is only a heuristic support selector. The return
    value 'not_certified' must NEVER be interpreted as proof of feasibility.
    """
    rows = exact_network_rows(windows, max_tau, eps, curvature, sensor_indices)
    dimension = len(rows[0].coefficients)
    count = len(rows)

    try:
        equations = np.array([
            [float(row.coefficients[k]) for row in rows]
            for k in range(dimension)
        ] + [[1.0] * count])
        objective = np.array([float(row.rhs) for row in rows])
        candidate = linprog(
            objective, A_eq=equations,
            b_eq=np.r_[np.zeros(dimension), 1.0],
            bounds=(0, None), method="highs-ds"
        )
    except (OverflowError, ValueError) as exc:
        return NetworkAudit("not_certified", None,
                            "numerical candidate construction failed: " + str(exc),
                            rows)

    if candidate.status != 0 or candidate.x is None:
        return NetworkAudit(
            "not_certified", None,
            "no usable floating dual candidate: " + candidate.message, rows
        )
    if candidate.fun >= -1e-8:
        return NetworkAudit("not_certified", None,
                            "no robust negative floating dual candidate", rows)

    for threshold in (1e-7, 1e-9, 1e-11, 1e-13):
        support = np.flatnonzero(candidate.x > threshold).tolist()
        if not support:
            continue
        matrix = [
            [rows[j].coefficients[k] for j in support]
            for k in range(dimension)
        ] + [[Fraction(1)] * len(support)]
        target = [Fraction(0)] * dimension + [Fraction(1)]
        exact_weights = _rational_linear_solve(matrix, target)
        if exact_weights is None or any(w < 0 for w in exact_weights):
            continue

        weights = [Fraction(0)] * count
        for j, value in zip(support, exact_weights):
            weights[j] = value
        rhs = sum(
            (weight * row.rhs for weight, row in zip(weights, rows)),
            Fraction(0)
        )
        cert = NetworkFarkasCertificate(tuple(weights), rhs)
        if cert.check(rows):
            return NetworkAudit("certified_infeasible", cert,
                                f"exact rational Farkas proof using {len(support)} candidate rows",
                                rows)
    return NetworkAudit(
        "not_certified", None,
        "floating dual support could not be reconstructed as an exact proof",
        rows
    )


def summarize_network_certificate(audit):
    """Export checkable nonzero witnesses; nothing is claimed if not proven."""
    if audit.status != "certified_infeasible" or audit.certificate is None:
        return {"status": "not_certified", "reason": audit.message}
    if not audit.certificate.check(audit.rows):
        raise ArithmeticError("Farkas witness failed exact verification")
    return {
        "status": audit.status,
        "contradiction": f"0 <= {audit.certificate.weighted_rhs}",
        "nonzero_weights": [
            {"row": row.label, "lambda": str(weight)}
            for row, weight in zip(audit.rows, audit.certificate.weights)
            if weight > 0
        ],
    }
