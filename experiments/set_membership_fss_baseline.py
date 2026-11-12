"""Exact 2D feasible-parameter-set intersection baseline for CIFC sensor pairs.

Conceptual reference: classical invariant-parameter set-membership diagnosis,
NOT a reproduction of the 2005 paper's parallelotope approximation algorithm.
The parameter polygon is clipped by every observation's affine half-plane.
Does not import the Fourier--Motzkin verifier or SciPy.
"""
from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
import math


def q(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, (str, Integral)):
        return Fraction(value)
    f = float(value)
    if not math.isfinite(f):
        raise ValueError("nonfinite input")
    return Fraction.from_float(f)


def bound_vector(value, n):
    if isinstance(value, (list, tuple)) or (
            hasattr(value, "ndim") and value.ndim != 0):
        result = tuple(q(v) for v in value)
        if len(result) != n:
            raise ValueError("invalid bound length")
    else:
        result = (q(value),) * n
    if any(v < 0 for v in result):
        raise ValueError("bounds must be nonnegative")
    return result


def clip(polygon, a, b, c):
    """Exact rational polygon intersection with a*x + b*y <= c."""
    if not polygon:
        return ()
    out = []
    previous = polygon[-1]
    previous_distance = a*previous[0] + b*previous[1] - c
    for point in polygon:
        distance = a*point[0] + b*point[1] - c
        was_in = previous_distance <= 0
        is_in = distance <= 0
        if was_in != is_in:
            ratio = previous_distance / (previous_distance-distance)
            intersection = (
                previous[0] + ratio*(point[0]-previous[0]),
                previous[1] + ratio*(point[1]-previous[1]),
            )
            out.append(intersection)
        if is_in:
            out.append(point)
        previous, previous_distance = point, distance
    deduped = []
    for point in out:
        if not deduped or point != deduped[-1]:
            deduped.append(point)
    if len(deduped) > 1 and deduped[0] == deduped[-1]:
        deduped.pop()
    return tuple(deduped)


@dataclass(frozen=True)
class BaselineVerdict:
    status: str
    vertices: tuple
    first_inconsistent_window: int | None


def intersect_parameter_strips(windows, pair, max_tau, eps, curvature):
    """Intersect all constraints independently by exact polygon clipping.

    'windows' contains ExactWindow-like summaries of the *recorded* numeric
    values interpreted as rational numbers. This is a feasible solution set
    for the relaxed pairwise model, not the full sensor ODE. It assumes each
    sensor has one constant response-time parameter across every window.
    """
    if not windows:
        raise ValueError("at least one window required")
    n = len(windows[0].integrals)
    if (len(pair) != 2 or pair[0] == pair[1]
            or any(not isinstance(i, Integral) or i < 0 or i >= n
                   for i in pair)):
        raise ValueError("invalid sensor pair")
    if any(len(w.integrals) != n or len(w.increments) != n
           or q(w.length) <= 0 or q(w.max_gap) <= 0
           or q(w.max_gap) > q(w.length) for w in windows):
        raise ValueError("invalid window")
    M = q(max_tau)
    if M < 0:
        raise ValueError("response-time bound must be nonnegative")
    errors, curves = bound_vector(eps, n), bound_vector(curvature, n)
    i, j = pair
    polygon = ((Fraction(0), Fraction(0)), (M, Fraction(0)),
               (M, M), (Fraction(0), M))
    for window_index, w in enumerate(windows):
        Ii, Ij = q(w.integrals[i]), q(w.integrals[j])
        di, dj = q(w.increments[i]), q(w.increments[j])
        L, h = q(w.length), q(w.max_gap)
        D = L*(errors[i]+errors[j]) + L*h*h*(curves[i]+curves[j])/12
        # Require overlap of each sensor's possible shared-input integral.
        # These conditions independently reconstruct the per-window strip.
        polygon = clip(polygon, di-2*errors[i], -(dj+2*errors[j]),
                       Ij-Ii+D)
        polygon = clip(polygon, -(di+2*errors[i]), dj-2*errors[j],
                       Ii-Ij+D)
        if not polygon:
            return BaselineVerdict("exact_infeasible", (), window_index)
    return BaselineVerdict("feasible_relaxation", polygon, None)
