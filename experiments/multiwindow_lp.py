"""Exploratory multi-window *relaxed* consistency check using linear programming.

This solves a necessary outer relaxation of the common-input first-order model;
it does NOT prove healthy sensors when feasible, identify defective hardware, or
use certified interval/exact arithmetic to prove infeasibility numerically.
Requires numpy and scipy (see requirements-research.txt).
"""

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.optimize import linprog


@dataclass(frozen=True)
class WindowSummary:
    """One window's observed trapezoid integrals and endpoint increments."""
    integrals: np.ndarray
    increments: np.ndarray
    length: float
    step: float


@dataclass(frozen=True)
class LPFinding:
    """Relaxation status; 'solver_infeasible' is NOT an exact certificate."""
    status: str
    indices: tuple[int, ...]
    time_constants: np.ndarray | None
    integrated_inputs: np.ndarray | None
    max_constraint_violation: float | None
    message: str


def summarize_window(t, z):
    """Build window summary; z shape (sensors, samples), uniform time grid."""
    t = np.asarray(t, dtype=float)
    z = np.asarray(z, dtype=float)
    if t.ndim != 1 or len(t) < 2:
        raise ValueError("t must contain at least two time samples")
    if z.ndim != 2 or z.shape[1] != len(t) or z.shape[0] < 1:
        raise ValueError("z must have shape (n_sensors, len(t))")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(z)):
        raise ValueError("samples must be finite")
    diffs = np.diff(t)
    if np.any(diffs <= 0):
        raise ValueError("timestamps must be strictly increasing")
    h = float(diffs[0])
    if not np.allclose(diffs, h, rtol=1e-10, atol=1e-12):
        raise ValueError("timestamps must be uniformly spaced")
    return WindowSummary(
        integrals=np.asarray(np.trapezoid(z, t, axis=1), dtype=float),
        increments=np.asarray(z[:, -1] - z[:, 0], dtype=float),
        length=float(t[-1] - t[0]),
        step=h,
    )


def _bounds(values, n, name):
    a = np.asarray(values, dtype=float)
    if a.ndim == 0:
        a = np.full(n, float(a))
    if a.shape != (n,) or not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError(name + " must be a finite nonnegative scalar or length-n array")
    return a


def _check(windows, indices, max_tau, eps, curvature):
    if not windows:
        raise ValueError("at least one observation window is required")
    n = len(windows[0].integrals)
    if n < 1:
        raise ValueError("at least one sensor is required")
    for w in windows:
        if not isinstance(w, WindowSummary):
            raise TypeError("expected WindowSummary; use summarize_window")
        if (w.integrals.shape != (n,) or w.increments.shape != (n,)
                or not np.all(np.isfinite(w.integrals))
                or not np.all(np.isfinite(w.increments))
                or not np.isfinite(w.length) or w.length <= 0
                or not np.isfinite(w.step) or w.step <= 0
                or w.step > w.length):
            raise ValueError("invalid or inconsistent window summaries")
    if not np.isfinite(max_tau) or max_tau < 0:
        raise ValueError("max_tau must be finite and nonnegative")
    e = _bounds(eps, n, "eps")
    k = _bounds(curvature, n, "curvature")
    chosen = tuple(range(n)) if indices is None else tuple(indices)
    if not chosen or len(chosen) != len(set(chosen)):
        raise ValueError("choose nonempty, distinct sensor indices")
    if any(not isinstance(i, (int, np.integer)) or i < 0 or i >= n for i in chosen):
        raise ValueError("indices out of range")
    return chosen, e, k


def joint_relaxation(windows, max_tau, eps, curvature, indices=None):
    """LP for one fixed tau_i per chosen sensor across every observation window.

    Each healthy sensor must admit J_w = integral(T) satisfying
       I_iw + tau_i (delta_iw - 2 eps_i) - D_iw <= J_w
       J_w <= I_iw + tau_i (delta_iw + 2 eps_i) + D_iw,
       D_iw = L_w eps_i + L_w h_w^2 K_i/12.

    Feasibility is NOT sufficiency for the differential model. Numerical solver
    'infeasible' is a candidate contradiction; certify away from tolerance or
    verify via exact/rational arithmetic before treating it as a proof.
    """
    chosen, eps, curvature = _check(windows, indices, max_tau, eps, curvature)
    s, w_count = len(chosen), len(windows)
    A, b = [], []
    for w_idx, w in enumerate(windows):
        for local, i in enumerate(chosen):
            I, d = float(w.integrals[i]), float(w.increments[i])
            D = w.length * eps[i] + w.length * w.step**2 * curvature[i] / 12.0
            # I + tau * (d - 2eps) - D <= J
            row = np.zeros(s + w_count)
            row[local] = d - 2 * eps[i]
            row[s + w_idx] = -1.0
            A.append(row)
            b.append(D - I)
            # J <= I + tau * (d + 2eps) + D
            row = np.zeros(s + w_count)
            row[local] = -(d + 2 * eps[i])
            row[s + w_idx] = 1.0
            A.append(row)
            b.append(I + D)
    A, b = np.asarray(A), np.asarray(b)
    result = linprog(
        np.zeros(s + w_count), A_ub=A, b_ub=b,
        bounds=[(0.0, float(max_tau))] * s + [(None, None)] * w_count,
        method="highs",
    )
    if result.status == 0 and result.x is not None:
        residual = A @ result.x - b
        return LPFinding("feasible_relaxation", chosen, result.x[:s],
                         result.x[s:], float(np.max(residual)), result.message)
    if result.status == 2:
        return LPFinding("solver_infeasible", chosen, None, None, None, result.message)
    return LPFinding("inconclusive", chosen, None, None, None,
                     "HiGHS status %s: %s" % (result.status, result.message))


def signaware_window_passes(w, i, j, max_tau, eps, curvature):
    """Exact-form single-window sign-aware condition (float arithmetic)."""
    n = len(w.integrals)
    if i == j or min(i, j) < 0 or max(i, j) >= n:
        raise ValueError("two different valid sensor indices are required")
    e, k = _bounds(eps, n, "eps"), _bounds(curvature, n, "curvature")
    if not np.isfinite(max_tau) or max_tau < 0:
        raise ValueError("max_tau must be nonnegative")
    d_i, d_j = w.increments[i], w.increments[j]
    D = w.length * (e[i] + e[j]) + w.length * w.step**2 * (k[i] + k[j]) / 12.0
    lo = max_tau * min(0.0, d_i - 2 * e[i]) - max_tau * max(0.0, d_j + 2 * e[j]) - D
    hi = max_tau * max(0.0, d_i + 2 * e[i]) - max_tau * min(0.0, d_j - 2 * e[j]) + D
    Q = w.integrals[j] - w.integrals[i]
    return bool(lo <= Q <= hi)


def pairwise_results(windows, max_tau, eps, curvature):
    """Pair-specific LP findings. Non-edges do not prove a healthy pair."""
    n = len(windows[0].integrals)
    return {(i, j): joint_relaxation(windows, max_tau, eps, curvature, (i, j))
            for i, j in combinations(range(n), 2)}
