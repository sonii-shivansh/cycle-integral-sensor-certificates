"""Demonstrate float64 cancellation in a mathematically healthy sensor ensemble.

This script is an implementation precision probe, not a physical benchmark.
The naive uncentered solver status can depend on the numerical library/version.
"""

from dataclasses import replace

import numpy as np
from scipy.optimize import linprog

from experiments.multiwindow_lp import summarize_window, joint_relaxation


def compare():
    t = np.arange(201, dtype=float)
    tau = np.array([0.0, 1.0, 2.0])
    # Exact model: T(t) = 2**52 + t, y_i(t) = T(t) - tau_i.
    # Every recorded sample is exactly representable in float64.
    z = 2.0**52 + t[None, :] - tau[:, None]
    centered = []
    uncentered = []
    for a, b in [(0, 50), (50, 120), (120, 200)]:
        t_w = t[a:b + 1]
        z_w = z[:, a:b + 1]
        summary = summarize_window(t_w, z_w)
        centered.append(summary)
        # Intentionally emulate the former numerically unstable algorithm.
        naive = np.trapezoid(z_w, t_w, axis=1)
        uncentered.append(replace(summary, integrals=naive))

    naive = joint_relaxation(uncentered, 3.0, 0.0, 0.0)
    fixed = joint_relaxation(centered, 3.0, 0.0, 0.0)
    print("Naive absolute integral LP status:", naive.status)
    print("Reference-centered integral LP status:", fixed.status)
    print("First-window naive differences:", np.diff(uncentered[0].integrals))
    print("First-window centered differences:", np.diff(centered[0].integrals))
    print("Centered parameter estimates:", fixed.time_constants)
    if fixed.status != "feasible_relaxation":
        raise RuntimeError("Healthy exactly representable case failed after centering")
    return naive.status, fixed.status


if __name__ == "__main__":
    compare()
