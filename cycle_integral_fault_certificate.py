"""Cycle-Integral Fault Certificate (CIFC), a proof-of-concept.

Requires NumPy. Run: python cycle_integral_fault_certificate.py

Soundness assumptions for any two 'healthy' sensors:
  * A shared continuous temperature T(t),
  * tau_i * dy_i/dt + y_i = T(t), with 0 <= tau_i <= TAU_MAX,
  * sampled observations z_i(k) with |z_i(k)-y_i(t_k)| <= EPS,
  * |d^2 y_i/dt^2| <= CURVATURE_MAX between samples.
Uniform sampling on [0, PERIOD].
The trapezoid integration error is accounted for explicitly.

This algorithm detects inconsistencies with the model, not guaranteed
hardware faults. No alarms means 'not disproved', not 'healthy'.
"""

import itertools
import math
import numpy as np

PERIOD = 20.0
DT = 0.05
AMBIENT_BASE = 25.0
AMBIENT_AMP = 8.0
W = 2 * math.pi / PERIOD
TAU_MAX = 15.0
EPS = 0.02
TAUS = np.array([1.0, 2.0, 3.0, 5.0, 8.0, 13.0])
# The unknown healthy physical response is bounded by this curvature.
CURVATURE_MAX = AMBIENT_AMP * W**2


def readings(seed=123, faults=None):
    """Make exactly solvable first-order periodic sensor signals plus bounded sample noise."""
    t = np.linspace(0.0, PERIOD, round(PERIOD / DT) + 1)
    x = W * t
    traces = np.array([
        AMBIENT_BASE + AMBIENT_AMP * (np.sin(x) - (W * tau) * np.cos(x))
        / (1.0 + (W * tau)**2)
        for tau in TAUS
    ])
    for index, offset in (faults or {}).items():
        traces[index] += offset
    rng = np.random.default_rng(seed)
    z = traces + rng.uniform(-EPS, EPS, size=traces.shape)
    return t, z


def certificates(t, z):
    """Pair violations are sound: no pair of healthy sensors can violate the bound."""
    n = z.shape[0]
    length = t[-1] - t[0]
    h = float(t[1] - t[0])
    assert np.allclose(np.diff(t), h)
    flags = []
    all_tests = []
    for i, j in itertools.combinations(range(n), 2):
        integral = float(np.trapezoid(z[j] - z[i], t))
        delta_i = float(z[i, -1] - z[i, 0])
        delta_j = float(z[j, -1] - z[j, 0])
        quadrature = length * h*h * (2 * CURVATURE_MAX) / 12.0
        limit = (TAU_MAX * (abs(delta_i) + abs(delta_j))
                 + (2 * EPS) * (length + 2 * TAU_MAX)
                 + quadrature)
        broken = abs(integral) > limit + 1e-10
        all_tests.append((i + 1, j + 1, integral, limit, broken))
        if broken:
            flags.append((i, j))
    return flags, all_tests



def directional_certificates(t, z):
    """Stronger sign-aware interval certificate under the SAME physical assumptions.

    The unknown nonnegative response times and noisy endpoint changes give
    asymmetric feasible bounds on the signed integral. This test is sound,
    but does not guarantee detection of every faulty sensor.

    Returns (violating_edges, all_tests), where each test includes the
    observed integral, lower/upper feasible interval, and violation flag.
    """
    n = z.shape[0]
    length = float(t[-1] - t[0])
    h = float(t[1] - t[0])
    assert np.allclose(np.diff(t), h)
    error = (2 * EPS) * length + length * h*h * (2 * CURVATURE_MAX) / 12.0
    flags, all_tests = [], []
    for i, j in itertools.combinations(range(n), 2):
        integral = float(np.trapezoid(z[j] - z[i], t))
        delta_i = float(z[i, -1] - z[i, 0])
        delta_j = float(z[j, -1] - z[j, 0])
        lower_i = min(0.0, TAU_MAX * (delta_i - 2 * EPS))
        upper_i = max(0.0, TAU_MAX * (delta_i + 2 * EPS))
        lower_j = min(0.0, TAU_MAX * (delta_j - 2 * EPS))
        upper_j = max(0.0, TAU_MAX * (delta_j + 2 * EPS))
        lower = lower_i - upper_j - error
        upper = upper_i - lower_j + error
        violated = integral < lower - 1e-10 or integral > upper + 1e-10
        all_tests.append((i + 1, j + 1, integral, lower, upper, violated))
        if violated:
            flags.append((i, j))
    return flags, all_tests


def maximum_matching_size(n, edges):
    """Exhaustive DP for the tiny demonstration; use blossom at scale."""
    edge_set = {tuple(sorted(e)) for e in edges}
    from functools import lru_cache

    @lru_cache(None)
    def solve(available):
        if not available:
            return 0
        i = available[0]
        rest = available[1:]
        best = solve(rest)
        for j in rest:
            if (i, j) in edge_set:
                best = max(best, 1 + solve(tuple(x for x in rest if x != j)))
        return best

    return solve(tuple(range(n)))


def minimum_vertex_covers(n, edges):
    """Exhaustive, exact for six sensors; use ILP for large graphs."""
    for k in range(n + 1):
        covers = []
        for subset in itertools.combinations(range(n), k):
            selected = set(subset)
            if all(i in selected or j in selected for i, j in edges):
                covers.append(tuple(i + 1 for i in subset))
        if covers:
            return k, covers
    raise AssertionError('Every finite graph has a vertex cover')


def run_demo():
    t, z = readings(seed=123, faults={2: 1.2, 5: -1.5})
    edges, tests = certificates(t, z)
    n = len(TAUS)
    lower = maximum_matching_size(n, edges)
    min_size, min_covers = minimum_vertex_covers(n, edges)
    print('Contradictory pairs (1-indexed):', [(i + 1, j + 1) for i, j in edges])
    print('Contradiction edges:', len(edges), 'of', n*(n-1)//2)
    print('Maximum matching size (provable fault-count lower bound):', lower)
    print('Minimum vertex-cover size:', min_size)
    print('All minimum covers:', min_covers)
    print('Healthy-vs-healthy false alarms:', sum(
        1 for i, j in edges if i not in (2, 5) and j not in (2, 5)))
    # Require the claimed example to check out.
    assert len(edges) == 9, edges
    assert lower == 2
    assert min_size == 2 and min_covers == [(3, 6)], min_covers
    assert not any(i not in (2, 5) and j not in (2, 5) for i, j in edges)
    # Numerical soundness check for 100 separately generated noise streams.
    for seed in range(100):
        ts, healthy = readings(seed=seed)
        e, _ = certificates(ts, healthy)
        assert e == [], (seed, e)
    print('Additional synthetic healthy-array trials:', 100, '(all passed)')
    return tests


if __name__ == '__main__':
    run_demo()