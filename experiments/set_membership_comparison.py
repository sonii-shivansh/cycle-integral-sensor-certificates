"""Fair, exact set-membership baseline comparison on the same six-sensor samples.

The geometric solver independently implements *classical invariant-parameter
feasible-set intersection*. It is NOT the parallelotope approximation
algorithm from Ingimundarson et al. (CDC/ECC 2005). The two implementations
solve the SAME relaxed parameter-feasibility problem, so equal decisions
are expected. This benchmark cannot establish an advantage for CIFC.

Observation values are treated as exact binary64 rationals for the geometric
test. The existing LP uses SciPy/HiGHS in floating point.
"""

from itertools import combinations

import numpy as np

from experiments.exact_pairwise_certificate import exact_window
from experiments.multiwindow_case_study import evaluate, generate, TAUS, EPS, K, M
from experiments.set_membership_fss_baseline import intersect_parameter_strips


def recorded_windows(seed, drift, faults):
    """Independently reconstruct exact recorded windows from the benchmark."""
    t = np.linspace(0.0, 13.5, 271)
    ys = []
    for tau in TAUS:
        p1, p2 = .25*t + .3, .4*t-.2
        ys.append(
            25 + drift*(t-tau)
            + 6*(np.sin(p1)-.25*tau*np.cos(p1))/(1+(.25*tau)**2)
            + 2*(np.cos(p2)+.4*tau*np.sin(p2))/(1+(.4*tau)**2)
        )
    z = np.asarray(ys)
    if faults:
        z[2] += 1.2
        z[5] -= 1.5
    z += np.random.default_rng(seed).uniform(-EPS, EPS, z.shape)
    return [exact_window(t[a:b+1], z[:,a:b+1])
            for a,b in [(0,60), (60,120), (120,180), (180,270)]]


def compare(trials=20):
    """Return (drift, biases, symmetric, sign-aware, fast-LP, FSS-baseline).

    Each count totals pairwise conflict *edges* across 'trials' seeds.
    FSS and fast-LP must agree on this constructed data; disagreement raises.
    """
    rows = []
    for drift in (0.0, 0.3, 0.6, 1.0):
        for faults in (False, True):
            counts = np.zeros(4, dtype=int)
            for seed in range(trials):
                standard_summaries = generate(seed, drift, faults)
                symmetric, signaware, fast_lp = evaluate(standard_summaries)
                raw_windows = recorded_windows(seed, drift, faults)
                fss = [
                    (i,j) for i,j in combinations(range(len(TAUS)), 2)
                    if intersect_parameter_strips(raw_windows, (i,j), M, EPS, K).status
                    == "exact_infeasible"
                ]
                if set(fss) != set(fast_lp):
                    raise AssertionError({
                        "seed":seed, "drift":drift, "faults":faults,
                        "fss_only":sorted(set(fss)-set(fast_lp)),
                        "fast_lp_only":sorted(set(fast_lp)-set(fss)),
                    })
                counts += [len(symmetric), len(signaware), len(fast_lp), len(fss)]
            rows.append((drift, faults, *map(int, counts)))
    return rows


if __name__ == "__main__":
    print("drift_degC_s, biases_injected, symmetric_edges, signaware_edges,"
          " multiwindow_lp_edges, classical_fss_edges")
    for row in compare():
        print(*row, sep=", ")
