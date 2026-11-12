"""Reproducible *illustrative* network-size benchmark for the CIFC LP relaxation.

This compares (a) ONE global LP detecting any network inconsistency with
(b) ALL n*(n-1)/2 pairwise LPs returning an entire contradiction graph.
These solve DIFFERENT reporting tasks; timings must not be interpreted as
a fair speedup for identical functionality. No formal worst-case timing claim.

The signals satisfy a shared-input first-order model with bounded sample
errors. Timings depend on Python, NumPy, SciPy, HiGHS and hardware.
"""

from itertools import combinations
from statistics import median
from time import perf_counter
import sys

import numpy as np
import scipy

from experiments.multiwindow_lp import summarize_window, joint_relaxation, pairwise_results


M = 15.0
EPS = 0.02
K = 8 * (2*np.pi/20)**2


def healthy_windows(n, seed=11):
    if not isinstance(n,int) or n < 2:
        raise ValueError("n must be an integer >= 2")
    t = np.linspace(0.0,13.5,271)
    tau_values = np.linspace(0.5,14.5,n)
    ys = []
    for tau in tau_values:
        p1, p2 = 0.25*t+0.3, 0.4*t-0.2
        ys.append(
            25 + .3*(t-tau)
            + 6*(np.sin(p1)-.25*tau*np.cos(p1))/(1+(.25*tau)**2)
            + 2*(np.cos(p2)+.4*tau*np.sin(p2))/(1+(.4*tau)**2)
        )
    z = np.asarray(ys) + np.random.default_rng(seed).uniform(-EPS,EPS,(n,len(t)))
    # Disjoint interiors; adjacent windows share the boundary node.
    return [summarize_window(t[a:b+1],z[:,a:b+1])
            for a,b in ((0,60),(60,120),(120,180),(180,270))]


def benchmark_one(n, repetitions=2):
    windows = healthy_windows(n)
    full_durations, pair_durations = [], []
    for _ in range(repetitions):
        start = perf_counter()
        network = joint_relaxation(windows,M,EPS,K)
        full_durations.append(perf_counter()-start)
        if network.status != "feasible_relaxation":
            raise AssertionError("valid healthy synthetic model rejected by global LP")

        start = perf_counter()
        pairs = pairwise_results(windows,M,EPS,K)
        pair_durations.append(perf_counter()-start)
        if any(v.status != "feasible_relaxation" for v in pairs.values()):
            raise AssertionError("valid healthy synthetic model rejected by pairwise LP")
    return {
        "n":n, "windows":len(windows),
        "pairwise_calls":len(pairs),
        "network_calls":1,
        "network_median_seconds":round(median(full_durations),6),
        "pairwise_all_median_seconds":round(median(pair_durations),6),
    }


def run(sizes=(6,12,24,48),repetitions=2):
    return [benchmark_one(int(n), repetitions=repetitions) for n in sizes]


if __name__ == "__main__":
    print("Python:",sys.version.split()[0],"NumPy:",np.__version__,
          "SciPy:",scipy.__version__)
    print("NOTE: global decision and pairwise edge graph are different output tasks.")
    print("n, W, network LP calls, pairwise LP calls, global seconds, all-pairs seconds")
    for result in run():
        print(result["n"],result["windows"],result["network_calls"],
              result["pairwise_calls"],result["network_median_seconds"],
              result["pairwise_all_median_seconds"],sep=", ")
