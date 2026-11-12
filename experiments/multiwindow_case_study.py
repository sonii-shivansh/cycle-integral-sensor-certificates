"""Fixed-parameter six-sensor comparison across identical observation windows.

Synthetic test only. Faults are injected biases, not physically defective
hardware; reported contradictions do not establish correct identification.
"""
from itertools import combinations

import numpy as np

from experiments.multiwindow_lp import (summarize_window,
                                        signaware_window_passes,
                                        pairwise_results)

TAUS = np.array([1., 2., 3., 5., 8., 13.])
EPS = 0.02
K = 8 * (2*np.pi/20)**2
M = 15.
STEP = .05


def generate(seed=0, drift=.3, faults=True):
    t = np.linspace(0, 13.5, 271)
    outputs = []
    for tau in TAUS:
        p1, p2 = .25*t + .3, .4*t-.2
        outputs.append(25 + drift * (t-tau)
                       + 6 * (np.sin(p1) - .25*tau*np.cos(p1))/(1+(.25*tau)**2)
                       + 2 * (np.cos(p2) + .4*tau*np.sin(p2))/(1+(.4*tau)**2))
    z = np.asarray(outputs)
    if faults:
        z[2] += 1.2
        z[5] -= 1.5
    z += np.random.default_rng(seed).uniform(-EPS, EPS, z.shape)
    windows = [(0,60), (60,120), (120,180), (180,270)]
    return [summarize_window(t[a:b+1], z[:,a:b+1]) for a,b in windows]


def symmetric_passes(w, i, j):
    q = w.integrals[j] - w.integrals[i]
    d = w.increments
    limit = M*(abs(d[i])+abs(d[j])) + 2*EPS*(w.length+2*M) + w.length*w.step**2*2*K/12
    return abs(q) <= limit


def evaluate(windows):
    sym, directional, lp = [], [], []
    all_lp = pairwise_results(windows, M, EPS, K)
    for i,j in combinations(range(len(TAUS)), 2):
        if any(not symmetric_passes(w,i,j) for w in windows):
            sym.append((i,j))
        if any(not signaware_window_passes(w,i,j,M,EPS,K) for w in windows):
            directional.append((i,j))
        if all_lp[(i,j)].status == 'solver_infeasible':
            lp.append((i,j))
        assert all_lp[(i,j)].status != 'inconclusive'
    assert set(sym) <= set(directional) <= set(lp)
    return sym, directional, lp


def run(trials=20):
    rows=[]
    for drift in (0., .3, .6, 1.):
        for faults in (False, True):
            totals=np.zeros(3, dtype=int)
            trial_counts=np.zeros(3, dtype=int)
            for seed in range(trials):
                edges=evaluate(generate(seed,drift,faults))
                totals += [len(e) for e in edges]
                trial_counts += [bool(e) for e in edges]
            rows.append((drift, faults, *map(int, totals), *map(int, trial_counts)))
    return rows


if __name__ == '__main__':
    print('drift, faults, sum_edges_original, sum_edges_directional, sum_edges_LP,'
          ' trials_original, trials_directional, trials_LP')
    for r in run():
        print(*r,sep=', ')
