"""Adversarial regression checks for exact two-variable CIFC feasibility.

These establish algebraic correctness for *represented rational LP inputs*.
They are not evidence of real-world sensor accuracy or novelty.
"""
import unittest
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog

from experiments.exact_pairwise_certificate import (
    FarkasCertificate, Inequality, _eliminate, exact_window,
    exact_pairwise_verdict,
)


class ExactPairwiseTests(unittest.TestCase):
    def test_contradiction_produces_exact_farkas_certificate(self):
        w = exact_window([0, 1], [[0, 1], [2, 2]])
        v = exact_pairwise_verdict([w], (0, 1), 0, 0, 0)
        self.assertEqual(v.status, "exact_infeasible")
        self.assertTrue(v.certificate.check(v.inequalities))
        self.assertLess(v.certificate.combined_rhs, 0)

    def test_large_common_baseline_healthy(self):
        t = np.arange(201, dtype=float)
        tau = np.array([0.0, 1.0, 2.0])
        z = 2.0**52 + t[None, :] - tau[:, None]
        windows = [exact_window(t[a:b+1], z[:, a:b+1])
                   for a, b in [(0, 50), (50, 120), (120, 200)]]
        self.assertEqual(windows[0].integrals, tuple(map(F, [1250, 1200, 1150])))
        for pair in [(0, 1), (1, 2), (0, 2)]:
            self.assertEqual(
                exact_pairwise_verdict(windows, pair, 3, 0, 0).status,
                "feasible_relaxation",
            )

    def test_extremely_nonuniform_grid_healthy(self):
        s = 2.0**-44
        x = np.array([0.0, 1.0, 2.0, 8.0])
        t = s*x
        A = 1/s
        y1 = A*(x**3-9*x**2+8*x)
        y2 = A*(x**3-12*x**2+32*x-32)
        w = exact_window(t, [y1, y2])
        self.assertEqual(w.max_gap, F.from_float(6*s))
        self.assertEqual(w.integrals[1]-w.integrals[0], F(-109))
        self.assertEqual(
            exact_pairwise_verdict([w], (0, 1), s, 0, 30/s**3).status,
            "feasible_relaxation",
        )

    def test_changing_time_constant_joint_conflict_exactly_proved(self):
        t1 = np.linspace(0, 0.5, 51)
        t2 = np.linspace(0.5, 1.5, 101)
        y1 = 10*np.exp(-t1/0.5)
        y2 = y1[-1]*np.exp(-(t2-0.5)/2.0)
        w1 = exact_window(t1, np.vstack([y1, np.zeros(len(t1))]))
        w2 = exact_window(t2, np.vstack([y2, np.zeros(len(t2))]))
        for w in (w1, w2):
            self.assertEqual(exact_pairwise_verdict([w], (0, 1), 3, 0, [40, 0]).status,
                             "feasible_relaxation")
        result = exact_pairwise_verdict([w1, w2], (0, 1), 3, 0, [40, 0])
        self.assertEqual(result.status, "exact_infeasible")
        self.assertTrue(result.certificate.check(result.inequalities))

    def test_exact_2d_fourier_motzkin_matches_100_small_lp_oracles(self):
        rng = np.random.default_rng(12345)
        for trial in range(100):
            with self.subTest(trial=trial):
                m = 7
                inequalities = [
                    Inequality(F(int(a)), F(int(b)), F(int(c)), f"r{k}")
                    for k, (a,b,c) in enumerate(rng.integers(-5, 6, (m,3)))
                ]
                inequalities.extend([
                    Inequality(F(-1), F(0), F(0), "x>=0"),
                    Inequality(F(1), F(0), F(4), "x<=4"),
                    Inequality(F(0), F(-1), F(0), "y>=0"),
                    Inequality(F(0), F(1), F(4), "y<=4"),
                ])
                count = len(inequalities)
                rows = [
                    ((p.a, p.b), p.rhs, tuple(F(int(k == j)) for k in range(count)))
                    for j, p in enumerate(inequalities)
                ]
                for variable in (0, 1):
                    rows = _eliminate(rows, variable)
                conflicts = [(c, rhs, weight) for c, rhs, weight in rows if rhs < 0]
                oracle = linprog(
                    [0, 0],
                    A_ub=[[float(p.a), float(p.b)] for p in inequalities],
                    b_ub=[float(p.rhs) for p in inequalities],
                    bounds=[(None, None), (None, None)],
                    method="highs",
                )
                if conflicts:
                    cert = FarkasCertificate(conflicts[0][2], conflicts[0][1])
                    self.assertTrue(cert.check(inequalities))
                    self.assertEqual(oracle.status, 2)
                else:
                    self.assertTrue(all(c == (0, 0) and rhs >= 0 for c, rhs, _ in rows))
                    self.assertEqual(oracle.status, 0)

    def test_six_sensor_example_exact_edge_counts_one_seed(self):
        # Independently regenerate raw recorded samples, not floating-point LP
        # summaries; test exactly the four data windows used in PR #4.
        from itertools import combinations
        taus = np.array([1., 2., 3., 5., 8., 13.])
        t = np.linspace(0, 13.5, 271)
        K = 8*(2*np.pi/20)**2
        for drift, expected in [(0., 9), (0.3, 9), (0.6, 9), (1., 7)]:
            with self.subTest(drift=drift):
                outputs = []
                for tau in taus:
                    p1,p2 = .25*t+.3, .4*t-.2
                    outputs.append(
                        25+drift*(t-tau)
                        + 6*(np.sin(p1)-.25*tau*np.cos(p1))/(1+(.25*tau)**2)
                        + 2*(np.cos(p2)+.4*tau*np.sin(p2))/(1+(.4*tau)**2)
                    )
                z = np.asarray(outputs)
                z[2] += 1.2
                z[5] -= 1.5
                z += np.random.default_rng(0).uniform(-0.02, 0.02, z.shape)
                windows = [exact_window(t[a:b+1], z[:,a:b+1])
                           for a,b in [(0,60),(60,120),(120,180),(180,270)]]
                flags = [
                    (i,j) for i,j in combinations(range(6), 2)
                    if exact_pairwise_verdict(windows,(i,j),15.,0.02,K).status
                       == "exact_infeasible"
                ]
                self.assertEqual(len(flags), expected)


if __name__ == "__main__":
    unittest.main()
