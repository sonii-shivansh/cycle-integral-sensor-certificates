"""Independent regression tests for general network Farkas proof extraction.

The exact verifier can prove SOME infeasible rational LPs but may return
'not_certified' even when they are mathematically infeasible.
"""
import unittest
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog

from experiments.exact_pairwise_certificate import (
    exact_window, exact_pairwise_verdict
)
from experiments.network_consistency_gap import observations
from experiments.network_exact_certificate import (
    certify_network_infeasibility, exact_network_rows,
    summarize_network_certificate,
)
from experiments.set_membership_comparison import recorded_windows
from experiments.multiwindow_case_study import K


class GeneralNetworkProofTests(unittest.TestCase):
    def test_three_way_inconsistency_with_every_pair_feasible(self):
        windows = [exact_window(t, z) for t, z in observations()]
        for i, j in [(0, 1), (0, 2), (1, 2)]:
            self.assertEqual(
                exact_pairwise_verdict(windows, (i, j), 2, 0, 0).status,
                "feasible_relaxation"
            )
        result = certify_network_infeasibility(windows, 2, 0, 0)
        self.assertEqual(result.status, "certified_infeasible", result.message)
        self.assertTrue(result.certificate.check(result.rows))
        self.assertLess(result.certificate.weighted_rhs, 0)
        summary = summarize_network_certificate(result)
        self.assertEqual(summary["status"], "certified_infeasible")
        self.assertTrue(summary["nonzero_weights"])
        # Altering a true multiplier invalidates the Farkas proof.
        from dataclasses import replace
        weights = list(result.certificate.weights)
        nonzero = next(j for j, value in enumerate(weights) if value > 0)
        weights[nonzero] = -1
        fake = replace(result.certificate, weights=tuple(weights))
        self.assertFalse(fake.check(result.rows))

    def test_high_offset_healthy_network_not_falsely_certified(self):
        t = np.arange(201, dtype=float)
        tau = np.array([0.0, 1.0, 2.0])
        z = 2.0**52 + t[None, :] - tau[:, None]
        windows = [exact_window(t[a:b+1], z[:, a:b+1])
                   for a, b in [(0, 50), (50, 120), (120, 200)]]
        result = certify_network_infeasibility(windows, 3, 0, 0)
        self.assertEqual(result.status, "not_certified", result.message)
        self.assertIsNone(result.certificate)

    def test_switching_time_constant_has_exact_network_proof(self):
        t1 = np.linspace(0, .5, 51)
        t2 = np.linspace(.5, 1.5, 101)
        y1 = 10 * np.exp(-t1/.5)
        y2 = y1[-1] * np.exp(-(t2-.5)/2.)
        windows = [
            exact_window(t1, np.vstack([y1, np.zeros(len(t1))])),
            exact_window(t2, np.vstack([y2, np.zeros(len(t2))]))
        ]
        for window in windows:
            result = certify_network_infeasibility([window], 3, 0, [40, 0])
            self.assertEqual(result.status, "not_certified")
        result = certify_network_infeasibility(windows, 3, 0, [40, 0])
        self.assertEqual(result.status, "certified_infeasible", result.message)
        self.assertTrue(result.certificate.check(result.rows))

    def test_six_sensor_bias_network_cases_are_certified(self):
        for drift in [0.0, 0.3, 0.6, 1.0]:
            with self.subTest(drift=drift):
                healthy = recorded_windows(0, drift, False)
                self.assertEqual(
                    certify_network_infeasibility(
                        healthy, 15, .02, K
                    ).status, "not_certified"
                )
                injected = recorded_windows(0, drift, True)
                result = certify_network_infeasibility(
                    injected, 15, .02, K
                )
                self.assertEqual(result.status, "certified_infeasible",
                                 result.message)
                self.assertTrue(result.certificate.check(result.rows))

    def test_35_random_small_integer_lp_against_independent_primal(self):
        rng = np.random.default_rng(20261112)
        for trial in range(35):
            with self.subTest(trial=trial):
                w1 = exact_window([0, 1], rng.integers(-6, 7, (4, 2)))
                w2 = exact_window([3, 4], rng.integers(-6, 7, (4, 2)))
                windows = [w1, w2]
                result = certify_network_infeasibility(
                    windows, 3, F(1, 10), 0
                )
                rows = exact_network_rows(windows, 3, F(1, 10), 0)
                dimension = len(rows[0].coefficients)
                oracle = linprog(
                    [0.0] * dimension,
                    A_ub=[[float(c) for c in r.coefficients] for r in rows],
                    b_ub=[float(r.rhs) for r in rows],
                    bounds=[(None, None)]*dimension, method="highs"
                )
                if result.status == "certified_infeasible":
                    self.assertEqual(oracle.status, 2)
                    self.assertTrue(result.certificate.check(result.rows))
                else:
                    # This is a conservative verifier: it may decline to
                    # prove infeasibility and must never announce feasibility.
                    self.assertEqual(result.status, "not_certified")

    def test_tiny_contradiction_is_not_falsely_claimed_certified(self):
        # The normalized floating dual objective is below the heuristic
        # threshold. The old exact pairwise solver can still prove this.
        windows = [exact_window([0, 1], [[0, 0], [1e-12, 1e-12]])]
        old = exact_pairwise_verdict(windows, (0, 1), 0, 0, 0)
        self.assertEqual(old.status, "exact_infeasible")
        new = certify_network_infeasibility(windows, 0, 0, 0)
        self.assertEqual(new.status, "not_certified")
        self.assertIsNone(new.certificate)

    def test_invalid_bounds_and_indices_rejected(self):
        windows = [exact_window([0, 1], [[0, 0], [1, 1]])]
        with self.assertRaises(ValueError):
            certify_network_infeasibility(windows, -1, 0, 0)
        with self.assertRaises(ValueError):
            certify_network_infeasibility(windows, 1, -.1, 0)
        with self.assertRaises(ValueError):
            certify_network_infeasibility(windows, 1, 0, 0, (0, 0))


if __name__ == "__main__":
    unittest.main()
