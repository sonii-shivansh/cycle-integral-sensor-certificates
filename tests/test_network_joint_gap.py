"""Exact triple-only inconsistency versus pairwise consistency.

This test is about the relaxed LP, not a new theorem or physical fault labels.
"""
import unittest
from itertools import combinations

from experiments.exact_pairwise_certificate import exact_pairwise_verdict
from experiments.network_consistency_gap import (
    exact_network_contradiction, exact_summaries, float_summaries, evaluate,
)
from experiments.network_conflicts import (
    network_consistency, inclusion_minimal_conflict,
)


class NetworkGapTests(unittest.TestCase):
    def test_every_pair_feasible_but_three_sensor_network_infeasible(self):
        exact = exact_summaries()
        for pair in combinations(range(3), 2):
            with self.subTest(pair=pair):
                self.assertEqual(
                    exact_pairwise_verdict(exact, pair, 2, 0, 0).status,
                    "feasible_relaxation",
                )
        z = float_summaries()
        self.assertEqual(network_consistency(z,2,0,0).status,
                         "solver_infeasible")
        self.assertEqual(inclusion_minimal_conflict(z,2,0,0),(0,1,2))
        result = evaluate()
        self.assertEqual(len(result["pairwise_tau_witnesses"]),3)

    def test_exact_farkas_network_contradiction(self):
        rows, (lhs, rhs) = exact_network_contradiction()
        self.assertEqual(lhs, (0,0,0,0,0))
        self.assertEqual(rhs, -1)
        self.assertEqual(len(rows), 4)
        # All multipliers are exactly +1, so the rows yield 0 <= -1.
        self.assertTrue(all(row_rhs.denominator == 2 or
                            row_rhs.denominator == 1 for _,row_rhs in rows))

    def test_nonconflicting_network_yields_no_minimal_conflict(self):
        from experiments.multiwindow_lp import summarize_window
        # Identical constant sensor responses obey the shared-input model.
        windows = [summarize_window([0,1], [[1,1],[1,1],[1,1]]),
                   summarize_window([3,4], [[2,2],[2,2],[2,2]])]
        self.assertEqual(network_consistency(windows,2,0,0).status,
                         "feasible_relaxation")
        self.assertEqual(inclusion_minimal_conflict(windows,2,0,0),())

    def test_pair_conflict_uses_only_two_sensors(self):
        from experiments.multiwindow_lp import summarize_window
        summary = [summarize_window([0,1], [[0,0],[5,5],[0,0]])]
        self.assertEqual(inclusion_minimal_conflict(summary,0,0,0),(1,2))


if __name__ == "__main__":
    unittest.main()
