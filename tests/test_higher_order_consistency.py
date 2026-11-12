"""A graph of feasible PAIRS need not imply feasibility of the whole network.

A constructive three-sensor, two-window counterexample to pairwise completeness.
This is a property of LP feasibility / set-intersection geometry, not a new
scientific theorem or hardware result.
"""

from itertools import combinations
import unittest

import numpy as np

from experiments.multiwindow_lp import summarize_window, joint_relaxation


class HigherOrderConsistencyTests(unittest.TestCase):
    def test_three_way_conflict_with_all_pairs_feasible(self):
        # Each observation window has linear (zero-curvature) traces.
        # The 10-second gap permits smooth interpolation between the windows;
        # the check itself uses only the two observed windows.
        #
        # Window A: I=(0,0,1), d=(1,0,-1).
        # Window B: I=(0,0,0), d=(0,1,1).
        # With eps=0, maximum tau=1 and K=1, the exact model relaxation
        # has pairwise feasible parameter choices but NO joint choice.
        windows = []
        for start in (0.0, 11.0):
            t = np.linspace(start, start + 1.0, 101)
            s = t - start
            if start == 0.0:
                z = np.vstack([s - 0.5, np.zeros_like(s), 1.5 - s])
            else:
                z = np.vstack([np.zeros_like(s), s - 0.5, s - 0.5])
            windows.append(summarize_window(t, z))

        # All three pairwise tests pass: the pairwise contradiction graph has
        # no edges and a vertex cover / matching lower bound of zero.
        for i, j in combinations(range(3), 2):
            with self.subTest(pair=(i, j)):
                r = joint_relaxation(windows, 1.0, 0.0, 1.0, indices=(i, j))
                self.assertEqual(r.status, "feasible_relaxation", r.message)

        # Nevertheless no single set of three tau values and common
        # integrated inputs satisfies all three sensors simultaneously.
        network = joint_relaxation(windows, 1.0, 0.0, 1.0)
        self.assertEqual(network.status, "solver_infeasible", network.message)


if __name__ == "__main__":
    unittest.main()
