"""Independent comparison against classical invariant-parameter feasible-set tests.

This is a scientific *negative result*: identical feasible-set constraints
must yield the same exact detection decisions, whichever solver is used.
"""
import unittest
from itertools import combinations
from fractions import Fraction as F

from experiments.exact_pairwise_certificate import exact_window, exact_pairwise_verdict
from experiments.set_membership_fss_baseline import intersect_parameter_strips, clip
from experiments.set_membership_comparison import recorded_windows, compare


class SetMembershipBaselineTests(unittest.TestCase):
    def test_full_20_seed_benchmark_matches_same_fss(self):
        self.assertEqual(compare(trials=20), [
            (0.0, False, 0, 0, 0, 0),
            (0.0, True, 16, 58, 180, 180),
            (0.3, False, 0, 0, 0, 0),
            (0.3, True, 35, 35, 180, 180),
            (0.6, False, 0, 0, 0, 0),
            (0.6, True, 0, 74, 180, 180),
            (1.0, False, 0, 0, 0, 0),
            (1.0, True, 0, 20, 133, 133),
        ])

    def test_geometry_independently_matches_exact_farkas_verifier(self):
        for drift in (0.0, 0.3, 1.0):
            for fault in (False, True):
                windows = recorded_windows(3, drift, fault)
                for pair in combinations(range(6), 2):
                    with self.subTest(drift=drift, fault=fault, pair=pair):
                        polygon = intersect_parameter_strips(
                            windows, pair, 15, .02, 8*(2*3.141592653589793/20)**2
                        )
                        farkas = exact_pairwise_verdict(
                            windows, pair, 15, .02, 8*(2*3.141592653589793/20)**2
                        )
                        self.assertEqual(polygon.status, farkas.status)
                        if polygon.status == "exact_infeasible":
                            self.assertTrue(farkas.certificate.check(farkas.inequalities))

    def test_healthy_high_offset_and_nonuniform_inputs(self):
        import numpy as np
        t = np.arange(201, dtype=float)
        z = 2.0**52 + t[None,:] - np.array([0,1,2])[:,None]
        windows = [exact_window(t[a:b+1], z[:,a:b+1])
                   for a,b in ((0,50),(50,120),(120,200))]
        for pair in ((0,1),(0,2),(1,2)):
            self.assertEqual(
                intersect_parameter_strips(windows, pair, 3, 0, 0).status,
                "feasible_relaxation"
            )
        s = 2.0**-44
        x = np.array([0., 1., 2., 8.])
        nonuniform = exact_window(s*x, [
            (x**3-9*x*x+8*x)/s,
            (x**3-12*x*x+32*x-32)/s
        ])
        self.assertEqual(
            intersect_parameter_strips([nonuniform], (0,1), s, 0, 30/s**3).status,
            "feasible_relaxation"
        )

    def test_convex_polygon_tangent_and_empty(self):
        square = ((F(0),F(0)),(F(2),F(0)),(F(2),F(2)),(F(0),F(2)))
        # A line that touches only the upper-right vertex remains feasible.
        tangent = clip(square, F(-1), F(-1), F(-4))
        self.assertEqual(tangent, ((F(2),F(2)),))
        # Excluding even the tangent point makes it empty.
        self.assertEqual(clip(square, F(-1), F(-1), F(-5)), ())

    def test_invalid_inputs_rejected(self):
        w = exact_window([0,1], [[0,0],[0,0]])
        for kwargs in [
            {"pair":(0,0), "max_tau":1, "eps":0, "curvature":0},
            {"pair":(0,2), "max_tau":1, "eps":0, "curvature":0},
            {"pair":(0,1), "max_tau":-1, "eps":0, "curvature":0},
            {"pair":(0,1), "max_tau":1, "eps":-1, "curvature":0},
        ]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                intersect_parameter_strips([w], **kwargs)


if __name__ == "__main__":
    unittest.main()
