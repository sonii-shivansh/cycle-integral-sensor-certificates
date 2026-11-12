"""Network scaling structure without fragile wall-clock assertions."""
import unittest

from experiments.network_scaling_benchmark import benchmark_one, healthy_windows
from experiments.multiwindow_lp import joint_relaxation


class NetworkScalingTests(unittest.TestCase):
    def test_healthy_network_24_sensor_global(self):
        windows = healthy_windows(24, seed=11)
        result = joint_relaxation(windows, 15, .02, 8*(2*3.141592653589793/20)**2)
        self.assertEqual(result.status,"feasible_relaxation")

    def test_count_and_deterministic_healthy_feasibility(self):
        for n in (3,6,12):
            with self.subTest(n=n):
                report=benchmark_one(n,repetitions=1)
                self.assertEqual(report["n"],n)
                self.assertEqual(report["pairwise_calls"],n*(n-1)//2)
                self.assertEqual(report["network_calls"],1)
                self.assertEqual(report["windows"],4)
                self.assertGreaterEqual(report["network_median_seconds"],0)
                self.assertGreaterEqual(report["pairwise_all_median_seconds"],0)

    def test_bad_sensor_count_rejected(self):
        with self.assertRaises(ValueError):
            healthy_windows(1)


if __name__=="__main__":
    unittest.main()
