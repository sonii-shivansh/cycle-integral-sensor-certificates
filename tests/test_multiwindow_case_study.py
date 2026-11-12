"""Small fixed-seed regression slice of the larger experimental comparison."""
import unittest
from experiments.multiwindow_case_study import run


class MultiWindowBenchmarkTests(unittest.TestCase):
    def test_matched_windows_three_methods(self):
        self.assertEqual(run(trials=5), [
            (0.0, False, 0, 0, 0, 0, 0, 0),
            (0.0, True, 3, 14, 45, 3, 5, 5),
            (0.3, False, 0, 0, 0, 0, 0, 0),
            (0.3, True, 9, 9, 45, 5, 5, 5),
            (0.6, False, 0, 0, 0, 0, 0, 0),
            (0.6, True, 0, 18, 45, 0, 5, 5),
            (1.0, False, 0, 0, 0, 0, 0, 0),
            (1.0, True, 0, 5, 32, 0, 5, 5),
        ])


if __name__ == '__main__':
    unittest.main()
