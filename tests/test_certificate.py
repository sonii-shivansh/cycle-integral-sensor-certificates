"""Small deterministic regression tests for the CIFC demonstration.

These tests exercise the implementation; they are not experimental validation
of hardware performance or a substitute for independent mathematical review.
"""

import unittest

import numpy as np

import cycle_integral_fault_certificate as cifc


class CertificateTests(unittest.TestCase):
    def test_healthy_synthetic_arrays_have_no_conflicts(self):
        for seed in (0, 1, 42, 99):
            with self.subTest(seed=seed):
                t, z = cifc.readings(seed=seed)
                edges, _ = cifc.certificates(t, z)
                self.assertEqual(edges, [])

    def test_original_injected_bias_example(self):
        t, z = cifc.readings(seed=123, faults={2: 1.2, 5: -1.5})
        edges, reports = cifc.certificates(t, z)
        self.assertEqual(len(reports), 15)
        self.assertEqual(len(edges), 9)
        self.assertEqual(cifc.maximum_matching_size(6, edges), 2)
        self.assertEqual(cifc.minimum_vertex_covers(6, edges), (2, [(3, 6)]))
        self.assertTrue(all(2 in e or 5 in e for e in edges))

    def test_common_offset_is_undetectable_by_pairwise_differences(self):
        t, baseline = cifc.readings(seed=7)
        offsets = {i: 1.0 for i in range(len(cifc.TAUS))}
        _, shifted = cifc.readings(seed=7, faults=offsets)
        np.testing.assert_allclose(shifted - baseline, 1.0, atol=1e-12)
        edges_baseline, _ = cifc.certificates(t, baseline)
        edges_shifted, _ = cifc.certificates(t, shifted)
        self.assertEqual(edges_baseline, edges_shifted)

    def test_simple_conflict_graph_lower_bound(self):
        edges = [(0, 1), (1, 2), (2, 3)]
        self.assertEqual(cifc.maximum_matching_size(4, edges), 2)
        size, covers = cifc.minimum_vertex_covers(4, edges)
        self.assertEqual(size, 2)
        self.assertIn((1, 3), covers)


if __name__ == "__main__":
    unittest.main()
