"""Independent synthetic checks for the directional consistency certificate.

These are mathematical/software regression tests, not hardware validation.
"""

import unittest

import cycle_integral_fault_certificate as cifc
from experiments import nonperiodic_benchmark as benchmark


class NonperiodicCertificateTests(unittest.TestCase):
    def test_exact_response_quadrature_stays_within_bound(self):
        worst, bound = benchmark.validate_integral_identity()
        self.assertLessEqual(worst, bound)

    def test_directional_bound_never_loses_original_contradictions(self):
        for drift in (0.0, 0.3, 0.6, 1.0):
            for seed in (0, 7, 42, 99):
                for faults in (None, {2: 1.2, 5: -1.5}):
                    t, z = benchmark.generate(seed, offsets=faults, drift=drift)
                    original, _ = cifc.certificates(t, z)
                    new, _ = cifc.directional_certificates(t, z)
                    self.assertLessEqual(set(original), set(new))

    def test_nonperiodic_100_seed_benchmark(self):
        rows = benchmark.compute_results(num_trials=100)
        self.assertEqual(rows, [
            (0.0, 0, 0, 200, 500, 100, 100),
            (0.3, 0, 0, 0, 218, 0, 100),
            (0.6, 0, 0, 0, 0, 0, 0),
            (1.0, 0, 0, 0, 0, 0, 0),
        ])

    def test_graph_explanation_can_understate_true_injected_bias_count(self):
        # With 0.3 degC/s drift, both sensors 3 and 6 receive offsets, but
        # the seed-123 contradiction graph has a size-one minimum cover.
        t, z = benchmark.generate(123, offsets={2: 1.2, 5: -1.5}, drift=0.3)
        edges, _ = cifc.directional_certificates(t, z)
        self.assertEqual(cifc.maximum_matching_size(6, edges), 1)
        self.assertEqual(cifc.minimum_vertex_covers(6, edges), (1, [(6,)]))


if __name__ == "__main__":
    unittest.main()
