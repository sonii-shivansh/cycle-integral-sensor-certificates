"""Randomized algebraic regression checks for the sign-aware certificate.

Randomized simulation is not a substitute for a proof. Test signals are exact
closed-form first-order solutions and have a guaranteed curvature envelope.
"""

import unittest

import numpy as np

import cycle_integral_fault_certificate as cifc


class DirectionalPropertyTests(unittest.TestCase):
    def test_250_randomized_healthy_ensembles_including_boundary_time_constants(self):
        rng = np.random.default_rng(20261111)
        t = np.linspace(0.0, 17.0, 341)
        for trial in range(250):
            with self.subTest(trial=trial):
                taus = rng.uniform(0.0, cifc.TAU_MAX, size=6)
                taus[0] = 0.0
                taus[-1] = cifc.TAU_MAX
                w1, w2 = rng.uniform(0.1, 0.4, size=2)
                a1, a2 = rng.uniform(0.2, 4.0, size=2)
                curvature = a1 * w1*w1 + a2 * w2*w2
                if curvature > 0.7:
                    a1 *= 0.7 / curvature
                    a2 *= 0.7 / curvature
                self.assertLessEqual(
                    a1 * w1*w1 + a2 * w2*w2, cifc.CURVATURE_MAX
                )
                p1, p2 = rng.uniform(-3.0, 3.0, size=2)
                drift = rng.uniform(-1.0, 1.0)
                y = []
                for tau in taus:
                    one, two = w1*t+p1, w2*t+p2
                    y.append(
                        25.0 + drift*(t-tau)
                        + a1*(np.sin(one)-w1*tau*np.cos(one))/(1.0+(w1*tau)**2)
                        + a2*(np.cos(two)+w2*tau*np.sin(two))/(1.0+(w2*tau)**2)
                    )
                y = np.asarray(y)
                # Adversarial-sign bounded sample errors (the noise bound is tight).
                z = y + rng.choice([-cifc.EPS, cifc.EPS], size=y.shape)
                old, _ = cifc.certificates(t, z)
                new, _ = cifc.directional_certificates(t, z)
                self.assertEqual(old, [])
                self.assertEqual(new, [])


if __name__ == "__main__":
    unittest.main()
