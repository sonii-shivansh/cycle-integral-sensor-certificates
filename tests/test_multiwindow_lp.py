"""Tests of the relaxed multi-window consistency LP (not full physics proof)."""
import unittest
import numpy as np

from experiments.multiwindow_lp import (summarize_window, joint_relaxation,
                                        signaware_window_passes)


class MultiWindowTests(unittest.TestCase):
    def test_two_windows_pass_individually_but_inconsistent_constant_tau(self):
        # Exact first-order responses with T=0. Sensor 1 has a time-constant
        # switch from 0.5s to 2s; sensor 2 is an ideal instantaneous sensor.
        h = 0.01
        t1 = np.linspace(0, 0.5, 51)
        t2 = np.linspace(0.5, 1.5, 101)
        y1 = 10 * np.exp(-t1 / 0.5)
        y2 = y1[-1] * np.exp(-(t2 - 0.5) / 2.0)
        w1 = summarize_window(t1, np.vstack([y1, np.zeros(len(t1))]))
        w2 = summarize_window(t2, np.vstack([y2, np.zeros(len(t2))]))
        for w in (w1, w2):
            self.assertTrue(signaware_window_passes(w, 0, 1, 3.0, 0.0, [40, 0]))
            self.assertEqual(joint_relaxation([w], 3.0, 0.0, [40, 0]).status,
                             'feasible_relaxation')
        two = joint_relaxation([w1, w2], 3.0, 0.0, [40, 0])
        self.assertEqual(two.status, 'solver_infeasible')

    def test_healthy_constant_tau_and_zero_common_input_feasible(self):
        t = np.linspace(0, 1.5, 151)
        taus = [0.5, 1.2]
        z = np.array([10*np.exp(-t/tau) for tau in taus])
        summary = [summarize_window(t[a:b+1], z[:, a:b+1])
                   for a,b in [(0,50),(50,100),(100,150)]]
        r = joint_relaxation(summary, 2.0, 0.0, [40, 7])
        self.assertEqual(r.status, 'feasible_relaxation')
        self.assertLessEqual(r.max_constraint_violation, 1e-6)

    def test_large_exact_common_offset_does_not_create_false_conflicts(self):
        # Exact healthy response: T(t) = 2**52 + t, tau = (0,1,2).
        # All sample values are exactly representable binary64 integers;
        # eps = K = 0, so the physical ODE and bounds hold exactly.
        # Integrating each ~1e17 raw signal separately can lose increments
        # to floating-point rounding and falsely label the LP infeasible.
        t = np.arange(201, dtype=float)
        tau = np.array([0.0, 1.0, 2.0])
        y = 2.0**52 + t[None, :] - tau[:, None]
        windows = [summarize_window(t[a:b+1], y[:, a:b+1])
                   for a, b in [(0, 50), (50, 120), (120, 200)]]
        np.testing.assert_array_equal(windows[0].integrals,
                                      [1250.0, 1200.0, 1150.0])
        np.testing.assert_array_equal(windows[1].integrals,
                                      [2450.0, 2380.0, 2310.0])
        np.testing.assert_array_equal(windows[2].integrals,
                                      [3200.0, 3120.0, 3040.0])
        r = joint_relaxation(windows, 3.0, 0.0, 0.0)
        self.assertEqual(r.status, 'feasible_relaxation', r.message)
        np.testing.assert_allclose(r.time_constants, tau, atol=1e-8)

    def test_nonuniform_microsteps_use_valid_quadrature_bound(self):
        # Exact healthy first-order pair with a shared polynomial input.
        # tau_1=0, tau_2=s, and y_2+s*y_2' = y_1 = T.
        # An old atol=1e-12 allclose check wrongly accepted this very
        # nonuniform grid and used the *first* gap in the curvature bound.
        s = 2.0**-44
        x = np.array([0.0, 1.0, 2.0, 8.0])
        t = s * x
        amplitude = 1.0 / s
        z = amplitude * np.array([
            x**3 - 9*x**2 + 8*x,
            x**3 - 12*x**2 + 32*x - 32
        ])
        w = summarize_window(t, z)
        self.assertEqual(w.step, 6*s)
        self.assertAlmostEqual(w.integrals[1] - w.integrals[0], -109.0)
        # The true curvature bound is max(K_1,K_2)=30*amplitude/s**2.
        K = 30 * amplitude / s**2
        # A first-gap bound would be 40, falsely smaller than |Q|=109.
        self.assertTrue(signaware_window_passes(w, 0, 1, s, 0.0, K))
        self.assertEqual(joint_relaxation([w], s, 0.0, K).status,
                         'feasible_relaxation')

    def test_invalid_times_and_shapes_rejected(self):
        with self.assertRaises(ValueError):
            summarize_window([0,0,1], np.ones((2,3)))
        # Nonuniform increasing samples are allowed with a conservative
        # maximum-step trapezoid bound.
        w_irregular = summarize_window([0, 1, 3], np.ones((2, 3)))
        self.assertEqual(w_irregular.step, 2.0)
        with self.assertRaises(ValueError):
            summarize_window([0,1], [[1, np.nan]])
        w = summarize_window([0,1], [[0,1],[0,0]])
        with self.assertRaises(ValueError):
            joint_relaxation([w],2,-0.1,0)
        with self.assertRaises(ValueError):
            joint_relaxation([w],2,0,-0.1)
        with self.assertRaises(ValueError):
            joint_relaxation([w],2,0,0,indices=(1,1))


if __name__ == '__main__':
    unittest.main()
