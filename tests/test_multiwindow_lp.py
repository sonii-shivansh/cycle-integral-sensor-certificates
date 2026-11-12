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

    def test_invalid_times_and_shapes_rejected(self):
        with self.assertRaises(ValueError):
            summarize_window([0,0,1], np.ones((2,3)))
        with self.assertRaises(ValueError):
            summarize_window([0,1,3], np.ones((2,3)))
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
