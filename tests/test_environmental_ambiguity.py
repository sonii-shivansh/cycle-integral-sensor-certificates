"""Regression tests for exact sensor-bias/environmental-input ambiguity.

Results concern model identifiability and conditional consistency only, not
classification accuracy or confirmed hardware failures.
"""
import unittest
from fractions import Fraction as F

import numpy as np

from experiments.environmental_ambiguity_study import (
    OFFSETS, EPS, K, M, indistinguishable_worlds,
    summarize_observations, evaluate_trial, run
)
from experiments.multiwindow_lp import (
    joint_relaxation, pairwise_results, signaware_window_passes
)


class EnvironmentalAmbiguityTests(unittest.TestCase):
    def test_exact_rational_two_worlds_are_identical(self):
        # Exact proof on a ramp: T=10+2t, tau=(0,1/2,2).
        # Observed additive b is either instrument bias or a shift in local T.
        taus=(F(0),F(1,2),F(2))
        shifts=(F(0),F(3,2),F(-1))
        for t in map(F,(0,1,2,7)):
            common=F(10)+2*t
            for tau,b in zip(taus,shifts):
                response=common-2*tau
                # World A: shared input, output bias b
                observed_a=response+b
                # World B: local input common+b, corresponding response.
                local_input=common+b
                response_b=local_input-2*tau
                self.assertEqual(observed_a,response_b)
                self.assertEqual(tau*F(2)+response_b,local_input)

    def test_exact_analytic_harmonic_ode_two_worlds(self):
        for drift in (0.0,0.3,1.0):
            for seed in (0,7,42):
                x=indistinguishable_worlds(seed,drift)
                self.assertTrue(np.array_equal(
                    x["z_output_bias"],x["z_local_temperature"]))
                self.assertLess(x["common_ode_residual_max"],1e-12)
                self.assertLess(x["local_ode_residual_max"],1e-12)
                self.assertLessEqual(np.max(np.abs(x["noise"])),EPS)

    def test_bounded_spatial_temperature_difference_restores_consistency(self):
        # Every environmental-world sensor satisfies tau*y'+y=T+b_i.
        # A physical bound gamma_i >= |b_i| must therefore prevent exact
        # model inconsistency (although numerical solvers still need care).
        for drift in (0.0,0.3,0.6,1.0):
            for seed in (0,3,9):
                result=evaluate_trial(seed,drift)
                self.assertEqual(result["bounded_spread_edges"],0)
                self.assertEqual(result["bounded_spread_global"],
                                 "feasible_relaxation")

    def test_more_spatial_allowance_never_adds_contradictions_in_stress_cases(self):
        data=indistinguishable_worlds(1,0.6)
        w=summarize_observations(data["t"],data["z_output_bias"])
        edges=[]
        for spread in (0.0,np.abs(OFFSETS),2*np.abs(OFFSETS)):
            results=pairwise_results(w,M,EPS,K,input_spread=spread)
            self.assertTrue(all(
                v.status!="inconclusive" for v in results.values()))
            edges.append({k for k,v in results.items()
                          if v.status=="solver_infeasible"})
        self.assertGreater(len(edges[0]),0)
        self.assertTrue(edges[2] <= edges[1] <= edges[0])
        self.assertEqual(edges[1],set())
        for summary in w:
            for i in range(len(OFFSETS)):
                for j in range(i+1,len(OFFSETS)):
                    self.assertTrue(signaware_window_passes(
                        summary,i,j,M,EPS,K,input_spread=np.abs(OFFSETS)))

    def test_omitting_input_spread_keeps_old_lp(self):
        data=indistinguishable_worlds(0,0.3)
        windows=summarize_observations(data["t"],data["z_output_bias"])
        old=joint_relaxation(windows,M,EPS,K)
        explicit=joint_relaxation(windows,M,EPS,K,input_spread=0.0)
        self.assertEqual(old.status,explicit.status)
        self.assertEqual(
            {k:v.status for k,v in pairwise_results(windows,M,EPS,K).items()},
            {k:v.status for k,v in pairwise_results(
                windows,M,EPS,K,input_spread=0.0).items()}
        )

    def test_invalid_environmental_bounds_rejected(self):
        data=indistinguishable_worlds(0,.3)
        windows=summarize_observations(data["t"],data["z_output_bias"])
        for invalid in (-.1,np.nan,np.inf,[0.0]*5,
                        [0.0,0.0,0.0,-.1,0.0,0.0]):
            with self.subTest(invalid=str(invalid)):
                with self.assertRaises(ValueError):
                    joint_relaxation(windows,M,EPS,K,
                                     input_spread=invalid)
                with self.assertRaises(ValueError):
                    pairwise_results(windows,M,EPS,K,
                                     input_spread=invalid)

    def test_full_20_seed_stress_table_has_no_relaxed_spatial_conflicts(self):
        rows=run(trials=20)
        self.assertEqual(len(rows),4)
        for row in rows:
            self.assertEqual(row[2],0)
            self.assertEqual(row[4],0)
            self.assertEqual(row[3],20)


if __name__=="__main__":
    unittest.main()
