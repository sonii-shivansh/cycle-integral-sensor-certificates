"""Proof-carrying conflict hyperedge tests; no physical sensor fault claims."""
import unittest
from dataclasses import replace
from fractions import Fraction as F

import numpy as np

from experiments.exact_pairwise_certificate import (
    exact_window, exact_pairwise_verdict,
)
from experiments.network_consistency_gap import observations
from experiments.network_exact_certificate import exact_network_rows
from experiments.certified_conflict_hyperedges import (
    RationalPrimal, certified_conflict_hyperedge, certify_primal_feasible,
    explain_conflict,
)
from experiments.set_membership_comparison import recorded_windows
from experiments.multiwindow_case_study import K


class CertifiedHyperedgeTests(unittest.TestCase):
    def test_network_only_triple_is_proved_inclusion_minimal(self):
        windows = [exact_window(t,z) for t,z in observations()]
        result = certified_conflict_hyperedge(windows, 2, 0, 0)
        self.assertEqual(result.status, "certified_minimal", result.message)
        self.assertEqual(result.sensors,(0,1,2))
        self.assertEqual(len(result.deletion_witnesses),3)
        self.assertTrue(result.verify(windows,2,0,0))
        for removed, witness in result.deletion_witnesses:
            remaining=tuple(i for i in result.sensors if i!=removed)
            self.assertTrue(witness.check(exact_network_rows(
                windows,2,0,0,remaining
            )))
            self.assertEqual(exact_pairwise_verdict(
                windows,remaining,2,0,0
            ).status,"feasible_relaxation")
        self.assertEqual(explain_conflict(result)["minimality"],
                         "inclusion-minimal with exact deletion witnesses")

    def test_pairwise_conflict_has_singleton_deletion_witnesses(self):
        win = exact_window([0,1],[[0,0],[2,2]])
        result = certified_conflict_hyperedge([win],0,0,0)
        self.assertEqual(result.status,"certified_minimal",result.message)
        self.assertEqual(result.sensors,(0,1))
        self.assertTrue(result.verify([win],0,0,0))
        self.assertEqual(
            set(i for i,_ in result.deletion_witnesses),{0,1})

    def test_four_sensor_superset_reduces_to_three(self):
        data=[]
        for t,readings in observations():
            data.append(exact_window(t,readings+[readings[0]]))
        result=certified_conflict_hyperedge(data,2,0,0)
        self.assertEqual(result.status,"certified_minimal",result.message)
        self.assertEqual(len(result.sensors),3)
        self.assertTrue(result.verify(data,2,0,0))
        self.assertTrue({1,2}.issubset(set(result.sensors)))

    def test_healthy_network_is_not_falsely_certified_inconsistent(self):
        t=np.arange(201,dtype=float)
        tau=np.array([0.0,1.0,2.0])
        z=2.0**52+t[None,:]-tau[:,None]
        windows=[exact_window(t[a:b+1],z[:,a:b+1])
                 for a,b in [(0,50),(50,120),(120,200)]]
        result=certified_conflict_hyperedge(windows,3,0,0)
        self.assertEqual(result.status,"not_certified")
        self.assertIsNone(result.contradiction)

    def test_tiny_inconsistency_is_reported_not_certified(self):
        windows=[exact_window([0,1],[[0,0],[1e-12,1e-12]])]
        result=certified_conflict_hyperedge(windows,0,0,0)
        self.assertEqual(result.status,"not_certified")

    def test_tampered_deletion_witness_fails_verification(self):
        windows=[exact_window(t,z) for t,z in observations()]
        result=certified_conflict_hyperedge(windows,2,0,0)
        self.assertEqual(result.status,"certified_minimal")
        removed, witness=result.deletion_witnesses[0]
        values=list(witness.values)
        values[0]=F(-1)
        forged=replace(
            result,
            deletion_witnesses=((removed,RationalPrimal(tuple(values))),)
                               + result.deletion_witnesses[1:]
        )
        self.assertFalse(forged.verify(windows,2,0,0))

    def test_duplicate_deletion_witness_cannot_fake_minimality(self):
        windows=[exact_window(t,z) for t,z in observations()]
        good=certified_conflict_hyperedge(windows,2,0,0)
        self.assertEqual(good.status,"certified_minimal")
        # A full list of proofs is essential: repeating a valid proof for
        # one removed sensor must not replace a missing proof for another.
        forged=replace(
            good,
            deletion_witnesses=(
                good.deletion_witnesses[0],
                good.deletion_witnesses[0],
                good.deletion_witnesses[2]
            )
        )
        self.assertFalse(forged.verify(windows,2,0,0))

    def test_25_random_feasible_exact_models_prove_primal_when_possible(self):
        rng=np.random.default_rng(24011)
        successes=0
        for seed in range(25):
            # All sampled readings identical across sensors and integer
            # affine in each window. A common T and tau=0 explain each.
            a=int(rng.integers(-10,11))
            b=int(rng.integers(-4,5))
            records=np.array([[a,a+b],[a,a+b],[a,a+b]],dtype=int)
            w=exact_window([0,1],records)
            rows=exact_network_rows([w],2,0,0)
            result=certify_primal_feasible(rows)
            if result.status=="certified_feasible":
                successes+=1
                self.assertTrue(result.witness.check(rows))
        self.assertEqual(successes,25)

    def test_six_sensor_biased_example_has_checked_inconsistency(self):
        windows=recorded_windows(0,0.3,True)
        result=certified_conflict_hyperedge(windows,15,.02,K)
        self.assertIn(result.status,
                      ("certified_minimal","certified_inconsistent"))
        self.assertTrue(result.verify(windows,15,.02,K))
        self.assertLessEqual(len(result.sensors),6)
        if result.status=="certified_minimal":
            self.assertEqual(len(result.deletion_witnesses),
                             len(result.sensors))


if __name__=="__main__":
    unittest.main()
