"""Network-wide consistency can reject when every pairwise relaxed LP passes.

This is a small exact-arithmetic example of a familiar global-vs-pairwise
consistency gap, NOT a new theorem or evidence of hardware fault detection.

Two disjoint 1-second observation windows. Each observed sensor is linear
on each observed window, so the claimed K=0 window curvature bounds hold.
The intervening unobserved gap is deliberately not modeled.
"""
from fractions import Fraction as F
from itertools import combinations
import numpy as np

from experiments.exact_pairwise_certificate import exact_window, exact_pairwise_verdict
from experiments.multiwindow_lp import summarize_window, joint_relaxation


def observations():
    """Return two windows: (times, 3x2 recorded readings)."""
    return [
        ([0, 1], [[F(-1,2), F(1,2)],
                  [F(1), F(1)],
                  [F(5,2), F(3,2)]]),
        ([3, 4], [[F(0), F(0)],
                  [F(-1,2), F(1,2)],
                  [F(-1,2), F(1,2)]])
    ]


def exact_summaries():
    return [exact_window(t,z) for t,z in observations()]


def float_summaries():
    return [summarize_window(np.array(t,dtype=float),
                             np.array(z,dtype=float))
            for t,z in observations()]


def _model_equality_row(w, i, sign):
    """A signed half of J_w - tau_i*d_iw = I_iw, as rational Ax <= b."""
    windows = exact_summaries()
    target = windows[w]
    a = [F(0)] * 5  # tau0,tau1,tau2,J0,J1
    a[i] = -sign * target.increments[i]
    a[3 + w] = sign
    rhs = sign * target.integrals[i]
    return tuple(a), rhs


def exact_network_contradiction():
    """Direct Farkas witness for the three-sensor/noiseless example.

    Four valid linear inequalities from the joint relaxed LP add to
    0 <= -1. All weights are +1. This witness uses *only* exact integers
    and fractions, no external optimization solver.
    """
    terms = [
        (0, 1, 1),   # +E(sensor2, window1) <= 0
        (0, 2, -1),  # -E(sensor3, window1) <= 0
        (1, 2, -1),  # -E(sensor3, window2) <= 0
        (1, 0, 1),   # +E(sensor1, window2) <= 0
    ]
    rows = [_model_equality_row(*term) for term in terms]
    lhs = tuple(sum((row[0][k] for row in rows),F(0)) for k in range(5))
    rhs = sum((row[1] for row in rows),F(0))
    if lhs != (F(0),)*5 or rhs != F(-1):
        raise AssertionError("exact Farkas witness verification failed")
    return rows, (lhs,rhs)


def evaluate():
    exact = exact_summaries()
    numeric = float_summaries()
    pairs = {}
    for pair in combinations(range(3),2):
        pair_exact = exact_pairwise_verdict(exact,pair,2,0,0)
        pair_float = joint_relaxation(numeric,2,0,0,indices=pair)
        if pair_exact.status != "feasible_relaxation":
            raise AssertionError("a pairwise relaxed model was rejected")
        if pair_float.status != "feasible_relaxation":
            raise AssertionError("numerical solver rejected a feasible pair")
        pairs[pair] = tuple(map(float,pair_float.time_constants))
    global_status = joint_relaxation(numeric,2,0,0)
    if global_status.status != "solver_infeasible":
        raise AssertionError("network-wide contradiction not detected")
    rows, proof = exact_network_contradiction()
    return {
        "pairwise_tau_witnesses": pairs,
        "network_numerical_status": global_status.status,
        "network_exact_witness_rhs": proof[1],
        "network_exact_combined_coefficients": proof[0],
        "network_exact_witness_rows": rows,
    }


if __name__ == "__main__":
    result = evaluate()
    print("Pairwise feasible response time constants:")
    for pair,taus in result["pairwise_tau_witnesses"].items():
        print(pair,taus)
    print("Network LP:",result["network_numerical_status"])
    print("Exact global contradiction: 0 <=",result["network_exact_witness_rhs"])
