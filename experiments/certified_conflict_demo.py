"""Display a checkable inclusion-minimal network conflict proof.

This example concerns the finite relaxed LP, not physical hardware.
"""
from experiments.exact_pairwise_certificate import exact_window
from experiments.network_consistency_gap import observations
from experiments.certified_conflict_hyperedges import (
    certified_conflict_hyperedge, explain_conflict,
)


def main():
    windows = [exact_window(t,z) for t,z in observations()]
    result = certified_conflict_hyperedge(windows,2,0,0)
    if result.status != "certified_minimal":
        raise RuntimeError("Expected proof-carrying minimality; found "+
                           result.status + ": " + result.message)
    if not result.verify(windows,2,0,0):
        raise RuntimeError("Independent proof verification failed")
    print("Result:",result.status)
    print("Sensors (1-indexed):",explain_conflict(result)["sensors_one_indexed"])
    print("Exact Farkas contradiction: 0 <=",result.contradiction.weighted_rhs)
    print("Verified feasible one-sensor deletions:")
    for sensor,proof in result.deletion_witnesses:
        print("  remove sensor %d; exact feasible tau+J=%s"%
              (sensor+1,tuple(map(str,proof.values))))
    return result


if __name__ == "__main__":
    main()
