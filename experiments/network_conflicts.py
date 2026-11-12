"""Experimental global-model consistency and inclusion-minimal conflict extraction.

A network-wide LP with a shared unknown integrated input J_w may reject a
network even when each pairwise relaxed LP is individually feasible.

Inclusion-minimal conflicts are NOT minimum-cardinality fault explanations;
no observed contradiction automatically establishes a broken sensor.
Uses numerical HiGHS, so 'solver_infeasible' is not exact proof in general.
"""
from experiments.multiwindow_lp import joint_relaxation


def network_consistency(windows, max_tau, eps, curvature):
    """Run the entire common-input relaxed model at once."""
    return joint_relaxation(windows, max_tau, eps, curvature)


def inclusion_minimal_conflict(windows, max_tau, eps, curvature):
    """Return one inclusion-minimal *candidate* sensor conflict.

    Repeatedly remove a sensor while the remaining network LP is still
    solver-infeasible. This is a standard deletion-based MUS algorithm,
    NOT a novel method and NOT guaranteed minimum-cardinality.

    For rigorous inconsistency, validate a Farkas certificate independently.
    For feasible or inconclusive initial statuses, return an empty tuple or
    raise; do not pretend that a conflict was found.
    """
    n = len(windows[0].integrals)
    indices = list(range(n))
    initial = joint_relaxation(windows, max_tau, eps, curvature,indices=indices)
    if initial.status == "feasible_relaxation":
        return ()
    if initial.status != "solver_infeasible":
        raise RuntimeError("global LP inconclusive: " + initial.message)
    for candidate in tuple(indices):
        reduced = [i for i in indices if i != candidate]
        if not reduced:
            continue
        test = joint_relaxation(windows, max_tau, eps, curvature,indices=reduced)
        if test.status == "inconclusive":
            raise RuntimeError("conflict minimization inconclusive: " + test.message)
        if test.status == "solver_infeasible":
            indices = reduced
    return tuple(indices)
