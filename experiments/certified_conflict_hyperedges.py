"""Proof-carrying inclusion-minimal sensor-model conflict extraction.

NOT a physical hardware fault identifier. Every returned certified conflict
contains a checked exact-rational Farkas proof, and every immediate deletion
has a checked exact-rational primal feasible assignment.

This is a conservative, *incomplete* procedure. LP solvers are used ONLY to
propose candidate witnesses; all positive proof claims are independently
checked with Python's Fraction exact arithmetic.

'not_certified' never means healthy. 'certified_inconsistent' does not mean
minimal. 'certified_minimal' means inclusion-minimal, NOT minimum cardinality.

The model is the finite multi-window CIFC *outer relaxation* formed from
recorded numeric input values interpreted as exact rational numbers.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from numbers import Integral
import math

import numpy as np
from scipy.optimize import linprog

from experiments.network_exact_certificate import (
    NetworkFarkasCertificate,
    certify_network_infeasibility,
    exact_network_rows,
)


def _f(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    x = float(value)
    if not math.isfinite(x):
        raise ValueError("nonfinite numeric value")
    return Fraction.from_float(x)


@dataclass(frozen=True)
class RationalPrimal:
    """Exact feasible vector for an already specified Ax <= b."""
    values: tuple[Fraction, ...]

    def check(self, rows):
        return (
            bool(rows)
            and len(self.values) == len(rows[0].coefficients)
            and all(
                sum((a*x for a, x in zip(row.coefficients, self.values)),
                    Fraction(0)) <= row.rhs
                for row in rows
            )
        )


@dataclass(frozen=True)
class FeasibilityAudit:
    status: str  # certified_feasible | not_certified
    witness: RationalPrimal | None
    message: str


def _affine_repair(rows, active, seed):
    """Solve selected active equalities exactly; use seed for free variables.

    Any inconsistent set of guessed active rows is rejected, never interpreted
    as evidence of infeasibility.
    """
    d = len(seed)
    m = len(active)
    if m == 0:
        return tuple(seed)
    table = [
        list(rows[j].coefficients) + [rows[j].rhs]
        for j in active
    ]
    pivot_columns = []
    rank = 0
    for col in range(d):
        pivot = next((r for r in range(rank, m) if table[r][col] != 0), None)
        if pivot is None:
            continue
        table[rank], table[pivot] = table[pivot], table[rank]
        scale = table[rank][col]
        table[rank] = [u/scale for u in table[rank]]
        for r in range(m):
            if r != rank and table[r][col]:
                factor = table[r][col]
                table[r] = [
                    u-factor*v for u, v in zip(table[r], table[rank])
                ]
        pivot_columns.append((col, rank))
        rank += 1
        if rank == m:
            break
    if any(all(x == 0 for x in row[:d]) and row[d] != 0
           for row in table):
        return None
    result = list(seed)
    for col, row in reversed(pivot_columns):
        result[col] = table[row][d] - sum(
            (table[row][j]*result[j] for j in range(d) if j != col),
            Fraction(0)
        )
    return tuple(result)


def certify_primal_feasible(rows):
    """Find and *verify* an exact rational primal witness, or return unsure.

    Fast numerical candidate generation does not establish feasibility by
    itself. The final predicate uses exact rational sums and comparisons.
    """
    if not rows:
        raise ValueError("expected at least one inequality")
    dimension = len(rows[0].coefficients)
    if dimension == 0 or any(len(r.coefficients) != dimension for r in rows):
        raise ValueError("invalid row dimensions")
    try:
        A = np.asarray([[float(c) for c in r.coefficients] for r in rows])
        b = np.asarray([float(r.rhs) for r in rows])
        if not (np.all(np.isfinite(A)) and np.all(np.isfinite(b))):
            return FeasibilityAudit("not_certified", None,
                                    "coefficients too large for numeric witness search")
        solution = linprog(
            np.zeros(dimension), A_ub=A, b_ub=b,
            bounds=[(None, None)] * dimension, method="highs"
        )
    except (OverflowError, ValueError) as exc:
        return FeasibilityAudit("not_certified", None,
                                "numerical primal construction failed: "+str(exc))

    if solution.status != 0 or solution.x is None:
        return FeasibilityAudit("not_certified", None,
                                "no usable numerical primal candidate")

    # The raw binary64 candidate and several low-complexity rationalizations.
    seeds = [tuple(Fraction.from_float(float(v)) for v in solution.x)]
    for denominator in (1000, 1000000, 10**9, 10**12):
        seeds.append(tuple(s.limit_denominator(denominator) for s in seeds[0]))
    for seed in seeds:
        proof = RationalPrimal(seed)
        if proof.check(rows):
            return FeasibilityAudit("certified_feasible", proof,
                                    "exact rational primal witness")

    residuals = A @ solution.x - b
    scales = 1.0 + abs(b) + np.sum(abs(A) * (abs(solution.x)[None, :] + 1), axis=1)
    normalized = np.abs(residuals) / scales
    for threshold in (1e-12, 1e-10, 1e-8, 1e-6):
        close = [int(i) for i in np.argsort(normalized)
                 if normalized[i] <= threshold]
        # Inexact active constraints may be mutually inconsistent as exact
        # equations. Add only rows compatible with the current active system.
        active = []
        for j in close:
            if len(active) >= dimension:
                break
            if _affine_repair(rows, active+[j], seeds[0]) is not None:
                active.append(j)
        if not active:
            continue
        for seed in seeds:
            repaired = _affine_repair(rows, active, seed)
            if repaired is not None:
                proof = RationalPrimal(repaired)
                if proof.check(rows):
                    return FeasibilityAudit(
                        "certified_feasible", proof,
                        "exact rational repaired active-set witness"
                    )
            # A reduced active set can avoid accidentally forcing a near-tight
            # but nonbinding numeric inequality.
            for count in range(len(active)-1, 0, -1):
                repaired = _affine_repair(rows, active[:count], seed)
                if repaired is not None:
                    proof = RationalPrimal(repaired)
                    if proof.check(rows):
                        return FeasibilityAudit(
                            "certified_feasible", proof,
                            "exact rational reduced-active-set witness"
                        )

    return FeasibilityAudit("not_certified", None,
                            "candidate feasible point failed exact verification")


@dataclass(frozen=True)
class CertifiedConflict:
    status: str  # certified_minimal | certified_inconsistent | not_certified
    sensors: tuple[int, ...]
    contradiction: NetworkFarkasCertificate | None
    contradiction_rows: tuple
    deletion_witnesses: tuple[tuple[int, RationalPrimal], ...]
    message: str

    def verify(self, windows, max_tau, eps, curvature):
        """Independently check all recorded witnesses against rebuilt exact LPs."""
        if self.status == "not_certified":
            return self.contradiction is None
        if self.contradiction is None or not self.contradiction.check(
            exact_network_rows(windows, max_tau, eps, curvature, self.sensors)
        ):
            return False
        if self.status == "certified_minimal":
            if len(self.deletion_witnesses) != len(self.sensors):
                return False
            for removed, witness in self.deletion_witnesses:
                if removed not in self.sensors:
                    return False
                subset = tuple(i for i in self.sensors if i != removed)
                if not subset:
                    return False
                if not witness.check(
                    exact_network_rows(windows, max_tau, eps, curvature, subset)
                ):
                    return False
        return True


def certified_conflict_hyperedge(windows, max_tau, eps, curvature,
                                 sensor_indices=None):
    """Certify one inclusion-minimal conflict when all needed proofs exist.

    A checked infeasibility certificate for the full candidate is mandatory.
    Greedy deletions occur only when the smaller set is *also* certified
    infeasible. To declare inclusion-minimality, EVERY final one-sensor
    deletion must have a checked exact primal feasible witness.

    Incomplete candidate search never upgrades to a mathematical claim.
    """
    if not windows:
        raise ValueError("need observation windows")
    n = len(windows[0].integrals)
    members = (tuple(range(n)) if sensor_indices is None
               else tuple(sensor_indices))
    if (not members or len(set(members)) != len(members)
            or any(not isinstance(i, Integral) or i < 0 or i >= n
                   for i in members)):
        raise ValueError("invalid sensor selection")

    whole = certify_network_infeasibility(
        windows, max_tau, eps, curvature, members
    )
    if whole.status != "certified_infeasible":
        return CertifiedConflict("not_certified", members, None, (), (),
                                 "no exact whole-subset contradiction")
    final = whole
    for candidate in tuple(members):
        reduced = tuple(i for i in members if i != candidate)
        if not reduced:
            continue
        test = certify_network_infeasibility(
            windows, max_tau, eps, curvature, reduced
        )
        if test.status == "certified_infeasible":
            members = reduced
            final = test

    # Rebuild/check the final full-set certificate. It might have changed in
    # the reduction pass.
    final = certify_network_infeasibility(
        windows, max_tau, eps, curvature, members
    )
    if final.status != "certified_infeasible":
        return CertifiedConflict("not_certified", members, None, (), (),
                                 "final subset contradiction could not be reverified")

    proofs = []
    for removed in members:
        reduced = tuple(i for i in members if i != removed)
        # A singleton contradiction cannot occur under correctly formed
        # single-sensor integrated-input relaxation with nonnegative bounds,
        # but keep the algorithm explicit.
        if not reduced:
            return CertifiedConflict(
                "certified_inconsistent", members, final.certificate,
                final.rows, tuple(proofs), "empty deletion not implemented"
            )
        child_rows = exact_network_rows(windows, max_tau, eps, curvature,
                                       reduced)
        child = certify_primal_feasible(child_rows)
        if child.status != "certified_feasible":
            return CertifiedConflict(
                "certified_inconsistent", members, final.certificate,
                final.rows, tuple(proofs),
                f"deletion of sensor {removed} lacks an exact feasible witness"
            )
        proofs.append((removed, child.witness))

    result = CertifiedConflict(
        "certified_minimal", members, final.certificate, final.rows,
        tuple(proofs), "exact Farkas contradiction and every exact deletion witness"
    )
    if not result.verify(windows, max_tau, eps, curvature):
        raise ArithmeticError("certified conflict failed independent recheck")
    return result


def explain_conflict(result):
    """Serializable summary; do not mislabel inconclusive cases."""
    if result.status == "not_certified":
        return {"status":result.status,"message":result.message}
    output = {
        "status":result.status,
        "sensors_zero_indexed":list(result.sensors),
        "sensors_one_indexed":[i+1 for i in result.sensors],
        "farkas_contradiction":str(result.contradiction.weighted_rhs),
        "minimality":(
            "inclusion-minimal with exact deletion witnesses"
            if result.status == "certified_minimal"
            else "unresolved"
        ),
        "feasible_deletions":[i+1 for i,_ in result.deletion_witnesses],
    }
    return output
