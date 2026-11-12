# Proof-carrying inclusion-minimal conflict hyperedges

**Date:** 12 November 2026. **Status:** mathematical/software research prototype. Not a new theorem, independent human peer review, or verified physical sensor fault detection. Published as part of the unmerged research PR #4; original GitHub \`main\` and the fixed Zenodo \`v0.1.0\` archive remain unchanged.

## What problem does this solve?

An exact infeasibility certificate for a network of sensor models does **not** prove that every sensor in that inconsistent group is necessary for the contradiction. A valid smaller conflict may exist. Calling an arbitrary incompatible group a *minimal* conflict would be mathematically unjustified.

The problem is amplified by global-but-not-pairwise constraints: our [three-sensor, two-window counterexample](network_wide_helly_consistency_2026-11-12.md) has every pairwise LP feasible while the full three-sensor LP is infeasible. An ordinary pairwise graph has **no edges** for that example. We therefore work with **conflict hyperedges**, not necessarily edges.

A subset \(S\) of sensor-model assumptions is **inclusion-minimal inconsistent** for a fixed constraint system if:

\[
F(S)=\varnothing
\quad\text{and}\quad
F(S\setminus\{i\})\neq\varnothing\ \text{for every } i\in S.
\]

Here \(F(U)\) denotes the feasible set of the **same finite multi-window outer-relaxation LP** restricted to sensors in \(U\). This is not necessarily the conflict subset of **smallest cardinality**. Because dropping constraints cannot destroy feasibility, checking all \(|S|\) immediate single-sensor deletions implies that *every proper subset* is feasible.

## Proof-carrying two-sided verification

The new [certified-conflict implementation](../experiments/certified_conflict_hyperedges.py) uses two separate classes of exact rational witnesses.

### 1. Exact Farkas contradiction for the full subset

The network [exact certificate module](../experiments/network_exact_certificate.py) rebuilds \(Ax\le b\) from **rational representations of the recorded input samples, timestamps, and numerical assumptions**. A SciPy/HiGHS dual solver proposes a sparse multiplier support, but a result is accepted only if the returned rational vector \(\lambda\) independently satisfies

\[
\lambda\ge0,\qquad \lambda^\mathsf T A=0,\qquad
\lambda^\mathsf T b<0.
\]

Then \(0\le\lambda^\mathsf Tb<0\) is an exact impossibility. All such verified proofs establish **model-relaxation inconsistency conditional on the assumed bounds**; they do not identify defective hardware.

### 2. Exact primal feasibility proof for *each* deletion

For each \(i\in S\), the verifier must construct a rational vector \(x_i\) satisfying **every** inequality of the reduced exact LP:

\[
A_{S\setminus\{i\}} x_i\le b_{S\setminus\{i\}}.
\]

A floating-point LP first proposes a point. The program rationalizes its coordinates and, where necessary, tries to reconstruct approximate active-constraint equalities by **exact Gaussian elimination**. **No candidate becomes a proof until every rational inequality is rechecked exactly.**

If any one deletion lacks a verified rational feasible point, the program **does not assert minimality**. It reports an unresolved proof obligation, even if the floating solver says that deletion is feasible.

### 3. Conservative outcome taxonomy

| Status | What is established | What is NOT established |
| --- | --- | --- |
| \`certified_minimal\` | Exact contradiction for the displayed sensor subset **and** exact feasibility after removing any one sensor | Minimum cardinality among *all* conflicts; hardware faults; physical model correctness |
| \`certified_inconsistent\` | Exact contradiction for the displayed subset | Minimality |
| \`not_certified\` | No exact contradiction was extracted | It does **not** mean model consistency or health |

This distinction is essential. A finite-precision solver can supply *candidates*, but it cannot substitute for the exact checker.

## The exact triple-only demonstration

For two one-second windows, three sensors, and zero sample/curvature uncertainty, the prior example yields the possible integrated shared-input regions

\[
C_1=\{(u,0):0\le u\le2\},\quad
C_2=\{(1,v):0\le v\le2\},\quad
C_3=\{(2-w,w):0\le w\le2\}.
\]

Every pair intersects:

- \(C_1\cap C_2=\{(1,0)\}\)
- \(C_1\cap C_3=\{(2,0)\}\)
- \(C_2\cap C_3=\{(1,1)\}\)

The triple intersection is empty and has an exact Farkas contradiction \(0\le-1\). **All three single-sensor deletions have explicitly checkable feasible rational parameter and integrated-input assignments**, so \(\{1,2,3\}\) is inclusion-minimal for the relaxed model. This example proves that an empty pairwise conflict graph does not imply global consistency for two or more windows.

The fact that this familiar geometric/constraint-programming phenomenon arises in the model is useful for software diagnosis, but **not evidence that it was first discovered here**.

## Validation and reproducibility

The tests in [test_certified_conflict_hyperedges.py](../tests/test_certified_conflict_hyperedges.py) cover:

- Exact minimal triple-only conflict and every exact feasible deletion;
- A two-sensor minimal conflict with two feasible one-sensor deletions;
- A redundant four-sensor superset reduced to a triple;
- A mathematically healthy high-common-offset example that is **not** falsely labeled inconsistent;
- Tiny real inconsistencies that remain **not certified** rather than triggering fabricated proofs;
- Deliberate tampering with a proposed feasible witness, which must fail validation;
- Randomized exact feasible models and a six-sensor synthetic injected-bias case.

Reproduce from repository root:

~~~bash
python -m pip install -r requirements-research.txt
python -m unittest discover -s tests -p 'test_certified_conflict_hyperedges.py' -v
python -m unittest discover -s tests -v
~~~

Actual CI run and exact test count must be checked after the final branch commit. Passing tests show executed correctness checks, **not** an exhaustive proof of every possible input.

## Relation to established research

Finding inclusion-minimal inconsistent subsets/unsatisfiable cores is **established model-based diagnosis and constraint-solving practice**, related to the classical conflict and hitting-set work of Reiter (1987) and de Kleer & Williams (1987). The exact Farkas and primal-witness checks use established linear programming/convex optimization.

This prototype contributes **transparent proof-carrying implementation**, not a validated new general algorithm.

## Open research gaps

1. **Incomplete candidate search:** The exact verifier may fail to find a certificate that exists; the primal feasible-witness generator may fail on difficult degenerate cases. These failures must remain inconclusive.
2. **Complexity and minimum-cardinality diagnoses:** Greedy deletions need multiple LP checks. Inclusion-minimality is weaker than minimum cardinality; many minimal conflicts may exist. Enumerating all or minimizing hitting-set diagnoses is not implemented.
3. **Hardware/model ambiguity:** Shared-temperature, synchronization, noise, curvature, calibration and constant-time-constant assumptions can fail globally. A hyperedge says the chosen hypotheses are incompatible **conditional on the model**, not how many devices have failed.
4. **Exactness versus measurement:** Arithmetic is exact for the **stored rational interpretation** of inputs, but physical quantities may have rounding and unmodeled uncertainty.
5. **Real data and peer review:** There is no calibrated physical dataset, published baseline hardware comparison, or independent human proof review in this experiment.
6. **Scientific novelty:** No new theorem or demonstrated practical advantage over prior diagnosis and set-membership literature has been established.

**Scientific decision:** Keep draft PR #4 unmerged. Review proof semantics and exact input handling before making claims beyond finite LP model incompatibility.
