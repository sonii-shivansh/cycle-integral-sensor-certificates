# Beyond pairwise graphs: an exact three-sensor network contradiction

**Research status:** 12 November 2026. An independently constructed mathematical example and a consequence of **classical Helly's theorem**. **Not a new theorem, not human peer review, and not evidence of faulty hardware**. Code and tests are in the draft branch \`research/multiwindow-lp-validation\`, not the Zenodo v0.1.0 software archive.

## Research question

Does a network-wide, shared-input LP capture any contradiction that is **invisible to all pairwise** multi-window LP tests?

**Answer: yes, if at least two observation windows are used.** We give a fully explicit counterexample for three sensors and two windows in exact rational arithmetic, with a Farkas witness \(0\le-1\). For one window, pairwise feasibility **does** imply global feasibility of the same relaxed LP (interval intersection). These are classical convex-feasibility facts, not unreported scientific discoveries.

## 1. Model and sensor-specific convex sets

For each healthy sensor \(i\), the original common-input first-order model is
\[
\tau_i y_i'(t)+y_i(t)=T(t),\qquad 0\le\tau_i\le M.
\]
For each observation window \(w\), the necessary integrated *outer relaxation* uses observed trapezoid \(I_{iw}\), endpoint increment \(d_{iw}\), and bounds \(\epsilon_i,K_i\):
\[
I_{iw}+\tau_i(d_{iw}-2\epsilon_i)-D_{iw}\le J_w\le
I_{iw}+\tau_i(d_{iw}+2\epsilon_i)+D_{iw},
\]
where \(D_{iw}=L_w\epsilon_i+L_wh_{\max,w}^2K_i/12\) and \(J_w=\int_{a_w}^{b_w}T(t)\,dt\). Applying the *same reference subtraction* to every sensor integral and \(J_w\) in a window preserves feasibility.

Write the possible **shared-input integral vectors** for sensor \(i\):
\[
C_i=\left\{(J_1,\ldots,J_W)\in\mathbb R^W:
\begin{array}{l}
\exists\,\tau_i\in[0,M]\text{ such that}\\
\forall w,\ I_{iw}+\tau_i(d_{iw}-2\epsilon_i)-D_{iw}\le J_w\\
\hspace{36mm}\le I_{iw}+\tau_i(d_{iw}+2\epsilon_i)+D_{iw}
\end{array}\right\}.
\]
Every \(C_i\) is **nonempty, compact and convex**: the joint set in \((\tau_i,J)\) is an intersection of affine halfspaces, with \(\tau_i\) in a closed bounded interval; it is bounded also in all \(J_w\), and projection preserves convexity and compactness.

The network LP is feasible **if and only if** \(\bigcap_i C_i\ne\emptyset\). A pairwise LP checks only \(C_i\cap C_j\ne\emptyset\).

## 2. Exact three-sensor example

Take \(M=2\), zero sample noise \(\epsilon_i=0\), zero curvature \(K_i=0\) *within each observed window*, and two **nonoverlapping** one-second windows, \([0,1]\) and \([3,4]\). The observed values are linear within each window. Endpoints:

| Sensor | Window 1 readings at (0,1) | Window 2 readings at (3,4) |
| --- | --- | --- |
| 1 | \((-1/2,1/2)\) | \((0,0)\) |
| 2 | \((1,1)\) | \((-1/2,1/2)\) |
| 3 | \((5/2,3/2)\) | \((-1/2,1/2)\) |

The exact composite-trapezoid integrals \(I_{iw}\) and endpoint changes \(d_{iw}\) are:

\[
I_1=(0,0),\ d_1=(1,0);
\qquad I_2=(1,0),\ d_2=(0,1);
\qquad I_3=(2,0),\ d_3=(-1,1).
\]

With zero observation and quadrature uncertainty, the relaxed per-sensor feasible shared-input sets are the line segments in \(\mathbb R^2\):
\[
C_1=\{(\tau_1,0):0\le\tau_1\le2\},
\quad C_2=\{(1,\tau_2):0\le\tau_2\le2\},
\quad C_3=\{(2-\tau_3,\tau_3):0\le\tau_3\le2\}.
\]

Every pair intersects, but at **three different** points:

| Pair | Shared \(J=(J_1,J_2)\) | Feasible time constants |
| --- | --- | --- |
| (1,2) | \((1,0)\) | \(\tau_1=1,\ \tau_2=0\) |
| (1,3) | \((2,0)\) | \(\tau_1=2,\ \tau_3=0\) |
| (2,3) | \((1,1)\) | \(\tau_2=1,\ \tau_3=1\) |

To satisfy **all three**, \(J_2=0\) from sensor 1. Sensor 2 then forces \(\tau_2=0\) and \(J_1=1\). Sensor 3 forces \(\tau_3=0\) and \(J_1=2\). This is a contradiction. No single triple-wide parameter choice exists, even though every pair was feasible with *its own* parameter choices.

### Explicit, exact Farkas certificate

Let \(E_{iw}=J_w-I_{iw}-\tau_i d_{iw}\). In the zero-error model, healthy observations require \(E_{iw}=0\). The following **four admissible inequalities**, each taken from the corresponding equality, are:

\[
E_{2,1}\le0,\qquad
-E_{3,1}\le0,\qquad
-E_{3,2}\le0,\qquad
E_{1,2}\le0.
\]

Adding all four inequalities with nonnegative unit multipliers **cancels every unknown \(J_1,J_2,\tau_1,\tau_2,\tau_3\)**, and the right-hand sides produce
\[
\boxed{0\le -1},
\]
an impossible condition. This certifies infeasibility of the **three-sensor LP in exact rational arithmetic**.

The [executable example](../experiments/network_consistency_gap.py) reconstructs the readings, verifies all three pairwise cases with the exact-rational verifier, checks global numerical infeasibility, and reconstructs the four-row Farkas witness using \`fractions.Fraction\`. No floating-point decision is needed for the exact proof.

This example is *not* a set of three simultaneously healthy physical sensors under a shared input. It is deliberately incompatible, and the joint LP correctly detects that. Measurements are only specified on the observed windows; no physical assumptions are made about unobserved gaps.

## 3. Conditional theorem: the Helly bound

**Proposition (classical Helly consequence).** For \(n\) sensors and \(W\ge1\) observed windows under the finite-dimensional LP outer relaxation above, if the full network is infeasible, there exists an infeasible subset of **at most \(\min(n,W+1)\) sensors**.

**Proof.** Each sensor's \(\,C_i\subseteq\mathbb R^W\) is convex. By [Helly's theorem](https://mathworld.wolfram.com/HellysTheorem.html), if every subfamily of \(W+1\) convex sets has nonempty intersection, all have nonempty intersection (for \(n\ge W+1\)). The contrapositive proves that an empty global intersection has an empty subfamily of at most \(W+1\). If \(n<W+1\), the full family is already such a subset. QED.

**Corollaries.**
- \(W=1\): pairwise feasibility is sufficient for *global LP feasibility* (the \(C_i\) are intervals on the real line).
- \(W=2\): at most three sensors suffice to witness any *global LP inconsistency*. The exact example shows the bound **three** can be necessary.
- For arbitrary \(W\), inclusion-minimal contradiction subsets cannot exceed \(W+1\) sensors in exact arithmetic, but **this does not imply** at least that many physically defective sensors.

This is a useful **classical theorem application**, not a new Helly theorem or a newly invented fault diagnosis algorithm. Helly-type reasoning appears throughout convex feasibility and constraint programming.

## 4. An experimental conflict hypergraph

The existing CIFC pairwise contradiction graph can be **empty even when a network-wide contradiction exists**, as in the three-sensor example. A richer representation is a hypergraph whose hyperedges are *minimal inconsistent sensor subsets* (potentially of size three or more). This is standard conflict diagnosis / unsatisfiable-core reasoning. We added:

- [Network-wide LP and deletion-based inclusion-minimal conflict search](../experiments/network_conflicts.py)
- [Regression tests](../tests/test_network_joint_gap.py)

The deletion procedure gives an **inclusion-minimal** *candidate* inconsistent subset, **not** necessarily the subset of minimum size; it relies on floating-point LP statuses, which need independent exact verification before being called mathematical certificates. For the exact constructed example, the candidate is the entire triple \((1,2,3)\), backed by the explicit rational certificate above.

One inconsistent hyperedge means **the shared-input/model assumptions cannot all hold simultaneously for its members**, not that there are three broken sensors. The same caution applies to the original pairwise graph's vertex-cover conclusions.

## 5. Computational scaling: clearly distinct tasks

Run the [network benchmark](../experiments/network_scaling_benchmark.py):

~~~bash
python -m pip install -r requirements-research.txt
python -m experiments.network_consistency_gap
python -m experiments.network_scaling_benchmark
python -m unittest discover -s tests -v
~~~

It times, on healthy analytical multi-sensor datasets:
- **One global LP**, which gives a yes/no consistency answer for the entire network;
- **All \(\binom{n}{2}\) pairwise LPs**, which return a complete pairwise contradiction graph.

The global LP has \(n+W\) unknowns and \(2nW\) interval inequalities plus \(n\) box constraints. Computing *all* pairwise LPs takes \(\binom{n}{2}\) small problems. The asymptotic number of LP calls is respectively 1 versus \(\Theta(n^2)\), **but these answer different questions**, so this difference is not a proved runtime superiority or an apples-to-apples comparison of identical output functionality. Timing is environment-dependent; no measured speedup or maximum-supported-network-size claim is made.

## 6. Scientific conclusion

1. **Mathematically established:** the current pairwise graph can miss global LP inconsistency; the three-sensor example has an independently checkable exact rational network witness.
2. **Classical theory, not new invention:** global vs pairwise consistency and Helly's \(W+1\) bound come from well-known convex geometry and set-membership diagnosis.
3. **Experimental:** a small candidate conflict hypergraph routine and size-scaling script have been added for future benchmarking.
4. **Still unresolved:** real sensor calibration, physical shared-input verification, certified *general* network-wide LP infeasibility, fair runtime comparison against a published *implementation*, and a distinctive new scientific advantage.

**Research decision:** Keep PR #4 as draft, and preserve the original published Zenodo archive. Avoid describing the classical Helly application as a new theorem.
