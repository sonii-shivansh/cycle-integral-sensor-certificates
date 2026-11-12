# Independent verification and scientific corrections — 12 November 2026

**Status:** AI-assisted internal mathematical/numerical audit. Not human peer review, not a novelty determination, not a physical calibration study.

## Version and reproducibility

The audit considered draft PR #4, \`research/multiwindow-lp-validation\`. Its original benchmark source was reconstructed in a **separate Python implementation** that used a *two-variable pairwise LP with the shared integral eliminated*, rather than calling \`joint_relaxation\`. The independent environment was Python **3.13.5**, NumPy **2.3.5**, SciPy **1.17.0**. The audit was updated after the PR changed to add reference-centered trapezoidal integrals; the benchmark numbers remained unchanged at ordinary sensor magnitudes.

Check the [pull request](https://github.com/sonii-shivansh/cycle-integral-sensor-certificates/pull/4) for the exact latest SHA and CI status. The fixed Zenodo v0.1.0 archive is untouched.

## Mathematical verdict

**Correct as necessary conditions in exact real arithmetic:**

1. The integrated identity \(\int(y_j-y_i)=\tau_i\Delta y_i-\tau_j\Delta y_j\) follows by subtracting two first-order response equations with a common \(T(t)\) and constant \(\tau_i,\tau_j\).
2. The symmetric and sign-aware single-window inequalities follow from bounded endpoint errors, nonnegative response times, positive-weight trapezoidal integration and the standard curvature error bound. They are **conservative**, not sufficient descriptions of the ODE model.
3. The sign-aware feasible interval is contained in the original symmetric interval when evaluated in exact arithmetic.
4. The multi-window LP is a **necessary outer relaxation** for a set of healthy sensors sharing an input and keeping one response-time constant across all windows. If it is **exactly infeasible**, the joint assumptions cannot all be correct. Numerical HiGHS \`solver_infeasible\` must not be labelled a formal mathematical certificate without a verified infeasibility witness.

The LP unknowns are one \(\tau_i\) per selected sensor and one integrated input \(J_w\) per window. Its inequalities are linear in those unknowns. They do **not** enforce the complete continuous-time ODE, cross-window endpoint-noise correlations or complete trajectory coupling. Passing the LP does not establish that a physical explanation exists.

## Independent simulation results

The independent implementation reconstructed the six-sensor nonperiodic-window generator using its stated constants and 20 NumPy seeds per condition. It produced exactly the following **total contradictory pair counts**:

| Drift °C/s | Injected biases? | Original | Sign-aware | Pairwise multi-window LP |
|---:|:---:|---:|---:|---:|
| 0.0 | No | 0 | 0 | 0 |
| 0.0 | Yes | 16 | 58 | 180 |
| 0.3 | No | 0 | 0 | 0 |
| 0.3 | Yes | 35 | 35 | 180 |
| 0.6 | No | 0 | 0 | 0 |
| 0.6 | Yes | 0 | 74 | 180 |
| 1.0 | No | 0 | 0 | 0 |
| 1.0 | Yes | 0 | 20 | 133 |

A separate 250-seed randomized **healthy** model experiment included \(\tau=0\), \(\tau=M\), multiple harmonics, bounded worst-case sample noise and random signed drift, producing **3,750 pair LP feasibility checks with no solver-infeasible result**. This is useful regression evidence, not proof of numerical safety in every regime or an empirical field false-alarm estimate.

In the 20-seed benchmark, minimum-vertex-cover explanations selected both injected-bias sensors uniquely in **20/20** trials at each drift 0.0, 0.3, 0.6 and in **14/20** trials at drift 1.0. The other six high-drift trials did **not** uniquely identify that exact pair. None of these counts establish a guaranteed identification rate for physical hardware.

## Numerical conditioning: common reference centering matters

Construct three mathematically healthy ramp sensors:
\[
T(t)=2^{52}+t,\qquad y_i(t)=T(t)-\tau_i,\qquad\tau=(0,1,2)\,\mathrm{s},
\]
observed at integer times \(0,1,\dots,200\), with \(\epsilon=K=0\). Every observed value is exactly representable in binary64.

Integrating each \(\sim 10^{15}\)-scale trace independently using floating-point trapezoids and forming the LP yields a **false numerical infeasible** status, even though the true \(\tau\) vector satisfies the physical model.

Subtracting the *same* reference value from every sensor trace **within each window** is an exact mathematical transformation: both \(I_{i,w}\) and \(J_w\) are shifted by \(L_w\,c_w\), leaving feasibility unchanged. The latest draft branch's \`summarize_window\` performs this centering. In an independent numerical reconstruction the centered LP became feasible and recovered \(\tau=(0,1,2)\). This is an important demonstrated mitigation, not a universal guarantee of floating-point soundness.

## New counterexample: every pair passes, the network fails

**A pairwise contradiction graph is not a complete test of multi-sensor consistency.**

Take three sensors, two observation windows and \(0\le\tau_i\le1\). In the zero-noise, zero-quadrature-error idealization, let measured integral and endpoint increment summaries be

\[
\begin{array}{c|cc}
&\text{window A}&\text{window B}\\
i=1&(I,d)=(0,1)&(0,0)\\
i=2&(0,0)&(0,1)\\
i=3&(1,-1)&(0,1)
\end{array}
\]

Then any sensor's model requires a common integrated-input vector \((J_A,J_B)\). The three constraints individually trace **line segments**

\[
S_1=\{(t,0):0\le t\le1\},\quad
S_2=\{(0,t):0\le t\le1\},\quad
S_3=\{(1-t,t):0\le t\le1\}.
\]

They intersect **pairwise** at \((0,0)\), \((1,0)\), and \((0,1)\) respectively, but their **triple intersection is empty**. Thus *every pairwise LP is feasible* while the network LP is infeasible.

This is not an impossible observation dataset: the windows can be [0,1] and [11,12] seconds. Use a linear trace on each window, with sensor readings:
- Window A: \(z_1=s-\tfrac12\), \(z_2=0\), \(z_3=\tfrac32-s\), \(s=t\in[0,1]\).
- Window B: \(z_1=0\), \(z_2=s-\tfrac12\), \(z_3=s-\tfrac12\), \(s=t-11\in[0,1]\).

The separate window traces can be connected by **quintic Hermite interpolation** over the ten-second gap, matching the values, first derivatives and zero second derivatives at both ends. An explicit rational-polynomial calculation bounds the maximum \(|y''|\) on the three bridges by approximately **0.423**, so a global \(C^2\) interpolation with \(K=1\) exists. Setting a global \(K=1\), sample step \(h=0.01\), and \(\epsilon=0\) gives the small, conservative quadrature radius \(D=1(0.01)^2/12\) per sensor per window. An independent SciPy LP still found all three pairs feasible and the three-sensor joint problem infeasible. See [tests/test_higher_order_consistency.py](../tests/test_higher_order_consistency.py).

**Implication:** The pairwise graph can have **no edges** and hence matching/vertex-cover lower bound zero despite a genuine higher-order joint contradiction. If network-level localization is required, investigate **minimal infeasible subsets / conflict hypergraphs**, not only pairwise edges and vertex covers. This idea is closely related to classical consistency-based diagnosis and is **not claimed to be an original mathematical theorem**.

## Fundamental identifiability limitation (unresolved)

If the true common input is \(T(t)=T_0+ct\), then a sensor with steady ramp response \(y_i=T-c\tau_i\) and constant bias \(B_i\) has the same observed trace as a purportedly healthy sensor with \(\widetilde{\tau}_i=\tau_i-B_i/c\), whenever \(\widetilde{\tau}_i\in[0,M]\). With identical slope \(c\) across all windows, **even the multi-window LP cannot distinguish these cases**. A verified four-window zero-noise example with actual \(\tau_i=2\), \(B_i=1.2\), \(c=1\) is feasible as a model with \(\widetilde{\tau}_i=0.8\).

## Noise treatment: correct but conservative

For fixed nonnegative \(\tau_i\), the exact extreme contribution of bounded discrete sample error to the combined endpoint and trapezoid term is

\[
\epsilon_i\left(L-h+\max(2\tau_i,h)\right).
\]

This expression was checked by enumerating all error-sign vertices of small sampling grids. It is generally tighter than the original separate \(L\epsilon_i+2\tau_i\epsilon_i\) bound. The published inequalities remain **sound**, because they overbound the uncertainty. Applying the tighter factor to the LP without compromising linearity and solver soundness needs a separate, formally justified design.

## Prior art and originality

The literature includes:
- de Kleer & Williams (1987), *Diagnosing multiple faults*, DOI \`10.1016/0004-3702(87)90063-4\`: conflict/hitting-set reasoning.
- Loiez & Taillibert (IJCAI 1997), *Polynomial Temporal Band Sequences for Analog Diagnosis*: integration-based analog model consistency under bounded uncertainty. https://ijcai.org/proceedings/1997/70
- Tagawa & Ohta (1997), *Two-thermocouple probe for fluctuating temperature measurement in combustion*, DOI \`10.1016/S0010-2180(97)00044-8\`: shared-input thermocouple dynamics.
- Ingimundarson, Bravo, Puig & Alamo (CDC/ECC 2005), *Robust Fault Diagnosis using Parallelotope-Based Set-Membership Consistency Tests*: feasible invariant parameter sets for robust fault diagnosis. https://skoge.folk.ntnu.no/prost/proceedings/cdc-ecc05/pdffiles/papers/2056.pdf

**No world-first or patentability claim can be supported** by this audit. The exact closed-form packaging might be distinctive, but there is significant related mathematics.

## Unfinished requirements before merging / publishing claims

1. Obtain external *human* mathematical and control-theory review.
2. Validate outward-rounded or rational LP infeasibility certificates where a provable result is required.
3. Test sensitivity to nonuniform grids, time synchronization, parameter drift, unknown \(K\) and shared temperature failure.
4. Compare against a properly implemented, published set-membership fault-diagnosis baseline.
5. Quantify higher-order conflict localization; pairwise graph explanations cannot capture all joint inconsistencies.
6. Acquire calibrated sensor measurements and document absolute measurement bounds.
7. Treat CI successes and synthetic observations as software validation, not as peer review or novelty proof.
