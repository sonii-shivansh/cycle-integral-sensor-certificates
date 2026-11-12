# Independent classical set-membership baseline: scientific comparison

**Date:** 12 November 2026. **Research branch:** \`research/multiwindow-lp-validation\`. **Status:** computational cross-check of an established method, not evidence of new mathematics, superior field performance, or a physical fault detector.

## 1. Direct primary-source comparison

Ingimundarson, A., Bravo, J. M., Puig, V., and Alamo, T. (2005), *Robust Fault Diagnosis using Parallelotope-Based Set-Membership Consistency Tests*, Proceedings of IEEE CDC/ECC 2005, pp. 993–998, DOI [10.1109/CDC.2005.1582287](https://doi.org/10.1109/CDC.2005.1582287), [open proceedings manuscript](https://skoge.folk.ntnu.no/prost/proceedings/cdc-ecc05/pdffiles/papers/2056.pdf).

The paper's Section II, equations (6)–(11), defines feasible parameter sets given deterministic output-error bounds. Section IV-A discusses the case of **time-invariant parameters**. It explicitly identifies an LP containing constraints from **all previous observations** as the least conservative test in that model, and studies approximate feasible-set tracking (parallelotopes) to reduce computational burden.

**Important distinction:** the original paper uses a *general regression model* and proposes a *parallelotope-based approximation* of the feasible solution set. Our implementation below is **NOT a reproduction of its parallelotope algorithm or four-tank experiment**. We implement a simpler **exact two-dimensional feasible-set intersection**, a standard approach explicitly described by the authors as a conceptual reference. The sensor-specific integrated inequalities are the CIFC specialization.

This prior art establishes that “use the same unknown parameters across observations and check feasible-set emptiness” is **not a new research idea**.

Further related primary sources:
- Adrot & Ploix, *Fault Detection Based on Set-Membership Inversion* (2006), DOI [10.3182/20060829-4-CN-2909.00095](https://doi.org/10.3182/20060829-4-CN-2909.00095): set membership with uncertain analytic-redundancy relations.
- Puig (2010), *Fault diagnosis and fault tolerant control using set-membership approaches: Application to real case studies*, DOI [10.2478/v10006-010-0046-y](https://doi.org/10.2478/v10006-010-0046-y): review of consistency checks with bounded uncertainty.

## 2. Why the baseline and CIFC LP should agree

For sensor \(i\) and window \(w\), the sound *necessary* integrated model constraint is

\[
I_{iw} + \tau_i(d_{iw}-2\epsilon_i)-D_{iw}
\ \le J_w\le\
I_{iw}+\tau_i(d_{iw}+2\epsilon_i)+D_{iw}
\]

where \(J_w\) is the unknown integrated common input and \(D_{iw}\) incorporates deterministic sample noise and curvature/quadrature error.

For two sensors \(i,j\), a common \(J_w\) exists **if and only if** their two closed intervals intersect. This yields exactly two affine halfplanes in the **constant** unknown parameters \((\tau_i,\tau_j)\):

\[
(d_{iw}-2\epsilon_i)\tau_i-(d_{jw}+2\epsilon_j)\tau_j
\le I_{jw}-I_{iw}+D_{iw}+D_{jw},
\]
\[
-(d_{iw}+2\epsilon_i)\tau_i+(d_{jw}-2\epsilon_j)\tau_j
\le I_{iw}-I_{jw}+D_{iw}+D_{jw}.
\]

The feasible solution set for **all** windows is the intersection of these halfplanes with the response-time box \([0,M]^2\).

**Exact equivalence proposition.** Pairwise multi-window LP feasibility, exact polygon clipping of these halfplanes, and exact Fourier–Motzkin elimination applied to the same halfplanes have **the same mathematical yes/no answer**. Proof: each is a different exact decision procedure for the *identical finite system of linear inequalities* in \((\tau_i,\tau_j)\). Polygon clipping starts at the box and intersects halfplanes one at a time; Fourier–Motzkin projects away the variables using nonnegative combinations; linear programming searches the same feasible set. Numerical HiGHS outcomes can differ due to floating-point error, so that comparison is empirical rather than a floating-point identity.

This result is a standard convex-geometry/linear-programming fact, **not an original theorem**.

## 3. Independent reference implementation

- [Exact geometric set-membership baseline](../experiments/set_membership_fss_baseline.py): independently implemented Sutherland–Hodgman clipping with rational numbers. It never calls SciPy or imports the existing Fourier–Motzkin checker. It processes window constraints incrementally and, on contradiction, records the first inconsistent window.
- [Matched six-sensor benchmark](../experiments/set_membership_comparison.py): independently reconstructs raw readings and exact trapezoid integrals from the same synthetic generator; compares exact feasible-set clipping against the existing SciPy multi-window LP for **identical** windows, noise, curvature and response-time bounds.
- [Regression tests](../tests/test_set_membership_fss.py): cross-check the geometric solver with exact Farkas elimination, healthy extreme numerical cases, tangent/degenerate parameter sets, input validation and the full 20-seed comparison.

The exact baseline and the exact Farkas solver both interpret **actual binary64 values as exact rational numbers**. This avoids introducing numerical solver error into the geometric baseline but cannot recover precision lost by the original physical recording/clock.

## 4. Reproduced comparison

Run at repository root:

~~~bash
python -m pip install -r requirements-research.txt
python -m experiments.set_membership_comparison
python -m unittest discover -s tests -v
~~~

20 deterministic random seeds per drift per healthy/bias condition; six sensors and four identical windows per trial. Counts are **total pairwise contradiction edges** (15 possible pairs per trial, 300 possible across 20 trials).

| Drift °C/s | Biases injected | Symmetric | Sign-aware | CIFC multi-window HiGHS | Classical exact FSS clipping |
|---:|:---:|---:|---:|---:|---:|
| 0.0 | No | 0 | 0 | 0 | 0 |
| 0.0 | Yes | 16 | 58 | 180 | 180 |
| 0.3 | No | 0 | 0 | 0 | 0 |
| 0.3 | Yes | 35 | 35 | 180 | 180 |
| 0.6 | No | 0 | 0 | 0 | 0 |
| 0.6 | Yes | 0 | 74 | 180 | 180 |
| 1.0 | No | 0 | 0 | 0 | 0 |
| 1.0 | Yes | 0 | 20 | 133 | 133 |

**Important negative result:** The classical geometric feasible-set baseline and CIFC's multi-window LP make **identical decisions on every tested sensor pair** under these shared constraints. That is consistent with the mathematical equivalence, and shows **no detection advantage** over this standard invariant-parameter set-membership formulation.

This conclusion does **not** imply that the original authors' *parallelotope approximation algorithm* necessarily gives identical decisions or has identical speed; it has not been implemented or benchmarked here. No claim of runtime superiority, novel theory, or field performance is made.

## 5. Scientific implications and remaining tasks

1. **Novelty:** Intersecting bounded-error feasible parameter regions across windows is established. CIFC's novelty cannot rest on this alone. Any publishable innovation needs a new guarantee, a demonstrably better numerical/real-world tradeoff, or another independently validated contribution.
2. **Correctness:** Three computational paths now cross-check the same pairwise relaxed LP: SciPy/HiGHS, exact rational Fourier–Motzkin with verifiable Farkas multipliers, and an independent exact rational polygon-clipping baseline. Agreement helps find implementation bugs but is **not peer review** or proof of the entire physical model.
3. **Model assumptions:** A conflict certifies model inconsistency *conditional* on shared input, bounded noise/curvature, and invariant response times. It cannot prove defective hardware.
4. **Baseline scope:** To compare an actual published implementation, separately implement the paper's stated parallelotope procedure or locate an authorized reference implementation; do not label this generic set-membership baseline as the authors' software.
5. **Practical validation:** The next discriminating experiment should use actual calibrated co-located temperature sensors, synchronized sampling and justified deterministic uncertainty bounds, rather than additional hand-picked synthetic detection counts.

**Decision:** Keep PR #4 draft. Reframe its main contribution as a reproducible, inspectable specialization of established set-membership methods. Do not claim that multi-window LP alone is a research discovery.
