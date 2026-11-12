# Scientific integrity and adversarial verification audit

**Date:** 12 November 2026. **Scope:** CIFC original and sign-aware mathematical claims, draft PR #4 multi-window LP, the author's Gemini/Codex AI-generated referee reports, and independent numerical reconstructions.

**Pinned comparison base (main):** \`c09d5367715dd0dd7fc16dee3faf860536d87587\`. **Validated research-branch implementation prior to adding this note:** \`3008078e42a08a34011a69f0a0c0250496062d52\`. **GitHub Actions:** run [43159443323](https://github.com/sonii-shivansh/cycle-integral-sensor-certificates/actions/runs/43159443323), completed successfully, including regression tests and the benchmark. These results are not human peer review.

## Evidence taxonomy

- **Mathematical deduction:** the shared-input integral identity, original conservative bound, sign-aware conservative interval, and multi-window LP *necessary* constraints follow from a constant first-order sensor model and valid deterministic uncertainty limits.
- **Independent numerical reconstruction:** a separate Python/NumPy/SciPy implementation produced exactly the 20-seed, four-drift, healthy/fault synthetic contradiction counts shown below.
- **Confirmed code-specific problems:** two adversarial constructions expose a catastrophic-cancellation false-infeasibility result and an incorrect quadrature bound for a nonuniform grid accepted by the old validator. Both have corrections and regression tests *on this draft branch*.
- **Not established:** field false-alarm rates, statistical detection performance, identifying defective hardware, superiority over state-of-the-art methods, originality, fully verified floating-point infeasibility, or peer-review acceptance.

## 1. Proof-level audit

For healthy sensor \(i\), let \(0\le\tau_i\le M\), \(\tau_i y'_i+y_i=T\), bounded noise \(|z_i-y_i|\le\epsilon_i\), and bounded curvature \(|y_i''|\le K_i\). The exact integrated identity is

\[
\int_a^b(y_j-y_i)\,dt=\tau_i\Delta y_i-\tau_j\Delta y_j.
\]

The original symmetric and sign-aware inequalities follow by interval arithmetic, trapezoid-error bounds and the triangle inequality; they are **necessary, not sufficient**, and their endpoint/integral noise terms are **correlated**. The sign-aware accepted interval is a subset of the original symmetric accepted interval in exact arithmetic. Careful outward rounding would still be needed for a mathematically certified float implementation.

For any healthy subset \(S\), the exact per-window integrated input is \(J_w=\int T\), so each sensor must satisfy

\[
 I_{iw}+\tau_i(d_{iw}-2\epsilon_i)-D_{iw}
 \le J_w\le
 I_{iw}+\tau_i(d_{iw}+2\epsilon_i)+D_{iw},
\quad
 D_{iw}=L_w\epsilon_i+\frac{L_w h_{\max,w}^2K_i}{12}.
\]

For the *reference-centered* implementation, replace \(I_{iw}\) and \(J_w\) by their common-reference-shifted values. Exact infeasibility of these constraints rejects the joint model assumptions; numerical HiGHS status 2 **alone is not a formal contradiction certificate**. A feasible solution is not sufficient to establish physically valid trajectories.

A matching or vertex cover counts necessary removals of *model-labeled sensor hypotheses under globally valid assumptions*. A spatially varying temperature field can generate many conflicts despite zero defective devices. The graph's minimum cover is an explanation, not a verified hardware diagnosis.

## 2. Reproduced matched-window numerical result

20 seeded trials per drift and per sensor state; same sample times, four observation windows, six sensor models, and deterministic uncertainty limits. Total contradiction **edges summed across 20 trials**, not classification accuracy:

| Drift (°C/s) | Faults injected? | Symmetric | Sign-aware | Multi-window LP |
| ---: | :---: | ---: | ---: | ---: |
| 0.0 | No | 0 | 0 | 0 |
| 0.0 | Yes | 16 | 58 | 180 |
| 0.3 | No | 0 | 0 | 0 |
| 0.3 | Yes | 35 | 35 | 180 |
| 0.6 | No | 0 | 0 | 0 |
| 0.6 | Yes | 0 | 74 | 180 |
| 1.0 | No | 0 | 0 | 0 |
| 1.0 | Yes | 0 | 20 | 133 |

This reconstruction matches [the branch's benchmark](../experiments/multiwindow_case_study.py). The sign-aware and LP methods are not compared against the best available published set-membership implementation. Fault identities are not established by these edge counts.

## 3. Confirmed issue A: catastrophic loss of significance

Consider exact healthy sensors with \(T(t)=2^{52}+t\), \(\tau_i=(0,1,2)\), \(t=0,1,\dots,200\), no noise and zero curvature; their exact outputs are \(y_i(t)=2^{52}+t-\tau_i\), and each sample is exactly representable in binary64.

With three windows [0,50], [50,120], [120,200], the original uncentered computation integrated each \(\sim10^{15}\) sensor reading separately. The first window's exact integrated intersensor differences are 50 and 50 but the float calculation returned **64 and 64**. An independently reconstructed HiGHS LP reported **infeasible (status 2)** even though \((0,1,2)\) is a valid parameter assignment.

The draft-branch correction subtracts the **same reference sensor sample from every sensor within a window before integration**, a change of coordinate that preserves the model's exact feasibility. It yields first-window integrals [1250,1200,1150], and HiGHS returns feasible. The precision-probe script and regression test preserve this failure mode.

This improvement is not universally certified numerical behavior at arbitrary scales.

## 4. Confirmed issue B: old time-grid check underestimated quadrature error

The previous \`np.allclose(..., atol=1e-12)\` accepted the nonuniform grid \(t=s[0,1,2,8]\), \(s=2^{-44}\), whose largest gap is six times the first. The original bound used \(h=s\) even though the actual composite trapezoid has a larger error.

Let \(x=t/s\), \(A=1/s\), \(\tau=(0,s)\), \(T(t)=A(x^3-9x^2+8x)\), \(y_1=T\), and
\[
y_2(t)=A(x^3-12x^2+32x-32).
\]
Then \(s\,y_2'+y_2=T\) exactly. Both endpoint increments vanish; \(|y_i''|\le K=30/s^3\). The observed intersensor trapezoid difference is **−109**. The accepted-but-incorrect first-gap formula gives a bound of **40**, falsely rejecting a healthy pair under the stated physical model.

**Correction:** the code uses \(h_{\max}=\max_k(t_{k+1}-t_k)\), so
\[
 \Bigl|\mathcal T(z)-\int z\Bigr|
 \le\frac{K}{12}\sum_k(t_{k+1}-t_k)^3
 \le\frac{KLh_{\max}^2}{12}.
\]
The corrected pairwise bound is **1440** and does not reject the healthy example. The experiment uses extreme mathematical scales to expose a logical input-validation defect, not to model realistic thermometer data.

## 5. Previous AI-referee reports

The Gemini and Codex reports support the algebraic *necessary-condition* statements, warn against hardware-failure interpretation, and identify substantial antecedents. Gemini explicitly could not fetch GitHub files and ran different parameters; Codex did not pin a SHA and the pasted report was partial. They are **AI-assisted review inputs, not human peer review**. See [AI synthesis](ai_review_synthesis_2026-11-12.md).

## 6. Prior-art scope

This is a compact *integration* of established methods, not a proven novel theory:

- de Kleer & Williams (1987): conflict-based multiple-fault diagnosis.
- Loiez & Taillibert (1997): bounded analog diagnosis via integrated differential constraints.
- Tagawa & Ohta (1997): common-input thermocouple dynamics.
- Ingimundarson, Bravo, Puig & Alamo (2005): set-membership consistency tests for systems with invariant parameters and bounded uncertainty, a particularly direct antecedent to the shared-parameter multi-window extension.

The last paper's archived full text is available at https://skoge.folk.ntnu.no/prost/proceedings/cdc-ecc05/pdffiles/papers/2056.pdf . Its existence undermines claims that merely intersecting feasible parameter sets across time is an unprecedented innovation.

## 7. Remaining blockers

- No formally checked Farkas infeasibility certificate or directed-rounding interval arithmetic.
- No human independent proof and implementation review.
- No physical sensor experiment with verified homogeneous input and uncertainty calibration.
- No quantitative comparison against a properly implemented published set-membership method.
- No general fault identifiability or reliability guarantees.
- No test of all possible signals, faults, time grids, floating-point environments or model mismatch.

**Decision:** Keep PR #4 as a draft and do not merge into \`main\` or produce another Zenodo release until these blockers are addressed. The current synthetic improvement is scientifically meaningful as an engineering exercise, but its novelty and practical advantage are not established.
