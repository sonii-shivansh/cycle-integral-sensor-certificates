# Multi-window LP validation: exploratory research branch

**Date:** 12 November 2026. **Status:** Mathematical outer-relaxation prototype plus *synthetic* numerical checks; **not peer reviewed**, **not hardware validated**, and **not established as novel**. This document does not replace the fixed Zenodo v0.1.0 snapshot or the original single-window implementation.

## Motivation and scientific context

AI-generated referee reports (one Codex report pasted by the repository owner, and an attached 24-page Gemini audit) independently support the original integral identity and conservative sign-aware inequality as *necessary* conditions, but emphasize severe practical and identifiability limitations. They are **AI checks, not independent human peer review**. Gemini did not directly inspect the repository and ran a different six-sensor simulation. Its proposed multi-window shared-parameter LP is related to longstanding set-membership and bounded-error diagnosis literature; **it must not be presented as a new paradigm**.

The core issue is that separately testing each window permits a different implicit response time on each interval even though the physical model assumes that every healthy sensor has one constant parameter across all intervals. This branch explores the standard joint-feasibility strategy.

## Derived necessary linear constraints

For each sensor \(i\) and observation window \(w\) of length \(L_w\) and sample step \(h_w\), assume

\[
\tau_i \dot y_i+y_i=T(t), \quad 0\le \tau_i \le M;
\quad |z_i(t_k)-y_i(t_k)|\le\epsilon_i,\quad |y_i''(t)|\le K_i.
\]

The **same** \(\tau_i\) must apply on every window, and healthy sensors must share \(T(t)\). Let \(I_{i,w}\) be the composite trapezoid integral of measured \(z_i\), \(d_{i,w}=z_i(b_w)-z_i(a_w)\), and \(J_w=\int_{a_w}^{b_w}T(t)\,dt\).

Define a conservative uncertainty radius:

\[
D_{i,w}=L_w\epsilon_i+\frac{L_w h_w^2K_i}{12}.
\]

A necessary condition for a healthy sensor is

\[
\boxed{I_{i,w}+\tau_i(d_{i,w}-2\epsilon_i)-D_{i,w}
\;\le\; J_w \;\le\;
I_{i,w}+\tau_i(d_{i,w}+2\epsilon_i)+D_{i,w}.}
\]

This is **linear in** the unknowns \(\tau_i\) and \(J_w\), because \(d_{i,w}\), \(D_{i,w}\) and \(\epsilon_i\) are fixed observed/bounded inputs. The LP has \(n+W\) unknowns and \(2nW\) inequality constraints, plus \(n\) bounds \(0\le\tau_i\le M\).

### Why infeasibility is sound in exact arithmetic

For truly healthy sensors, choose their actual constant time constants \(\tau_i^*\) and actual integrated common inputs \(J_w^*=\int T\). The integrated first-order ODE gives

\[
J_w^* = \int y_i+\tau_i^* (y_i(b_w)-y_i(a_w)).
\]

Bounded endpoint errors place the true endpoint increment inside \([d_{i,w}-2\epsilon_i,\ d_{i,w}+2\epsilon_i]\). Bounded sample errors and the composite-trapezoid estimate place the true integral within \([I_{i,w}-D_{i,w},\ I_{i,w}+D_{i,w}]\). Since \(\tau_i^*\ge0\), the exact variables satisfy *all* the displayed LP inequalities. Therefore **an exactly infeasible LP excludes the simultaneous healthy-model assumptions**.

**Critical qualification:** This is an **outer relaxation**, not an exact characterization of whether measured traces satisfy the original ODEs. Endpoint and integral errors, overlapping windows, and cross-sensor trajectory dynamics remain correlated. A feasible solution does not prove sensor health. SciPy's HiGHS reports numerical infeasibility, which is **not automatically a verified mathematical certificate** near solver tolerances. The code reports \`solver_infeasible\`, not \`proven_defective\`.

### Why multi-window can be stronger than single-window

A joint feasible parameter must satisfy the constraints from every window simultaneously, so its feasible set is the intersection of per-window feasible sets. When a pair is contradictory in any single window, the joint pairwise LP must also be infeasible (in exact arithmetic). The reverse need not hold: windows may individually admit different response times while their intersection is empty.

This is a standard set-membership idea. Stronger detection relative to *our own* weak single-window test does **not** establish superiority over published set-membership, observer, or interval-diagnosis methods.

## Experiment 1 — Changing time constant (two windows)

Use the common input \(T(t)=0\) and two sensors, with sensor 2 instantaneous and constant at zero. Sensor 1 starts at \(y_1(0)=10\) and evolves exponentially, with \(\tau_1=0.5\) seconds on \(t\in[0,0.5]\) and \(\tau_1=2.0\) seconds on \(t\in[0.5,1.5]\). Its response is continuous but changes dynamics at \(t=0.5\), violating the assumed *globally constant* \(\tau_1\). Each window separately has a valid first-order response. Set \(M=3\), \(\epsilon_1=\epsilon_2=0\), \(K_1=40\), \(K_2=0\), and \(h=0.01\) seconds; curvature bounds hold **within** each window.

For sensor 1, because sensor 2 has \(I_{2,w}=d_{2,w}=D_{2,w}=0\), the LP reduces to \(|I_{1,w}+\tau_1d_{1,w}|\le D_{1,w}\). The implied feasible intervals for \(\tau_1\) from observed trapezoid data are:

- Window 1: approximately **[0.49999030, 0.50004303]** seconds.
- Window 2: approximately **[1.99977388, 2.00023445]** seconds.

Both windows pass the individual sign-aware test. Those intervals are disjoint, so the joint feasibility relaxation is rejected. These numbers use finite-precision samples, but the separation exceeds numerical quadrature and solver error by orders of magnitude.

**This demonstrates sensitivity to a time-constant change**, not detection of a uniquely defective device. The joint constant-parameter model is intentionally violated.

## Experiment 2 — Constant time constants stay feasible

Two exponential responses to a common zero input, with constant \(\tau=(0.5,1.2)\), are observed in three windows. The LP remains feasible with curvature bounds \(K=(40,7)\), \(M=2\), zero sample noise and \(h=0.01\). Input-validation tests reject nonfinite data, duplicate/nonuniform timestamps, and negative uncertainty bounds.

## Experiment 3 — Six sensors, matched windows, 20 seeded trials per setting

Using the same synthetic six-sensor dynamics as the original nonperiodic experiment, split the **same** 13.5-second trace into four windows: 0–3, 3–6, 6–9, and 9–13.5 seconds. Every method receives the same sample times, data, bounds (\(M=15,\epsilon=0.02,K=8(2\pi/20)^2\)), and windows. Add fixed biases +1.2°C (sensor 3) and −1.5°C (sensor 6) for the injected-bias conditions. Each drift condition uses 20 fixed seeds, with *matched noise realizations* in the healthy and injected-bias data sets (the same seed is reused for each pair of conditions).

A pair is flagged by the original or sign-aware baseline if **any** window violates that inequality. For the LP, a pair is flagged if the **same pair's** response-time variables cannot satisfy all four windows (HiGHS status \`solver_infeasible\`).

| Drift °C/s | Biases injected? | Original: total conflicting edges | Sign-aware: total edges | Multi-window LP: total edges |
|---:|:---:|---:|---:|---:|
| 0.0 | No | 0 | 0 | 0 |
| 0.0 | Yes | 16 | 58 | 180 |
| 0.3 | No | 0 | 0 | 0 |
| 0.3 | Yes | 35 | 35 | 180 |
| 0.6 | No | 0 | 0 | 0 |
| 0.6 | Yes | 0 | 74 | 180 |
| 1.0 | No | 0 | 0 | 0 |
| 1.0 | Yes | 0 | 20 | 133 |

Totals sum edges over 20 trials: one trial has 15 possible sensor-pair edges. The LP found at least one contradictory pair in **all 20 injected-bias trials at each drift**, even at 1.0°C/s. None of the 80 healthy synthetic trials produced a contradiction.

These are **descriptive simulated results**, not estimated real-world probabilities. The LP has more information than a one-window standalone method, but the two single-window baselines also use all four windows independently. That makes this a fair comparison for the *specific assumptions*, though it is **not** a comparison with state-of-the-art published set-membership tools.

Importantly, contradictions are **not equivalent to correct fault localization**; no claim of uniquely identifying both biased sensors is made.

## Reproduce locally

From the repository root:

~~~bash
python -m pip install -r requirements-research.txt
python -m unittest discover -s tests -v
python -m experiments.multiwindow_case_study
~~~

An additional fixed 5-seed regression test runs in CI; the full command computes 20-seed results. Use a fresh environment to record Python, NumPy, SciPy and HiGHS versions.

## Next review criteria

1. Obtain a human expert assessment of the exact LP soundness argument, uncertainty assumptions, and reproducibility.
2. Construct verified infeasibility certificates using rational arithmetic or certified LP bounds for cases close to numerical tolerance.
3. Study overlapping-window coupling, endpoint/integral error correlation, and parameter-identifiability limitations.
4. Compare with a published interval/set-membership fault-detection implementation under matched information.
5. Test calibrated hardware with justified noise and curvature bounds; a common temperature input is essential.
6. Preserve all failures and limitations in the permanent record.

**Nothing in this branch changes the 2026 Zenodo v0.1.0 archive.**
