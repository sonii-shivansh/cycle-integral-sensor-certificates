# Identifiability barrier: sensor bias versus local temperature differences

**Research date:** 12 November 2026. **Status:** exact model-level indistinguishability proof, controlled simulations, and conditional uncertainty extension. **Not a novel theorem or validated hardware classifier.** All code belongs to draft PR #4; the original Zenodo v0.1.0 archive and GitHub main remain unchanged.

## 1. The scientific question

CIFC certificates assume healthy first-order sensors experience exactly the same temperature \(T(t)\). But physical sensors need not be perfectly co-located, thermally equilibrated, or synchronized. A pairwise or network model contradiction may represent:

- A measurement-system offset;
- A genuine *local temperature difference* between sensing locations;
- A wrong response-time, curvature, noise, or timing assumption;
- Some combination of the above.

Can the recorded trajectories by themselves discriminate constant sensor bias from a constant local environmental temperature difference?

**Answer: no, not under the original first-order observation model. The two mechanisms can produce exactly the same data for all times.**

## 2. Exact observational-equivalence proposition

Assume an ideal response satisfies

\[
\tau_i y'_i(t)+y_i(t)=T(t)
\]

for an unknown shared input \(T\), a constant \(\tau_i\ge0\), and a sensor-specific constant \(b_i\). Consider two models:

**World A — biased reported measurements.** All sensors are exposed to the **same** \(T(t)\), but the recorded output is
\[
z_i^A(t)=y_i(t)+b_i+e_i(t),
\]
where \(e_i\) is the remaining bounded observational noise.

**World B — physically different local inputs.** Each sensor is correctly calibrated, but it experiences \(T_i(t)=T(t)+b_i\). The corresponding physical response is \(y_i^B(t)=y_i(t)+b_i\) (its initial condition is shifted by \(b_i\)), and it records
\[
z_i^B(t)=y_i^B(t)+e_i(t).
\]

**Proposition.** These two worlds produce exactly identical observations, for **every** \(t\) and every choice of observation noise \(e_i\) shared between the two hypothetical worlds:
\[
\boxed{z_i^A(t)=z_i^B(t)\quad\text{for all }t.}
\]

**Proof.** A constant offset has derivative zero. Therefore

\[
\tau_i\frac{d}{dt}(y_i+b_i)+(y_i+b_i)
=\tau_i y'_i+y_i+b_i
=T(t)+b_i=T_i(t).
\]

Thus the shifted response obeys the local-input model. By construction \(z_i^B=y_i+b_i+e_i=z_i^A\). QED.

The equivalence is even stronger than an equality of summary statistics: **the complete sampled (or continuous) readings are identical**. No algorithm whose input consists only of those readings, regardless of runtime, number of windows, or mathematical sophistication, can universally and correctly choose between World A and World B. An additional measurement, independently established environmental homogeneity, or another external assumption is necessary.

The constant offsets must be present in the physical response/initial condition in World B; abruptly changing the local environment without allowing the response to settle is a different experiment and need not be observationally identical. The theorem assumes the stated first-order model and constant offsets.

**Scientific prior art caution:** separation of sensor faults and process/environment changes is a long-recognized fault-diagnosis issue. Relevant background includes [Kroschel & Wernz, 1990, sensor fault detection and localization](https://doi.org/10.1016/0924-4247(90)87007-6) and [NASA TM-2012-217586 on identifiability of bias-type sensor/actuator faults](https://ntrs.nasa.gov/citations/20120011947). These sources address the broader identifiability problem; we do **not** claim they present exactly the same shifted-temperature example word-for-word.

## 3. A conservative local-input mismatch model

Generalize each sensor's physical input to

\[
\tau_i y'_i+y_i=T(t)+u_i(t),
\quad |u_i(t)|\le\gamma_i
\quad\text{for every }t\text{ in every window}.
\]

Here \(\gamma_i\ge0\) is a **pre-calibrated environmental bound**, in temperature units, on the local difference from a shared hypothetical reference temperature \(T(t)\). It is **not** the measurement error \(\epsilon_i\). The constants \(\gamma_i\) must be supplied from **independent environmental evidence**: they cannot be chosen after looking at fault labels just to make all discrepancies disappear.

Let \(I_{i,w}\) be the observed (common-reference-centered) composite-trapezoid integral, \(d_{i,w}\) the measured endpoint increment, \(L_w\) the window length, \(h_{\max,w}\) its largest sampling gap, and \(|y_i''|\le K_i\). Define

\[
D_{i,w}^{\mathrm{environment}}=
 L_w(\epsilon_i+\gamma_i)
 +\frac{L_wh_{\max,w}^2K_i}{12}.
\]

**Sound necessary LP bounds.** For every truly healthy sensor complying with this more permissive physical model, there must be constants \(\tau_i\in[0,M]\) and common window integrals \(J_w=\int T\) such that

\[
\boxed{
I_{i,w}+\tau_i(d_{i,w}-2\epsilon_i)-D_{i,w}^{\mathrm{environment}}
\le J_w\le
I_{i,w}+\tau_i(d_{i,w}+2\epsilon_i)+D_{i,w}^{\mathrm{environment}}.
}
\]

**Proof.** Integrating the local-input ODE gives
\(\int y_i+\tau_i\Delta y_i=J_w+\int u_i\), hence
\(J_w=\int y_i+\tau_i\Delta y_i-\int u_i\).
The usual trapezoidal and bounded-measurement-error estimates contribute at most \(L_w\epsilon_i+L_wh_{\max,w}^2K_i/12\) to the integral. Endpoint error changes contribute at most \(2\epsilon_i\) to \(d_{i,w}\); because \(\tau_i\ge0\), they enter the displayed LP as signed linear terms. Finally, \(|\int u_i|\le L_w\gamma_i\). Combining bounds proves both inequalities.

The sign-aware pairwise certificate has the analogous extra additive radius \(L_w(\gamma_i+\gamma_j)\).

### Monotonicity and cost of more environmental tolerance

If \(\gamma_i'\ge\gamma_i\) for every sensor, every accepted feasible assignment at \(\gamma\) is still feasible at \(\gamma'\). Thus increasing allowed environmental mismatch can only **relax** the exact LP model. No additional *mathematically sound* contradictions can appear from increasing \(\gamma\); some previous contradictions may disappear. Independently constraining sensor environments is essential for meaningful fault detection.

This is a **conservative uncertainty-model extension** using standard interval/set-membership reasoning, not a new physical law. Feasibility remains insufficient to prove physical model correctness or hardware health.

## 4. Reproducible two-world stress test

[The experiment](../experiments/environmental_ambiguity_study.py) uses the same six-sensor, 13.5-second nonperiodic-window response family as the earlier CIFC benchmark, with time constants \((1,2,3,5,8,13)\) seconds and drift \(v=0,0.3,0.6,1.0\) °C/s. The ideal shared input combines a temperature ramp and two harmonics. The code analytically verifies the underlying first-order ODE residual.

The observed readings for each world are **exactly equal at every sampled time**. We use the same bounded random noise array in both constructed worlds. The constant offsets, in °C, are:

\[
b=(0,\;0,\;+1.2,\;0,\;0,\;-1.5).
\]

- **World A:** those are measurement biases on sensors 3 and 6.
- **World B:** those are genuine constant *local temperature differences*, and all six sensors are physically healthy in their corresponding local inputs.

A strict common-input CIFC model has \(\gamma=0\). The relaxed environmental model uses \(\gamma_i=|b_i|\), deliberately *chosen from known simulation truth* to demonstrate the effect of physical uncertainty.

**Critical experimental limitation:** In a real deployment \(|b_i|\) would **not be known from the sensor readings themselves**. Using those exact values here is a constructed comparison, **not a validated method for estimating environmental bounds**.

Run:

~~~bash
python -m pip install -r requirements-research.txt
python -m unittest discover -s tests -p 'test_environmental_ambiguity.py' -v
python -m experiments.environmental_ambiguity_study
~~~

The tests also contain an **exact Fraction-rational ramp construction** (not dependent on floating-point simulation), explicit local-input ODE residual checks, monotonicity tests and independent rational feasible witnesses for the relaxed bound.

A strict-model alarm in World B indicates that the common-input assumption was violated. It is **not** evidence of broken hardware. Conversely, increasing \(\gamma\) may make an actual measurement bias in World A no longer detectable. This is the unavoidable **sensitivity versus modeling-robustness tradeoff**.

## 5. What would count as real physical validation?

A useful controlled data collection protocol should **predefine** these independently observable quantities:

1. **Sensor transfer-model calibration:** estimate each response-time interval \([0,M]\) or sensor-specific range, curvature bound \(K_i\), measurement uncertainty \(\epsilon_i\), and timing/synchronization uncertainty *without using the same fault labels later scored for detection*.
2. **Environmental uniformity:** measure, using independent reference thermometry, whether all sensing locations see sufficiently similar temperatures. If not, calibrate conservative local-input spread bounds \(\gamma_i\) *before* evaluating test records.
3. **Separate experimental conditions:** record (a) co-located healthy sensors, (b) healthy sensors at known different local temperatures, and (c) known injected electronic/calibration biases while separately monitoring environmental homogeneity.
4. **Meaningful evaluation:** report model-conflict detection **separately** from correct fault classification, false alarms when assumptions hold, sensitivity to bound estimation, CPU time, and inconclusive exact proofs.
5. **Independent cause labels:** electronic self-tests, an external reference thermometer, or independently controlled spatial conditions are necessary to distinguish the two observationally identical causes in this note.

No such calibrated physical dataset has been supplied or verified here. Synthetic "ground truth" does not substitute for independent hardware measurements.

## 6. Implementation scope and remaining research

- [\`joint_relaxation(..., input_spread=...)\`](../experiments/multiwindow_lp.py) and the sign-aware helper now use the \(\gamma_i\) margin; the default is **0**, retaining previous results and assumptions.
- [\`exact_network_rows(..., input_spread=...)\`](../experiments/network_exact_certificate.py) also supports the model, so exact-rational LP constraints can be audited with environmental allowance.
- [Tests](../tests/test_environmental_ambiguity.py) check exact observational equivalence, sound retained feasibility, invalid assumptions, and the original-versus-relaxed comparison.
- **No demonstrated novelty:** the addition is a straightforward known robust set-membership technique. It offers a way to **express** environmental uncertainty, not to **discover** its correct magnitude from raw observations.
- The physical and calibration gaps still block claims of real-world fault identity, industrial reliability, or publication-worthy originality.

**Decision:** Keep PR #4 draft and the archived release unchanged. The strongest immediate research priority is independent environmental/measurement calibration and a real dataset that includes both causes, rather than additional hand-selected synthetic fault simulations.
