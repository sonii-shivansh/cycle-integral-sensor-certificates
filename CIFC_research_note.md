# Cycle-Integral Fault Certificate (CIFC): a first-principles research candidate

**Status (11 November 2026):** Mathematical statement and proof under explicit assumptions, plus synthetic numerical verification. **Not verified as original in the worldwide literature and not experimentally validated.** The physics, mathematical tools, and parts of the approach have well-established precedents.

## 1. Problem

Several thermometers monitor the same time-varying temperature, but their response times differ and the true temperature is not observed directly. Can we **certify** that at least some sensors are incompatible with the shared-temperature model, without estimating the temperature or any particular response time?

## 2. Model

A healthy sensor `i` has a (possibly unknown) *constant* first-order lag `tau_i`:

`tau_i * (dy_i/dt)(t) + y_i(t) = T(t)`, with `0 <= tau_i <= M`.

Here `T` is the **same** continuously varying ambient temperature for all healthy sensors, and `M` is a trustworthy known upper bound. We observe synchronous samples `z_i(t_k)` at uniform spacing `h`, on an interval `[a,b]` of length `L`. All healthy measurements obey `|z_i(t_k)-y_i(t_k)| <= eps_i`. The underlying **healthy** response obeys `|y_i''(t)| <= K_i` for all times on the interval. Both the noise and curvature bounds must be physically justified.

No condition of periodicity is needed for the soundness theorem below.

## 3. Integral invariant

For healthy sensors `i,j`, subtract their differential equations and integrate:

`I_ij := integral_a^b (y_j(t)-y_i(t)) dt = tau_i [y_i(b)-y_i(a)] - tau_j [y_j(b)-y_j(a)]`.

The unknown shared input `T(t)` cancels **exactly**.

Let `Q_ij` be the composite-trapezoidal approximation to the integral of the *observed sampled difference* `z_j-z_i`. Let `d_i=z_i(b)-z_i(a)` and `d_j=z_j(b)-z_j(a)`. Put

`q_ij = (L*h^2/12)*(K_i+K_j)`.

## 4. Soundness theorem

**Theorem 1.** If both sensors `i,j` satisfy the model and all bounds, then

`|Q_ij| <= M (|d_i|+|d_j|) + (eps_i+eps_j)(L+2M) + q_ij`.

Consequently, **if the observed left-hand side exceeds the right-hand side, at least one sensor violates one or more assumptions of the healthy-sensor model**. In particular, no two genuinely healthy sensors can produce such a contradiction.

**Proof.** Let `E=eps_i+eps_j`. The absolute difference between the trapezoidal sum applied to noisy samples and the trapezoidal sum applied to the true outputs is at most `EL`, because the trapezoidal weights are nonnegative and sum to `L`. The classical composite-trapezoid bound gives error at most `q_ij` for the integral of the true difference. Each endpoint difference differs from its noiseless counterpart by at most `2 eps_i` or `2 eps_j`. Apply the exact integral invariant, then triangle inequality and `0<=tau_i,tau_j<=M`. QED.

This theorem is a **one-sided certificate**: a violation proves model inconsistency; the absence of a violation does not prove that a sensor is healthy.

## 5. Fault-count theorem

For `n` sensors, build a graph `G` with one vertex per sensor. Put an edge between `i,j` exactly when Theorem 1 produces a violation on **at least one** tested interval.

**Theorem 2.** If the graph contains a matching of `r` mutually vertex-disjoint edges, at least `r` distinct sensors violate the healthy-model assumptions.

**Proof.** Every contradiction edge has at least one unhealthy endpoint by Theorem 1. A matching consists of edges with no shared endpoints; therefore one sensor cannot explain violations on two matched edges. At least one distinct unhealthy sensor is required for each of the `r` edges. QED.

A sharper purely combinatorial bound is the size of a *minimum vertex cover* of `G`: the actual set of unhealthy sensors must be a vertex cover. Computing the minimum is combinatorially hard in general. For six sensors it can be found by brute force; for larger networks, integer programming or graph algorithms can be used. A unique minimum cover is a *smallest explanatory hypothesis*, **not** automatically a list of certainly faulty devices unless an independent upper bound on the number of faults is available.

## 6. Detection guarantee in a repeating environment

Suppose the underlying healthy outputs are exactly periodic over the chosen window, so `y_i(b)=y_i(a)` and `integral_a^b (y_j-y_i) = 0`. Suppose the observed outputs have fixed unknown biases `beta_i,beta_j`, plus sample noise with the stated bounds.

**Theorem 3.** The pair is guaranteed to be flagged for every allowed noise realization if

`|beta_j - beta_i| > 2E + (4 M E + 2 q_ij)/L`, where `E=eps_i+eps_j`.

**Proof.** The trapezoidal integral of the smooth unbiased difference differs from zero by at most `q_ij`. A constant relative bias contributes exactly `L (beta_j-beta_i)` to the trapezoidal sum. Sample noise changes that sum by at most `EL`, giving `|Q_ij| >= L|beta_j-beta_i|-q_ij-EL`. The observed endpoint changes are at most `2eps_i` and `2eps_j`. Theorem 1's right side is therefore at most `4ME+LE+q_ij`. The strict inequality above makes the lower bound for `|Q_ij|` exceed this upper bound. QED.

## 7. Exact synthetic physics and numerical validation

Choose `T(t)=25+8 sin(2*pi*t/20)` degrees Celsius. A steady-state sensor of time constant `tau` has the exact output

`y_tau(t) = 25 + 8 [sin(omega*t) - omega*tau*cos(omega*t)]/[1+(omega*tau)^2]`, with `omega=2*pi/20`.

Six physical time constants are `[1,2,3,5,8,13]` seconds; all are below `M=15` seconds. We sample at `h=0.05` seconds for a full `L=20` seconds, with independently drawn uniform noise bounded by `eps_i=0.02` degrees Celsius. Bias sensor **3** by `+1.2` degrees and sensor **6** by `-1.5` degrees. A conservative curvature cap for every healthy signal is `K=8*omega^2 = 0.7895683521` degrees Celsius per second squared.

For a pair, `q <= 0.0065797363` degree-seconds. Theorem 3 guarantees detection of any fixed *relative* bias strictly larger than `0.2006579736` degrees Celsius, in the full-cycle setting with these bounds.

The runnable script `cycle_integral_fault_certificate.py` generated the following result (seed 123):

- 9 of 15 pairs contradicted the healthy model; exactly the pairs with sensor 3 or sensor 6.
- No healthy-vs-healthy pair was rejected.
- Maximum matching size: **2**, certifying at least **two** incompatible sensors.
- Unique minimum vertex cover: **{3,6}** (sensor numbers start at 1).
- In another 100 separately generated all-healthy, bounded-noise synthetic trials, no contradiction was raised.

These numerical results confirm that the code behaves as predicted on this constructed test. **The proofs, not the simulations, establish the mathematical guarantees.** There has been no physical hardware experiment.

## 8. An impossibility result

An identical bias added to *all* sensors cannot be identified from these traces alone. Given any healthy solution `(y_i,T)`, define `y_i^*(t)=y_i(t)+c` and `T^*(t)=T(t)+c`. Then

`tau_i (dy_i^*/dt)+y_i^*=T^*`.

Thus a world with a universally offset set of sensors is observationally equivalent to a world whose ambient temperature is truly shifted by that same amount. **Without an external absolute reference, the two worlds cannot be distinguished by any algorithm relying only on the shared-response model.** This is an information-theoretic limitation, not an implementation defect.

## 9. Real-world limitations

- Must know justified error bounds, a trustworthy `M`, and credible curvature bounds `K_i`. Violating these assumptions may create false alarms relative to real hardware health.
- Healthy sensors must measure essentially the same `T(t)` and share a synchronized clock; spatial gradients, clocks, conduction/radiation effects, and changing time constants can invalidate the model.
- The first-order constant-lag approximation may be unsuitable for some sensors or environments.
- The method detects violations of the specified model; it does not automatically distinguish a defective instrument from an incorrect environmental model.
- Common-mode errors are impossible to detect from these traces alone. Slowly varying or small differential faults may remain undetected.
- No alarms do not prove absence of faults. A unique minimum vertex cover is not automatically proof of which individual sensors failed.

## 10. Novelty audit, not a claim of worldwide priority

The mathematical derivation and proof of the specific test here were constructed for this research proposal, but an exhaustive novelty claim is **not** justified. Essential nearby areas are established:

- Tagawa and Ohta, 1997, *Two-thermocouple probe for fluctuating temperature measurement*, https://www.sciencedirect.com/science/article/pii/S0010218097000448 — simultaneous in-situ time-constant and temperature reconstruction.
- Li et al., 2018, *Blind system identification of two-thermocouple sensor based on cross-relation method*, https://pubmed.ncbi.nlm.nih.gov/29604754/ — common-input blind identification.
- Bravo, Alamo and Camacho, 2006, *Bounded Error Identification of Systems With Time-Varying Parameters*, https://idus.us.es/items/6bf3a897-a6bf-459d-a771-6396a5f0a16a — guaranteed bounded-error parameter identification.
- *Symmetric PMC model of diagnosis, b-matchings in graphs and fault identification*, 2021, https://doi.org/10.1016/j.tcs.2021.08.025 — graph matching and vertex-cover concepts in fault diagnosis.
- *Spectral-Bounded Consistency Certification for Heterogeneous Wireless Sensor Networks*, 2026, https://doi.org/10.1109/jiot.2026.3685152 — related consistency certificate research.

**Conclusion:** a coherent, useful *candidate engineering algorithm*, with soundness, fault-count, and differential-bias sensitivity guarantees under explicit assumptions. Its precise novelty requires systematic full-text prior-art and patent review, and its practical usefulness requires real-sensor testing and comparison with established fault-detection algorithms. This document does **not** declare it to be an unprecedented scientific discovery.