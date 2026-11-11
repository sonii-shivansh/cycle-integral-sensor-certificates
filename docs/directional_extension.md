# Sign-aware interval extension and nonperiodic stress test

**Status: 11 November 2026.** This note documents a *follow-up mathematical refinement and reproducible synthetic benchmark*, not a peer-reviewed paper, an independent validation, or an originality claim. Its ingredients are basic interval reasoning and the original first-order sensor model; close relationships to published set-membership/interval methods remain to be examined.

## Why revisit the original bound?

The original CIFC inequality uses an absolute-value bound on two endpoint changes. This throws away their **signs**, even though the response times are nonnegative. Retaining signs yields a narrower feasible interval for the observed integrated difference, without changing the underlying physical assumptions.

## Claim: a sign-aware necessary interval

Let each healthy sensor satisfy

\[
\tau_i y_i'(t)+y_i(t)=T(t),\quad 0\le\tau_i\le M,
\]

with *the same* unknown input \(T(t)\). On a uniform grid \(t_k\in[a,b]\) of step \(h\), let the observed samples obey \(|z_i(t_k)-y_i(t_k)|\le\epsilon_i\) and the true responses satisfy \(|y_i''(t)|\le K_i\). Set \(L=b-a\), \(d_i=z_i(b)-z_i(a)\), \(E=\epsilon_i+\epsilon_j\), and let \(Q_{ij}\) be the composite-trapezoid sum over samples of \(z_j-z_i\).

Define
\[
  \ell_i=M\min(0,d_i-2\epsilon_i),\quad
  u_i=M\max(0,d_i+2\epsilon_i),
\]
and likewise for sensor \(j\). Put
\[
D=LE+\frac{Lh^2(K_i+K_j)}{12}.
\]

**Proposition (one-sided, sign-aware certificate).** If both sensors satisfy the joint healthy-model assumptions, then

\[
\boxed{\ell_i-u_j-D \le Q_{ij}\le u_i-\ell_j+D}.
\]

**Proof.** Because endpoint measurement errors are bounded, the true endpoint change \(\Delta y_i\) lies in \([d_i-2\epsilon_i,d_i+2\epsilon_i]\). Multiplying an unknown member of this interval by a nonnegative \(\tau_i\in[0,M]\) produces a number within \([\ell_i,u_i]\). The exact integrated identity
\[
\int_a^b(y_j-y_i)\,dt=\tau_i\Delta y_i-\tau_j\Delta y_j
\]
therefore lies in \([\ell_i-u_j,u_i-\ell_j]\). The trapezoidal quadrature error plus bounded noise changes the observed trapezoid sum by at most \(D\). Expanding the interval by \(D\) proves the claim.

**Strength comparison.** For the *same bounds*, this feasible interval is contained in the symmetric interval permitted by the original CIFC inequality
\[
|Q_{ij}|\le M(|d_i|+|d_j|)+2ME+LE+q,\quad q=\frac{Lh^2(K_i+K_j)}{12}.
\]
Indeed, \(\ell_i\ge -M(|d_i|+2\epsilon_i)\), \(u_i\le M(|d_i|+2\epsilon_i)\), and similarly for \(j\). Thus **every pair rejected by the original test is also rejected by the sign-aware test**, apart from numerical-tolerance edge cases. It can reject additional pairs; this does **not** guarantee correct physical fault identification or detection of all bias values.

The repository's [Python implementation](../cycle_integral_fault_certificate.py) uses a common \(M,\epsilon,K\) across sensors. This theorem allows differing sensor-specific error and curvature bounds, but that generality is not implemented.

## Reproducible nonperiodic-window experiment

Run from the repository root:

~~~bash
python -m pip install -r requirements.txt
python -m experiments.nonperiodic_benchmark
python -m unittest discover -s tests -v
~~~

The experiment generates **exact analytical first-order responses** (not approximations to the ODE) to a shared input

\[
T(t)=25+v t+6\sin(0.25t+0.3)+2\cos(0.40t-0.2).
\]

Observations cover \(0\le t\le13.5\) seconds at \(h=0.05\). All six healthy time constants are \([1,2,3,5,8,13]\) seconds and satisfy \(M=15\). The curvature bound \(6(0.25)^2+2(0.40)^2=0.695\) is below the original script's allowed \(K\). For each condition we generate **100 deterministic seeded bounded-noise arrays**, with both all-healthy measurements and measurements containing biases **+1.2°C at sensor 3** and **−1.5°C at sensor 6**. The \(v=0\) signal is periodic globally, but the *window is not one full period*; positive \(v\) also makes the input genuinely nonperiodic.

| Drift \(v\) (°C/s) | Healthy contradictions, original/new (total across 100 trials) | Injected-bias contradictions, original | Injected-bias contradictions, sign-aware | Trials with ≥1 alarm, original/new |
| --- | ---: | ---: | ---: | ---: |
| 0.0 | 0 / 0 | 200 | 500 | 100 / 100 |
| 0.3 | 0 / 0 | 0 | 218 | 0 / 100 |
| 0.6 | 0 / 0 | 0 | 0 | 0 / 0 |
| 1.0 | 0 / 0 | 0 | 0 | 0 / 0 |

These are **counts of contradiction edges**, not accuracy percentages or estimates of a false-alarm probability. For \(v=0.3\), the seed-123 sign-aware graph has minimum cover \(\{6\}\), even though two sensors were deliberately biased—again illustrating that minimum covers are explanations, **not verified fault identities**.

As a separate numerical check, the largest discrepancy between composite-trapezoid evaluation of the noiseless sensor integral and the exact ODE identity across 15 sensor pairs was **0.000395147°C·s**, below the theoretical quadrature bound **0.00444132°C·s**.

## What we learned

- **The new inequality is mathematically at least as sensitive as the original bound** for the same set of assumptions.
- In this controlled nonperiodic-window experiment it detects additional modeled inconsistencies.
- It still **misses both injected fixed biases at higher drift**, revealing a severe sensitivity limitation for this chosen set of bounds.
- The naive full-period-average comparator from our first experiment is not valid for this partial, drifting interval. It cannot be used as a fair general baseline without further assumptions.
- The nonperiodic demonstration is **not** evidence of superiority to existing interval-consistency, analytic-redundancy or set-membership methods, which may already subsume this refinement.

## Next scientific checks

1. Verify both interval propositions independently, including degenerate \(\tau=0\), adversarial bounded errors, and non-smooth sample noise.
2. Compare against **published** bounded-error and parity-relation methods, not only the original CIFC test.
3. Explore tighter *joint* feasible-set tests retaining cross-time correlations and parameter consistency, while carefully measuring computational cost.
4. Test calibrated real measurements with justified sampling, curvature and noise bounds.
5. Record negative findings and correct papers/notes as needed; do not imply a new scientific discovery on the strength of algebra alone.
