# Exact-rational verification of pairwise multi-window LP contradictions

**Research date:** 12 November 2026. **Status:** independent numerical-verification prototype; **not peer reviewed, not physically validated, and not claimed to be mathematically novel**.

## Why this exists

Our earlier SciPy/HiGHS-based multi-window implementation could return \`solver_infeasible\` when mathematically healthy sensors generated extremely large absolute measurements. Reference centering and conservative nonuniform-grid quadrature estimates repaired two tested failure modes, but **no finite set of regression tests can make arbitrary floating-point LP infeasibility a proof**.

This follow-up adds a *separate*, simple two-variable exact-rational feasibility checker at [experiments/exact_pairwise_certificate.py](../experiments/exact_pairwise_certificate.py). It reconstructs the necessary LP constraints **from raw recorded samples and timestamps**, performs Fourier–Motzkin elimination using Python's \`fractions.Fraction\`, and produces independently checkable **Farkas multipliers** if the relaxed pairwise LP is infeasible.

The checker neither calls SciPy nor relies on HiGHS decisions. It is deliberately specialized to **two sensors**, for no more than 24 windows, and is not efficient for large problems.

## What exactly is being proved?

Let the supplied sensor sample values \(z_i(t_k)\), timestamps \(t_k\), bounds \(\epsilon_i,K_i,M\) be interpreted as exact rational numbers. For inputs specified by binary64 floating-point numbers, this means treating **their actual stored binary64 values as exact**, not pretending a rounded number equals an intended decimal or a more precise unknown physical measurement.

Let \(I_{i,w}\) be the **exact** composite trapezoid of recorded sensor samples for window \(w\), after subtracting one **common** arbitrary reference value from *all sensors* on that window. Let \(d_{i,w}=z_i(b_w)-z_i(a_w)\), \(L_w=b_w-a_w\), \(h_{\max,w}\) the largest sampling gap, and
\[
D_{i,w}=L_w\epsilon_i+\frac{L_w h_{\max,w}^2K_i}{12}.
\]
The common-offset subtraction only changes both the observed integrals \(I_{i,w}\) and the latent integrated input \(J_w\) by the same amount. It does not change feasibility.

Under the assumed first-order model, any healthy pair \(i,j\) must admit *constant* \(\tau_i,\tau_j\in[0,M]\) on all windows, and, for each window, some common integrated input \(J_w\) satisfying
\[
I_{i,w}+\tau_i(d_{i,w}-2\epsilon_i)-D_{i,w}\le J_w
\le I_{i,w}+\tau_i(d_{i,w}+2\epsilon_i)+D_{i,w}
\]
and the analogous inequalities for \(j\). Eliminating \(J_w\) is exact **for these relaxed inequalities** and gives the two pairwise constraints
\[
(d_{i,w}-2\epsilon_i)\tau_i-(d_{j,w}+2\epsilon_j)\tau_j
\le I_{j,w}-I_{i,w}+D_{i,w}+D_{j,w},
\]
\[
-(d_{i,w}+2\epsilon_i)\tau_i+(d_{j,w}-2\epsilon_j)\tau_j
\le I_{i,w}-I_{j,w}+D_{i,w}+D_{j,w}.
\]
Add the four box constraints \(0\le\tau_i,\tau_j\le M\). All coefficients are rational. There are \(2W+4\) inequalities in two unknowns.

**A mathematical LP infeasibility certificate** is a nonnegative vector \(\lambda\) (one multiplier per inequality) satisfying:
\[
\sum_k\lambda_k a_k=0,\quad
\sum_k\lambda_k b_k=0,\quad
\boxed{\sum_k\lambda_k c_k<0},
\]
where inequality \(k\) has the form \(a_k\tau_i+b_k\tau_j\le c_k\).

Multiplying all constraints by their \(\lambda_k\) and adding forces \(0\le\sum_k\lambda_k c_k<0\), an exact contradiction. This is the elementary Farkas lemma. Every returned \`FarkasCertificate\` is independently rechecked by the implementation with rational arithmetic.

### How the checker finds such certificates

Fourier–Motzkin elimination combines an inequality with a positive coefficient and one with a negative coefficient using **nonnegative** multipliers to cancel that variable. The checker retains the combination weights back to its original inequalities. Eliminating both \(\tau_i\) and \(\tau_j\) leaves inequalities \(0\le c\). If any such \(c\) is negative, the stored weights constitute the Farkas certificate above.

Unlike simply checking SciPy's numerical status, this establishes mathematical infeasibility **of the rational LP represented by the input values**.

A \`feasible_relaxation\` output only says no exact LP contradiction exists; it does **not** certify sensor health or even full differential-equation feasibility.

## A one-window exact contradiction example

For \(t=(0,1)\), sensor observations \(z_i=(0,1)\), \(z_j=(2,2)\), and bounds \(M=\epsilon_i=\epsilon_j=K_i=K_j=0\), the integrated model allows no lag and no uncertainty. The measured difference is incompatible with a shared input. The checker produces a rational combination equivalent to
\[
\boxed{0\le -2}.
\]
It can be independently verified by \`certificate.check(inequalities)\`.

## Healthy numerical stress tests

1. **Huge common offset:** exactly representable binary64 \(T(t)=2^{52}+t\), \(\tau=(0,1,2)\), with three windows and zero uncertainty. The *old floating implementation* reported infeasible due to catastrophic cancellation. The exact rational checker correctly finds a feasible relaxation for every pair.
2. **Nonuniform tiny time grid:** \(t=2^{-44}[0,1,2,8]\), exact healthy polynomial sensors with \(\tau=(0,2^{-44})\). The invalid first-step curvature bound would cause a false contradiction. The exact checker uses the maximum gap and does not reject the healthy pair.
3. **Deliberate time-constant switch:** an exponentially decaying signal changes response time from 0.5 to 2 seconds between two windows. Each window's LP is feasible individually, but the combined constant-\(\tau\) two-window LP produces a checkable exact-rational contradiction.
4. **Small rational algorithm audit:** 100 seeded random two-variable LPs with box bounds were compared with an independent SciPy/HiGHS feasibility oracle. Every exact-rational contradiction returned valid Farkas weights and agreed with HiGHS on these nondegenerate integer-coefficient tests.

## Matched-window benchmark (one independent seeded case per drift)

The exact checker rebuilt all sample integrals directly from the recorded raw float64 readings (interpreted as exact binary rationals), not from rounded \`WindowSummary.integrals\` produced by NumPy. One example seed (\`seed=0\`) yielded the following numbers of exact pairwise LP contradictions among 15 sensor pairs for the four-window experimental data:

| Drift (°C/s) | All healthy | Biases +1.2°C and −1.5°C |
| ---: | ---: | ---: |
| 0.0 | 0 | 9 |
| 0.3 | 0 | 9 |
| 0.6 | 0 | 9 |
| 1.0 | 0 | 7 |

For each of the exact contradictions, nonnegative rational Farkas multipliers can be extracted and checked. These results agree with the corresponding one-seed portion of the existing synthetic benchmark. They do **not** quantify real-world detection or prove unique faulty-sensor identities.

## Reproduction

From the GitHub repository root:

~~~bash
python -m pip install -r requirements-research.txt
python -m unittest discover -s tests -p 'test_exact_pairwise_certificate.py' -v
python -m unittest discover -s tests -v
~~~

The module is intentionally free of SciPy; SciPy is used only as an **independent test oracle** in one test.

## Important limits

- The theorem establishes a contradiction only for the **specified deterministic assumptions**, not a hardware failure.
- The checker decides a **relaxed LP**, not the full continuous-time ODE sample-feasibility problem.
- The proof treats supplied float64 values as exact binary rationals. If timestamps or sensor values were rounded during acquisition, a true physical soundness claim also requires bounds accounting for **measurement representation, timing, model, and calibration errors**.
- This does not provide a certified fast float LP solver. Use the exact verifier to **audit** particular flagged pairwise contradictions. Do not silently convert a HiGHS infeasibility result into a formal physical claim.
- The exact checker is **pairwise only**. A network-wide LP involving shared parameters among three or more sensors can detect additional contradictions that pairwise tests miss. General exact LP certificates and scalable verification are not implemented.
- Fourier–Motzkin elimination can have severe intermediate-expression and time growth. The 24-window cap is a practical guard, not a complexity guarantee.
- The method is a standard application of exact linear feasibility and Farkas' lemma, **not an original mathematical invention**.

**Next scientific task:** independently verify the exact-certificate module against additional LP oracles and published interval/set-membership implementations; investigate certified network-wide infeasibility and calibrated sensor data.
