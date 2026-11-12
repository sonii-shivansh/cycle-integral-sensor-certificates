# Network-wide exact-rational Farkas certificates

**Research status: 12 November 2026.** A *verification-engineering milestone*, not a new mathematical theorem, a physical fault diagnosis, or independent human peer review.

## Why extend the pairwise exact verifier?

The earlier [two-sensor rational checker](../experiments/exact_pairwise_certificate.py) verifies selected pairwise contradictions. However, a previously constructed [three-sensor, two-window example](network_wide_helly_consistency_2026-11-12.md) has **all three pairwise LPs feasible** while the joint three-sensor relaxation is infeasible. A pairwise verifier cannot certify that particular full-network contradiction.

This branch adds a separate [general network proof extractor](../experiments/network_exact_certificate.py). It reconstructs the linear model exactly from rational recorded-sample summaries and validates **Farkas multipliers** for a chosen sensor subset or the whole sensor network.

## Exact finite LP being checked

Healthy sensor \(i\) is assumed to follow \(y_i+\tau_i y_i'=T(t)\) with one *constant* \(\tau_i\in[0,M]\) across all windows, bounded sample errors \(\epsilon_i\), curvature bounds \(K_i\), and a genuinely shared input \(T\). For recorded window \(w\), let \(I_{i,w}\) be the **exact composite-trapezoid integral** of the represented sample values after a *common* reference shift, \(d_{i,w}\) its endpoint difference, \(L_w\) its length, and \(h_{\max,w}\) its largest sampling gap. The numerical uncertainty bound is

\[
D_{i,w}=L_w\epsilon_i + \frac{L_w h_{\max,w}^2K_i}{12}.
\]

The finite LP uses unknown constants \(\tau_i\) and integrated shared inputs \(J_w\):

\[
I_{i,w} + \tau_i(d_{i,w}-2\epsilon_i)-D_{i,w}\le J_w\le
I_{i,w}+\tau_i(d_{i,w}+2\epsilon_i)+D_{i,w},
\qquad 0\le\tau_i\le M.
\]

Write all such inequalities as \(Ax\le b\), with \(n+W\) freely represented variables \((\tau_1,\ldots,\tau_n,J_1,\ldots,J_W)\) and explicit box constraints on \(\tau_i\).

**Soundness proposition (conditional).** If any chosen subset of sensors truly follows the stated shared-input model, with all specified deterministic error/time/curvature limits satisfied, then the corresponding *exact* LP relaxation is feasible.

**Proof.** Choose the actual constant response times \(\tau_i^*\) and exact shared-input integrals \(J_w^*=\int T\), with the same reference centering. Integration of the first-order differential equation gives \(J_w^*=\int y_i+\tau_i^*\Delta y_i\). Bounded observation and quadrature errors place the exact integral within \(D_{i,w}\) of the recorded trapezoid, while endpoint differences differ by at most \(2\epsilon_i\). Thus \((\tau_i^*,J_w^*)\) satisfies all displayed inequalities. Therefore exact LP infeasibility excludes the joint modeling assumptions. QED.

This theorem does **not** establish a physical fault, prove that a feasible LP corresponds to a valid continuous-time response, or establish mathematical novelty.

## How the witness is constructed

For any rational LP \(Ax\le b\) with unconstrained \(x\) (all variable bounds are explicitly included among the inequalities), a Farkas witness is a rational vector \(\lambda\ge0\) with

\[
\boxed{\lambda^\mathsf TA=0,\qquad \lambda^\mathsf Tb<0.}
\]

Multiplying by \(\lambda\) gives the impossible statement \(0\le\lambda^\mathsf Tb<0\). This is **classical Farkas' lemma** (see, e.g. [CGAL's exact LP certificate documentation](https://doc.cgal.org/latest/QP_solver/index.html)).

The new algorithm separates **candidate search** from **proof checking**:

1. A floating-point SciPy/HiGHS dual problem searches for a nonnegative multiplier vector satisfying \(A^\mathsf T\lambda=0\) and the normalization \(\sum_k\lambda_k=1\), minimizing \(b^\mathsf T\lambda\). This is only a *heuristic support selector*.
2. For candidate nonzero multiplier indices, the checker re-solves the equality system with Python \`fractions.Fraction\`, **not floating-point rounding or rationalizing each floating multiplier independently**. Candidate free variables are set to zero.
3. The resulting \(\lambda\) is independently checked with exact arithmetic for every equality, nonnegativity, and a *strictly negative* exact value of \(b^\mathsf T\lambda\).
4. A witness is released **only** if the exact verification passes. Otherwise the result is **\`not_certified\`**, which explicitly means neither feasible nor proved infeasible.

The result is a *sound but incomplete certificate-search procedure*. It cannot promise a proof for every infeasible rational LP, nor runtime scalability. Failure to find a witness is **not evidence of sensor health**.

## Reproduction and verified evidence

~~~bash
python -m pip install -r requirements-research.txt
python -m unittest discover -s tests -p 'test_network_exact_certificate.py' -v
python -m unittest discover -s tests -v
~~~

- [Core implementation](../experiments/network_exact_certificate.py)
- [Independent regression tests](../tests/test_network_exact_certificate.py)
- [Original three-sensor exact contradiction](../experiments/network_consistency_gap.py)

[GitHub Actions run #43226515287](https://github.com/sonii-shivansh/cycle-integral-sensor-certificates/actions/runs/43226515287) **passed 41 tests**, including the seven new network-exact tests.

| Check | Finding | What it means |
| --- | --- | --- |
| Three-sensor/two-window example | Network LP **certified infeasible**, while exact pairwise checks remain feasible | A higher-order conflict can be formally proved |
| Huge exactly represented common offset | **No certificate** | Did not falsely reject a healthy exact model; not a proof of feasibility |
| Sensor response-time change between windows | **Certified infeasible** | Incompatible with the constant-time-constant assumption |
| Six-sensor synthetic biased observations, one seed for each of four drifts | **Exact network contradictions** at all four drift settings | Reproducible conditional model contradictions, not hardware fault identification |
| 35 random small integer-valued network LP inputs | All returned exact witnesses verified; contradictory cases agreed with a separate SciPy primal feasibility oracle | Additional implementation cross-check |
| Very small genuine contradiction (approximately \(10^{-12}\)) | **Not certified**, despite exact pairwise infeasibility | Demonstrates deliberate incompleteness of the floating-supported candidate search |
| Input validation and tampered multipliers | Rejected or failed exact verification | Distinguishes a valid Farkas certificate from an arbitrary solver output |

## Limitations and future experiments

1. **Exact numerical input meaning:** a float64 value is interpreted as its exact *stored binary fraction*, not as an ideal physical quantity. If a sensor sample or timestamp was rounded, the stated measurement-error/time uncertainties must cover the representation effects.
2. **Numerical candidate limitations:** candidate support, truncation thresholds, conditioning, and missing support variables can prevent finding a witness; the algorithm then returns \`not_certified\`. It does not silently claim a proof.
3. **Computational scaling:** rational Gaussian elimination can create huge numerators and denominators; there are no verified large-network runtime or memory bounds here.
4. **LP relaxation vs physics:** exact LP feasibility is only a necessary condition for the underlying differential model, not sufficient.
5. **Conflict interpretation:** a contradiction could arise from shared environmental/model assumptions rather than physical sensor defects.
6. **Prior art:** rational Farkas certificates, LP duality, set-membership feasibility, and global conflict diagnosis are established mathematics. This project currently demonstrates *transparent reproducibility*, **not a novel discovery**.
7. **Incomplete comparison:** the established parallelotope implementation of Ingimundarson et al. has not been reconstructed or fairly benchmarked; calibrated physical-sensor data are absent.

**Recommended next scientific step:** use verified Farkas witnesses to extract and rigorously justify conflict *hyperedges*, test solver robustness at different scales, and obtain independent human scientific review before promoting the software to a validated diagnostic tool.

**Publication decision:** keep PR #4 in draft; do not alter \`main\`, the existing Zenodo archive, or submit to arXiv/JOSS merely because the new software tests pass.
