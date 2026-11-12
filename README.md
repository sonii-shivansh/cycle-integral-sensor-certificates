# Cycle-Integral Sensor Certificates (CIFC)

**Research prototype · mathematical consistency checks · synthetic validation · open for scientific review**

> **Status (11 November 2026):** Preliminary technical research note and reproducible Python simulation. **Not peer reviewed, not hardware validated, and not established as a novel or world-first algorithm.** Closely related approaches to integrated consistency checks, thermocouple response models, and conflict-based diagnosis predate this work.

## Overview

CIFC explores a deterministic consistency test for several sensors measuring a common, unknown, time-varying temperature. Each sensor may have an *unknown but bounded* first-order response time. Instead of reconstructing the true temperature, the method eliminates that shared input algebraically.

Under explicit physical and measurement-error assumptions, a violation of the pairwise bound **certifies that at least one assumption for one of the two sensors is false**. It does **not** by itself establish hardware failure.

The prototype combines:
- An integrated cross-sensor relation for first-order systems;
- A conservative threshold accounting for bounded sample noise, time constants, and quadrature error;
- A contradiction graph and a matching-based lower bound on incompatible sensors;
- An exact minimum-vertex-cover calculation for a small demonstration network.

These ingredients are related to established prior research; **no originality or patentability claim is made**.

## Mathematical idea

For sensor \(i\), assume

\[
\tau_i y_i'(t)+y_i(t)=T(t),\qquad 0\le\tau_i\le M.
\]

Two sensors observing the same \(T(t)\) obey

\[
\int_a^b (y_j-y_i)\,dt
=\tau_i[y_i(b)-y_i(a)]-\tau_j[y_j(b)-y_j(a)].
\]

Let \(Q_{ij}\) be the composite trapezoid estimate of the *observed* sampled difference, \(d_i=z_i(b)-z_i(a)\), \(L=b-a\), \(h\) the sampling step, \(\epsilon_i\) a hard bound on sample error, and \(K_i\) a hard bound on \(|y_i''|\). If both sensors satisfy the stated model, then

\[
|Q_{ij}| \le M(|d_i|+|d_j|)
+(\epsilon_i+\epsilon_j)(L+2M)
+\frac{Lh^2(K_i+K_j)}{12}.
\]

A strict violation is a **model-inconsistency certificate**. The graph matching bound is a standard consequence of pairwise conflicts; see the [research note](CIFC_research_note.md) for assumptions and proofs.

## Reproduce the simulation

Requires **Python 3.10+** and **NumPy 2.x**.

~~~bash
python -m pip install -r requirements.txt
python cycle_integral_fault_certificate.py
~~~

The constructed demonstration uses six sensors and known injected biases on sensors 3 and 6. It reports:

- **9 of 15** pairwise contradictions;
- A maximum matching of size **2**, implying at least two sensors incompatible with the assumed model;
- A **unique minimum vertex-cover hypothesis** consisting of sensors **3 and 6**;
- No contradictions in **100** additional healthy-array synthetic noise trials.

These results are limited to the specified simulation; they do not demonstrate field accuracy or superiority over existing methods.

## Repository guide

| File | Content |
| --- | --- |
| [cycle_integral_fault_certificate.py](cycle_integral_fault_certificate.py) | Executable NumPy proof of concept |
| [CIFC_research_note.md](CIFC_research_note.md) | Model, integral identity, formal proofs, and synthetic experiment |
| [CIFC_novelty_audit_2026-11-11.md](CIFC_novelty_audit_2026-11-11.md) | Prior-art investigation and a simpler baseline comparison |
| [REFERENCES.md](REFERENCES.md) | Core related academic work and links |
| [LIMITATIONS.md](LIMITATIONS.md) | What the method does and does not certify |
| [CONTRIBUTING.md](CONTRIBUTING.md) | How to propose corrections, replication studies, and independent review |

## Citation and release archiving

A machine-readable [CITATION.cff](CITATION.cff) is included so GitHub can display **Cite this repository** and Zenodo can read the metadata for an archived software release. The citation is prepared for version **0.1.0**. **No DOI has been assigned or verified yet**; do not cite a placeholder DOI.

The creator name matches the public copyright notice in [LICENSE](LICENSE). When a Zenodo DOI is available, use the **version-specific DOI** to cite a particular archived release.

For the first official GitHub release and Zenodo archive, follow [PUBLISHING.md](PUBLISHING.md). Enable GitHub integration in Zenodo **before** publishing the release.

## Regression checks

In addition to the example above, run:

~~~bash
python -m unittest discover -s tests -v
~~~

GitHub Actions runs this test suite and the example script on pushes and pull requests. These are software regression checks, not independent experimental validation.

## Follow-up mathematical check: sign-aware interval certificate

A follow-up [derivation and controlled stress test](docs/directional_extension.md) keeps the signs of the sensor endpoint changes instead of immediately taking their absolute values. Under the **same** first-order model and deterministic bounds, the sign-aware feasible interval is contained in the original one, so it can only add contradictions, not remove them. This is an elementary **interval refinement**, **not a verified new research discovery**.

For a 13.5-second partial-period/nonperiodic-window synthetic benchmark, two fixed sensor biases were added while temperature drift varied. Across **100 seeded trials per drift setting**:

| Temperature drift (°C/s) | Original test: trials with any contradiction | Sign-aware test: trials with any contradiction |
| --- | ---: | ---: |
| 0.0 | 100/100 | 100/100 |
| 0.3 | 0/100 | 100/100 |
| 0.6 | 0/100 | 0/100 |
| 1.0 | 0/100 | 0/100 |

On corresponding all-healthy synthetic arrays, **neither method raised a contradiction**. Another **250 randomized healthy ensembles** also passed. These counts are specific to exactly described simulated physics, chosen bounds and injected biases—not hardware false-alarm or detection rates. At larger drift **both tests miss the injected biases**; there is no demonstrated advantage over published interval or set-membership approaches.

Reproduce locally:

~~~bash
python -m experiments.nonperiodic_benchmark
python -m unittest discover -s tests -v
~~~

The original implementation is preserved in \`certificates()\`, and the new one is \`directional_certificates()\`. See [nonperiodic experiment source](experiments/nonperiodic_benchmark.py) and [formal extension with limitations](docs/directional_extension.md). Existing Zenodo release \`v0.1.0\` is a fixed historical snapshot and **does not contain this follow-up work**.

## Experimental multi-window LP (proposed, not yet merged)

A separate experimental implementation is under review in [experiments/multiwindow_lp.py](experiments/multiwindow_lp.py). It uses one unknown but **constant** response-time parameter per sensor across several observation windows, with an LP outer relaxation of the same first-order sensor model. It **does not certify physical hardware failures**, and a numerical solver result is not a formal proof of infeasibility near roundoff tolerances.

- [Derivation, caveats, and matched-window results](docs/multiwindow_validation_2026-11-12.md)
- [Source-comparison case study](experiments/multiwindow_case_study.py)
- [Attribution and limitations of the Gemini/Codex AI-generated referee reports](docs/ai_review_synthesis_2026-11-12.md)

Run the research-only comparison (requires SciPy):

~~~bash
python -m pip install -r requirements-research.txt
python -m unittest discover -s tests -v
python -m experiments.multiwindow_case_study
~~~

**This extension is experimental; the Zenodo v0.1.0 archive has not been changed.** No novelty, peer-review, hardware-accuracy or superiority-to-published-methods claim is made.

## Exact-rational contradiction verification (experimental branch)

The multi-window research branch now includes an **independent two-sensor LP verifier** that recomputes trapezoid integrals from recorded samples as exact rational numbers, uses Fourier–Motzkin elimination, and produces a checkable **Farkas infeasibility certificate** when the relaxed pairwise constraints conflict. It is separate from the fast, floating-point SciPy/HiGHS solver.

- [Proof, scope and reproducibility guide](docs/exact_lp_certification_2026-11-12.md)
- [Pure-Python rational certificate checker](experiments/exact_pairwise_certificate.py)
- [Adversarial regression tests](tests/test_exact_pairwise_certificate.py)
- [Independent scientific integrity audit](docs/verification_audit_2026-11-12.md)

This verifies the **rational LP built from the exact represented numeric inputs**, not sensor calibration, scientific novelty, or full differential-equation feasibility. It does **not** justify converting any numerical solver's rejection into a certified physical hardware failure. No source files from Zenodo \`v0.1.0\` have been altered.

## Independent set-membership prior-art comparison (negative result)

We independently reconstructed a classical **invariant-parameter feasible-solution-set (FSS)** baseline using exact-rational convex-polygon clipping, following the conceptual framework discussed in Ingimundarson et al. (2005). This is **not a replication of their parallelotope approximation code**. For the **same affine constraints and matched sensor data**, this independent baseline and CIFC's multi-window LP found **exactly the same sensor-pair contradictions** across all 20 seeds in every tested condition. This is consistent with their proven mathematical equivalence, and does **not** demonstrate any advantage for CIFC over the standard set-membership formulation.

- [Prior-art comparison, full table and mathematical equivalence](docs/set_membership_prior_art_comparison_2026-11-12.md)
- [Independent rational geometric baseline](experiments/set_membership_fss_baseline.py)
- [Fair matched-window comparison script](experiments/set_membership_comparison.py)
- [Adversarial baseline agreement tests](tests/test_set_membership_fss.py)

Run:

~~~bash
python -m pip install -r requirements-research.txt
python -m experiments.set_membership_comparison
~~~

**Research integrity:** All figures are synthetic contradiction-edge counts, not sensor identification rates, and the project is neither novel nor validated merely because different implementations agree. The fixed Zenodo v0.1.0 remains unchanged.

## Exact network-only inconsistency and Helly bound (research draft)

The complete network LP can be infeasible **even though all pairwise multi-window LPs are feasible**. A [three-sensor, two-window example](docs/network_wide_helly_consistency_2026-11-12.md) constructs exact rational observations where every pair separately admits compatible response-time parameters but the network jointly yields a checkable **Farkas contradiction \`0 <= -1\`**.

Each sensor's relaxed feasible shared-input integral set is convex in \(\mathbb R^W\). By **classical Helly's theorem**, any infeasible network relaxation has an inconsistent subset of at most \(\min(n,W+1)\) sensors. With **one window**, pairwise checks suffice for the relaxed network LP; with **two windows**, three sensors may be necessary. This is an **application of established convex geometry**, not a new theorem or evidence of broken hardware.

- [Exact example and theorem proof](docs/network_wide_helly_consistency_2026-11-12.md)
- [Three-sensor rational proof script](experiments/network_consistency_gap.py)
- [Network LP and inclusion-minimal candidate conflict extraction](experiments/network_conflicts.py)
- [Reproducible network scaling benchmark](experiments/network_scaling_benchmark.py)

A [GitHub Actions run](https://github.com/sonii-shivansh/cycle-integral-sensor-certificates/actions/runs/43223985881) passed 34 tests and printed numerical timings for up to 48 sensors. **One network-consistency decision and a full pairwise contradiction graph are different outputs**; no apples-to-apples speed advantage has been established. The fixed Zenodo \`v0.1.0\` record is unchanged.

## Certified multi-sensor network contradictions (research-only)

The research branch now contains a **general multi-sensor Farkas proof extractor**. A floating-point dual solver proposes a candidate contradiction, but no result is called certified unless it is reconstructed and checked with **exact rational arithmetic**. The output is strictly \`certified_infeasible\` or \`not_certified\`; **not certified does not mean healthy or feasible**.

- [Formal explanation and tests](docs/network_exact_farkas_2026-11-12.md)
- [Network-wide rational certificate verifier](experiments/network_exact_certificate.py)
- [Adversarial and randomized tests](tests/test_network_exact_certificate.py)

The new checker verifies *finite relaxed LP infeasibility* for the represented data, **not physical sensor failure**. Exact certificates are established linear-programming mathematics, not a scientific novelty claim. A [GitHub Actions run](https://github.com/sonii-shivansh/cycle-integral-sensor-certificates/actions/runs/43226515287) passed 41 tests including these checks. It does not change the archived Zenodo \`v0.1.0\`.

## Important findings and constraints

- The certificate is **one-sided**: an alarm implies model inconsistency under the stated bounds, but no alarm cannot establish health.
- Sensors must share the same temperature input and synchronized sampling; thermal gradients and clock errors can create apparent contradictions.
- The bound relies on trustworthy maximum response times, sample-noise bounds, and curvature bounds.
- A common offset affecting all sensors is not identifiable using only this model.
- On the initial *periodic* synthetic benchmark, a simpler period-average test detects the **same nine pairs**, with a tighter guaranteed bias threshold for that special setting. The current demonstration therefore shows **no performance advantage** for CIFC over that baseline.
- A unique minimum vertex cover is a *smallest explanation*, **not a proven identification of the failed devices**.

See [LIMITATIONS.md](LIMITATIONS.md) before interpreting any result.

## Scientific prior art and attribution

This work is **not** presented as a new physical law. In particular, Loiez & Taillibert (1997) used integrated differential constraints and bounded observation representations for diagnosis; de Kleer & Williams (1987) developed conflict-based multiple-fault diagnosis; and existing thermocouple literature studies first-order response dynamics. See [REFERENCES.md](REFERENCES.md) and the [novelty audit](CIFC_novelty_audit_2026-11-11.md).

## Academic publication preparation (not submitted)

The Zenodo software release is an archive, **not** a peer-reviewed research paper. We are preparing separate, preliminary manuscript drafts:

- **[arXiv draft and preparation notes](papers/arxiv/README.md)** — [LaTeX manuscript source](papers/arxiv/main.tex); not submitted and requires independent scientific review.
- **[JOSS preparation notes](papers/joss/README.md)** — [JOSS-formatted manuscript draft](papers/joss/paper.md) and [bibliography](papers/joss/paper.bib); **not eligible yet** under the journal's current sustained public-development, completeness, and impact requirements.

Both drafts disclose AI assistance and substantial prior work. Neither represents an arXiv endorsement, JOSS submission, journal acceptance, or a claim of scientific novelty. Check each venue's latest policies before submission.

## Open invitation for scientific review

**Researchers and engineers are invited to challenge, reproduce, or improve the work.** Especially valuable are:
1. Counterexamples to any theorem or hidden violations of assumptions;
2. Corrections to the numerical bounds, derivation, or implementation;
3. Earlier publications or patents that subsume the exact method;
4. Reproducible comparisons against interval-consistency, parity-relation, and model-based diagnosis baselines;
5. Real-world sensor data and experiments involving nonperiodic inputs.

Please [open an issue](https://github.com/sonii-shivansh/cycle-integral-sensor-certificates/issues) or submit a pull request, citing the theorem, file, or line concerned. Disagreement, replication failures, and negative results are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Research and authorship transparency

This initial research note, derivation, accompanying code, and documentation were developed with assistance from OpenAI's ChatGPT and have **not** undergone independent expert review. The repository owner is publishing them for scrutiny and correction. No new discovery, peer-review status, experimental validation, or specific author credentials are implied.

**License:** This repository includes the [MIT License](LICENSE), with copyright attributed to Shivansh Soni (2026). See the license file for its terms.
