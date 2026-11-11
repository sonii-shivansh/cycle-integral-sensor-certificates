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

**Licensing:** No reuse license is included in this initial publication. The repository owner can choose and add an appropriate license separately.
