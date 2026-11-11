# Prior-Art and Novelty Audit — Cycle-Integral Fault Certificate (CIFC)

**Audit date:** 11 November 2026  
**Research candidate:** `CIFC_research_note.md` and `cycle_integral_fault_certificate.py`  
**Disposition:** **Do not claim a novel foundational discovery, world-first algorithm, or patentability based on the current record.** The exact formula and assembled implementation have not been conclusively matched to a single pre-existing publication, but central techniques are clearly established.

## Scope and limits

This was a broad, targeted public-web search of scientific publications, publishers, conference proceedings, academic repositories, and Google Patents, including equivalent mathematical descriptions rather than only the label “CIFC”. It is **not** a complete search of all websites, subscription databases, unindexed materials, private reports, unpublished filings, or every national patent database. Absence of an exact-title or exact-formula result cannot prove worldwide novelty. Patent search is not legal clearance or an infringement opinion.

## What CIFC claims

Given first-order lag measurements of the same unknown input, with `tau_i y_i' + y_i = T` and bounded `tau_i`, measurement errors, and curvature, CIFC removes the unknown input by subtracting two differential equations and integrating. It then derives a conservative deterministic threshold for the observed trapezoidal integral of their difference. An exceeding pair becomes an edge in a conflict graph. A matching proves a lower bound on the number of incompatible sensors; a minimum vertex cover identifies a smallest explanation. A periodic-signal case gives a sufficient detectable relative-bias size.

## Source-based prior-art comparison

| CIFC ingredient | Documented earlier work | Assessment |
|---|---|---|
| Shared first-order thermocouple dynamics and unknown common temperature | Tagawa & Ohta, 1997, *Two-thermocouple probe for fluctuating temperature measurement*, https://doi.org/10.1016/S0010-2180(97)00044-8 ; Li et al., 2018, *Blind system identification of two-thermocouple sensor based on cross-relation method*, https://doi.org/10.1063/1.5019965 | **Known underlying model and common-input cross-relation concept**; do not represent elimination of the shared temperature as a new physical principle. |
| Integrate differential model constraints to avoid numerical differentiation and test bounded observation consistency with no false conflicts | Loiez & Taillibert, 1997, *Polynomial Temporal Band Sequences for Analog Diagnosis*, IJCAI, https://ijcai.org/proceedings/1997/70 | **Strong conceptual match.** Their implementation uses temporal polynomial bands instead of CIFC’s trapezoidal curvature-error envelope. The exact expressions are not established as identical. |
| Fault detection using analytical redundancy, bounded uncertainty, and interval-consistency thresholds | *Model-based approach for fault diagnosis using set-membership formulation* (published in *Engineering Applications of Artificial Intelligence*), https://www.sciencedirect.com/science/article/abs/pii/S0952197616301324 ; 2017 interval redundancy study https://www.sciencedirect.com/science/article/abs/pii/S0888327017300869 | **Known methodology.** CIFC’s explicit closed-form inequality may still be a distinct implementation. |
| Pairwise sensor conflicts and compute minimum hitting sets (equivalent to vertex cover for two-element conflicts) | de Kleer & Williams, 1987, *Diagnosing multiple faults*, https://doi.org/10.1016/0004-3702(87)90063-4 ; Jonas Biteus, 2005, *Distributed Diagnosis and Simulation Based Residual Generators*, https://www.vehicular.isy.liu.se/Publications/Lic/05_LIC_1176_JB.pdf | **Direct pre-existing logic.** The 2005 thesis specifically uses two sensors measuring the same temperature as an example of pairwise conflict and describes diagnosis as a minimum hitting-set problem. |
| More recent sensor-network inconsistency certificates and graph-based fault localization | Baskar, George & Roberts, 2026, *Spectral-Bounded Consistency Certification for Heterogeneous Wireless Sensor Networks*, https://doi.org/10.1109/JIOT.2026.3685152 | **Related ongoing field**, not verified identical to CIFC. |
| Prior patented redundant-sensor fault detection based on uncertainties | Analog Devices, US12399044B2, *Diagnostic coverage improvements in redundant sensor element fault detection architectures*, https://patents.google.com/patent/US12399044B2/en ; US7841771B2, *Self-validating thermocouple*, https://patents.google.com/patent/US7841771B2/en | **Adjacent patent landscape**, not evidence of an identical claimed invention. A patent attorney would need to read claims and relevant family members to assess patentability or freedom to operate. |

## Independent numerical baseline check

The archived CIFC code was executed successfully. In its six-sensor, one-period sinusoidal demonstration with fixed biases on sensors 3 and 6:

- CIFC found **9/15** contradictory pairs; a minimum conflict cover of `{3, 6}` and a maximum matching of size 2.
- A **much simpler** comparator flagged a pair when the absolute trapezoidal integrated difference exceeded `L*(epsilon_i+epsilon_j)`. This comparator is justified for this particular *exactly periodic sinusoidal, steady-state, full-cycle* setup because healthy signals have zero integrated relative difference; sampled sinusoidal trapezoidal quadrature is exact over the full cycle.
- On the **same generated signals**, this simple comparator also found **the same 9/15 edges**.
- Under bounded adversarial sample noise with two epsilons of `0.02°C`, for these sinusoidal traces the simple period-average comparator guarantees detection of fixed relative biases strictly greater than **0.08°C**, while the CIFC note's general inequality produces a bound of about **0.20066°C**. Thus the specific synthetic showcase does **not** establish CIFC’s superiority or optimal detection sensitivity.
- The simple period-average rule does not work for arbitrary nonperiodic intervals. CIFC’s general test has broader assumptions than that simple baseline. A fair benchmark should test nonperiodic inputs as well as several stronger published analytical-redundancy and set-membership methods.

## Mathematically sound does not mean mathematically new

The CIFC bound is a valid sufficient contradiction test under the paper's stated first-order model, synchronized bounded sample noise, upper-bounded lag time, and curvature assumptions. The conflict-graph matching result is an immediate graph-theoretic consequence of the established hitting-set diagnosis framework. The guaranteed bias threshold is algebra following from the same elementary error bounds. A mathematically correct new *presentation* or exact numerical coefficient does not itself establish a substantial novel result.

The method checks *model inconsistency*, not guaranteed physical component failure. Incorrect shared-temperature assumptions, time-varying lag, neglected radiation or conductivity, and clock skew can create apparent conflicts despite healthy hardware. The all-sensor common-offset unobservability is correct for the stated model, but should not be marketed as a novel discovery.

## Decision

1. **Foundational physics discovery? No evidence.** The physics and common-input cancellation belong to established first-order sensor theory.
2. **Brand-new mathematical algorithm? Not established.** Integrated residuals with bounded consistency checks and graph-based diagnosis have direct precedents, dating to at least the 1980s and 1990s.
3. **Exact CIFC formula already published word for word? Not identified in this targeted search; absence is not proof of novelty.**
4. **Potential practical engineering contribution? Possibly**, if a new theorem, reduced computational cost, calibration-free field robustness, or reliable advantage over published methods is demonstrated with rigorous baselines and real sensor data.
5. **Ready to claim patentability or publish as a world-first result? No.**

**Next research steps:** full-text mathematical comparison to Loiez & Taillibert (1997), thermal cross-relation papers, interval parity relations, and de Kleer–Williams; independent theorem review and sensitivity analysis; hardware experiments; benchmarking on nonperiodic transients; comprehensive patent search using CPC/IPC classes and claim-by-claim legal review if commercialization is planned.