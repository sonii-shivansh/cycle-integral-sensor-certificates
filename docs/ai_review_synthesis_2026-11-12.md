# AI-assisted scientific referee synthesis (preliminary)

**Recorded:** 12 November 2026. **Code snapshot reviewed in follow-up:** \`main\` at commit \`c09d5367715dd0dd7fc16dee3faf860536d87587\` (the last identified repository commit before creating the experimental research branch). New computational evidence is in branch \`research/multiwindow-lp-validation\`.

## Provenance and scope

The repository owner provided two model-generated review documents:

1. A **Codex-labelled referee report pasted in chat**, including mathematical proofs, a critique of graph fault identification, and a substantial prior-art list. The pasted text was incomplete (ended mid-sentence); the reviewer could access repository text but reported inability to pin a commit SHA.
2. A **Gemini-labelled 24-page scientific audit PDF**, which explicitly states it could **not fetch the repository source tree**. It independently formulated proofs and ran its **own synthetic experiments with different parameters**, rather than reproducing the exact repository Python implementation.

These are **AI-generated analyses, not human peer review**, and their reported execution counts, literature claims, and URLs should not be treated as independently verified merely because they appear in the reports. The uploaded PDF and chat text are not reproduced wholesale in this repository; this file preserves a careful summary for future human reviewers.

## Convergent findings

- The integrated common-input identity is correct under the stated constant first-order model and suitable regularity.
- The original symmetric and sign-aware inequalities are sound *necessary* conditions with deterministic noise, curvature and time-constant bounds.
- The sign-aware acceptance interval is contained within the original symmetric interval **in exact arithmetic**.
- Both inequalities are conservative outer bounds, not sufficient characterizations of physical sensor health.
- A matching or minimum vertex cover explains pairwise contradictions *conditional on* a shared-input/localizable-model interpretation. It is **not** an unconditional hardware fault count.
- A shared spatial temperature gradient or other global model failure can produce many conflict edges without any sensor being defective.
- Substantial related literature predates the project: integrated analog diagnosis, common-input sensor dynamics, set-membership parameter estimation, and conflict/hitting-set diagnosis.

## Specific findings needing additional work

1. **Endpoint-noise correlation:** The two reports observe that separately bounding endpoint changes and trapezoidal noise can be unnecessarily loose. A refined exact single-sensor noise factor for a *fixed* \(\tau\) is \((L-h+\max(2\tau,h))\epsilon\), before accounting for curvature. This does not make the current certificate incorrect.
2. **Floating-point rounding:** Gemini reports a one-ULP-scale discrepancy in interval containment for some synthetic parameter values. The existing Python predicates have a fixed \(10^{-10}\) comparison tolerance, so that example is not automatically a decision-level failure. A scale-aware, directed-rounding analysis remains necessary.
3. **Bias/lag aliasing:** Under linear input drift with slope \(c\neq0\), a constant bias can imitate a shifted valid time constant \(\tau-B/c\). A one-window method cannot universally distinguish these cases.
4. **Multi-window shared-parameter LP:** Gemini proposes preserving the same \(\tau_i\) over windows. Its constraints are a valid **linear outer relaxation**, not an exact full dynamic-model feasibility test. The new research branch independently implements and tests this direction using SciPy HiGHS; numerical infeasibility must still be scrutinized near tolerances.
5. **Novelty:** None of the AI audits establishes a new theoretical result. Relevant set-membership feasibility methods are classical; the next comparison must use a published implementation rather than only the original CIFC code.

## Current action and remaining questions

The accompanying branch contains a separate LP prototype, unit tests, a two-window example with changing time constant, and a matched-window 6-sensor benchmark. No earlier code is replaced and the Zenodo v0.1.0 archive is unchanged.

Before treating any of this as a new research contribution:
- Independently verify the multi-window proof with a control-theory expert.
- Find and compare published multi-window set-membership methods.
- Introduce certified numerical infeasibility checks and tests of error correlations.
- Test actual calibrated hardware under a genuinely shared temperature input.
- Preserve counterexamples and negative outcomes.

**Source caution:** This synthesis records claims from the user's provided AI reports; it does not attest that their citations, simulations or access assertions have all been independently checked.
