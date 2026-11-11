# Contributing and scientific review

We welcome **independent, critical review** of this preliminary research prototype. Corrections and negative findings are as important as confirmations.

## Request for review

Please open a GitHub issue (or a pull request with a reproducible correction) in one of these areas:

1. **Mathematical rigor:** check the shared-input cancellation, composite trapezoid error, endpoint noise terms, soundness claim, matching argument, and periodic bias threshold.
2. **Prior art:** identify an earlier equivalent integral bound, interval parity check, fault-certificate graph, or patent disclosure. Provide DOI/URL and the exact analogous result.
3. **Counterexamples:** provide a mathematical counterexample *within the published assumptions*, or clearly document which physical assumption fails.
4. **Baseline benchmarking:** compare with established model-based residuals, set-membership diagnosis, and a periodic average when appropriate. State identical noise regimes and datasets.
5. **Physical validation:** propose experiments using documented, calibrated temperature sensors and synchronized logging. Never treat this prototype as a safety controller.
6. **Software correctness:** tests, portability, documentation, scientific reproducibility, and clearer diagnostics.

## Useful issue format

- **Summary / claim being reviewed:**
- **Relevant theorem, equation, source file, or lines:**
- **Assumptions / environment / versions:**
- **Evidence or minimal reproduction:**
- **Expected result vs. observed result:**
- **Earlier references, if applicable:**

## Reproduction

~~~bash
python -m pip install -r requirements.txt
python cycle_integral_fault_certificate.py
~~~

Please include interpreter and NumPy versions, complete command lines, random seeds, and raw numerical outputs with submitted experiments. A pull request should explain the changed assumption or guarantee and add a reproducible check where relevant.

## Tone and scientific integrity

Respectful, evidence-based criticism is explicitly encouraged. Do not claim originality, hardware accuracy, statistical significance, or peer review without supporting evidence. When a claim is disproved, document the correction rather than hiding negative results.

## Authorship

The initial content was developed with AI assistance. Contributors should identify their own contributions accurately; collaboration on code or review does not automatically confer academic authorship. Bibliographic and acknowledgment decisions should follow relevant journal and institutional norms.
