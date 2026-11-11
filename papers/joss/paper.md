---
title: "Cycle-Integral Sensor Certificates: An Exploratory Python Prototype for Bounded Sensor Consistency Tests"
tags:
  - Python
  - sensor networks
  - model-based diagnosis
  - control systems
authors:
  - name: Shivansh Soni
    affiliation: 1
affiliations:
  - name: "Affiliation to be verified by author before submission"
    index: 1
bibliography: paper.bib
---

**PREPARATION DRAFT — NOT SUBMITTED TO JOSS.** The current software does not yet meet the journal's public-development-history, feature-completeness, or demonstrated-impact screening requirements. This draft intentionally identifies missing evidence; it is not a claim of eligibility or acceptance.

# Summary

Cycle-Integral Sensor Certificates (CIFC) is a small Python research prototype for detecting inconsistencies among sensors responding to the same unknown, time-varying physical quantity. Its initial example concerns temperature sensors with different response speeds. When measurements disagree, a difference might indicate a sensor fault, or it might reflect the expected delay of a slower sensor. CIFC uses an assumed first-order response model to form tests that do not require knowing the actual temperature.

For pairs of sensors, the implementation approximates an integrated difference in observed outputs, compares it to a deterministic bound based on endpoint changes, measurement-noise limits, response-time limits, and numerical-integration error, and constructs a graph of contradictory sensor pairs. Standard graph algorithms provide a lower bound on the number of incompatible sensor-model assumptions and a smallest set of sensors that could explain the detected conflicts. A flagged contradiction means the stated assumptions cannot all hold; it does not automatically identify broken equipment.

# Statement of need

Researchers evaluating redundancy in sensor networks need ways to reason about disagreement when sensor dynamics differ and the true input is unknown. A straightforward comparison of instantaneous readings may produce misleading alarms: two sensors can report different values while obeying the same physical input because their response times differ. A model-consistency formulation provides an alternative framework for asking whether observed differences are *compatible* with assumed dynamics and bounded measurement error.

CIFC's current purpose is educational and investigative: to provide a compact executable demonstration of this formulation, not to replace calibrated industrial fault diagnosis. The intended users are researchers and students working in system identification, model-based diagnosis, or sensor reliability. The package exposes a synthetic signal generator, pairwise consistency checks, an exhaustive maximum-matching demonstration, and a minimum-vertex-cover search for small graphs. These functions offer a transparent point of departure for controlled benchmark studies.

# State of the field

The foundational ideas are established. Sensor-dynamics compensation and unknown-input thermal reconstruction have long been studied, including two-thermocouple methods [@tagawa1997]. Integrated differential constraints and bounded observation representations were developed for analog diagnosis before this project [@loiez1997]. Conflict sets and minimal multiple-fault explanations are also classical in model-based diagnosis [@dekleer1987]. CIFC makes **no claim** to have invented these concepts. A rigorous full-text comparison against modern interval methods, parity relations, and available software packages is still required.

An important open question for this project is whether its particular combination offers any advantage in a meaningful regime. In the initial periodic synthetic benchmark, a simpler full-period average comparison finds the same nine contradictory pairs and can achieve a tighter sufficient detection bound. Therefore current results do **not** demonstrate an improvement over a basic special-case comparator, much less over state-of-the-art fault-diagnosis packages. This deficiency must be addressed before any future journal submission.

# Software design

The prototype is written in Python using NumPy. Each healthy sensor is assumed to satisfy a first-order differential equation with a shared input and an unknown but bounded, nonnegative time constant. Integrating the difference of two such equations removes the shared input. The implementation then computes a composite-trapezoid approximation from synchronized samples.

The inequality uses conservative deterministic assumptions: a maximum sensor time constant, a hard observation-error bound, and a curvature bound controlling quadrature error. A strict violation is one-sided evidence of inconsistency. Passing the test neither establishes health nor proves that a fault is absent. An identical offset added to all sensor readings cannot be distinguished from an equivalent change in shared temperature under this model.

Pairwise contradictions form edges in a graph whose vertices represent sensors. A maximum matching gives a standard lower bound on the number of incompatible vertices; a minimum vertex cover gives a smallest explanation for the conflicts. The current graph algorithms deliberately use exhaustive search to keep a six-sensor illustration understandable. They are **not scalable** implementations suitable for large research networks. Future maintainable software would need a proper library interface, more efficient combinatorial algorithms, input validation, extensible models, versioned APIs, and independent benchmark tests.

# Research impact statement

At the time of this draft, CIFC has one six-sensor, synthetic periodic example and deterministic regression checks. The constructed case includes deliberately offset sensors 3 and 6; nine of fifteen pairs are flagged and the unique smallest graph explanation selects those two indices. These tests support reproducibility of the example, **not** real-world accuracy or originality. No independent users, external deployments, research papers enabled by the library, or hardware benchmarks have been documented. Accordingly, the current repository has **not demonstrated** the research impact or sustained open development that JOSS requires. This section must be updated with verified evidence, not aspirational claims.

# AI usage disclosure

OpenAI ChatGPT (GPT-6, this drafting session) assisted with documentation, manuscript drafting, regression-test scaffolding, and presentation of mathematical arguments. Earlier project artifacts also state ChatGPT assistance; their exact model versions should be verified by the author before submission. Human independent review, editing, and validation of all AI-assisted scientific claims are **not yet established**. This is a **draft disclosure**, not the mandatory pre-submission human-verification attestation. A future author must accurately identify all AI tools and versions, verify every result, make and document substantive design decisions, and personally handle any editorial/reviewer communication as required by JOSS policy.

# Acknowledgements

AI assistance is disclosed above. Funding sources and any additional human contributors must be confirmed by the submitting author before publication.

# References
