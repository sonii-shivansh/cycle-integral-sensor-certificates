# References and related work

**Scope:** The references below identify relevant existing methods. Their inclusion does **not** establish that CIFC is novel, patented, peer reviewed, or distinct from every earlier formulation. See the [novelty audit](CIFC_novelty_audit_2026-11-11.md) for an ingredient-by-ingredient comparison.

## Prior scientific work

1. **de Kleer, J., & Williams, B. C. (1987).** *Diagnosing multiple faults.* Artificial Intelligence, 32(1), 97–130. https://doi.org/10.1016/0004-3702(87)90063-4  
   Established conflict-based reasoning and multiple-fault diagnosis, relevant to minimum hitting sets and the graph explanation step.

2. **Loiez, E., & Taillibert, P. (1997).** *Polynomial Temporal Band Sequences for Analog Diagnosis.* Proceedings of IJCAI-97. https://www.ijcai.org/Proceedings/97-1/Papers/070.pdf  
   Closely related use of integrated differential constraints and bounded observations for analog-model diagnosis.

3. **Tagawa, M., & Ohta, Y. (1997).** *Two-thermocouple probe for fluctuating temperature measurement.* Combustion and Flame. https://doi.org/10.1016/S0010-2180(97)00044-8  
   Established thermocouple dynamics, coupled sensor response, and shared temperature measurement.

4. **Li and coauthors (2018).** *Blind system identification of two-thermocouple sensor based on cross-relation method.* Review of Scientific Instruments. https://doi.org/10.1063/1.5019965  
   Related unknown-common-input and cross-relation approaches.

5. **Biteus, J. (2005).** *Distributed Diagnosis and Simulation Based Residual Generators.* Linköping University. https://www.vehicular.isy.liu.se/Publications/Lic/05_LIC_1176_JB.pdf  
   Fault-diagnosis context, including pairwise sensor conflicts and candidate faulty-component explanations.

## Mathematical foundations

6. **Composite trapezoidal quadrature error:** the standard error estimate for twice-differentiable functions, used in the observation bound. See a standard numerical analysis text for precise regularity conditions.
7. **Graph matching and vertex cover:** established graph theory; the lower-bound argument follows from a vertex-disjoint collection of conflicting pairs and is **not** claimed to be a new graph-theoretic theorem.

## Additional nearby work / patent landscape

The [novelty audit](CIFC_novelty_audit_2026-11-11.md) includes additional bounded-error identification, set-membership consistency, sensor-network research, and adjacent patent links. **Some sources have only been screened through abstracts or available online material, not systematically evaluated theorem-by-theorem.** Patent links are informational, not legal clearance.

## Suggestions?

Please submit a GitHub issue with **full bibliographic information, link/DOI, and a short explanation of which equation, theorem, or experiment the source anticipates.** Earlier counterexamples and negative replication studies are welcome.
