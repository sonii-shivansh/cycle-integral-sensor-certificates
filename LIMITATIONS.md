# Limitations, assumptions, and scope

CIFC is an **exploratory, one-sided model-inconsistency certificate**. It has not been independently peer reviewed or validated on physical sensors.

## Exact assumptions behind the bound

For each healthy sensor \(i\):
- A **common** continuous temperature input \(T(t)\) exists for every healthy sensor.
- A **constant first-order response time** \(\tau_i\in[0,M]\) satisfies \(\tau_i y_i'(t)+y_i(t)=T(t)\).
- Time samples are **synchronous and uniformly spaced** across the tested interval.
- Every measured sample satisfies an **absolute deterministic error bound** \(|z_i(t_k)-y_i(t_k)|\le\epsilon_i\).
- Every true healthy sensor trace satisfies a trusted curvature bound \(|y_i''(t)|\le K_i\) over the entire interval.
- Parameters \(M,\epsilon_i,K_i\) are **physical bounds**, not estimated optimistically after seeing the data.

The provided script uses a common noise limit and curvature bound for all sensors, and a single synthetic steady-state sinusoidal example.

## What a detected contradiction means

If a pair strictly violates the derived inequality, **both sensor models cannot simultaneously satisfy all assumptions**. At least one sensor must be *incompatible with the model*. A hardware device may still be working correctly if the common-input assumption, synchronization, or model itself is invalid.

## What a non-detection means

A passing pair is **not certified healthy**. Small relative biases, changes that cancel over the observation interval, shared deviations, and other faults may be invisible.

## What the graph calculation means

- A matching of size \(r\) proves that at least \(r\) distinct sensor-model assumptions must be incompatible when all tested assumptions apply.
- A minimum vertex cover is a **smallest explanation for the conflicts**. It cannot automatically identify failed physical components.
- The exhaustive matching and vertex-cover implementation is intended only for the six-sensor demonstration; do not scale it naively to large networks.

## Fundamentally unobservable common-mode errors

An identical offset applied to every sensor can be absorbed into a shifted shared temperature input. An external absolute reference is needed to distinguish those worlds.

## Practical constraints

- Temperature gradients, sensor placement differences, clock drift, radiation/conduction effects, nonconstant dynamics, or sensor nonlinearities may invalidate assumptions.
- Stated numerical guarantees require **worst-case noise and curvature bounds**; an informal standard-deviation estimate is not interchangeable with a deterministic limit.
- This is **not a deployed safety system**. Do not use its warnings or absence of warnings as the sole basis for critical control, industrial safety, or maintenance decisions.

## Benchmark limitations

Our published test uses one steady-state sinusoidal period with deliberately injected fixed offsets. A simple full-period average consistency check finds the same nine conflicting pairs and has a **tighter** sufficient guaranteed threshold in that special example. CIFC has not demonstrated improved sensitivity, speed, or fault localization against established competing methods.

The 100 healthy synthetic trials are regression checks, **not** empirical field false-alarm rates. Larger evaluation requires independent hardware data, nonperiodic inputs, sensitivity studies, and comparisons with existing interval/parity-based diagnosis literature.

## Novelty and citation limits

Substantial prior art exists. The particular closed-form inequality and its packaging have not been proven unique. Bibliographic searches cannot cover all unpublished or private research. Consult [REFERENCES.md](REFERENCES.md) and the [novelty audit](CIFC_novelty_audit_2026-11-11.md).
