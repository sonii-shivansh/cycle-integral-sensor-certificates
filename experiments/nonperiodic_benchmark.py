"""Reproducible nonperiodic synthetic benchmark of two necessary tests.

Healthy traces are exact first-order solutions of an unknown shared input:
T(t) = 25 + drift*t + 6*sin(0.25*t+0.3) + 2*cos(0.4*t-0.2).

The input is not periodic over the 13.5-second observation interval. The
underlying response curvature obeys the original implementation's bound.

The example is an illustrative stress test, not real-world validation or
proof of scientific priority. Both certificates can fail to detect biases.
"""

import numpy as np

import cycle_integral_fault_certificate as cifc

DURATION = 13.5
STEP = 0.05
W1 = 0.25
W2 = 0.4
A1 = 6.0
A2 = 2.0
CURVATURE_BOUND = A1 * W1**2 + A2 * W2**2
assert CURVATURE_BOUND <= cifc.CURVATURE_MAX


def traces(t, drift):
    """Exact responses to a linear drift plus two harmonic inputs."""
    outputs = []
    for tau in cifc.TAUS:
        p1 = W1 * t + 0.3
        p2 = W2 * t - 0.2
        outputs.append(
            25.0 + drift * (t - tau)
            + A1 * (np.sin(p1) - W1 * tau * np.cos(p1)) / (1.0 + (W1 * tau)**2)
            + A2 * (np.cos(p2) + W2 * tau * np.sin(p2)) / (1.0 + (W2 * tau)**2)
        )
    return np.asarray(outputs)


def generate(seed=123, offsets=None, drift=0.3):
    t = np.linspace(0.0, DURATION, round(DURATION / STEP) + 1)
    ys = traces(t, drift)
    for index, bias in (offsets or {}).items():
        ys[index] += bias
    rng = np.random.default_rng(seed)
    observations = ys + rng.uniform(-cifc.EPS, cifc.EPS, size=ys.shape)
    return t, observations


def validate_integral_identity():
    """Check numerical quadrature against the exact ODE identity."""
    t = np.linspace(0.0, DURATION, round(DURATION / STEP) + 1)
    y = traces(t, 0.3)
    qbound = DURATION * STEP**2 * 2 * cifc.CURVATURE_MAX / 12.0
    errors = []
    for i in range(len(cifc.TAUS)):
        for j in range(i + 1, len(cifc.TAUS)):
            observed = np.trapezoid(y[j] - y[i], t)
            exact = (cifc.TAUS[i] * (y[i, -1] - y[i, 0])
                     - cifc.TAUS[j] * (y[j, -1] - y[j, 0]))
            errors.append(abs(observed - exact))
    assert max(errors) <= qbound + 1e-10, (max(errors), qbound)
    return max(errors), qbound


def compute_results(num_trials=100):
    """Counts summed over seeds; no metrics interpreted as field error rates."""
    rows = []
    for drift in [0.0, 0.3, 0.6, 1.0]:
        old_healthy_edges = new_healthy_edges = 0
        old_fault_edges = new_fault_edges = 0
        old_alarm_trials = new_alarm_trials = 0
        for seed in range(num_trials):
            t, healthy = generate(seed, drift=drift)
            old_he, _ = cifc.certificates(t, healthy)
            new_he, _ = cifc.directional_certificates(t, healthy)
            old_healthy_edges += len(old_he)
            new_healthy_edges += len(new_he)
            assert set(old_he).issubset(set(new_he)), (drift, seed, 'healthy')
            t, injected = generate(seed, offsets={2: 1.2, 5: -1.5}, drift=drift)
            old_fe, _ = cifc.certificates(t, injected)
            new_fe, _ = cifc.directional_certificates(t, injected)
            old_fault_edges += len(old_fe)
            new_fault_edges += len(new_fe)
            old_alarm_trials += bool(old_fe)
            new_alarm_trials += bool(new_fe)
            assert set(old_fe).issubset(set(new_fe)), (drift, seed, 'fault')
        rows.append((drift, old_healthy_edges, new_healthy_edges,
                     old_fault_edges, new_fault_edges, old_alarm_trials,
                     new_alarm_trials))
    return rows


if __name__ == '__main__':
    worst, q = validate_integral_identity()
    print('Maximum trapezoidal identity error: %.6g; bound: %.6g' % (worst, q))
    print('drift | healthy-old | healthy-new | fault-old | fault-new '
          '| trials-with-old-alarm | trials-with-new-alarm')
    for row in compute_results():
        print(' | '.join(str(x) for x in row))
