"""An exact physical ambiguity: measurement bias vs local temperature offset.

For first-order sensors tau_i*y_i' + y_i = T, constant observation biases
b_i are observationally indistinguishable from perfectly calibrated sensors
whose LOCAL inputs are T_i(t)=T(t)+b_i, with correspondingly shifted initial
conditions. This is an algebraic NON-IDENTIFIABILITY result, not fault detection.

Demonstration includes matched synthetic nonperiodic data and a *sensitivity*
study: using an independently justified input_spread bound can prevent model
contradictions due solely to temperature heterogeneity. In this simulation
the environmental allowance is deliberately chosen equal to the injected
offset magnitudes, which would NOT be known from readings alone.
"""

import numpy as np

from experiments.multiwindow_lp import (
    summarize_window, joint_relaxation, pairwise_results
)
from experiments.multiwindow_case_study import TAUS, EPS, K, M


OFFSETS = np.array([0.0, 0.0, 1.2, 0.0, 0.0, -1.5])
WINDOWS = ((0, 60), (60, 120), (120, 180), (180, 270))


def analytic_first_order_sensors(t, drift):
    """Exactly solve the common-input first-order ODE for the chosen harmonics."""
    t = np.asarray(t, dtype=float)
    p1, p2 = .25 * t + .3, .4 * t - .2
    signals, derivatives = [], []
    for tau in TAUS:
        a, b = .25*tau, .4*tau
        response = (
            25.0 + drift * (t-tau)
            + 6.0 * (np.sin(p1) - a*np.cos(p1))/(1+a*a)
            + 2.0 * (np.cos(p2) + b*np.sin(p2))/(1+b*b)
        )
        derivative = (
            drift + 6.0*.25*(np.cos(p1)+a*np.sin(p1))/(1+a*a)
            + 2.0*.4*(-np.sin(p2)+b*np.cos(p2))/(1+b*b)
        )
        signals.append(response)
        derivatives.append(derivative)
    T = 25.0 + drift*t + 6.0*np.sin(p1) + 2.0*np.cos(p2)
    return T, np.asarray(signals), np.asarray(derivatives)


def indistinguishable_worlds(seed=0, drift=.3):
    """Return two *identical* observed arrays from different physical causes.

    World A: all sensors see the same T, but sensors 3 and 6 have additive
             measurement biases outside the bounded regular noise.
    World B: all sensors are perfectly calibrated apart from small regular
             observation noise, but local physical temperatures are T+b_i,
             with shifted sensor trajectories and initial conditions.

    Both worlds use the same bounded noise array, so indistinguishability
    holds sample-for-sample, not merely as a statistical observation.
    """
    t = np.linspace(0.0,13.5,271)
    common_input, healthy, derivatives = analytic_first_order_sensors(t, drift)
    noise = np.random.default_rng(seed).uniform(
        -EPS,EPS,size=healthy.shape
    )
    # World A: faults appear at the OUTPUT/measurement stage.
    z_bias = (healthy + OFFSETS[:,None]) + noise
    # World B: shifted LOCAL input and corresponding exact physical response.
    local_inputs = common_input[None,:] + OFFSETS[:,None]
    local_responses = healthy + OFFSETS[:,None]
    z_environment = local_responses + noise
    # Direct differential identity checks against each world's own inputs.
    residual_common = TAUS[:,None]*derivatives + healthy - common_input
    residual_local = TAUS[:,None]*derivatives + local_responses - local_inputs
    return {
        "t":t,
        "common_input":common_input,
        "healthy_responses":healthy,
        "noise":noise,
        "z_output_bias":z_bias,
        "z_local_temperature":z_environment,
        "common_ode_residual_max":float(np.max(np.abs(residual_common))),
        "local_ode_residual_max":float(np.max(np.abs(residual_local))),
    }


def summarize_observations(t,z):
    return [summarize_window(t[a:b+1],z[:,a:b+1]) for a,b in WINDOWS]


def evaluate_trial(seed=0, drift=.3):
    data = indistinguishable_worlds(seed,drift)
    # The arrays are identical for ALL individual observations: a diagnostic
    # receiving only z can never determine which physical world generated it.
    if not np.array_equal(
        data["z_output_bias"],data["z_local_temperature"]
    ):
        raise AssertionError("two worlds should be observation-identical")
    w = summarize_observations(data["t"],data["z_output_bias"])
    baseline = pairwise_results(w,M,EPS,K)
    protected = pairwise_results(w,M,EPS,K,input_spread=np.abs(OFFSETS))
    baseline_bad = {pair for pair,report in baseline.items()
                    if report.status=="solver_infeasible"}
    protected_bad = {pair for pair,report in protected.items()
                     if report.status=="solver_infeasible"}
    if any(v.status=="inconclusive" for v in
           (*baseline.values(),*protected.values())):
        raise RuntimeError("floating LP result inconclusive")
    nominal_global = joint_relaxation(w,M,EPS,K)
    protected_global = joint_relaxation(
        w,M,EPS,K,input_spread=np.abs(OFFSETS)
    )
    if nominal_global.status=="inconclusive" or protected_global.status=="inconclusive":
        raise RuntimeError("floating global LP result inconclusive")
    if protected_bad or protected_global.status!="feasible_relaxation":
        raise AssertionError("valid temperature-heterogeneity model falsely rejected")
    return {
        "drift":drift,
        "seed":seed,
        "nominal_edges":len(baseline_bad),
        "bounded_spread_edges":len(protected_bad),
        "nominal_global":nominal_global.status,
        "bounded_spread_global":protected_global.status,
    }


def run(trials=20):
    rows=[]
    for drift in (0.0,0.3,0.6,1.0):
        totals=[0,0]
        global_rejections=[0,0]
        for seed in range(trials):
            outcome=evaluate_trial(seed,drift)
            totals[0]+=outcome["nominal_edges"]
            totals[1]+=outcome["bounded_spread_edges"]
            global_rejections[0]+=outcome["nominal_global"]=="solver_infeasible"
            global_rejections[1]+=outcome["bounded_spread_global"]=="solver_infeasible"
        rows.append((drift,*totals,*global_rejections))
    return rows


if __name__=="__main__":
    print("drift,edges_common_input,edges_with_spatial_allowance,"
          "trials_global_inconsistent_common,with_spatial_allowance")
    for row in run():
        print(*row,sep=",")
    print("Two causes produce identical z samples for every tested seed.")
