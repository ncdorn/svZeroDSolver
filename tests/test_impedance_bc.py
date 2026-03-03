import copy
import math

import numpy as np
import pytest

import pysvzerod


def _base_impedance_config():
    period = 1.0
    num_period_steps = 8
    times = np.linspace(0.0, period, num_period_steps + 1)
    flows = [1.0 + 0.4 * math.sin(2.0 * math.pi * t / period) for t in times]

    return {
        "simulation_parameters": {
            "number_of_cardiac_cycles": 4,
            "number_of_time_pts_per_cardiac_cycle": num_period_steps + 1,
            "output_all_cycles": True,
            "steady_initial": False,
            "absolute_tolerance": 1e-10,
        },
        "boundary_conditions": [
            {
                "bc_name": "INFLOW",
                "bc_type": "FLOW",
                "bc_values": {"Q": flows, "t": list(times)},
            },
            {
                "bc_name": "OUT",
                "bc_type": "IMPEDANCE",
                "bc_values": {
                    "period": period,
                    "Pd": 2.0,
                    "z": [4.0 * (0.25**i) for i in range(num_period_steps)],
                },
            },
        ],
        "vessels": [
            {
                "boundary_conditions": {"inlet": "INFLOW", "outlet": "OUT"},
                "vessel_id": 0,
                "vessel_length": 10.0,
                "vessel_name": "branch0_seg0",
                "zero_d_element_type": "BloodVessel",
                "zero_d_element_values": {"R_poiseuille": 10.0},
            }
        ],
    }


def _run(config):
    return pysvzerod.simulate(config)


def test_impedance_exact_matches_discrete_reference():
    config = _base_impedance_config()
    result = _run(config)

    pd = config["boundary_conditions"][1]["bc_values"]["Pd"]
    kernel = config["boundary_conditions"][1]["bc_values"]["z"]
    q = result["flow_out"].to_numpy()
    p = result["pressure_out"].to_numpy()
    rho_infty = config["simulation_parameters"].get("rho_infty", 0.5)
    alpha_f = 1.0 / (1.0 + rho_infty)

    num_period_steps = len(kernel)
    head = num_period_steps - 1
    accepted_steps = 0
    history = np.zeros(num_period_steps)
    p_ref = np.zeros_like(p)
    p_ref[0] = p[0]

    for n in range(1, len(q)):
        conv_sum = 0.0
        if accepted_steps > num_period_steps:
            for m in range(1, num_period_steps):
                idx = (head - (m - 1)) % num_period_steps
                conv_sum += kernel[m] * history[idx]

        # Algebraic equations are enforced at generalized-alpha midpoint:
        # y_af = y_n + alpha_f * (y_{n+1} - y_n).
        rhs = pd + conv_sum
        s_prev = p_ref[n - 1] - kernel[0] * q[n - 1]
        s_new = (rhs - (1.0 - alpha_f) * s_prev) / alpha_f
        p_ref[n] = s_new + kernel[0] * q[n]

        head = (head + 1) % num_period_steps
        history[head] = q[n]
        accepted_steps += 1

    assert np.max(np.abs(p[1:] - p_ref[1:])) < 1e-8


def test_impedance_fails_on_kernel_period_mismatch():
    config = _base_impedance_config()
    config["boundary_conditions"][1]["bc_values"]["z"] = [1.0] * 7

    with pytest.raises(RuntimeError, match="kernel length|Expected"):
        _run(config)


def test_impedance_fails_on_cardiac_period_mismatch():
    config = _base_impedance_config()
    config["simulation_parameters"]["cardiac_period"] = 1.2

    with pytest.raises(RuntimeError, match="Inconsistent cardiac cycle period"):
        _run(config)


def test_impedance_truncated_full_terms_matches_exact():
    exact_config = _base_impedance_config()
    trunc_config = copy.deepcopy(exact_config)

    z = trunc_config["boundary_conditions"][1]["bc_values"]["z"]
    trunc_config["boundary_conditions"][1]["bc_values"]["convolution_mode"] = "truncated"
    trunc_config["boundary_conditions"][1]["bc_values"]["num_kernel_terms"] = len(z)

    exact = _run(exact_config)
    trunc = _run(trunc_config)

    assert np.allclose(
        exact["pressure_out"].to_numpy(),
        trunc["pressure_out"].to_numpy(),
        atol=1e-10,
        rtol=1e-10,
    )


def test_impedance_truncated_has_bounded_error():
    exact_config = _base_impedance_config()
    trunc_config = copy.deepcopy(exact_config)

    trunc_config["boundary_conditions"][1]["bc_values"]["convolution_mode"] = "truncated"
    trunc_config["boundary_conditions"][1]["bc_values"]["num_kernel_terms"] = 3

    exact = _run(exact_config)
    trunc = _run(trunc_config)

    max_abs_diff = np.max(
        np.abs(
            exact["pressure_out"].to_numpy() - trunc["pressure_out"].to_numpy()
        )
    )
    assert max_abs_diff < 2e-1
