import os
import json
import numpy as np
import pandas as pd
import pytest
import pysvzerod

import sys
sys.path.append(os.path.dirname(__file__))

from .utils import run_with_reference, RTOL_PRES, RTOL_FLOW

EXPECTED_FAILURES = {
    'closedLoopHeart_singleVessel_mistmatchPeriod.json',
    'pulsatileFlow_R_RCR_mismatchPeriod.json'
}

@pytest.mark.parametrize("testfile", ['steadyFlow_R_R.json', 
                                      'pulsatileFlow_R_coronary_cycle_error.json', 
                                      'pulsatileFlow_R_coronary.json', 
                                      'coupledBlock_closedLoopHeart_withCoronaries.json', 
                                      'pulsatileFlow_R_RCR_derivative.json', 
                                      'steadyFlow_bifurcationR_R1.json', 
                                      'closedLoopHeart_singleVessel.json', 
                                      'pulsatileFlow_R_RCR_variable.json', 
                                      'coupledBlock_closedLoopHeart_singleVessel.json', 
                                      'steadyFlow_R_steadyPressure.json', 
                                      'steadyFlow_bifurcationR_R1_blockNames.json', 
                                      'pulsatileFlow_R_RCR.json', 
                                      'pulsatileFlow_R_impedance.json',
                                      'closedLoopHeart_withCoronaries.json', 
                                      'steadyFlow_RL_R.json', 
                                      'steadyFlow_bifurcationR_R2.json', 
                                      'pulsatileFlow_R_RCR_derivative_variable.json', 
                                      'steadyFlow_confluenceR_R.json', 
                                      'steadyFlow_R_RCR.json', 
                                      'steadyFlow_stenosis_R.json',
                                      'pulsatileFlow_R_RCR_mean.json', 
                                      'pulsatileFlow_R_RCR_mean_derivative.json', 
                                      'pulsatileFlow_CStenosis_steadyPressure.json', 
                                      'pulsatileFlow_R_RCR_mean_variable.json', 
                                      'steadyFlow_RC_R.json', 
                                      'steadyFlow_R_coronary.json', 
                                      'steadyFlow_RLC_R.json', 
                                      'steadyFlow_blood_vessel_junction.json', 
                                      'valve_tanh.json', 
                                      'pulsatileFlow_bifurcationR_RCR_cycle_error.json', 
                                      'pulsatileFlow_R_RCR_mean_derivative_variable.json',
                                      'closedLoopHeart_singleVessel_mistmatchPeriod.json',
                                      'pulsatileFlow_R_RCR_mismatchPeriod.json',
                                      'pulsatileFlow_CStenosis_steadyPressure_definedPeriod.json',
                                      'chamber_sphere.json',
                                      'piecewise_Chamber_and_Valve.json',
                                      'closed_loop_two_hill.json',
                                      'pulsatileFlow_CRL.json',
                                      'pulsatileFlow_R_coronary_varres.json'
                                      ])
def test_solver(testfile):
    '''
    run all test cases and compare against stored reference solution
    '''

    # default tolerances
    rtol_pres = RTOL_PRES
    rtol_flow = RTOL_FLOW
    if 'coupledBlock_closedLoopHeart_withCoronaries.json' in testfile:
        rtol_pres = 2.0e-1
        rtol_flow = 2.0e-1

    this_file_dir = os.path.abspath(os.path.dirname(__file__))

    results_dir = os.path.join(this_file_dir, 'cases', 'results')

    if testfile in EXPECTED_FAILURES:
        pytest.xfail(reason=f"Known failure for test case: {testfile}")

    ref = pd.read_json(os.path.join(results_dir, f'result_{testfile}'))

    run_with_reference(ref, os.path.join(this_file_dir, 'cases', testfile), rtol_pres, rtol_flow)


def test_solver_rejects_deprecated_chamber_elastance_inductor():
    config = {
        "simulation_parameters": {
            "number_of_cardiac_cycles": 1,
            "number_of_time_pts_per_cardiac_cycle": 2
        },
        "boundary_conditions": [],
        "chambers": [{
            "type": "ChamberElastanceInductor",
            "name": "legacy_chamber",
            "values": {
                "Emax": 1.0,
                "Emin": 0.1,
                "Vrd": 10.0,
                "Vrs": 5.0,
                "Impedance": 0.01
            },
            "activation_function": {
                "type": "half_cosine",
                "t_active": 0.2,
                "t_twitch": 0.3
            }
        }]
    }

    with pytest.raises(RuntimeError, match="ChamberElastanceInductor has been removed"):
        pysvzerod.simulate(config)


def test_solver_accepts_canonical_coupled_impedance_config_without_points_per_cycle():
    config = {
        "simulation_parameters": {
            "coupled_simulation": True,
            "number_of_time_pts": 2,
            "external_step_size": 0.001,
            "cardiac_period": 0.004,
            "output_all_cycles": True,
            "steady_initial": False,
            "absolute_tolerance": 1e-10,
        },
        "boundary_conditions": [
            {
                "bc_name": "IMP",
                "bc_type": "IMPEDANCE",
                "bc_values": {
                    "Pd": 5.0,
                    "z": [400.0, 120.0, 40.0, 10.0],
                },
            }
        ],
        "external_solver_coupling_blocks": [
            {
                "name": "FLOW_COUPLING",
                "type": "FLOW",
                "location": "inlet",
                "connected_block": "IMP",
                "periodic": False,
                "values": {
                    "t": [0.0, 0.004],
                    "Q": [1.0, 1.0],
                },
            }
        ],
        "junctions": [],
        "vessels": [],
    }

    result = pysvzerod.simulate(config)

    assert not result.empty
