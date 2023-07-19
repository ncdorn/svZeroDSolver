from svzerodsolver import runner
import csv
from pathlib import Path
import numpy as np
import json
from struct_tree_utils import *
from scipy.optimize import minimize, Bounds
from svzerodsolver.model.structuredtreebc import StructuredTreeOutlet



def optimize_preop_0d(clinical_targets: csv, input_file, log_file, unsteady=False, change_to_R=False):
    '''

    :param clinical_targets: clinical targets input csv
    :param input_file: 0d solver json input file name string
    :param output_file: 0d solver json output file name string
    :return: preop simulation with optimized BCs
    '''
    # get the clinical target values
    with open(log_file, "w") as log:
        log.write("Getting clinical target values... \n")
    bsa = float(get_value_from_csv(clinical_targets, 'bsa'))
    cardiac_index = float(get_value_from_csv(clinical_targets, 'cardiac index'))
    q = bsa * cardiac_index # cardiac output in L/min
    mpa_pressures = get_value_from_csv(clinical_targets, 'mpa pressures') # mmHg
    mpa_sys_p_target = int(mpa_pressures[0:2])
    mpa_dia_p_target = int(mpa_pressures[3:5])
    mpa_mean_p_target = int(get_value_from_csv(clinical_targets, 'mpa mean pressure'))
    target_ps = np.array([
        mpa_sys_p_target,
        mpa_dia_p_target,
        mpa_mean_p_target
    ])
    target_ps = target_ps * 1333.22 # convert to barye

    # load input json as a config dict
    with open(input_file) as ff:
        zeroD_config = json.load(ff)

    if not unsteady:
        make_inflow_steady(zeroD_config)
        with open(log_file, "a") as log:
            log.write("inlet BCs converted to steady \n")

    if change_to_R:
        Pd = convert_RCR_to_R(zeroD_config)
        with open(log_file, "a") as log:
            log.write("RCR BCs converted to R, Pd = " + str(Pd) + "\n")

    # get resistances from the zerod input file
    resistance = get_resistances(zeroD_config)
    # get the LPA and RPA branch numbers
    lpa_rpa_branch = ["V" + str(idx) for idx in zeroD_config["junctions"][0]["outlet_vessels"]]

    # scale the inflow

    # run zerod simulation to reach clinical targets
    def zerod_optimization_objective(r,
                                     input_config=zeroD_config,
                                     target_ps=None,
                                     unsteady=unsteady,
                                     lpa_rpa_branch=lpa_rpa_branch
                                     ):
        # r = abs(r)
        # r = [r, r]
        write_resistances(input_config, r)
        zerod_result = runner.run_from_config(input_config)
        mpa_pressures, mpa_sys_p, mpa_dia_p, mpa_mean_p  = get_mpa_pressure(zerod_result, branch_name='V0') # get mpa pressure

        # lpa_rpa_branch = ["V" + str(idx) for idx in input_config["junctions"][0]["outlet_vessels"]]

        q_MPA = get_df_data(zerod_result, branch_name='V0', data_name='flow_in')
        q_RPA = get_df_data(zerod_result, branch_name=lpa_rpa_branch[0], data_name='flow_in')
        q_LPA = get_df_data(zerod_result, branch_name=lpa_rpa_branch[1], data_name='flow_in')
        # print(q_RPA[-1], q_LPA[-1])
        if unsteady:
            pred_p = np.array([
                mpa_sys_p,
                mpa_dia_p,
                mpa_mean_p
            ])
            p_diff = np.sum(np.square(np.subtract(pred_p, target_ps)))
        else:
            p_diff = abs(target_ps[2] - mpa_mean_p) ** 2
        # SSE = np.sum(np.square(np.subtract(pred_p, target_p)))
        # MSE = np.square(np.subtract(pred_p, target_p)).mean()

        RPA_diff = abs((q_RPA[-1] - (0.52 * q_MPA[-1]))) ** 2
        min_obj = p_diff + RPA_diff
        return min_obj

    # write to log file for debugging
    with open(log_file, "a") as log:
        log.write("Optimizing preop outlet resistance... \n")
    # run the optimization algorithm
    result = minimize(zerod_optimization_objective, resistance, args=(zeroD_config, target_ps), method="CG", options={"disp": True})
    # write to log file for debugging
    with open(log_file, "a") as log:
        log.write("Outlet resistances optimized! " + str(result.x) +  "\n")

    R_final = result.x # get the array of optimized resistances
    write_resistances(zeroD_config, R_final)
    return zeroD_config, R_final


def construct_trees(config: dict, resistances=None, log_file=None):
    for vessel_config in config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:
                outlet_tree = StructuredTreeOutlet.from_outlet_vessel(vessel_config, config["simulation_parameters"])
                R = resistances[get_resistance_idx(vessel_config)]
                # write to log file for debugging
                with open(log_file, "a") as log:
                    log.write("** building tree for resistance: " + str(R) + " ** \n")
                # outlet_tree.optimize_tree_radius(R)
                outlet_tree.optimize_tree_radius(R, log_file)
                # write to log file for debugging
                with open(log_file, "a") as log:
                    log.write("     the number of vessels is " + str(len(outlet_tree.block_dict["vessels"])) + "\n")
                vessel_config["tree"] = outlet_tree.block_dict


def calculate_flow(config: dict, repair=None, repair_degree=1, log_file=None):
    preop_result = runner.run_from_config(config)
    lpa_rpa_branch = [idx for idx in config["junctions"][0]["outlet_vessels"]]
    repair_vessels=None
    # repair_vessels based on desired strategy: extensve, proximal or custom
    if repair == 'proximal':
        repair_vessels = lpa_rpa_branch
    elif repair == 'extensive':
        pass
    else:
        repair_vessels=repair

    repair_stenosis(config, repair_vessels, repair_degree, log_file=log_file)
    postop_result = runner.run_from_config(config)

    return preop_result, postop_result


def adapt_trees(config, preop_result, postop_result):
    preop_q = get_outlet_q(config, preop_result, steady=True)
    postop_q = get_outlet_q(config, postop_result, steady=True)
    # q_diff = [postop_q[i] - q_old for i, q_old in enumerate(preop_q)]
    # print("the change in q is: " + str(q_diff))
    adapted_config = config
    outlet_idx = 0 # index through outlets
    R_new = []

    for vessel_config in adapted_config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:
                outlet_tree = StructuredTreeOutlet.from_outlet_vessel(vessel_config, config["simulation_parameters"], tree_config=True)
                R_new.append(outlet_tree.adapt_constant_wss(preop_q[outlet_idx], postop_q[outlet_idx], disp=False))
                vessel_config["tree"] = outlet_tree.block_dict
                outlet_idx +=1

    write_resistances(adapted_config, R_new)
    return adapted_config


def run_final_flow(config, preop_result, postop_result, output_file, log_file):
    final_result = runner.run_from_config(config)
    with open(log_file, "a") as log:
        log.write("Writing result to file... \n")
    result = {'name': 'results!', 'data': final_result.to_dict('index')}
    with open(output_file, "w") as ff:
        json.dump(result, ff)

    preop_q = get_outlet_q(config, preop_result, steady=True)
    postop_q = get_outlet_q(config, postop_result, steady=True)
    final_q = get_outlet_q(config, final_result, steady=True)
    with open(log_file, "a") as log:
        log.write("** RESULT COMPARISON ** \n")
        log.write("     preop outlet flowrates: " + str(preop_q) + "\n")
        log.write("     pre-adaptation outlet flowrates: " + str(postop_q) + "\n")
        log.write("     post-adaptation outlet flowrates: " + str(final_q) + "\n \n")
        log.write("Simulation completed!")


def generate_results(result):
    pass


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    test_dir = Path("tree_tuning_test")
    dirname = 'LPA_RPA_0d_steady'

    input_file = test_dir / dirname / '{}.in'.format(dirname)
    log_file = test_dir / dirname / '{}.log'.format(dirname)
    output_file = test_dir / dirname / '{}.out'.format(dirname)

    config, R_final = optimize_preop_0d(test_dir / 'clinical_targets.csv',
                                     input_file, log_file, unsteady=False, change_to_R=True)

    construct_trees(config, R_final, log_file)

    preop_result, postop_result = calculate_flow(config, repair='proximal', log_file=log_file)

    adapted_config = adapt_trees(config, preop_result, postop_result)

    final_result = run_final_flow(adapted_config, preop_result, postop_result, output_file, log_file)


    # write the config to a file for observation
    with open(test_dir / dirname / "adapted_config.txt", "w") as ff:
        json.dump(adapted_config, ff)




# See PyCharm help at https://www.jetbrains.com/help/pycharm/
