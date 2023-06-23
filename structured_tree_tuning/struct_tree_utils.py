import csv
import numpy as np
import math
# utilities for working with structured trees

def get_branch_pressure(result_df, branch_name):
    # get the mpa pressure array for a given branch name from a zerod solver output file
    pressures = get_df_data(result_df, 'pressure_in', branch_name)
    systolic_p = np.min(pressures)
    diastolic_p = np.max(pressures)
    mean_p = np.mean(pressures)
    return pressures, systolic_p, diastolic_p, mean_p


def get_outlet_flowrate(config, result_df, steady=False):
    outlet_vessels = []
    for vessel_config in config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:
                outlet_vessels.append(vessel_config["vessel_id"])

    q_out = []
    for vessel in outlet_vessels:
        q_out.append(get_df_data(result_df, 'flow_out', vessel))
    if steady:
        return [q[-1] for q in q_out]
    else:
        return q_out


def get_df_data(result_df, data_name, branch_name):
    # get data from a dataframe based on a data name and branch name
    data = []
    for idx, name in enumerate(result_df['name']):
        if str(branch_name) in name:
            data.append(result_df.loc[idx, data_name])
    return data


def get_resistances(config):
    # get the resistances from an input config file
    # return a list of resistances
    resistance = []
    for bc_config in config["boundary_conditions"]:
        if bc_config["bc_type"] == 'RESISTANCE':
            resistance.append(bc_config['bc_values'].get('R'))
    np.array(resistance)
    return resistance


def write_resistances(config, resistances):
    idx = 0
    for bc_config in config["boundary_conditions"]:
        if bc_config["bc_type"] == 'RESISTANCE':
            bc_config['bc_values']['R'] = resistances[idx]
            idx += 1


def get_value_from_csv(csv_file, name):
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if name.lower() in row[0].lower():
                return row[1]  # Return the value in the same row

    return None  # Return None if the name is not found


def get_resistance_idx(vessel_config):
    name = vessel_config["boundary_conditions"]["outlet"]
    str_idx = 10
    idx = name[str_idx:]
    while not idx.isdigit():
        str_idx += 1
        idx = name[str_idx:]

    return int(idx)


def repair_stenosis(config, vessels=None, repair_all=False, log_file=None):
    if repair_all:
        with open(log_file, "a") as log:
            log.write("repairing all stenoses")
        for vessel_config in config["vessels"]:
            if vessel_config["zero_d_element_values"]["stenosis_coefficient"] > 0:
                vessel_config["zero_d_element_values"]["stenosis_coefficient"] = 0
                with open(log_file, "a") as log:
                    log.write("     vessel " + str(vessel_config["vessel_id"]) + " has been repaired \n")
    else:
        print("** repairing stenoses in vessels " + str(vessels) + " **")
        for vessel_config in config["vessels"]:
            if vessel_config["vessel_id"] in vessels:
                vessel_config["zero_d_element_values"]["stenosis_coefficient"] = 0
                with open(log_file, "a") as log:
                    log.write("     vessel " + str(vessel_config["vessel_id"]) + " has been repaired \n")


# this is some random code that may need to be used in the non-steady state case
'''
def scale_inflow(config, Q_target):
    for bc_config in config['boundary_conditions']:
        if bc_config["bc_name"] == 'inflow':

    Q_avg = np.mean(inflow_wave)
    scale_ratio = Q_target / Q_avg
    scaled_inflow = scale_ratio * inflow_wave
    return scaled_inflow
'''