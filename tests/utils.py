import json
import os
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

# global boolean to perform coverage testing
# (run executables instead of Python interface, much slower)
from pytest import coverage

import pysvzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

RTOL_PRES = 1.0e-7
RTOL_FLOW = 1.0e-7


def execute_pysvzerod(testfile, mode):
    """Execute pysvzerod (via Python interface or executable).

    Args:
        testfile: Path to the input file.
        mode: svZeroDSolver application (solver or calibrator).
    """
    assert mode in ["solver", "calibrator"], "unknown mode: " + mode

    # read configuration
    with open(testfile) as ff:
        config = json.load(ff)

    if coverage:
        # run via executable (slow)
        with TemporaryDirectory() as tempdir:
            out_name = os.path.join(tempdir, "out")
            exe = os.path.join(this_file_dir, "..", "Release", "svzerod")
            subprocess.run([exe + mode, testfile, out_name])
            if mode == "solver":
                result = pd.read_csv(out_name)
            elif mode == "calibrator":
                with open(out_name) as ff:
                    result = json.load(ff)
    else:
        # run via Python binding (fast)
        if mode == "solver":
            result = pysvzerod.simulate(config)
        elif mode == "calibrator":
            result = pysvzerod.calibrate(config)

    return result, config


def run_with_reference(
        ref,
        test_config
        ):


    res = pysvzerod.simulate(test_config)

    if res.shape[1] >= 6:
        # we have a result with fields [name, time, p_in, p_out, q_in, q_out] SOME HAVE GREATTER LENGTH NEED TO ADDRESS
        pressure_cols = ["pressures_in", "pressure_out"]
        flow_cols = ["flow_in", "flow_out"]

        bool_df = pd.DataFrame(index=res.index, columns=res.columns, dtype=bool)

        for col in res.columns:
            tol = RTOL_PRES if col in pressure_cols else RTOL_FLOW if col in flow_cols else None
            if tol is not None:
                rel_diff = np.abs(res[col] - ref[col]) - tol - tol * np.abs(ref[col])
                check = rel_diff <= 0.0
                bool_df[col] = check  # True if difference is below tolerance

        if not bool_df.all().all():  # Check if any value exceeded tolerance
            # Extract only differing rows/columns for a cleaner error message
            differing_locs = np.where(~bool_df)  # Shows only the out-of-tolerance values
            differing_indices = res.index[differing_locs[0]]
            differing_columns = res.columns[differing_locs[1]]

            if differing_indices.any():  # If any differences are found
                diff_locations = list(zip(differing_indices, differing_columns))
                # find the value where the difference exceeds tolerance
                for i, col in diff_locations:
                    if i in res.index and col in res.columns:
                        # Get the actual values from both result and reference
                        diff_values = (res.at[i, col], ref.at[i, col])
                        print(f"Test failed in field {col}: {diff_values[0]} (result) vs {diff_values[1]} (reference)")
                diff_values = [(res.at[i, col], ref.at[i, col]) for i, col in diff_locations]

                raise AssertionError(f"Differences exceed tolerance.")

    else:
        # we have a result with fields [name, time, y] and the result must be compared based on the name field. name is of format [flow:vessel:outlet]
        # we will compare the average of each branch

        # Merge both datasets on 'Object' and 't' to align values
        res_merged = ref.merge(res, on=["name", "time"], suffixes=("_expected", "_actual"))

        # Compute relative difference
        def compute_relative_difference(row):
            tol = RTOL_FLOW if "flow" in row["name"] else RTOL_PRES
            return abs(row["y_actual"] - row["y_expected"]) - tol - (tol * abs(row["y_expected"]))
        res_merged["Relative_Difference"] = res_merged.apply(compute_relative_difference, axis=1)

        # Apply object-specific tolerance
        res_merged["Within_Tolerance"] = res_merged.apply(
            lambda row: row["Relative_Difference"] <= 0.0, axis=1
        )

        if not res_merged["Within_Tolerance"].all():
            # Extract only differing rows for a cleaner error message
            differing_rows = res_merged[~res_merged["Within_Tolerance"]]
            if not differing_rows.empty:
                diff_info = differing_rows[["name", "time", "y_expected", "y_actual", "Relative_Difference"]]
                print("Test failed in the following rows:\n", diff_info.to_string(index=False))
                raise AssertionError("Differences exceed tolerance.")
            else:
                raise AssertionError("Differences exceed tolerance but no specific rows found.")

def run_test_case_by_name(name, output_variable_based=False, folder="."):
    """Run a test case by its case name.

    Args:
        name: Name of the test case.
        testdir: Directory for performing the simulation.
    """
    # file name of test case
    testfile = os.path.join(this_file_dir, "cases", name + ".json")

    # run test
    result, config = execute_pysvzerod(testfile, "solver")

    if not output_variable_based:
        output = {
            "pressure_in": {},
            "pressure_out": {},
            "flow_in": {},
            "flow_out": {},
        }

        for vessel in config["vessels"]:
            name = vessel["vessel_name"]
            branch_id, seg_id = name.split("_")
            branch_id, seg_id = int(branch_id[6:]), int(seg_id[3:])

            for f in ["pressure", "flow"]:
                for l in ["in"] * (seg_id == 0) + ["out"]:
                    n = f + "_" + l
                    ids = result.name == name
                    output[n][branch_id] = np.array(result[ids][n])

    else:
        output = result

    return output


def get_result(result_array, field, branch, time_step):
    """ "Get results at specific field, branch, branch_node and time step."""
    # extract result
    return result_array[field][branch][time_step]
