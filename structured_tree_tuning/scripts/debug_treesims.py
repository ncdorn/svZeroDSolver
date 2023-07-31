import pandas as pd
import numpy as np
import json
from pathlib import Path
from svzerodsolver import runner
from struct_tree_utils import *
from structured_tree_simulation import optimize_preop_0d

def return_dataframe(input_file):

    with open(input_file) as ff:
        zeroD_config = json.load(ff)
    zeroD_result = runner.run_from_config(zeroD_config)

    return zeroD_result

def get_config(input_file):
    with open(input_file) as ff:
        zeroD_config = json.load(ff)

    return zeroD_config

def test_get_resistance_idx(config):
    idxs = []
    for vessel_config in config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:
                idxs.append(get_resistance_idx(vessel_config))

    return idxs


if __name__ == '__main__':
    test_dir = Path("tree_tuning_test")
    config = get_config(test_dir / "LPA_RPA_0d_steady.in")
    idxs = test_get_resistance_idx(config)
    print(idxs)
