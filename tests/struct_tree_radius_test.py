import json
from svzerodsolver.model.structuredtreebc import StructuredTreeOutlet
from pathlib import Path


def test_tree_radius(config):
    simparams = config["simulation_parameters"]
    for vessel_config in config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:
                outlet_stree = StructuredTreeOutlet.from_outlet_vessel(vessel_config, simparams)
                optimal_radius = outlet_stree.optimize_tree_radius(Resistance=5.0)
                print("the optimal radius is r = " + str(optimal_radius))


def run_from_file(input_file, output_file):
    """Run the svZeroDSolver from file.

    Args:
        input_file: Input file with configuration.
        output_file: Output file with configuration.
    """
    with open(input_file) as ff:
        config = json.load(ff)
        test_tree_radius(config)



if __name__ == '__main__':
    test_dir = Path("../../zerodsolvertests/struct_tree_test")
    run_from_file(test_dir / "struct_tree_test.in", test_dir / "struct_tree_test.out")