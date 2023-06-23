import json
from svzerodsolver.model.structuredtreebc import StructuredTreeOutlet
from pathlib import Path

def test_tree_adapt(config):
    simparams = config["simulation_parameters"]
    test_dict = {}
    for vessel_config in config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:
                outlet_stree = StructuredTreeOutlet.from_outlet_vessel(vessel_config, simparams)
                # outlet_stree.build_tree_olufsen()
                #R_mat = outlet_stree.calculate_resistance()
                outlet_stree.optimize_tree_radius(Resistance=5.0)
                # R_1 = outlet_stree.calculate_resistance()
                # outlet_stree.adapt_constant_wss(2, 1, disp=True)

    return outlet_stree.block_dict


def run_from_file(input_file, output_file):
    """Run the svZeroDSolver from file.

    Args:
        input_file: Input file with configuration.
        output_file: Output file with configuration.
    """
    with open(input_file) as ff:
        config = json.load(ff)
        result = test_tree_adapt(config)
    # with open(output_file, "w") as ff:
    #     json.dump(result, ff)
    return result


if __name__ == '__main__':
    test_dir = Path("../../zerodsolvertests/struct_tree_test")
    result = run_from_file(test_dir / "struct_tree_test.in", test_dir / "struct_tree_test.out")