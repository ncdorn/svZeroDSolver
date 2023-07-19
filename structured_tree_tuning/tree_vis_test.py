from st_post_process import *
from pathlib import Path

if __name__ == '__main__':
    test_dir = Path("tree_tuning_test")
    dirname = 'LPA_RPA_0d_steady'
    filepath = test_dir / dirname / "adapted_config.txt"

    visualize_trees(filepath)
