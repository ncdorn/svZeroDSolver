from structured_tree_simulation import *
from struct_tree_utils import *


if __name__ == '__main__':
    # name of model to analyze
    model_name = 'LPA_RPA_0d_steady'
    # path to experiment file (json to be loaded as dict)
    expname = 'stenosis_degree_test_7.28.23.txt'
    '''
        Requirements of experiment params json file:
        "name": name of experiment
        "repair type": extensive, proximal or a list of vessels to repair
        "repair degrees": list of ints with degrees of stenosis repair
        '''
    run_simulation(model_name, expname, optimized=True, vis_trees=True)
