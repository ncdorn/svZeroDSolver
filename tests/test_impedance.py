import pysvzerod
import json
from svzerodtrees.utils import *

def test_impedance():
    '''
    test the impedance boundary condition
    '''

    with open('cases/pulsatileFlow_IMPEDANCE.json') as f:
        config = json.load(f)
    
    result = pysvzerod.simulate(config)

    with open('cases/results/result_pulsatileFlow_IMPEDANCE.json', 'w') as f:
        f.write(result.to_json())

    fig, axs = plt.subplots(2, 2)

    axs[0, 0].plot(result['time'], result['pressure_out'] / 1333.2, label='Pressure')
    axs[0, 0].set_title('Pressure')

    axs[0, 1].plot(result['time'], result['flow_out'], label='Flow')
    axs[0, 1].set_title('Flow')

    axs[1, 0].plot(result['flow_out'], result['pressure_out'] / 1333.2, label='Pressure-Flow')


    plt.show()



if __name__ == '__main__':
    test_impedance()




