import numpy as np
from .rh_mynardvalve import RH_Mynard_Valve
from .rvchambermodel import RV_Chamber_Model
from .pressurereferencebc import PressureReferenceBC

def create_RH_LPN(connections, block_dict, RH_constants, RH_params, period):
    """
    creates the RH LPN and adds it to the connections and block_dict
    Args:
        connections: list of connections in create_blocks from utils.py
        block_dict: block_dict from create_blocks in utils.py
        RH_constants: config file containing constants for the RH LPN
        RH_params: config dict containing tunable RH parameters

    Returns:
        updated connections list and block_dict with the mynard valve and rv chamber class instances

    """
    # right atrium
    ra_params = dict(
        time=[0.0, 1.0],
        P=[RH_constants['P_RA']] * 2,
    )
    RA = PressureReferenceBC(name='RA', params=ra_params)
    # tricuspid valve
    TV = RH_Mynard_Valve(name='TV', rho=RH_constants['rho'], l_v=RH_constants['l_tv'], S_v=RH_constants['S_TV'], K_vo=RH_constants['K_voT'], K_vc=RH_constants['K_vcT'])
    # right ventricle
    RV = RV_Chamber_Model(name='RV', Tc=period, v0=RH_constants['V_rv0'], RH_params=RH_params)
    # pulmonary valve
    PV = RH_Mynard_Valve(name='PV', rho=RH_constants['rho'], l_v=RH_constants['l_pv'], S_v=RH_constants['S_PV'], K_vo=RH_constants['K_voP'], K_vc=RH_constants['K_vcP'])

    # add th connections to the connections list
    connections.append(('RA', 'TV'))
    connections.append(('TV', 'RV'))
    connections.append(('RV', 'PV'))
    # add the blocks and names to the block dict
    block_dict['RA'] = RA
    block_dict['TV'] = TV
    block_dict['RV'] = RV
    block_dict['PV'] = PV