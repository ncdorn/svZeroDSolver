import numpy as np
from .block import Block

class RH_Mynard_Valve(Block):
    """
        Mynard valve for Right heart modeling
        models a heart valve in 0D according to Mynard et al. (2011)
    """

    _NUM_EQUATIONS = 3
    _NUM_INTERNAL_VARS = 1

    def __init__(self, params: dict = None, name: str = None, rho=None, l_v=None, S_v=None, K_vo=None, K_vc=None):
        """ create a new Mynard valve instance

            Args:
                Configuration parameters of the block. mostly constants.
            Name:
                name of the instance
        """
        super().__init__(params=params, name=name)

        # initialize parameters (maybe take from config file later
        self.rho = rho  # blood density
        self.l_v = l_v  # effective length
        self.S_v = S_v  # effective area (fully opened)
        self.K_vo = K_vo  # valve opening rate constant
        self.K_vc = K_vc  # valve closing rate constant
        self.beta = 0.5 * self.rho / (self.S_v ** 2)
        self.eps = 1.0e-7  # offset factor (prevent O_v = 0)
        # the ordering of the solution variables is : (P_in, Q_in, P_out, Q_out, O_v)

        #  initialize system of equations
        self._mat["E"] = np.zeros((3, 5), dtype=float)
        self._mat["dE"] = np.zeros((3, 5), dtype=float)
        self._mat["F"] = np.array([
            [1.0,    0, -1.0,    0,    0],
            [  0,    0,    0,    0,    0],
            [  0,  1.0,    0, -1.0,    0]
        ])
        self._vec["C"] = np.array([0.0, 0.0, 0.0], dtype=float)



    @classmethod
    def from_config(cls, config: dict):
        """ Create block from the boundary condition config dictionary
            Returns the block instance"""
        # need to think about this a bit more - how is the mynard valve generated from the config file
        # will have to generate a custom config file to test the RH model generation
        pass

    def valve_resistance(self, q_in, O_v):
        """Compute resistance imposed by the valve orifice."""
        # return (0.5 * self.rho * abs(q_in) * q_in ) / ( (self.S_v ** 2) * (O_v ** 2 + self.eps) )
        return (self.beta * abs(q_in) * q_in ) / (O_v ** 2 + self.eps)

    def d_res_d_q(self, q_in, O_v, eps=1e-3):
        # Central difference
        dq = max(abs(q_in)*eps, eps**2)
        return ( self.valve_resistance(q_in + dq, O_v) - self.valve_resistance(q_in - dq, O_v) ) / (2. * dq)

    def d_res_d_Ov(self, q_in, O_v, eps=1e-3):
        # Forward difference
        dOv = max(abs(O_v)*eps, eps**2)
        return ( self.valve_resistance(q_in, O_v + dOv) - self.valve_resistance(q_in, O_v) ) / dOv

    def update_solution(self, y: np.ndarray, ydot: np.ndarray) -> None:
        """Update solution dependent element contributions.

        Args:
            y: Current solution.
            ydot: current solution time derivative
        """
        q_in = np.abs(y[self.inflow_nodes[0].flow_dof])
        q_in_dot = np.abs(ydot[self.inflow_nodes[0].flow_dof])
        p_in = np.abs(y[self.inflow_nodes[0].pres_dof])
        p_out = np.abs(y[self.outflow_nodes[0].pres_dof])

        O_v = y[self._flat_col_ids[4]] # opening/closing internal variable

        # calculate blood inertance and make E matrix
        L = self.rho * self.l_v / (self.S_v * (O_v + self.eps))  # inertance of blood
        d_L_dOv = -self.rho * self.l_v / (self.S_v * (O_v + self.eps) ** 2.0)  # dL/dzeta

        # Constants from linearization of Eqn (1)
        c_lin_res = -self.valve_resistance(q_in, O_v) + \
                    q_in * self.d_res_d_q(q_in, O_v) + \
                    O_v * self.d_res_d_Ov(q_in, O_v)

        # Eqn (2)
        delta = (1.0 - O_v) * self.K_vo * (p_in - p_out) if p_in >= p_out else (O_v * self.K_vc * (p_in - p_out))
        d_del_d_pin = (1.0 - O_v) * self.K_vo if p_in >= p_out else (O_v * self.K_vc)
        d_del_d_pout = -d_del_d_pin
        d_del_d_Ov = -self.K_vo * (p_in - p_out) if p_in >= p_out else (self.K_vc * (p_in - p_out))
        c_lin_state = delta - p_in * d_del_d_pin - p_out * d_del_d_pout - O_v * d_del_d_Ov


        # add terms to E matrix
        self._mat["E"][0, 1] = -L
        self._mat["E"][1, 4] = -1.0

        # add terms to F matrix
        self._mat["F"][0, 1] = -self.d_res_d_q(q_in, O_v)
        self._mat["F"][0, 4] = -self.d_res_d_Ov(q_in, O_v)
        self._mat["F"][1, 0] = d_del_d_pin
        self._mat["F"][1, 2] = d_del_d_pout
        self._mat["F"][1, 4] = d_del_d_Ov

        # add terms to cvec
        self._vec["C"][0] = c_lin_res
        self._vec["C"][1] = c_lin_state

        # add terms to dE matrix
        self._mat["dE"][0, 4] = -d_L_dOv * q_in_dot












