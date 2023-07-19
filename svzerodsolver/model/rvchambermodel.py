import numpy as np
from .block import Block

class RV_Chamber_Model(Block):
    """
        Model of the right ventricle chamber
        Chamber Capacitance -- with direct prescription of pressure
    """
    _NUM_EQUATIONS = 3
    _NUM_INTERNAL_VARS = 1 # internal variable is volume

    def __init__(self, params: dict = None, name: str = None, Tc = None, v0 = None,
                 RH_params: dict = None):
        """ create a new Mynard valve instance

                    Args:
                        Configuration parameters of the block. mostly constants.
                        theta, c1, c2, c3, c4 are the tunable constants
                    Name:
                        name of the instance
                """
        super().__init__(params=params, name=name)

        self.Tc = Tc
        # RH_params0 = np.array([Ts, Tr, theta, c1, c2, c3, c4])
        self.Ts = RH_params["Ts"]
        self.Tr = RH_params["Tr"]
        self.theta = RH_params["theta"]
        self.c1 = RH_params["c1"]
        self.c2 = RH_params["c2"]
        self.c3 = RH_params["c3"]
        self.c4 = RH_params["c4"]
        self.V_rv0 = v0
        self.CONVERSION = 1333.2
        self.count = 0

        # initialize E matrix
        self._mat["E"] = np.array([
            [0, 0, 0, 0, 1.0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0]
        ])
        # initialize F matrix
        self._mat["F"] = np.array([
            [   0, -1.0,    0,  1.0,    0],
            [-1.0,    0,  1.0,    0,    0],
            [ 1.0,    0,    0,    0,    0]
        ])
        # initialize C vec
        self._vec["C"] = np.array([0, 0, 0])
        # initalize dC mat
        self._mat["dC"] = np.array([
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ])
        self.dcmxcoe = np.array([(0,) * 5,
                                 (0,) * 5,
                                 (0, 0, 0, 0, 0)])


    @classmethod
    def from_config(cls, config: dict):
        """ Create block from the boundary condition config dictionary
            Returns the block instance"""
        # need to think about this a bit more - how is the mynard valve generated from the config file
        # will have to generate a custom config file to test the RH model generation
        pass

    def linearize_func(self, v, func, eps_dv=1e-4):
        # Input  : func(t,v), typically t = t_n+1 ie t_n + dt
        # Return : d/dv (func) at (t,v) using central differences

        dv = eps_dv * v
        return (func(v + dv) - func(v - dv)) / (2. * dv)

    def linearize_func_t(self, t, v, func, dt=1e-3):
        """
        Args:
            t: time
            v: volume
            func: func(t,v)
            dt: timestep

        Returns: linearized function

        """
        return (func(t + dt, v) - func(t - dt, v)) / (2. * dt)

    def initialize_volume(self, sol_vec, scale_fact=1.5):
        """
        initialize volume variable
        Args:
            sol_vec: y
            Vu: initial volume
            scale_fact: scaling factor

        Returns:
            initial volume value
        """
        sol_vec[self._flat_col_ids[4]] = scale_fact * self.V_rv0

    def update_time(self, time: float) -> None:
        """
        update time dependent quantities
        Args:
            time: time vector

        Returns:

        """

        t_m = time % self.Tc
        k = 2.0 / (np.pi + 2.0)

        if t_m < (self.Ts / 2.0):
            self.a = 2.0 * (1.0 - k) * t_m / self.Ts
        elif t_m < self.Ts:
            self.a = 1.0 - k + k * np.sin(np.pi * (t_m / self.Ts - 0.5))
        elif t_m < (self.Ts + self.Tr):
            self.a = 1.0 - np.sin(np.pi * (t_m - self.Ts) / (2.0 * self.Tr))
        else:
            self.a = 0.0


    def update_solution(self, y: np.ndarray, ydot: np.ndarray) -> None:
        if self.count == 0:
            y[self._flat_col_ids[4]] = self.V_rv0

        V = y[self._flat_col_ids[4]] # not sure where I got this indexing from

        if V <= 0:
            print('Zero chamber volume detected in chamber ' + self.name + '\n')


        active_fcn = lambda V: self.CONVERSION * (self.a * self.theta * (self.c1 * V + self.c2))
        passive_fcn = lambda V: self.CONVERSION * ((1.0 - self.a) * (self.c3 * V + self.c4))

        # linearize above equations
        a_lin = self.linearize_func(V, active_fcn)
        p_lin = self.linearize_func(V, passive_fcn)

        # update c vec
        self._vec["C"][2] = -active_fcn(V) - passive_fcn(V)


        # update dC mat
        self._mat["dC"][2, 4] = -a_lin - p_lin

        self.count += 1 # bootleg initialization count





