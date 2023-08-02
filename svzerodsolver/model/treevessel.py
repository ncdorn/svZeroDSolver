import numpy as np

# class for recursive generation of a structured tree

class TreeVessel:
    # class to construct binary tree of blood vessels using recursion
    def __init__(self, info: dict, name: str = None):
        self.name = name # name of the, tree only has a str in the root node
        self.info = info # vessel info dict
        self.d = self.info["vessel_D"] # diameter
        self.l = self.info["vessel_length"]
        self._left = None
        self._right = None
        self.collapsed = False
        self.id = self.info["vessel_id"]
        self.gen = self.info["generation"]
        self.eta = self.info["viscosity"]
        self._R = self.info["zero_d_element_values"].get("R_poiseulle")
        self._R_eq = self._R # this is the initial value, not dependent on left and right



    @classmethod
    def create_vessel(cls, id, gen, diameter, eta):
        if diameter > .3:
            viscosity = eta
        else:
            viscosity = eta # implement fahraeus lindkvist later
        R, C, L, l = cls.calc_zero_d_values(cls, diameter, viscosity)
        # print(R, C, L, l)
        name = " "  # to implement later

        # generate essentially a config file for the BloodVessel instances
        vessel_info = {"vessel_id": id,  # mimic input json file
                       "vessel_length": l,
                       "vessel_D": diameter,
                       "vessel_name": name,
                       "generation": gen,
                       "viscosity": viscosity,
                       "zero_d_element_type": "BloodVessel",
                       "zero_d_element_values": {
                           "R_poiseulle": R,
                           "C": C,
                           "L": L
                       }}

        return cls(info=vessel_info)

    # property setters to dynamically update equivalent resistance
    @property
    def left(self):
        return self._left

    @left.setter # this method updates R_eq based on the left and right vessels
    def left(self, new_left):
        self._left = new_left
        if self.right is not None:
            self._update_R_eq()

    @property
    def right(self):
        return self._right

    @right.setter
    def right(self, new_right):
        self._right = new_right
        self._update_R_eq()

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, new_R):
        self._R = new_R
        if self.left is not None and self.right is not None:
            self._update_R_eq() # update R_eq if R is updated

    @property
    def R_eq(self):
        return self._R + (1 / (self._left._R_eq ** -1 + self._right._R_eq ** -1))

    def _update_R_eq(self):
        self._R_eq = self._R + (1 / (self._left._R_eq ** -1 + self._right._R_eq ** -1))



    def calc_zero_d_values(self, vesselD, eta):
        # calculate zero_d values based on an arbitrary vessel diameter
        r = vesselD / 2
        l = 12.4 * r ** 1.1  # from ingrid's paper, does this remain constant throughout adaptation?
        R = 8 * eta * l / (np.pi * r ** 4)
        C = 0  # to implement later
        L = 0  # to implement later

        return R, C, L, l

    def update_vessel_info(self):
        R, C, L, l = self.calc_zero_d_values(self.d, self.eta)
        self.info["vessel_length"] = l
        self.info["vessel_D"] = self.d
        self.info["zero_d_element_values"]["R_poiseulle"] = R
        self.info["zero_d_element_values"]["C"] = C
        self.info["zero_d_element_values"]["L"] = L
        self.R = R
        if not self.collapsed:
            self._update_R_eq()


