
import numpy as np
import random
from scipy.optimize import minimize, Bounds, LinearConstraint
import math
# from structured_tree_tuning.struct_tree_utils import create_vesselDlist

class StructuredTreeOutlet():
    """Structured tree microvascular adaptation to upstream changes
    input: R, length of a vessel with an outlet BC
    output: structured tree of BloodVessels and Junctions

    need to restructure this class, as it is not a block but rather a collection of blocks

    """
    def __init__(self, params: dict = None, name: str = None, config: dict = None):
        """Create a new structured tree instance
        Args:
            params: The configuration paramaters of the block. Mostly comprised
                of constants for element contribution calculation.
            name: Optional name of the block.
        """
        self.params = params
        # initial diameter of the vessel from which the tree starts
        self.initialD = ((128 * self.params["eta"] * self.params["l"]) / (np.pi * self.params["R"])) ** (1 / 4)
        # intialize resistance
        self.totalResistance = 0
        # set up empty block dict if not generated from pre-existing tree
        if config is None:
            self.name = name
            self.block_dict = {'name': name, 'origin_d': self.initialD, 'vessels': [], 'junctions': [], 'adaptations': 0}
            self.vesselDlist = []
        else:
            self.name = config["name"]
            self.block_dict = config
            self.vesselDlist = create_vesselDlist(config["vessels"])

    @classmethod
    def from_outlet_vessel(cls, config: dict, simparams: dict, tree_config=False) -> "StructuredTreeOutlet":
        """Creates instance from config dictionary of outlet vessel
            Args:
                config file of outlet vessel
                config file of simulation parameters to get viscosity
                tree config file, if the tree already exists
            Returns:
                instance of structured tree
        """
        params = dict(
            # need vessel length to determine vessel diameter
            l=config.get("vessel_length"),
            R=config["zero_d_element_values"].get("R_poiseuille"),
            # Probably don't need C and L, just getting them for the sake of due diligence I guess
            C=config["zero_d_element_values"].get("C", 0.0),
            L=config["zero_d_element_values"].get("L", 0.0),
            stenosis_coefficient=config["zero_d_element_values"].get(
                "stenosis_coefficient", 0.0
            ),
            eta=simparams.get("viscosity"),
        )
        if tree_config:
            return cls(params=params, config = config["tree"])
        else:
            return cls(params=params, name="OutletTree" + str(config["vessel_id"]))

    def reset_block_dict(self):
        """
        reset the block dict if you are generating many iterations of the structured tree to optimize the radius
        Returns: empty block_dict

        """
        self.block_dict["vessels"] = []
        self.block_dict["junctions"] = []
        self.vesselDlist = []


    def build_next_layer(self, layer):
        """generates the exponents for the next layer of blood vessels given an input layer
           Args:
               input layer (list of tuples)
           Returns:
               alpha, beta exponents of next layer (list of tuples)"""
        next_layer = []
        for vessel in layer:
            next_layer.append((vessel[0] + 1, vessel[1]))
            next_layer.append((vessel[0], vessel[1] + 1))

        return next_layer

    def build_tree(self, initial_r=None, r_min=0.049, optimized=False, alpha=None, beta=None):
        """ iteratively generates a tree with BloodVessel and (Junction?) instances
            using pre-determined alpha and beta parameters from Olufsen et al. (alpha=0.9, beta=0.6)
            Args:
                initial vessel radius, minimum radius
            Returns:
                updated self.blockdict

        """
        if optimized:
            self.reset_block_dict() # reset block dict if making the tree many times

        if initial_r is None: # set default value of initial_r
            initial_r = self.initialD / 2

        if alpha is None and beta is None:
            # initialize with Olufsen parameters
            alpha = 0.9 # initialize alpha
            beta = 0.6 #initialize beta

        else:
            # print("the chosen alpha and beta are: " + str([alpha, beta]))
            alpha = alpha
            beta = beta

        # layers = random.randint(4, 7) # initialize a random number of bifurcation layers
        # print(r_min, initial_r)
        if initial_r <= 0:
            raise Exception("initial_r is invalid, " + str(initial_r))
        if beta ==0:
            raise Exception("beta is zero")
        layers = math.ceil(math.log(r_min / initial_r) / math.log(beta)) # number of layers required to reach r_min
        if layers <=1:
            layers = 2 # make sure there are at least two layers
        elif layers >= 6:
            layers = 6 # prevent creation of too many layers during optimization

        tree = [[(0, 0)]] #initialize the tree
        # print(initial_r, layers)
        for i in range(layers):
            next_layer = self.build_next_layer(tree[i])
            tree.append(next_layer)

        for layer in tree:
            vessel_layer = []
            for exponents in layer:
                vessel_layer.append((alpha ** exponents[0]) * (beta ** exponents[1]) * initial_r * 2)
            self.vesselDlist.append(vessel_layer)

        # create the blood vessel instances
        vessel_id = 0
        for i, layer in enumerate(self.vesselDlist):
            for j, vesselD in enumerate(layer):
                R, C, L, l = self.calc_zero_d_values(vesselD)
                name = " " # to implement later
                # generate essentially a config file for the BloodVessel instances
                vessel_info = {"vessel_id": vessel_id,    # mimic input json file
                               "vessel_length": l,
                               "vessel_D": vesselD,
                               "vessel_name": name,
                               "generation": i,
                               "zero_d_element_type": "BloodVessel",
                               "zero_d_element_values": {
                                   "R_poiseulle": R,
                                   "C": C,
                                   "L": L
                               }}
                self.block_dict["vessels"].append(vessel_info)

                # make a junction coming off of this vessel
                if i < len(self.vesselDlist[:-1]):
                    inlet_vessels = [vessel_id]
                    inlet_area = [np.pi * ((vesselD / 2) ** 2)] # in case this is needed for BloodVesselJunction
                    outlet_vessels = [2*vessel_id + 1,
                                      2*vessel_id + 2]
                    outlet_area = [np.pi * ((self.vesselDlist[i + 1][2 * j] / 2) ** 2), # in case this is needed for BloodVesselJunction
                                   np.pi * ((self.vesselDlist[i + 1][(2 * j) + 1] / 2) ** 2)]
                    # make dictionary for junction info
                    junction_info = {"inlet_vessels": inlet_vessels,
                                     "junction_name": "J" + str(j),
                                     "outlet_vessels": outlet_vessels}
                    self.block_dict["junctions"].append(junction_info)
                # track vessel id in the tree
                vessel_id += 1

    def calc_zero_d_values(self, vesselD):
        # calculate zero_d values based on an arbitrary vessel diameter
        r = vesselD / 2
        R = 8 * self.params["eta"] * self.params["l"] / (np.pi * r ** 4)
        C = 0  # to implement later
        L = 0  # to implement later
        l = 12.4 * r ** 1.1  # from ingrid's paper

        return R, C, L, l


    def update_zero_d_params(self):
        # update the zero_d parameters of the block based on a change in vessel diameter
        for vessel in self.block_dict["vessels"]:
            R, C, L, l = self.calc_zero_d_values(vessel.get("vessel_D"))
            vessel["zero_d_element_values"]["R_poiseulle"] = R
            vessel["zero_d_element_values"]["C"] = C
            vessel["zero_d_element_values"]["L"] = L
            vessel["vessel_length"] = l


    def calculate_resistance(self):
        R_mat = []
        for layer in self.vesselDlist:
            for vesselD in layer:
                # need to add in fahraeus-Lindqvist effect
                pass
            R_mat.append([8 * self.params["eta"] * self.params["l"] / (np.pi * (vesselD / 2) ** 4) for vesselD in layer])
        R_tot = []
        for i, layer in enumerate(reversed(R_mat[:-1])): # loop through all but the last layer in a reversed manner
            for j, r in enumerate(layer):

                r_add = 1 / ((1 / R_mat[-1 - i][2*j]) + (1 / R_mat[-1 - i][2 * j + 1]))
                layer[j] = r + r_add
            R_tot.append(layer)
        R_tot = R_tot[-1][0]
        ## need to figure out how to enumerate this correctly in some way and update the resistances for each vessel as you move up the tree
        # print(R_mat[-1])
        return R_tot


    def adapt_constant_wss(self, Q, Q_new, disp=False):
        # adapt the radius of the vessel based on the constant shear stress assumption
        R_old = self.calculate_resistance() # calculate pre-adaptation resistance

        for vessel in self.block_dict["vessels"]: # update the block_dict
            r = vessel.get("vessel_D") / 2
            r_new = (Q_new / Q)**(1/3) * r
            vessel["vessel_D"] = r_new * 2

        self.vesselDlist = create_vesselDlist(self.block_dict["vessels"]) # update the vesselDlist
        self.block_dict["adaptations"] += 1  # keep track of how many times the tree has been adapted
        self.update_zero_d_params()
        R_new = self.calculate_resistance() # calculate post-adaptation resistance

        if disp: # display change in resistance if necessary
            print("the change in resistance is "+ str(R_new - R_old))

        return R_new


    def adapt_pries_secomb(self):
        # adapt the tree based the pries and secomb model for diameter change

        pass


    def optimize_tree_radius(self, Resistance=5.0, log_file=None):
        """ use nelder-mead to optimize the radius of the tree vessels with respect to the desired resistance
            Args:
                desired tree resistance, radius guess
            Returns:
                optimized radius, total resistance
                """

        r_guess = self.initialD / 2

        def r_min_objective(radius):
            self.build_tree(radius[0], optimized=True)
            # print(self.block_dict)
            R = self.calculate_resistance()
            R_diff = (Resistance - R)**2
            return R_diff

        bounds = Bounds(lb=0.049, ub=self.initialD) # minimum is r_min
        r_final = minimize(r_min_objective,
                           r_guess,
                           options={"disp": False},
                           bounds=bounds) # Nelder mead doesn't seem to work here
        R_final = self.calculate_resistance()
        with open(log_file, "a") as log:
            log.write("     the optimized radius is " + str(r_final.x))
        return r_final.x, R_final

    def optimize_alpha_beta(self, Resistance=5.0, log_file=None):
        """ use constrained optimization to optimize the diameter and alpha and beta
            Args:
                desired tree resistance, radius guess
            Returns:
                optimized radius, total resistance
        """

        def r_min_objective(params):
            self.build_tree(params[0], optimized=True, alpha=params[1], beta=params[2])
            R = self.calculate_resistance()
            R_diff = (Resistance - R) ** 2
            return R_diff

        r_guess = self.initialD / 2
        params_guess = np.array([r_guess, 0.9, 0.6]) # r, alpha, beta initial guess
        param_constraints = LinearConstraint([[0, 0, 0], [0, 1, 1], [0, 1, -1.5]], [0.0, 1, 0], [np.inf, np.inf, 0])
        param_bounds = Bounds(lb=[0.049, 0, 0], ub=[np.inf, 1, 1], keep_feasible=True)
        r_final = minimize(r_min_objective,
                           params_guess,
                           options={"disp": True},
                           method='trust-constr',
                           constraints=param_constraints,
                           bounds=param_bounds)
        R_final = self.calculate_resistance()
        # write a log file of the optimization results
        with open(log_file, "a") as log:
            log.write("     Resistance after optimization is " + str(R_final) + "\n")
            log.write("     the optimized radius is " + str(r_final.x[0]) + "\n")
            log.write("     the optimized alpha value is " + str(r_final.x[1]) + "\n")
            log.write("     the optimized alpha value is " + str(r_final.x[2]) + "\n")

        return r_final.x[0], R_final

# method to generate vesselDlist if tree config already exists
class VesselD:
    # class to construct binary tree of vesselD if tree config already exists
    def __init__(self, val):
        self.val = val
        self.left = None
        self.right = None

def create_vesselDlist(vessels):
    # takes in the list of vessels in a tree config dictionary
    n_vessels = len(vessels)
    vesselDlist = []
    vesselDlayer = []
    root = VesselD(vessels[0].get("vessel_D"))
    current_layer = [root]
    next_layer = []
    i = 1
    while i < n_vessels:
        for vessel in current_layer:
            vesselDlayer.append(vessel.val)
            if i < n_vessels:
                vessel.left = VesselD(vessels[i].get("vessel_D"))
                next_layer.append(vessel.left)
            i += 1

            if i < n_vessels:
                vessel.right = VesselD(vessels[i].get("vessel_D"))
                next_layer.append(vessel.right)
            i += 1

        vesselDlist.append(vesselDlayer)
        vesselDlayer = []
        current_layer = next_layer
        next_layer = []

    return vesselDlist