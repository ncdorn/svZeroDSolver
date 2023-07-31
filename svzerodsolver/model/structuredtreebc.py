
import numpy as np
import random
from scipy.optimize import minimize, Bounds, LinearConstraint
from .treevessel import TreeVessel
import math

class StructuredTreeOutlet():
    """Structured tree microvascular adaptation to upstream changes
    input: R, length of a vessel with an outlet BC
    output: structured tree of BloodVessels and Junctions

    need to restructure this class, as it is not a block but rather a collection of blocks

    """
    def __init__(self, params: dict = None, name: str = None, config: dict = None, root: TreeVessel = None):
        """Create a new structured tree instance
        Args:
            params: The configuration paramaters of the block. Mostly comprised
                of constants for element contribution calculation.
            name: Optional name of the block.
        """
        self.params = params
        # initial diameter of the vessel from which the tree starts
        self.initialD = ((128 * self.params["eta"] * self.params["l"]) / (np.pi * self.params["R"])) ** (1 / 4)
        # set up empty block dict if not generated from pre-existing tree
        if config is None:
            self.name = name
            self.block_dict = {'name': name, 'origin_d': self.initialD, 'vessels': [], 'junctions': [], 'adaptations': 0}
            # initialize the root of the structured tree
            self.root = None
        else:
            self.name = config["name"]
            self.block_dict = config
            # initialize the root of the structured tree
            self.root = root




    @classmethod
    def from_outlet_vessel(cls, config: dict, simparams: dict, tree_exists=False, root: TreeVessel = None, ) -> "StructuredTreeOutlet":
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
        if tree_exists:
            return cls(params=params, config = config["tree"], root=root)
        else:
            return cls(params=params, name="OutletTree" + str(config["vessel_id"]))


    def reset_tree(self):
        """
        reset the block dict if you are generating many iterations of the structured tree to optimize the radius
        Returns: empty block_dict

        """
        self.root = None
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



    def build_tree(self, initial_r=None, d_min=0.049, optimizing=False, alpha=0.9, beta=0.6):

        if optimizing:
            self.reset_tree() # reset block dict if making the tree many times

        if initial_r is None: # set default value of initial_r
            initial_r = self.initialD / 2

        if d_min <= 0:
            raise ValueError("The min diameter must be greater than 0.")

        if initial_r <= 0:
            raise Exception("initial_r is invalid, " + str(initial_r))
        if beta ==0:
            raise Exception("beta is zero")
        # add r_min into the block dict
        self.block_dict["D_min"] = d_min
        initial_d = initial_r * 2
        # initialize counting values
        vessel_id = 0
        junc_id = 0
        # initialize the root vessel of the tree
        self.root = TreeVessel.create_vessel(0, 0, initial_d, self.params["eta"])
        self.root.name = self.name
        self.block_dict["vessels"].append(self.root.info)
        queue = [self.root]

        while len(queue) > 0:
            q_id = 0
            current_vessel = queue.pop(q_id)
            creating_vessels = True
            while current_vessel.collapsed:
                if len(queue) == 0:
                    creating_vessels = False
                    break
                current_vessel = queue.pop(q_id)
            if not creating_vessels:
                break

            if not current_vessel.collapsed:
                next_gen = current_vessel.gen + 1
                # create left vessel
                vessel_id += 1
                left_dia = alpha * current_vessel.d
                current_vessel.left = TreeVessel.create_vessel(vessel_id, next_gen, left_dia, self.params["eta"])
                # current_vessel.left.create_vessel_info(vessel_id,next_gen, left_dia)
                queue.append(current_vessel.left)
                self.block_dict["vessels"].append(current_vessel.left.info)
                if left_dia < d_min:
                    current_vessel.left.collapsed = True

                # create right vessel
                vessel_id += 1
                right_dia = beta * current_vessel.d
                current_vessel.right = TreeVessel.create_vessel(vessel_id, next_gen, right_dia, self.params["eta"])
                # current_vessel.right.info = self.create_vessel_info(vessel_id, next_gen, right_dia)
                queue.append(current_vessel.right)
                self.block_dict["vessels"].append(current_vessel.right.info)
                if right_dia < d_min:
                    current_vessel.right.collapsed = True

                # add a junction
                junction_info = {"inlet_vessels": [current_vessel.id],
                                 "junction_name": "J" + str(junc_id),
                                 "outlet_vessels": [current_vessel.left.id, current_vessel.left.id]}
                self.block_dict["junctions"].append(junction_info)
                junc_id += 1


    def update_zero_d_params(self):
        # update the zero_d parameters of the block based on a change in vessel diameter
        for vessel in self.block_dict["vessels"]:
            R, C, L, l = self.calc_zero_d_values(vessel.get("vessel_D"))
            vessel["zero_d_element_values"]["R_poiseulle"] = R
            vessel["zero_d_element_values"]["C"] = C
            vessel["zero_d_element_values"]["L"] = L
            vessel["vessel_length"] = l


    def adapt_constant_wss(self, Q, Q_new, disp=False):
        R_old = self.root.R_eq  # calculate pre-adaptation resistance

        def constant_wss(d, Q=Q, Q_new=Q_new):
            # adapt the radius of the vessel based on the constant shear stress assumption
            return (Q_new / Q) ** (1 / 3) * d

        def update_diameter(vessel, update_func):
            # preorder traversal to update the diameters of all the vessels in the tree
            if vessel:

                vessel.d = update_func(vessel.d)
                vessel.update_vessel_info()
                update_diameter(vessel.left, update_func)
                update_diameter(vessel.right, update_func)

        update_diameter(self.root, constant_wss)

        # self.root.update_vessel_info()

        R_new = self.root.R_eq
        if disp: # display change in resistance if necessary
            print("R_new = " + str(R_new) + ", R_old = " + str(R_old) + "\n")
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
            self.build_tree(radius[0], optimizing=True)
            R = self.root.R_eq
            R_diff = (Resistance - R)**2
            return R_diff

        bounds = Bounds(lb=0.005) # minimum is r_min
        r_final = minimize(r_min_objective,
                           r_guess,
                           options={"disp": True},
                           method='Nelder-Mead',
                           bounds=bounds) # Nelder mead doesn't seem to work here
        R_final = self.root.R_eq
        if log_file is not None:
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
            self.build_tree(params[0], optimizing=True, alpha=params[1], beta=params[2])
            R = self.root.R_eq
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
        R_final = self.root.R_eq
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
