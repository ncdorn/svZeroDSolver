

import numpy as np
from .block import Block


class MPA_BloodVesselJunction(Block):
    """MPA Blood vessel junction.

    to accurately calculate flow split at teh LPA/RPA split
    Attributes:
        name: Name of the block.
        inflow_nodes: Inflow nodes of the element.
        outflow_nodes: Outflow nodes of the element.
    """

    def __init__(self, params: dict = None, name: str = None, RPA_flow_split=None):

        super().__init__(params=params, name=name)

        self.RPA_flow = RPA_flow_split
        self.LPA_flow = 1 - RPA_flow_split

    def setup_dofs(self, dofhandler):
        """Setup degree of freedoms of the block.

        Registers equations and internal variables at a DOF handler.

        Args:
            dofhandler: The DOF handler to register the variables and equations
                at.
        """
        # Derive number of inlets and outlets
        num_inlets = len(self.inflow_nodes)
        num_outlets = len(self.outflow_nodes)

        # Set number of equations of a junction block based on number of
        # inlets/outlets. Must be set before calling parent constructor
        self._NUM_EQUATIONS = num_inlets + num_outlets
        super().setup_dofs(dofhandler)

        # Set some constant element element contributions that needed
        # _NUM_EQUATIONS
        self._mat["F"] = np.zeros(
            (self._NUM_EQUATIONS, 6)
        )
        # pressure continuity
        for i in range(2):
            self._mat["F"][i, [0, 2 * i + 2]] = [1.0, -1.0]
        # flow continuity
        self._mat["F"][2, 1] = 1.0
        self._mat["F"][2, 3] = -1.0
        self._mat["F"][2, 5] = -1.0
        # RPA flow split
        # self._mat["F"][3, 1] = self.RPA_flow / (1 - self.RPA_flow)
        # self._mat["F"][3, 3] = -1.0

        # print("this is the F matrix for " + str(self.name) + ": " + str(self._mat["F"]))

    @classmethod
    def from_config(cls, config: dict) -> "Block":
        """Create block from config dictionary.

        Args:
            config: The configuration dict for the block.

        Returns:
            The block instance.
        """
        name = config["junction_name"]

        if not name.startswith("J") and not name[1].isnumeric():
            raise ValueError(
                f"Invalid junction name {name}. Junction names must "
                "start with J following by a number."
            )

        # need to add functionality to add a custom flow split value later
        return cls(name=name, RPA_flow_split=0.52)