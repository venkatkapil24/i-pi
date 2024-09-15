""" An interface for the [MACE](https://github.com/ACEsuit/mace) calculator """

import sys
import json
from ase.io import read
import numpy as np
from .dummy import Dummy_driver
from ipi.utils.units import unit_to_internal, unit_to_user

try:
    from mace.calculators import MACECalculator
except Exception as e:
    MACECalculator = None
    print(f"Error loading mace bindings: {e}", file=sys.stderr)

__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"

ERROR_MSG = """
MACE driver requires specification of an ASE formatted initialization file (e.g. an extended xyz file) and a json formatted file with the keyword arguments for the MACE calculator.

Example: python driver.py -m mace -u -o template.xyz,mace_calculator.json

If the json file is not provided, the default file 'mace_calculator.json' will be checked for.
"""

class MACE_driver(Dummy_driver):
    def __init__(self, args=None, verbose=False):
        if MACECalculator is None:
            raise ImportError("Couldn't load mace bindings")

        super().__init__(args, verbose, ERROR_MSG)

    def check_arguments(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        # Check and parse self.template from self.args[0]
        if self.args[0] != '':
            self.template = self.args[0]
        else:
            raise ValueError("Error: No structure template provided." + '\n' + ERROR_MSG)

        self.template_ase = read(self.template)

        # Reads the JSON file (self.args[1] or default)
        try:
            # Check if the file path is provided
            if len(self.args) > 1:
                file_path = self.args[1]
            else:
                file_path = 'mace_calculator.json'
                print(f"No file provided, using default: {file_path}")

            # Try to open the file
            with open(file_path, "r") as file:
                calc_params = json.load(file)

        except FileNotFoundError:
            raise FileNotFoundError(f"Error: File {file_path} not found." + '\n' + ERROR_MSG)

        except json.JSONDecodeError as e:
            raise json.JSONDecodeError(f"Error: Failed to decode JSON from {file_path}." + '\n' + ERROR_MSG, pos=e.pos, doc=e.doc)

        # Proceed with initializing the calculator if the file was read successfully
        self.ase_calculator = MACECalculator(**calc_params)

    def __call__(self, cell, pos):
        """Get energies, forces, and stresses from the ASE calculator
        This routine assumes that the client will take positions
        in angstrom, and return energies in electronvolt, and forces
        in ev/ang.
        """

        # ASE calculators assume angstrom and eV units
        pos = unit_to_user("length", "angstrom", pos)
        # ASE expects cell-vectors-as-rows
        cell = unit_to_user("length", "angstrom", cell.T)
       
       # applies the cell and positions to the template
        structure = self.template_ase.copy()
        structure.positions[:] = pos
        structure.cell[:] = cell
        structure.calc = self.ase_calculator

        # Do the actual calculation
        self.ase_calculator.calculate(atoms=structure)
        properties = self.ase_calculator.results

        pot = properties["energy"]
        force = properties["forces"]
        stress = properties["stress"]
        if len(stress) == 6:
            # converts from voight notation
            stress = np.array(stress[[0, 5, 4, 5, 1, 3, 4, 3, 2]])

        # converts to internal quantities
        pot_ipi = np.asarray(
            unit_to_internal("energy", "electronvolt", pot), np.float64
        )
        force_ipi = np.asarray(unit_to_internal("force", "ev/ang", force), np.float64)
        vir_calc = -stress * structure.get_volume()
        vir_ipi = np.array(
            unit_to_internal("energy", "electronvolt", vir_calc.T), dtype=np.float64
        )
        extras = ""

        return pot_ipi, force_ipi, vir_ipi, extras

