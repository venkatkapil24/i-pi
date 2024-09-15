""" An interface for the [MACE](https://github.com/ACEsuit/mace) calculator """

import sys
from .mace import MACE_driver
from ipi.utils.units import unit_to_internal, unit_to_user

try:
    from mace.calculators import MACECalculator
except:
    MACECalculator = None

__DRIVER_NAME__ = "mace-dielectric"
__DRIVER_CLASS__ = "MACEDielectric_driver"

ERROR_MSG = """
MACEDielectric driver requires specification of an ASE formatted initialization file (e.g. an extended xyz file), a json formatted file with the keyword arguments for the MACE calculator, and a Boolean specifying if the derivatives of the dielectric response should be calculated.

Example: python driver.py -m mace -u -o template.xyz,mace_calculator.json,True

If the json file is not provided, the default file 'mace_calculator.json' will be checked for. A default value of False will be assumed if the Boolean is not specified.
"""


class MACEDielectric_driver(MACE_driver):

    def check_arguments(self):
        """
        Check the arguments requuired to run the driver
        This loads the potential and atoms template in MACE
        """

        super().check_arguments()

        # Check and parse self.derivative from self.args[2] (if available)
        if len(self.args) > 2:
            if self.args[2] == "True":
                self.derivative = True
            elif self.args[2] == "False":
                self.derivative = False
            else:
                raise ValueError(f"Invalid value for derivative: {self.args[2]}. Expected 'True' or 'False'." + '\n' + ERROR_MSG)
        else:
            self.derivative = False


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

        # Do the actual calculation
        #properties = structure.get_properties(["energy", "forces", "stress"])
        #pot = properties["energy"]
        #force = properties["forces"]
        #stress = properties["stress"]
        #if len(stress) == 6:
            # converts from voight notation
        #    stress = np.array(stress[[0, 5, 4, 5, 1, 3, 4, 3, 2]])

        # converts to internal quantities
        #pot_ipi = np.asarray(
        #    unit_to_internal("energy", "electronvolt", pot), np.float64
        #)
        #force_ipi = np.asarray(unit_to_internal("force", "ev/ang", force), np.float64)
        #vir_calc = -stress * structure.get_volume()
        #vir_ipi = np.array(
        #    unit_to_internal("energy", "electronvolt", vir_calc.T), dtype=np.float64
        #)
        extras = {}

        pot_ipi = 0
        force_ipi = pos * 0
        vir_ipi = cell * 0

        self.ase_calculator.calculate(atoms=structure)
        extras = self.ase_calculator.results

        if self.derivative:
            r = self.ase_calculator.get_dielectric_derivatives(atoms=structure)
            extras['dipole_der'] = r[0]
            extras['polarizability_der'] = r[1]

        extras = str(extras)

        return pot_ipi, force_ipi, vir_ipi, extras
