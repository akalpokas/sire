__all__ = ["emle"]

from ..legacy import Convert as _Convert

from .. import use_new_api as _use_new_api

_use_new_api()


class EMLEEngine(_Convert._SireOpenMM.PyQMEngine):
    """A class to enable use of EMLE as a QM engine."""

    def get_forces(self):
        """
        Get the OpenMM forces for this engine. The first force is the actual
        EMLE force, which uses a CustomCPPForceImpl to calculate the electrostatic
        embedding force. The second is a null CustomBondForce that can be used to
        add a "lambda_emle" global parameter to a context to allow the force to be
        scaled.

        Returns
        -------

        emle_force : openmm.Force
            The EMLE force object to compute the electrostatic embedding force.

        interpolation_force : openmm.CustomBondForce
            A null CustomBondForce object that can be used to add a "lambda_emle"
            global parameter to an OpenMM context. This allows the electrostatic
            embedding force to be scaled.
        """

        from copy import deepcopy as _deepcopy
        from openmm import CustomBondForce as _CustomBondForce

        # Create a dynamics object for the QM region.
        d = self._mols["property is_perturbable"].dynamics(
            timestep="1fs",
            constraint="none",
            platform="cpu",
            qm_engine=self,
        )

        # Get the OpenMM EMLE force.
        emle_force = _deepcopy(d._d._omm_mols.getSystem().getForce(0))

        # Create a null CustomBondForce to add the EMLE interpolation
        # parameter.
        interpolation_force = _CustomBondForce("")
        interpolation_force.addGlobalParameter("lambda_emle", 1.0)

        # Return the forces.
        return emle_force, interpolation_force


def emle(
    mols,
    qm_atoms,
    calculator,
    cutoff="7.5A",
    neighbourlist_update_frequency=20,
    redistribute_charge=False,
    map=None,
):
    """
    Create an EMLE engine object to allow QM/MM simulations using sire.mol.dynamics.

    Parameters
    ----------

    mols : sire.system.System
        The molecular system.

    qm_atoms : str, int, list, molecule view/collection etc.
        Any valid search string, atom index, list of atom indicies,
        or molecule view/container that can be used to select
        qm_atoms from 'mols'.

    calculator : emle.calculator.EMLECalculator
        The EMLECalculator object to use for elecotrostatic embedding calculations.

    cutoff : str or sire.legacy.Units.GeneralUnit, optional, default="7.5A"
        The cutoff to use for the QM/MM calculation.

    neighbourlist_update_frequency : int, optional, default=20
        The frequency with which to update the neighbourlist.

    redistribute_charge : bool
        Whether to redistribute charge of the QM atoms to ensure that the total
        charge of the QM region is an integer. Excess charge is redistributed
        over the non QM atoms within the residues involved in the QM region.

    Returns
    -------

    engine : sire.qm.EMLEEngine
        The EMLE engine object.
    """

    try:
        from emle.calculator import EMLECalculator as _EMLECalculator
    except:
        raise ImportError(
            "Could not import emle. Please install emle-engine and try again."
        )

    from ..base import create_map as _create_map
    from ..mol import selection_to_atoms as _selection_to_atoms
    from ..system import System as _System
    from ..legacy import Units as _Units
    from ..units import angstrom as _angstrom
    from .. import u as _u

    if not isinstance(mols, _System):
        raise TypeError("mols must be a of type 'sire.System'")

    # Clone the system.
    mols = mols.clone()

    try:
        qm_atoms = _selection_to_atoms(mols, qm_atoms)
    except:
        raise ValueError("Unable to select 'qm_atoms' from 'mols'")

    if not isinstance(calculator, _EMLECalculator):
        raise TypeError(
            "'calculator' must be a of type 'emle.calculator.EMLECalculator'"
        )

    if not isinstance(cutoff, (str, _Units.GeneralUnit)):
        raise TypeError(
            "cutoff must be of type 'str' or 'sire.legacy.Units.GeneralUnit'"
        )

    if isinstance(cutoff, str):
        try:
            cutoff = _u(cutoff)
        except:
            raise ValueError("Unable to parse cutoff as a GeneralUnit")

    if not cutoff.has_same_units(_angstrom):
        raise ValueError("'cutoff' must be in units of length")

    if not isinstance(neighbourlist_update_frequency, int):
        raise TypeError("'neighbourlist_update_frequency' must be of type 'int'")

    if neighbourlist_update_frequency < 0:
        raise ValueError("'neighbourlist_update_frequency' must be >= 0")

    if not isinstance(redistribute_charge, bool):
        raise TypeError("'redistribute_charge' must be of type 'bool'")

    if map is not None:
        if not isinstance(map, dict):
            raise TypeError("'map' must be of type 'dict'")
    map = _create_map(map)

    # Create the EMLE engine.
    engine = EMLEEngine(
        calculator,
        "_sire_callback",
        cutoff,
        neighbourlist_update_frequency,
    )

    from ._utils import (
        _check_charge,
        _create_qm_mol_to_atoms,
        _configure_engine,
        _create_merged_mols,
        _get_link_atoms,
    )

    # Check that the charge of the QM region is integer valued.
    _check_charge(mols, qm_atoms, map, redistribute_charge)

    # Get the mapping between molecule numbers and QM atoms.
    qm_mol_to_atoms = _create_qm_mol_to_atoms(qm_atoms)

    # Get link atom information.
    mm1_to_qm, mm1_to_mm2, bond_scale_factors, mm1_indices = _get_link_atoms(
        mols, qm_mol_to_atoms, map
    )

    # Configure the engine.
    engine = _configure_engine(
        engine, mols, qm_atoms, mm1_to_qm, mm1_to_mm2, bond_scale_factors, map
    )

    # Create the merged molecule.
    qm_mols = _create_merged_mols(qm_mol_to_atoms, mm1_indices, map)

    # Update the molecule in the system.
    mols.update(qm_mols)

    # Bind the system as a private attribute of the engine.
    engine._mols = mols

    return mols, engine
