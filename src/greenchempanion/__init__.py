"""Your interactive companion for Green Chemistry."""

from __future__ import annotations
from .functions import Atom_Count_With_H, Reaction, compute_PMI, canonicalize_smiles, compute_E , structural_assessment,atoms_assessment
from .assessments import get_solvent_info, waste_efficiency, PMI_assessment, Atom_ec_assessment, Atom_ec_m_assessment, logP_assessment_molecule

__version__ = "0.0.1"