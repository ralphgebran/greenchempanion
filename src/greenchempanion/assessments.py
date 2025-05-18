from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors
from typing import Dict
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import TanimotoSimilarity


def get_solvent_info(extras: Dict[Chem.Mol, float]) -> tuple[str, str]:
    """
    Given extras={Mol: mass_in_g_per_kg_product}, returns a short verdict + a color hex code (red/yellow/green).
    """
    Green      = {"O", "CCO", "CC(=O)OCC", "CC1COCC1", "O=C=O", "CC(O)C", "CO"}
    Bad        = {"ClCCl", "ClC(Cl)Cl", "c1ccccc1", "ClC(Cl)(Cl)Cl", "CCCCCC", "CCCCC"}
    
    seen = {"Green": 0, "Acceptable": 0, "Bad": 0}
    for mol in extras:
        smi = Chem.MolToSmiles(mol)
        if smi in Green:
            seen["Green"] += 1
        elif smi in Bad:
            seen["Bad"] += 1
        else:
            seen["Acceptable"] += 1

    if seen["Bad"] > 0:
        return "⚠️ Bad solvent detected!", "#F03335"
    elif seen["Acceptable"] > 0:
        return "⚠️ Adequate solvents detected.", "#F6DF7E"
    else:
        return "✅ All solvents are Green!", "#88DF66"


def waste_efficiency(E : float) -> tuple[str,str] :
    """ Returns if the E factor is within normal ranges, along with a color hex code """
    if E <= 1:
        return "Great Waste Efficiency, with stellar E factor ✅", "#4BAF24"
    elif 1 < E <= 10:
        return "Good Waste Efficiency ✅","#88DF66"
    elif 10 < E <= 50:
        return "Average Waste Efficiency, could be improved ⚠️","#F6DF7E"
    elif 50 < E <= 100:
        return "Bad Waste Efficiency, should be improved 🚨","#F68B8C"
    elif E > 100:
        return "🚨 Very Bad Waste Efficiency, must be improved 🚨","#F03335"
    
    
def PMI_assesment(PMI : float) -> tuple[str,str] :
    """ Returns if the PMI factor is within normal ranges, along with a color hex code """
    if  PMI <= 10:
        return "Great Process Mass Intensity ✅","#88DF66"
    elif 10 < PMI <= 50:
        return "Average Process Mass Intensity, could be improved ⚠️","#F6DF7E"
    elif 50 < PMI <= 100:
        return "Bad Process Mass Intensity, should be improved 🚨","#F68B8C"
    elif PMI > 100:
        return "🚨 Very Bad Process Mass Intensity, must be improved 🚨","#F03335"
    
    
def Atom_ec_assesment(ae: float) -> tuple[str, str]:
    """
    Returns a qualitative assessment of Atom Economy (atom-count based),
    Returns a short explanation and a color code.
    """
    if 89 < ae <= 100:
        return "✅ Excellent — most atoms end up in the product", "#4BAF24"
    elif 79 < ae <= 89:
        return "✅ Very good — minimal atoms lost as by-products", "#88DF66"
    elif 59 < ae <= 79:
        return "⚠️ Moderate — some atoms are lost as by-products", "#F6DF7E"
    elif 39 < ae <= 59:
        return "🚨 Poor — many atoms wasted as by-products", "#F68B8C"
    elif ae <= 39:
        return "🚨 Very poor — significant atom wastage", "#F03335"
    elif ae > 100:
        return "❓ Check input — AE > 100% not possible", "#AAAAAA"


def Atom_ec_m_assesment(ae: float) -> tuple[str, str]:
    """
    Returns a qualitative assessment of Atom Economy (mass-based),
    returns a short explanation and a color code.
    """
    if 89 < ae <= 100:
        return "✅ Excellent incorporation of reactant mass — near-zero waste", "#4BAF24"
    elif 79 < ae <= 89:
        return "✅ Very good incorporation of reactant mass — minimal waste", "#88DF66"
    elif 59 < ae <= 79:
        return "⚠️ Moderate incorporation of reactant mass — room for improvement", "#F6DF7E"
    elif 39 < ae <= 59:
        return "🚨 Poor incorporation of reactant mass — significant waste", "#F68B8C"
    elif ae <= 39:
        return "🚨 Very poor incorporation of reactant mass — substantial loss", "#F03335"
    elif ae > 100:
        return "❓ Double-check your inputs — AE > 100% is not physically possible", "#AAAAAA"

  
def logP_assessment_molecule(lo: float) -> tuple[str,str]:
    """Environmental friendliness based on MolLogP (basis as ideal is 2)."""
    if 1.5<= lo <= 2.5:
        return f"Excellent Log P of {lo:.2f}", "#4BAF24"
    if 1 <= lo < 4:
        return f"✅ LogP of {lo:.2f}, within acceptable range", "#88DF66"
    elif 0 <= lo < 1 or 4 <= lo <= 4.5:
        return f"⚠️ LogP of {lo:.2f}, caution", "#F6DF7E"
    else:
        return f"🚨 LogP of {lo:.2f}, outside of acceptable limits, high concern", "#F03335"

