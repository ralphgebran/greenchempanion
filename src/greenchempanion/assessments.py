from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors
from typing import Dict
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import TanimotoSimilarity
from functions import Reaction


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

def atoms_assessment(react: Reaction) -> tuple[str,str]:
    """
    Warn if any risky elements are present in one of the products.
    Returns (message, color_hex).
    """
    # Make the lookup O(1) by storing the elements in a set
    RISKY_ATOMS = {
        "F", "Cl", "Br", "I", "Li", "Ti", "Sn", "Pb", "Pd",
        "Hg", "Cd", "As", "Cr", "Ni", "Se", "Tl", "Pt", "Rh",
    }

    risky_atoms = {
        atom.GetSymbol()
        for mol in react.products.keys()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() in RISKY_ATOMS
    }

    if risky_atoms:
        return f"⚠️ Concerning atoms: {', '.join(sorted(risky_atoms))}", "#F6DF7E"
    else:
        return "✅ All atoms are green", "#88DF66"


def structural_assessment(react: Reaction) -> tuple[str,str]:
    bad_smarts = {
        "Carbon oxides":      ["O=C=O", "C#O"],
        "Nitro-":             ["[N+](=O)[O-]", "O[N](=O)[O-]", "O=N[O-]", "[NX3+](=O)[O-,O]"],
        "Azo-":               ["N=N"],
        "Dichloro-aromatic":  ["c([Cl,Br])c.*c([Cl,Br])c"],
    }

    heavy_chain_flag = False
    bad_groups      = set()

    for mol in react.products.keys():
        if sum(1 for a in mol.GetAtoms() if a.GetSymbol() != "H") > 10:
            heavy_chain_flag = True

        mol_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

        for name, patterns in bad_smarts.items():
            for sm in patterns:
                ref = Chem.MolFromSmarts(sm)
                if ref is None:
                    continue

                # guard 1: skip any SMARTS that can't be sanitized
                try:
                    Chem.SanitizeMol(ref)
                except Exception:
                    continue

                if mol.HasSubstructMatch(ref):
                    bad_groups.add(name)
                    break
                try:
                    ref_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(ref, radius=2, nBits=2048)
                    sim    = TanimotoSimilarity(ref_fp, mol_fp)
                except Exception:
                    continue

                if sim > 0.15:
                    bad_groups.add(f"Similarity to {name}")
                    break

    issues = []
    if heavy_chain_flag:
        issues.append("Presence of long heavy-atom chain(s)")
    if bad_groups:
        issues.append("Problematic groups: " + ", ".join(sorted(bad_groups)))

    if issues:
        return "⚠️ " + "; ".join(issues), "#F03335"
    else:
        return "✅ No structural red flags", "#88DF66"
