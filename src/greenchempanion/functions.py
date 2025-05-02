from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors
from typing import Dict
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import TanimotoSimilarity

#Atom count function that takes hydrogen atoms into account

def Atom_Count_With_H(mol: Mol) -> int:
    """
    Returns the total number of atoms for a given molecule, counting the H atoms
    
    Parameters:
    - mol: RDKit Mol object
    
    Returns:
    - int: Number of atoms
    """
    
    H_mol = Chem.AddHs(mol)
    return H_mol.GetNumAtoms()


# Canonicalize Smiles

def canonicalize_smiles(smiles : str) -> str:
    """
    Returns the canonicalized SMILES for a given SMILES string.

    Parameters:
    - smile : A valid SMILES string representing a molecule.

    Returns:
    - str: The canonicalized SMILES string, or None if input is invalid
    """
    
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol)

# Reaction Class

class Reaction:
    """
    A class to represent a chemical reaction.

    Attributes:
    - reactants (Dict[Mol, int]): Dictionary of RDKit Mol objects with their Stoichiometric Coefficients, representing the reactants.
    - products (Dict[Mol, int]): Dictionary of RDKit Mol objects with their Stoichiometric Coefficients, representing the products.
    - main_product_index (int): The index of the main product, in the products list (Set as 0 by default)

    Methods:
    - Atom_Economy_A(): Calculates and returns the atom economy of the reaction, based on numbers of atoms.
    - Atom_Economy_M(): Calculates and returns the atom economy of the reaction, based on molar masses.
    """
    
    def __init__(self, reactants: Dict[Mol, int], products: Dict[Mol, int], main_product_index: int = 0):
        self.reactants = reactants
        self.products = products
        self.main_product_index = main_product_index

        # Error handling
        if not reactants:
            raise ValueError("Reactants list cannot be empty.")
        if not products:
            raise ValueError("Products list cannot be empty.")
        if main_product_index < 0 or main_product_index >= len(products):
            raise IndexError("main_product_index is out of range.")
        
        self.main_product = list(products.keys())[main_product_index]

    # Atom Economy function using number of atoms
    def Atom_Economy_A(self) -> float:
        """
        Calculates the atom economy of the reaction, based on numbers of atoms.

        (Number of atoms in main product / Number of atoms in all reactants) * 100

        Returns:
        - float: Atom economy as a percentage (%)
        """
        reactants_atoms = sum(Atom_Count_With_H(mol) * coeff for mol,coeff in self.reactants.items())
        main_product_atoms = Atom_Count_With_H(self.main_product) * self.products[self.main_product]
        economy = (main_product_atoms / reactants_atoms) * 100
        if economy > 100 :
            raise ValueError("Atom economy cannot be greater than 100, recheck your reaction")
        else :
            return economy
    
    # Atom Economy function using molar masses
    def Atom_Economy_M(self) -> float:
        """
        Calculates the atom economy of the reaction, based on molar masses.

        (Molar mass of main product / Total molar mass of reactants) * 100

        Returns:
        - float: Atom economy as a percentage (%)
        """
        reactants_mass = sum(Descriptors.MolWt(mol) * coeff for mol,coeff in self.reactants.items())
        main_product_mass = Descriptors.MolWt(self.main_product) * self.products[self.main_product]
        return (main_product_mass / reactants_mass) * 100
        
    # Reaction Balance check
    def isBalanced(self) -> bool:
        """
        Checks whether the reaction is balanced or not.

        Returns:
        - bool: True if the reaction is balanced, False if it's not.
        """
        def atom_count_dict(mol_dict: Dict[Mol, int]) -> Dict:
            atom_counts = {}
            for mol, coeff in mol_dict.items():
                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    if symbol not in atom_counts:
                        atom_counts[symbol] = 0
                    atom_counts[symbol] += coeff
            return atom_counts
        
        reac_count = atom_count_dict(self.reactants)
        prod_count = atom_count_dict(self.products)
        return reac_count == prod_count
    


# PMI and E-Factor Functions

def compute_PMI(reaction: Reaction, extras: Dict[Mol, float], prod_yield: float) -> float:
    """
    Computes the Process Mass Intensity (PMI) per kg of main product.
    
    (Total mass of chemical inputs / 1kg of product) * 100

    Inputs:
    - reaction: GCP Reaction object (with reactants, products, and main product index)
    - extras: Dictionary containing solvents, extraction material and others, with their mass [in grams] per kg of product {Mol: float}
    - prod_yield: Yield of main product (0 < float <= 1)
    
    Returns:
    - PMI (float): total input mass (g) per 1 kg of product
    """
    from rdkit.Chem import Descriptors

    # Error Handling
    if not isinstance(reaction, Reaction):
        raise TypeError("Expected a GCP Reaction object for 'reaction'")
    if not isinstance(extras, dict):
        raise TypeError("Expected a dictionary for 'extras'")
    if not all(isinstance(k, Mol) and isinstance(v, (int, float)) for k, v in extras.items()):
        raise TypeError("Keys in 'extras' must be RDKit Mol objects and values must be numbers")
    if not isinstance(prod_yield, (int, float)) or not (0 < prod_yield <= 1):
        raise ValueError("Yield must be a float between 0 (exclusive) and 1 (inclusive)")

    # Main Product Molar Mass
    product_mols = list(reaction.products.keys())
    main_product = product_mols[reaction.main_product_index]
    mp_coeff = reaction.products[main_product]
    main_product_mw = Descriptors.MolWt(main_product)

    # Moles for 1kg of product
    mass_target = 1000.0  # g
    mols_product_needed = mass_target / (prod_yield * main_product_mw)

    # Reactants masses
    total_input_mass = 0.0
    for mol, coeff in reaction.reactants.items():
        mw = Descriptors.MolWt(mol)
        total_input_mass += coeff * mw * mols_product_needed/mp_coeff

    # Masses from Extras
    total_input_mass += sum(extras.values())

    # PMI Output
    PMI = total_input_mass / mass_target
    return PMI


def compute_E(reaction: Reaction, extras: Dict[Mol, float], prod_yield: float) -> float:
    """
    Computes the Environmental Factor (E-Factor).

    E-Factor = Waste per kilogram of product

    Inputs:
    - reaction: GCP Reaction object
    - extras: Dictionary {Mol: float} for solvents, etc.
    - prod_yield: Yield of main product (0 < float <= 1)
    
    Returns:
    - E-Factor (float)
    """
    from rdkit.Chem import Descriptors

    # Error Handling
    if not isinstance(reaction, Reaction):
        raise TypeError("Expected a GCP Reaction object for 'reaction'")
    if not isinstance(extras, dict):
        raise TypeError("Expected a dictionary for 'extras'")
    if not all(isinstance(k, Mol) and isinstance(v, (int, float)) for k, v in extras.items()):
        raise TypeError("Keys in 'extras' must be RDKit Mol objects and values must be numbers")
    if not isinstance(prod_yield, (int, float)) or not (0 < prod_yield <= 1):
        raise ValueError("Yield must be a float between 0 (exclusive) and 1 (inclusive)")

    # Main Product Molar Mass
    product_mols = list(reaction.products.keys())
    main_product = product_mols[reaction.main_product_index]
    mp_coeff = reaction.products[main_product]
    main_product_mw = Descriptors.MolWt(main_product)

    # Moles for 1kg of product
    mass_target = 1000.0  # g
    mols_product_needed = mass_target / (prod_yield * main_product_mw)

    # Side Products masses
    total_waste_mass = 0.0
    for mol, coeff in reaction.products.items():
        if mol == main_product :
            continue
        mw = Descriptors.MolWt(mol)
        total_waste_mass += coeff * mw * mols_product_needed/mp_coeff

    # Masses from Extras
    total_waste_mass += sum(extras.values())

    # PMI Output
    E = total_waste_mass / mass_target
    return E

def get_solvent_info(extras: Dict[Chem.Mol, float]) -> tuple[str, str]:
    """
    Given extras={Mol: mass_in_g_per_kg_product}, returns
    a short verdict + a color hex code (red/yellow/green).
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

    # verdict logic
    if seen["Bad"] > 0:
        return "⚠️ Bad solvent detected!", "#F03335"
    elif seen["Acceptable"] > 0:
        return "⚠️ Adequate solvents detected.", "#F6DF7E"
    else:
        return "✅ All solvents are Green!", "#88DF66"


def waste_efficiency(E : float) -> tuple[str,str] :
    """ 
    Returns if the E factor is within normal ranges, along with a color hex code """
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
    """ 
    Returns if the PMI factor is within normal ranges, along with a color hex code 
    """
    if  PMI <= 10:
        return "Great Process Mass Intensity ✅","#88DF66"
    elif 10 < PMI <= 50:
        return "Average Process Mass Intensity, could be improved ⚠️","#F6DF7E"
    elif 50 < PMI <= 100:
        return "Bad Process Mass Intensity, should be improved 🚨","#F68B8C"
    elif PMI > 100:
        return "🚨 Very Bad Process Mass Intensity, must be improved 🚨","#F03335"


def Atom_ec_assesment(ae : float)  -> tuple[str,str] :
    """ 
    Returns if the atom economy is within normal ranges, along with a color hex code 
    """
    if 89 < ae <= 100:
        return "Amazing atom economy ✅", "#4BAF24"
    elif  79 < ae <= 89:
        return "Great atom economy ✅","#88DF66"
    elif 59 < ae <= 79:
        return " ⚠️ Could be improved, Average atom economy ","#F6DF7E"
    elif 39 < ae <= 59:
        return " 🚨 Should be improved, Bad Atom economy","#F68B8C"
    elif ae <= 39:
        return "🚨 Must be improved, Very Bad Atom economy","#F03335"
    elif ae >100:
        return " Doublecheck your reaction, Impossible Atom economy "
    
  
    
def logP_assessment_molecule(lo: float) -> tuple[str,str]:
    """Environmental friendliness based on MolLogP (basis as ideal is 2)."""
    if 1.5<= lo <= 2.5:
        return f"Excellent Log P of {lo:.2f}", "#4BAF24"
    if 0 <= lo <= 4:
        return f"✅ LogP of {lo:.2f}, within acceptable range", "#88DF66"
    elif -1 <= lo < 0 or 4 < lo <= 6:
        return f"⚠️ LogP of {lo:.2f}, slightly outside ideal", "#F6DF7E"
    else:
        return f"🚨 LogP of {lo:.2f}, outside acceptable limits", "#F03335"
        
        

def atoms_assessment(react: Reaction) -> tuple[str,str]:
    """
    Warn if any risky elements are present in one of the products.
    Returns (message, color_hex).
    """
    green_scores = {
        "C":100, "H":100, "O":100, "N":100, "P":100, "S":100,
        "B":85,  "Si":85, "Mg":80,  "Fe":75, "Al":70, "Zn":70,
        "F":60,  "Cl":60, "Br":60,  "I":60,  "Li":60, "Ti":50,
        "Sn":30, "Pb":20, "Pd":20,  "Hg":20,  "Cd":20, "As":20,
        "Cr":20, "Ni":20, "Se":20,  "Tl":20,  "Pt":10, "Rh":10
    }

    risky_atoms = set()
    for mol in react.products.keys():  
        risky_atoms |= {
            a.GetSymbol()
            for a in mol.GetAtoms()
            if green_scores.get(a.GetSymbol(), 50) <= 60
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