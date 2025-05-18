from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors
from typing import Dict
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import TanimotoSimilarity


def Atom_Count_With_H(mol: Mol) -> int:
    """
    Returns the total number of atoms for a given molecule, counting the H atoms
    
    Arguments:
    - mol: RDKit Mol object
    
    Returns:
    - int: Number of atoms
    """
    H_mol = Chem.AddHs(mol)
    return H_mol.GetNumAtoms()


def canonicalize_smiles(smiles : str) -> str:
    """
    Returns the canonicalized SMILES for a given SMILES string.

    Arguments:
    - smiles : A valid SMILES string representing a molecule.

    Returns:
    - str: The canonicalized SMILES string, or None if input is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol)


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
    
    def Atom_Economy_M(self) -> float:
        """
        Calculates the atom economy of the reaction, based on molar masses.
        (Molar mass of main product / Total molar mass of reactants) * 100

        Returns:
        - float: Atom economy as a percentage (%)
        """
        reactants_mass = sum(Descriptors.MolWt(mol) * coeff for mol,coeff in self.reactants.items())
        main_product_mass = Descriptors.MolWt(self.main_product) * self.products[self.main_product]
        economy_m = (main_product_mass / reactants_mass) * 100
        if economy_m > 100 : 
            raise ValueError("Atom Economy cannot be greater than 100%, recheck your reaction")
        return economy_m

    def isBalanced(self) -> bool:
        """
        Checks whether the reaction is balanced or not.

        Returns:
        - bool: True if the reaction is balanced, False if it's not.
        """
        def atom_count_dict(mol_dict: Dict[Mol, int]) -> Dict:
            atom_counts = {}
            for mol, coeff in mol_dict.items():
                h_mol = Chem.AddHs(mol)
                for atom in h_mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    if symbol not in atom_counts:
                        atom_counts[symbol] = 0
                    atom_counts[symbol] += coeff
            return atom_counts
        
        reac_count = atom_count_dict(self.reactants)
        prod_count = atom_count_dict(self.products)
        return reac_count == prod_count
    

def compute_PMI(reaction: Reaction, extras: Dict[Mol, float], prod_yield: float) -> float:
    """
    Computes the Process Mass Intensity (PMI) per kg of main product.
    (Total mass of chemical inputs / 1kg of product) * 100

    Arguments:
    - reaction: GCP Reaction object (with reactants, products, and main product index)
    - extras: Dictionary containing solvents, extraction material and others, with their mass [in grams] per kg of product {Mol: float}
    - prod_yield: Yield of main product (0 < float <= 1)
    
    Returns:
    - PMI (float): total input mass (g) per 1 kg of product
    """
    
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

    Arguments:
    - reaction: GCP Reaction object
    - extras: Dictionary {Mol: float} for solvents, etc.
    - prod_yield: Yield of main product (0 < float <= 1)
    
    Returns:
    - E-Factor (float)
    """

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
