from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors
from typing import Dict



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



# TESTING Atom_Count_With_H FUNCTION

#ethanol = Chem.MolFromSmiles("CCO")
#print(ethanol.GetNumAtoms()) #Should be 3
#print(Atom_Count_With_H(ethanol)) #Should be 9

#benzene = Chem.MolFromSmiles("c1ccccc1")
#print(benzene.GetNumAtoms()) #Should be 6
#print(Atom_Count_With_H(benzene)) #Should be 12



# Reaction Class

class Reaction:
    """
    A class to represent a chemical reaction.

    Attributes:
    - reactants (Dict[Mol, int]): Dictionary of RDKit Mol objects with their Stoichiometric Coefficients, representing the reactants.
    - products (List[Mol]): Dictionary of RDKit Mol objects with their Stoichiometric Coefficients, representing the products.
    - main_product_index (int): The index of the main product, in the products list (Set as 0 by default)

    Methods:
    - Atom_Economy_A(): Calculates and returns the atom economy of the reaction, based on numbers of atoms.
    - Atom_Economy_M(): Calculates and returns the atom economy of the reaction, based on molar masses.
    """
    
    def __init__(self, reactants: Dict[Mol, int], products: Dict[Mol, int], main_product_index: int = 0):
        self.reactants = reactants
        self.products = products

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
        return (main_product_atoms / reactants_atoms) * 100
    
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

 
        
#TESTING Atom_Economy_X FUNCTIONS

#methane = Chem.MolFromSmiles("C")
#dioxygen = Chem.MolFromSmiles("O=O")
#carbon_dioxide = Chem.MolFromSmiles("O=C=O")
#water = Chem.MolFromSmiles("O")

#methane_comb = Reaction({methane: 1, dioxygen:2}, {carbon_dioxide:1, water:2})
#print(methane_comb.Atom_Economy_A()) # SHOULD BE 33.33%
#print(methane_comb.Atom_Economy_M()) # SHOULD BE 55%