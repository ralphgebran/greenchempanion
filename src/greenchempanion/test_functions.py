from .functions import Atom_Count_With_H, Reaction, compute_PMI, Canonicalize_Smiles

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Descriptors
from typing import Dict



# TESTING Atom_Count_With_H FUNCTION

ethanol = Chem.MolFromSmiles("CCO")
benzene = Chem.MolFromSmiles("c1ccccc1")


def test_atomcount_ethanol():
    assert Atom_Count_With_H(ethanol) == 9, f"test_atomcount_ethanol failed: Expected output is 9 ; Actual output is {Atom_Count_With_H(ethanol)}"
test_atomcount_ethanol()

def test_atomcount_benzene():
    assert Atom_Count_With_H(benzene) == 12, f"test_atomcount_benzene failed: Expected output is 12 ; Actual output is {Atom_Count_With_H(benzene)}"
test_atomcount_benzene()


# TESTING canonicalization

def test_canonicalization():
    assert Canonicalize_Smiles("C(C)O") == "CCO", f"test_canonicalization failed: Expected output is CCO ; Actual output is {Canonicalize_Smiles('C(C)O')} "
    assert Canonicalize_Smiles("OCC") == "CCO", f"test_canonicalization failed: Expected output is CCO ; Actual output is {Canonicalize_Smiles('OCC')} "


#TESTING Atom_Economy_X FUNCTIONS

methane = Chem.MolFromSmiles("C")
dioxygen = Chem.MolFromSmiles("O=O")
carbon_dioxide = Chem.MolFromSmiles("O=C=O")
water = Chem.MolFromSmiles("O")

methane_comb = Reaction({methane: 1, dioxygen:2}, {carbon_dioxide:1, water:2})


def test_atomecon_a_methanecomb():
    result = methane_comb.Atom_Economy_A()
    expected = 33.33
    assert abs(result - expected) < 1e-2, f"test_atomecon_a_methanecomb failed: Expected output is 33.33 ; Actual output is {methane_comb.Atom_Economy_A()}"
test_atomecon_a_methanecomb()

def test_atomecon_m_methanecomb():
    result = methane_comb.Atom_Economy_M()
    expected = 55
    assert abs(result - expected) < 1e-1, f"test_atomecon_m_methanecomb failed: Expected output is 55 ; Actual output is {methane_comb.Atom_Economy_M()}"
test_atomecon_m_methanecomb()



# TESTING THE PMI FUNCTION

#First Reaction
sodium_azide = Chem.MolFromSmiles("[N-]=[N+]=[N-].[Na+]")
two_iodobutane = Chem.MolFromSmiles("CCC(I)C")
two_azidobutane = Chem.MolFromSmiles("CCC([N]=[N+]=[N-])C")
sodium_iodide = Chem.MolFromSmiles("[Na+].[I-]")

acetone = Chem.MolFromSmiles("CC(=O)C")
ether = Chem.MolFromSmiles("CCOCC")

sn_reaction = Reaction({two_iodobutane:1, sodium_azide:1}, {two_azidobutane:1, sodium_iodide:1})
sn_reaction_pmi = compute_PMI(sn_reaction, {acetone:17470, water:112000, ether:81800}, 0.9)


def test_pmi_sn():
    result = sn_reaction_pmi
    expected = 214
    assert abs(result - expected) < 1e-1, f"test_pmi_sn failed: Expected output is 214 ; Actual output is {sn_reaction_pmi}"
test_pmi_sn()

#Second Reaction
butan_two_ol = Chem.MolFromSmiles("CCC(O)C")
cr_oxide = Chem.MolFromSmiles("O=[Cr](=O)=O")
butan_two_one = Chem.MolFromSmiles("CCC(=O)C")
cr_acid = Chem.MolFromSmiles("O[Cr](=O)O")

ox_reaction = Reaction({butan_two_ol:1, cr_oxide:1}, {butan_two_one:1, cr_acid:1})
ox_reaction_pmi = compute_PMI(ox_reaction, {ether:37500}, 0.9)


def test_pmi_ox():
    result = ox_reaction_pmi
    expected = 40.2
    assert abs(result - expected) < 1e-1, f"test_pmi_ox failed: Expected output is 40.2 ; Actual output is {ox_reaction_pmi}"
test_pmi_ox()
