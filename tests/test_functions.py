from greenchempanion.functions import Atom_Count_With_H, Reaction, compute_PMI, canonicalize_smiles, compute_E, atoms_assessment, structural_assessment
from greenchempanion.assessments import get_solvent_info, waste_efficiency, PMI_assesment, Atom_ec_assesment, Atom_ec_m_assesment, logP_assessment_molecule

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
    assert canonicalize_smiles("C(C)O") == "CCO", f"test_canonicalization failed: Expected output is CCO ; Actual output is {canonicalize_smiles('C(C)O')} "
    assert canonicalize_smiles("OCC") == "CCO", f"test_canonicalization failed: Expected output is CCO ; Actual output is {canonicalize_smiles('OCC')} "


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
sn_reaction_e= compute_E(sn_reaction, {acetone:17470, water:112000, ether:81800}, 0.9)

def test_pmi_sn():
    result = sn_reaction_pmi
    expected = 214
    assert abs(result - expected) < 1e-1, f"test_pmi_sn failed: Expected output is 214 ; Actual output is {sn_reaction_pmi}"
test_pmi_sn()

def test_e_sn():
    result = sn_reaction_e
    expected = 212.9
    assert abs(result - expected) < 1e-1, f"test_e_sn failed: Expected output is 212.9 ; Actual output is {sn_reaction_e}"
test_e_sn()

#Second Reaction
butan_two_ol = Chem.MolFromSmiles("CCC(O)C")
cr_oxide = Chem.MolFromSmiles("O=[Cr](=O)=O")
butan_two_one = Chem.MolFromSmiles("CCC(=O)C")
cr_acid = Chem.MolFromSmiles("O[Cr](=O)O")

ox_reaction = Reaction({butan_two_ol:1, cr_oxide:1}, {butan_two_one:1, cr_acid:1})
ox_reaction_pmi = compute_PMI(ox_reaction, {ether:37500}, 0.9)
ox_reaction_e = compute_E(ox_reaction, {ether:37500}, 0.9)

def test_pmi_ox():
    result = ox_reaction_pmi
    expected = 40.2
    assert abs(result - expected) < 1e-1, f"test_pmi_ox failed: Expected output is 40.2 ; Actual output is {ox_reaction_pmi}"
test_pmi_ox()

def test_e_ox():
    result = ox_reaction_e
    expected = 39.1
    assert abs(result - expected) < 1e-1, f"test_e_ox failed: Expected output is 39.1 ; Actual output is {ox_reaction_e}"
test_e_ox()


# TESTING THE REACTION BALANCE METHOD

#Iron Sulfate Reaction
iron_bromide = Chem.MolFromSmiles("[Fe+3].[Br-].[Br-].[Br-]")
sulfuric_acid = Chem.MolFromSmiles("OS(=O)(=O)O")
iron_sulfate = Chem.MolFromSmiles("[Fe+3].[Fe+3].[O-]S(=O)(=O)[O-].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O")
hbr = Chem.MolFromSmiles("Br")

unbalanced_reaction = Reaction({iron_bromide:1, sulfuric_acid:1}, {iron_sulfate:1, hbr:1})
balanced_reaction = Reaction({iron_bromide:2, sulfuric_acid:3}, {iron_sulfate:1, hbr:6})
nonsense_reaction = Reaction({acetone:18, water:2}, {dioxygen:3, cr_acid:7})
one_c_reaction = Reaction({Chem.MolFromSmiles("[C]"):2, dioxygen:1}, {carbon_dioxide:1})

def test_balance_iron_false():
    assert unbalanced_reaction.isBalanced() == False, f"test_balance_iron_false failed: Expected output is False ; Actual output is {unbalanced_reaction.isBalanced()}"
test_balance_iron_false()

def test_balance_iron_true():
    assert balanced_reaction.isBalanced() == True, f"test_balance_iron_true failed: Expected output is True ; Actual output is {balanced_reaction.isBalanced()}"
test_balance_iron_true()

def test_balance_methane():
    assert methane_comb.isBalanced() == True, f"test_balance_methane failed: Expected output is True ; Actual output is {methane_comb.isBalanced()}"
test_balance_methane()

def test_balance_nonsense():
    assert nonsense_reaction.isBalanced() == False, f"test_balance_nonsense failed: Expected output is False ; Actual output is {nonsense_reaction.isBalanced()}"
test_balance_nonsense()

def test_balance_one_c():
    assert one_c_reaction.isBalanced() == False, f"test_balance_one_c failed: Expected output is False ; Actual output is {nonsense_reaction.isBalanced()}"
test_balance_one_c()


def test_get_solvent_info_green():
    mol = Chem.MolFromSmiles("O")
    verdict, color = get_solvent_info({mol: 1.0})
    assert "Green" in verdict and color == "#88DF66"

def test_waste_efficiency_bad():
    verdict, color = waste_efficiency(75)
    assert "Bad" in verdict and color == "#F68B8C"

def test_pmi_assessment_very_bad():
    verdict, color = PMI_assesment(150)
    assert "Very Bad" in verdict and color == "#F03335"

def test_atom_ec_assessment_excellent():
    verdict, color = Atom_ec_assesment(95)
    assert "Excellent" in verdict and color == "#4BAF24"

def test_atom_ec_m_assessment_moderate():
    verdict, color = Atom_ec_m_assesment(65)
    assert "Moderate" in verdict and color == "#F6DF7E"

def test_logP_assessment_out_of_range():
    verdict, color = logP_assessment_molecule(5.5)
    assert "outside" in verdict and color == "#F03335"
    