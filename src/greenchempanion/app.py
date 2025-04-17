import streamlit as st
import numpy
from rdkit import Chem
from rdkit.Chem import Draw
from functions import Atom_Count_With_H, Reaction

st.set_page_config(page_title="GreenChemPanion", page_icon="🍃")
st.title("🍃 GCP: GreenChemPanion") #TITLE
st.write("Interactive Streamlit Applet showcasing the functions for GCP!") #TEXT

# SMILES TO MOLECULE IMAGE DISPLAY
st.header("😊 SMILES to Molecule Converter")
M2S_smiles = st.text_input("Enter a SMILES:", key ="Mol2Smiles_input")

if M2S_smiles:

    #Developers Easter Egg
    if M2S_smiles == "VRMT":
        st.write("Easter Egg :o")
        st.image("../../assets/banner.png")

    elif Chem.MolFromSmiles(M2S_smiles) is None:
        st.write("⚠️ Enter a valid SMILES format [e.g.: CCO]")

    else:
        input_mol = Chem.MolFromSmiles(M2S_smiles)
        img = Draw.MolToImage(input_mol)
        st.image(img)

        #Atom Counts
        st.subheader("ATOM COUNTING 1️⃣2️⃣3️⃣")
        st.write(f"Number of large atoms: {input_mol.GetNumAtoms()}")
        st.write(f"Number of atoms (counting Hydrogens): {Atom_Count_With_H(input_mol)}")

#SEPARATOR IN STREAMLIT :O
st.markdown("---")

# REACTION SECTION
st.header("💥 Enter a Reaction")

# Storage values
if "reactants" not in st.session_state:
    st.session_state.reactants = {}
if "products" not in st.session_state:
    st.session_state.products = {}
if "main_product_index" not in st.session_state:
    st.session_state.main_product_index = 0  # default to first product

st.subheader("🧪 Add a Molecule")

# SMILES input
R_smiles = st.text_input("Enter a SMILES:", key ="ReactionSmiles_input")

# Stoichiometry Coefficient input
stoich = st.number_input("Stoichiometric Coefficient:", min_value=1, value=1, step=1)

# Reaction or Product
role = st.radio("Add as:", ["Reactant", "Product"], horizontal=True)

# Confirm button
if st.button("➕ Add Molecule"):
    R_mol = Chem.MolFromSmiles(R_smiles)
    if R_mol:
        if role == "Reactant":
            st.session_state.reactants[R_mol] = stoich
            st.success(f"Added reactant: {R_smiles} (x{stoich})")
        else:
            st.session_state.products[R_mol] = stoich
            st.success(f"Added product: {R_smiles} (x{stoich})")
    else:
        st.error("⚠️ Enter a valid SMILES format [e.g.: CCO]")


st.subheader("📦 Stored Reaction")

# Reactants section
st.write("### ⚗️ Reactants")
if st.session_state.reactants:
    for mol in list(st.session_state.reactants):  # make a copy of keys to allow deletion during loop
        qty = st.session_state.reactants[mol]
        col1, col2, col3 = st.columns([1, 3, 1])
        with col1:
            st.image(Draw.MolToImage(mol, size=(60, 60)))
        with col2:
            st.markdown(f"**SMILES:** `{Chem.MolToSmiles(mol)}`  \n**Stoichiometric Coefficient:** {qty}")
        with col3:
            if st.button("❌", key=f"del_reactant_{Chem.MolToSmiles(mol)}"):
                del st.session_state.reactants[mol]
                st.rerun()
else:
    st.info("No reactants added.")

# Products section
st.write("### 🔬 Products")


product_mols = list(st.session_state.products.keys())

# Main Product selection
if product_mols:
    selected = st.radio(
        "Select main product:",
        options=range(len(product_mols)),
        format_func=lambda i: f"{Chem.MolToSmiles(product_mols[i])} (x{st.session_state.products[product_mols[i]]})",
        index=st.session_state.main_product_index,
        key="main_product_selector"
    )
    st.session_state.main_product_index = selected

    # Display each product with delete button and highlight if main
    for i, mol in enumerate(product_mols):  # use fixed list for safe deletion
        qty = st.session_state.products[mol]
        col1, col2, col3 = st.columns([1, 3, 1])
        with col1:
            st.image(Draw.MolToImage(mol, size=(60, 60)))
        with col2:
            st.markdown(f"**SMILES:** `{Chem.MolToSmiles(mol)}`  \n**Stoichiometric Coefficient:** {qty}")
            if i == st.session_state.main_product_index:
                st.markdown("⭐ **Main product**")
        with col3:
            if st.button("❌", key=f"del_product_{Chem.MolToSmiles(mol)}"):
                del st.session_state.products[mol]
                # Adjust main index if necessary
                if st.session_state.main_product_index >= len(st.session_state.products):
                    st.session_state.main_product_index = max(0, len(st.session_state.products) - 1)
                st.rerun()
else:
    st.info("No products added.")

st.write("### 🧮 Atom Economy")
# === Compute Atom Balance button only if both reactants and products are set ===
if st.session_state.reactants and st.session_state.products:
    if st.button("Compute Atom Economy"):
        # Retrieve from session
        R_reactants = st.session_state.reactants
        R_products = st.session_state.products
        R_main_index = st.session_state.main_product_index

        # Convert to list to get main product
        product_mols = list(R_products.keys())
        main_product = product_mols[R_main_index]

        # Create the Reaction object
        input_reaction = Reaction(reactants=R_reactants, products=R_products, main_product_index=R_main_index)

        # Compute metrics (replace with your real logic)
        atom_economy_m_result = input_reaction.Atom_Economy_M()
        atom_economy_a_result = input_reaction.Atom_Economy_A()

        # Display results
        st.success(f"⚖️ Atom Economy based on Molar Mass: **{atom_economy_m_result:.2f}%**")
        st.success(f"⚛️ Atom Economy based on Number of Atoms: **{atom_economy_a_result:.2f}%**")
else:
    st.info("Add at least one reactant and one product to compute atom economy.")