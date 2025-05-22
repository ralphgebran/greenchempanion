import streamlit as st
import pandas as pd
import math
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from streamlit_ketcher import st_ketcher
from functions import Atom_Count_With_H, Reaction, compute_PMI, canonicalize_smiles, compute_E, atoms_assessment, structural_assessment
from assessments import get_solvent_info, waste_efficiency, PMI_assessment, Atom_ec_assessment, Atom_ec_m_assessment, logP_assessment_molecule
from app_utilities import (inject_base_css, show_info_box, html_box, tooltip_icon, missing_input_alert, dual_metric_box, title_with_icon, 
                           SOLVENT_TIP, LOGP_TIP, AA_TIP, SAA_TIP, E_FACTOR_TIP, PMI_TIP, AE_M_TIP, AE_A_TIP)


st.set_page_config(page_title="GreenChemPanion", page_icon="../../assets/logo.ico", layout= "wide")

col1, col2 = st.columns([2, 13])

with col1:
    st.image("../../assets/logo_icon.png", width=1000)

with col2:
    st.markdown("## GCP: GreenChemPanion\nInteractive Streamlit Applet showcasing the functions for GCP!")
    
inject_base_css()
st.markdown("---")

st.header("💱 SMILES Converters", anchor="Converters")


# SMILES TO MOLECULE IMAGE DISPLAY

with st.expander("😊 SMILES to Molecule Converter"):
    M2S_smiles = st.text_input("Enter a SMILES:", key ="Mol2Smiles_input")

    if M2S_smiles:
        #Developers Easter Egg
        if M2S_smiles == "VRMT"and not st.session_state.get("easter_egg_shown", False):
            st.balloons()
            st.snow()
            st.write("Easter Egg :o")
            st.image("../../assets/vrmt.png")
            st.session_state["easter_egg_shown"] = True
        elif M2S_smiles != "VRMT":
            st.session_state["easter_egg_shown"] = False
            
        if Chem.MolFromSmiles(M2S_smiles) is None:
            st.write("⚠️ Enter a valid SMILES format [e.g.: CCO]")

        else:
            input_mol = Chem.MolFromSmiles(M2S_smiles)
            img = Draw.MolToImage(input_mol)
            st.image(img)
            #Atom Counts
            st.subheader("ATOM COUNTING 1️⃣2️⃣3️⃣", anchor=False)
            st.write(f"Number of large atoms: {input_mol.GetNumAtoms()}")
            st.write(f"Number of atoms (counting Hydrogens): {Atom_Count_With_H(input_mol)}")


# CONVERT A MOLEUCLE DRAWING INTO ITS SMILES 

def draw_and_get_smiles() -> str:
    """
    Streamlit function to draw a molecule using Ketcher and return the SMILES string.
    
    Returns:
        str: The canonical SMILES of the drawn molecule.
    """

    st.subheader("Draw Molecule And Get SMILES", anchor= False)
    molecule = ""
    if "ketcher_key" not in st.session_state:
        st.session_state.ketcher_key = 0

    ketcher_smiles = st_ketcher(
        molecule, 
        height=600,
        key=f"ketcher_{st.session_state.ketcher_key}")
    
    st.button("🔄 Reset drawing and SMILES", on_click=lambda: st.session_state.update({"ketcher_key": st.session_state.ketcher_key + 1}))

    while ketcher_smiles:
        try:
            mol = Chem.MolFromSmiles(ketcher_smiles)
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                st.subheader("SMILES from your drawing:", anchor= False)
                st.success(canonical_smiles)
                return canonical_smiles
            else:
                st.error("❌ Invalid molecule drawn! Please try again.")
                return ""
        except Exception as e:
            st.error(f"⚠️ Error parsing molecule: {e}")
            return ""
    if not ketcher_smiles :
        st.info("🧪 Draw a molecule above to generate its SMILES")
        return ""

with st.expander("✏️ Drawn Molecule to SMILES Converter"):
    smiles_result = draw_and_get_smiles()
    if smiles_result:
        st.write(f"The SMILES you generated is: {smiles_result}")
        #Atom Counts
        st.subheader("ATOM COUNTING 1️⃣2️⃣3️⃣", anchor=False)
        st.write(f"Number of large atoms: {Chem.MolFromSmiles(smiles_result).GetNumAtoms()}")
        st.write(f"Number of atoms (counting Hydrogens): {Atom_Count_With_H(Chem.MolFromSmiles(smiles_result))}")


# REACTION SECTION

st.header("💥 Enter a Reaction", anchor="enter-reaction")

if "reactants" not in st.session_state:
    st.session_state.reactants = {}
if "products" not in st.session_state:
    st.session_state.products = {}
if "main_product_index" not in st.session_state:
    st.session_state.main_product_index = 0  # default to first product
if "extras" not in st.session_state:
    st.session_state.extras = {}
if "prod_yield" not in st.session_state:
    st.session_state.prod_yield = 1 # default to 100%

# Add a Molecule Expander

with st.expander("🧪 Add a Molecule"):

    # SMILES input
    R_smiles = st.text_input("Enter a SMILES:", key ="ReactionSmiles_input")

    # Stoichiometry Coefficient input
    stoich = st.number_input("Stoichiometric Coefficient:", min_value=1, value=1, step=1)

    # Reaction or Product
    role = st.radio("Add as:", ["Reactant", "Product"], horizontal=True)

    # Confirm button
    if st.button("➕ Add Molecule", key ="Molecule_Add"):
        R_mol = Chem.MolFromSmiles(R_smiles)
        smiles_reactants = [Chem.MolToSmiles(mol) for mol in st.session_state.reactants]
        smiles_products = [Chem.MolToSmiles(mol) for mol in st.session_state.products]
        if R_mol:
            R_smiles = canonicalize_smiles(R_smiles)
            if role == "Reactant":
                if R_smiles in smiles_reactants:
                    for mol in st.session_state.reactants:
                        if Chem.MolToSmiles(mol) == R_smiles:
                            corresponding_mol = mol
                    st.session_state.reactants[corresponding_mol] += stoich
                    st.success(f"Added {stoich} more of reactant {R_smiles} ")
                else:
                    st.session_state.reactants[R_mol] = stoich
                    st.success(f"Added reactant: {R_smiles} (x{stoich})")
            else:
                if R_smiles in smiles_products:
                    for mol in st.session_state.products:
                        if Chem.MolToSmiles(mol) == R_smiles:
                            corresponding_mol = mol
                    st.session_state.products[corresponding_mol] += stoich
                    st.success(f"Added {stoich} more of product {R_smiles} ")
                else:
                    st.session_state.products[R_mol] = stoich
                    st.success(f"Added a product: {R_smiles} (x{stoich})")
        else:
            st.error("⚠️ Enter a valid SMILES format [e.g.: CCO]", icon="🚨")

# Common solvent expander

with st.expander("💧 Add a Common Solvent (as a volume input)"):
    st.info("ℹ️ Densities are considered at NTP Conditions (20 °C, 1 atm)")

    #Creating data frame from csv
    solvent_df = pd.read_csv("../../data/solvent_data.csv", delimiter=";")
    solvent_list = solvent_df["Solvent"].tolist()

    #Solvent Selection
    solvent_select =st.selectbox("Choose a solvent:", solvent_list, index=None, placeholder="Solvent...")

    #Volume input
    solvent_vol = st.number_input("Enter the volume (L) used per kg of product:", min_value=0.0, step=0.01)

    # Confirm button
    if st.button("➕ Add Solvent", key ="Solvent_Add"):
        density = solvent_df.loc[solvent_df["Solvent"] == solvent_select, "Density"].values[0]
        solvent_smiles = solvent_df.loc[solvent_df["Solvent"] == solvent_select, "SMILES"].values[0]
        solvent_mass = float(round(density*solvent_vol*1000))

        S_mol = Chem.MolFromSmiles(solvent_smiles)
        smiles_solvent = [Chem.MolToSmiles(mol) for mol in st.session_state.extras]

        S_smiles = canonicalize_smiles(solvent_smiles)
        if S_smiles in smiles_solvent :
            for mol in st.session_state.extras:
                if Chem.MolToSmiles(mol) == S_smiles:
                    corresponding_mol= mol
            st.session_state.extras[corresponding_mol] += solvent_mass
            st.success(f"Added solvent : another {solvent_vol} L ({solvent_mass} g) of {solvent_select} ({S_smiles})")
        else :
            st.session_state.extras[S_mol] = solvent_mass
            st.success(f"Added solvent: {solvent_select} ({S_smiles}), {solvent_vol} L ({solvent_mass} g)")

# Solvents expander

with st.expander("🪣 Add extras (Yield, Solvents, Extraction material...)"):

    # SMILES input
    E_smiles = st.text_input("Enter a SMILES:", key ="ExtrasSmiles_input")

    # Mass input
    E_mass = st.number_input("Mass [in grams] per 1 kg of product:", min_value=0, value=1)

    # Confirm button
    if st.button("➕ Add Molecule", key ="Extras_Add"):
        E_mol = Chem.MolFromSmiles(E_smiles)
        smiles_extras = [Chem.MolToSmiles(mol) for mol in st.session_state.extras]
        if E_mol:
            E_smiles = canonicalize_smiles(E_smiles)
            if E_smiles in smiles_extras :
                for mol in st.session_state.extras:
                    if Chem.MolToSmiles(mol) == E_smiles:
                        corresponding_mol= mol
                st.session_state.extras[corresponding_mol] += E_mass
                st.success(f"Added extra : another ({E_mass} g ) of {E_smiles}")
            else :
                st.session_state.extras[E_mol] = E_mass
                st.success(f"Added extra: {E_smiles} ({E_mass} g)")
        else:
            st.error("⚠️ Enter a valid SMILES format [e.g.: CCO]")
    
    #YIELD input
    E_yield = st.number_input("Enter main product yield (%)", min_value=0.0, max_value=100.0, value =100.0, step = 0.01)
    st.session_state.prod_yield = E_yield/100


# DISPLAY OF REACTION AND COMPUTATION OF GREEN FACTORS

column1, column2 = st.columns(2)
# Stored reaction
with column1 :
    st.subheader("📦 Stored Reaction", anchor="stored-reaction")

    # Error message if Reaction is not balanced
    if st.session_state.reactants and st.session_state.products:
        test_reac = Reaction(st.session_state.reactants, st.session_state.products)
        if test_reac.isBalanced() == False:
            st.error("⚠️ The reaction entered is not balanced. Results shown may be incoherent.")

    # Reactants section
    st.write("##### ⚗️ Reactants")
    if st.session_state.reactants:
        for mol in list(st.session_state.reactants):  # make a copy of keys to allow deletion during loop
            qty = st.session_state.reactants[mol]
            col1, col2, col3 = st.columns([1, 3, 1])
            with col1:
                st.image(Draw.MolToImage(mol, size=(175, 175)))
            with col2:
                st.markdown(f"**SMILES:** `{Chem.MolToSmiles(mol)}`  \n**Stoichiometric Coefficient:** {qty}")
            with col3:
                if st.button("❌", key=f"del_reactant_{Chem.MolToSmiles(mol)}"):
                    del st.session_state.reactants[mol]
                    st.rerun()
    else:
        st.info("No reactants added.")

    # Products section
    st.write("##### 🔬 Products")

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
            img_col, txt_col, btn_col = st.columns([1, 3, 1])
            with img_col:
                st.image(Draw.MolToImage(mol, size=(175, 175)))
            with txt_col:
                st.markdown(f"**SMILES:** `{Chem.MolToSmiles(mol)}`  \n**Stoichiometric Coefficient:** {qty}")
                if i == st.session_state.main_product_index:
                    st.markdown("⭐ **Main product**")
            with btn_col:
                if st.button("❌", key=f"del_product_{Chem.MolToSmiles(mol)}"):
                    del st.session_state.products[mol]
                    # Adjust main index if necessary
                    if st.session_state.main_product_index >= len(st.session_state.products):
                        st.session_state.main_product_index = max(0, len(st.session_state.products) - 1)
                    st.rerun()
    else:
        st.info("No products added.")

        
    st.write(" ##### 🚰 Extras")
                
    st.markdown(f"Product yield: `{st.session_state.prod_yield * 100 } %` ")
    if st.session_state.extras:
        for mol in list(st.session_state.extras):  # make a copy of keys to allow deletion during loop
            qty = st.session_state.extras[mol]
            col1, col2, col3 = st.columns([1, 3, 1])
            with col1:
                st.image(Draw.MolToImage(mol, size=(175, 175)))
            with col2:
                st.markdown(f"**SMILES:** `{Chem.MolToSmiles(mol)}`  \n**Mass (/kg product):** {qty} g")
            with col3:
                if st.button("❌", key=f"del_extra_{Chem.MolToSmiles(mol)}"):
                    del st.session_state.extras[mol]
                    st.rerun()
    else:
        st.info("No extras added.")

# Green Chemistry Factors

from app_utilities import tooltip_icon

with column2:
    st.subheader("🧮 Compute Factors", anchor="compute-factors")

    # Compute E factor button only if both reactants and products are set
    if st.session_state.reactants and st.session_state.products:
        # Retrieve from session
        E_reactants = st.session_state.reactants
        E_products = st.session_state.products
        E_main_index = st.session_state.main_product_index

        E_extras = st.session_state.extras
        E_yield = st.session_state.prod_yield

        # Convert to list to get main product
        product_mols = list(E_products.keys())
        main_product = product_mols[E_main_index]

        # Create the Reaction object
        input_reaction = Reaction(reactants=E_reactants, products=E_products, main_product_index=E_main_index)

        # Compute metrics (replace with your real logic)
        E_result = compute_E(input_reaction, E_extras, E_yield)

        # Display results
        col_btn, col_tip = st.columns([8, 1])

        with col_btn:
            compute_e_clicked = st.button("Compute E factor", key="btn_e_factor")
 
        with col_tip:
            col_tip.markdown(
                tooltip_icon("ℹ️",
                    """
                    <b>E factor (Environmental factor)</b> is the ratio of the total mass of all waste generated 
                    (side-products, solvents, auxiliaries, etc.) to the mass of desired product obtained.
                    Lower values indicate a greener, more sustainable process,
                    as less waste is produced per kilogram of product.
                    """
                    ),
                unsafe_allow_html=True
            )
        if compute_e_clicked :
            st.success(f"🚮 E: **{E_result:.2f}**")
    else:
        missing_input_alert(
        "Add at least one reactant and one product to compute E factor.",
        """
        <b>E factor (Environmental factor)</b> is the ratio of the total mass of all waste generated 
        (side-products, solvents, auxiliaries, etc.) to the mass of desired product obtained.
        Lower values indicate a greener, more sustainable process,
        as less waste is produced per kilogram of product.
        """
        )

    # Compute PMI button only if both reactants and products are set
    if st.session_state.reactants and st.session_state.products:
        # Retrieve from session
        PMI_reactants = st.session_state.reactants
        PMI_products = st.session_state.products
        PMI_main_index = st.session_state.main_product_index

        PMI_extras = st.session_state.extras
        PMI_yield = st.session_state.prod_yield

        # Convert to list to get main product
        product_mols = list(PMI_products.keys())
        main_product = product_mols[PMI_main_index]

        # Create the Reaction object
        input_reaction = Reaction(reactants=PMI_reactants, products=PMI_products, main_product_index=PMI_main_index)

        # Compute metrics 
        PMI_result = compute_PMI(input_reaction, PMI_extras, PMI_yield)

        # Display results
        col_btn, col_tip = st.columns([8, 1])
        
        with col_btn:
            compute_clicked_2 = st.button("Compute PMI", key="btn_pmi")

        with col_tip:
            col_tip.markdown(
                tooltip_icon("ℹ️",
                    """
                    <b>Process Mass Intensity (PMI) </b> is the total mass of all materials
                    used in the process divided by the mass of product obtained.
                    Lower values (maximum&nbsp;100&nbsp;%) indicate a greener, more material-efficient process.
                    """
                    ),
                unsafe_allow_html=True
            )
        if compute_clicked_2:
            st.success(f"🅿️ PMI: **{PMI_result:.2f}**")

    else:
        missing_input_alert(
        "Add at least one reactant and one product to compute PMI.",
        """
        <b>Process Mass Intensity (PMI) </b> is the total mass of all materials
        used in the process divided by the mass of product obtained.
        Lower values (maximum&nbsp;100&nbsp;%) indicate a greener, more material-efficient process.
        """
        )

        
    # Compute Atom Balance button only if both reactants and products are set
    if st.session_state.reactants and st.session_state.products:
        # Retrieve from session
        R_reactants = st.session_state.reactants
        R_products = st.session_state.products
        R_main_index = st.session_state.main_product_index

        # Convert to list to get main product
        product_mols = list(R_products.keys())
        main_product = product_mols[R_main_index]

        # Create the Reaction object
        input_reaction = Reaction(reactants=R_reactants, products=R_products, main_product_index=R_main_index)

        col_btn, col_tip = st.columns([8, 1])

        with col_btn:
            compute_clicked = st.button("Compute Atom Economy", key="btn_atom_economy")

        with col_tip:
            col_tip.markdown(
                tooltip_icon("ℹ️", 
                        """<b>Atom Economy</b> is the percentage of the mass (or atoms) of reactants that ends up in the desired product. 
                        Higher values (max 100%) indicate a greener reaction."
                        """),
                unsafe_allow_html=True
            )
        if compute_clicked:
            try:
                atom_economy_m_result = input_reaction.Atom_Economy_M()
                atom_economy_a_result = input_reaction.Atom_Economy_A()
                st.success(f"⚖️ Atom Economy based on Molar Mass: **{atom_economy_m_result:.2f}%**")
                st.success(f"⚛️ Atom Economy based on Number of Atoms: **{atom_economy_a_result:.2f}%**")
            except ValueError as e:
                st.error(f"{e}", icon="🚨")
                    
    else:  
        missing_input_alert(
        "Add at least one reactant and one product to compute Atom Economy.",
        """
        <b>Atom Economy</b> is the percentage of the mass (or atoms) of reactants that ends up in the desired product. 
        Higher values (max 100%) indicate a greener reaction."
        """
        )


# GCP GREEN CHEMISTRY EVALUATION 

st.header("🧑‍🔬 GCP Green Chemistry Evaluation", anchor="GCP_evaluation")
from app_utilities import show_info_box
if st.session_state.reactants and st.session_state.products:
    col1, col2 = st.columns(2)

    with col1:
        if st.session_state.extras:
            solv_text, solv_color = get_solvent_info(st.session_state.extras)
            show_info_box(title_with_icon("🌱 Solvent Quality Check", SOLVENT_TIP), content=solv_text, color=solv_color)
        else:
            show_info_box("🌱 Solvent Quality Check",
              "Seems like no solvent was used. Nice Job!")
        Log_P = Descriptors.MolLogP(input_reaction.main_product)
        log_text, log_color = logP_assessment_molecule(Log_P)
        show_info_box(title_with_icon("🌍 Main Product LogP", LOGP_TIP), content=log_text, color=log_color)

        a_ass_text, a_ass_color = atoms_assessment(input_reaction)
        show_info_box(title_with_icon("🔎 Atom Assessment", AA_TIP), content=a_ass_text, color=a_ass_color)

        struct_text, struct_color = structural_assessment(input_reaction)
        show_info_box(title_with_icon("🧬 Structural Attributes Analysis", SAA_TIP), content=struct_text, color=struct_color)

    with col2:
        score_text, color = waste_efficiency(E_result)
        show_info_box(title_with_icon("🚮 E-Factor: Waste Efficiency", E_FACTOR_TIP), content=score_text, color=color)

        pmi_text, pmi_color = PMI_assessment(PMI_result)
        show_info_box(title_with_icon("🅿️ PMI",PMI_TIP), content=pmi_text, color=pmi_color)

        try:
            atom_economy_m_result = input_reaction.Atom_Economy_M()
            atom_economy_a_result = input_reaction.Atom_Economy_A()
            ae_m_text, ae_m_color = Atom_ec_m_assessment(atom_economy_m_result)
            ae_a_text, ae_a_color = Atom_ec_assessment(atom_economy_a_result)

            dual_metric_box(title_with_icon("⚖️ Atom Econ. (Molar Mass)", AE_M_TIP), ae_m_text, ae_m_color,
                title_with_icon("⚛️ Atom Econ. (Atom no.)", AE_A_TIP), ae_a_text, ae_a_color)

        except ValueError as e:
            st.error(f"{e}", icon="🚨")

else:
    col1, col2 = st.columns(2)

    with col1:
        if st.session_state.extras:
            solv_text, solv_color = get_solvent_info(st.session_state.extras)
            show_info_box(title_with_icon("🌱 Solvent Quality Check", SOLVENT_TIP), content=solv_text, color=solv_color)
        else:
            show_info_box(title_with_icon("🌱 Solvent Quality Check", SOLVENT_TIP), content="Seems like no solvent was used. Nice Job!")
            
        show_info_box(title_with_icon("🌍 Main Product LogP", LOGP_TIP), content="Add a Reaction!")
        show_info_box(title_with_icon("🔎 Atom Assessment", AA_TIP), content="Add a Reaction!")
        show_info_box(title_with_icon("🧬 Structural Attributes Analysis", SAA_TIP), content="Add a Reaction!")

    with col2:
        show_info_box(title_with_icon("🚮 E-Factor: Waste Efficiency", E_FACTOR_TIP), content="Add a Reaction!")
        show_info_box(title_with_icon("🅿️ PMI", PMI_TIP),content="Add a Reaction!")
        dual_metric_box(title_with_icon("⚖️ Atom Econ. (Molar Mass)", AE_M_TIP),"Add a Reaction!","#f9f9f9",
            title_with_icon("⚛️ Atom Econ. (Atom no.)", AE_A_TIP), "Add a Reaction!", "#f9f9f9", "#28a745","#28a745","#000")