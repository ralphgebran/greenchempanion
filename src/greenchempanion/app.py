import streamlit as st
import pandas as pd
import math
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from streamlit_ketcher import st_ketcher
from functions import Atom_Count_With_H, Reaction, compute_PMI, canonicalize_smiles, compute_E 
from functions import get_solvent_info, waste_efficiency, PMI_assesment, Atom_ec_assesment, logP_assessment_molecule, atoms_assessment, structural_assessment


st.set_page_config(page_title="GreenChemPanion", page_icon="../../assets/logo.ico", layout= "wide")

col1, col2 = st.columns([2, 13])

with col1:
    st.image("../../assets/logo_icon.png", width=1000)

with col2:
    st.markdown("## GCP: GreenChemPanion\nInteractive Streamlit Applet showcasing the functions for GCP!")

st.markdown(
    """
    <style>
        .tooltip      {display:inline-block; position:relative; font-size:15px;}
        .tooltiptext  {
            visibility:hidden; width:340px; background:#f9f9f9; color:#333;
            text-align:left; padding:6px 12px; border:2px solid #ccc;
            border-radius:6px; position:absolute; z-index:1;
            top:125%; left:50%; transform: translateX(-50%); box-shadow:0 2px 8px rgba(0,0,0,0.15);
            opacity:0; transition:opacity .25s;
        }
        .tooltip:hover .tooltiptext {visibility:visible; opacity:1;}
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <style>
    div.row-widget.stButton { 
        margin-bottom: 0rem
        margin-top: 0rem
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <style>
    /* Shrink the vertical gap that Streamlit puts after every columns row */
    div[data-testid="stHorizontalBlock"] {
        margin-bottom: 0.5rem !important;   /* 0.5 rem ≈ 8 px – adjust to taste */
    }
    </style>
    """,
    unsafe_allow_html=True,
)

#st.title("GCP: GreenChemPanion", anchor= False) #TITLE
#st.write("Interactive Streamlit Applet showcasing the functions for GCP!") 
st.markdown("---")

st.markdown(
    """
    <style>
      /* limit the main block container to 1100px instead of full width */
      .block-container {
        max-width: 1100px;
      }
    </style>
    """,
    unsafe_allow_html=True,
)

st.header("💱 SMILES Converters", anchor="Converters")

# SMILES TO MOLECULE IMAGE DISPLAY

with st.expander("😊 SMILES to Molecule Converter"):
    M2S_smiles = st.text_input("Enter a SMILES:", key ="Mol2Smiles_input")

    if M2S_smiles:
        #Developers Easter Egg
        if M2S_smiles == "VRMT":
            st.write("Easter Egg :o")
            st.image("../../assets/vrmt.png")


        elif Chem.MolFromSmiles(M2S_smiles) is None:
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

    st.subheader("Draw Molecule and Get SMILES", anchor= False)
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
        st.info("🧪 Draw a molecule above to generate SMILES automatically!")
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

#ADD MOLECULE EXPANDER

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


#COMMON SOLVENT EXPANDER

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


#SOLVENTS EXPANDER

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
with column2:
    st.subheader("🧮 Compute Factors", anchor="compute-factors")
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
            st.markdown(
                """
                <span class="tooltip">ℹ️
                    <span class="tooltiptext">
                        <b>Atom Economy</b> is the percentage of the mass (or atoms)
                        of reactants that ends up in the desired product.
                        Higher values (maximum&nbsp;100&nbsp;%) indicate a greener reaction.
                    </span>
                </span>
                """,
                unsafe_allow_html=True,
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
        st.markdown(
        """
        <style>
        .tooltip-wrapper {
            position: relative;
            display: inline-block;
            width: 100%;
        }
        
        .tooltip-wrapper .tooltip-bubble {
            display: none;
            position: absolute;
            top: 50%;
            left: 100%;
            transform: translateY(-50%);
            margin-left: 12px;
            background-color: white;
            color: black;
            padding: 10px;
            border-radius: 6px;
            font-size: 13px;
            line-height: 1.4;
            width: 260px;
            text-align: justify;
            box-shadow: 0 0 10px rgba(0,0,0,0.2);
            z-index: 10;
        }

        .tooltip-wrapper:hover .tooltip-bubble {
            display: block;
        }
        </style>

        <div class="tooltip-wrapper">
            <div style="
                background-color: rgba(188, 212, 180, 0.3);
                color: #28a745;
                padding: 0.75rem 1rem;
                border-radius: 0.25rem;
                font-size: 18px;
                margin-bottom: 40px;
                box-shadow:0 2px 8px rgba(0,0,0,.15);
            ">
                Add at least one reactant and one product to compute Atom Economy.
            </div>
            <div class="tooltip-bubble">
                <b>Atom Economy</b> is the percentage of the mass (or atoms) of reactants that ends up in the desired product. Higher values (max 100%) indicate a greener reaction.
            </div>
        </div>
        """,
        unsafe_allow_html=True,
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
            st.markdown(
                """
                <span class="tooltip">ℹ️
                    <span class="tooltiptext">
                        <b>Process Mass Intensity (PMI) </b> is the total mass of all materials
                        used in the process divided by the mass of product obtained.
                        Lower values (maximum&nbsp;100&nbsp;%) indicate a greener, more material-efficient process.
                    </span>
                </span>
                """,
                unsafe_allow_html=True,
            )
        
        if compute_clicked_2:
            st.success(f"🅿️ PMI: **{PMI_result:.2f}**")

    else:
        st.markdown(
        """
        <div class="tooltip-wrapper">
            <div style="
                background-color: rgba(188, 212, 180, 0.3);
                color: #28a745;
                padding: 0.75rem 1rem;
                border-radius: 0.25rem;
                font-size: 18px;
                margin-bottom: 40px;
                box-shadow:0 2px 8px rgba(0,0,0,.15);
            ">
                Add at least one reactant and one product to compute PMI.
            </div>
            <div class="tooltip-bubble">
                <b>Process Mass Intensity (PMI)</b> is the total mass of all materials used in the
                process divided by the mass of product obtained. Lower values (max 100%) indicate
                a greener, more material-efficient process.
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


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
            st.markdown(
                """
                <span class="tooltip">ℹ️
                    <span class="tooltiptext">
                        <b>E factor (Environmental factor)</b> is the ratio of the total mass of all waste generated 
                        (side-products, solvents, auxiliaries, etc.) to the mass of desired product obtained.
                        Lower values indicate a greener, more sustainable process,
                        as less waste is produced per kilogram of product.
                    </span>
                </span>
                """,
                unsafe_allow_html=True,
            )
        if compute_e_clicked :
            st.success(f"🚮 E: **{E_result:.2f}**")

    else:
        st.markdown(
        """
        <div class="tooltip-wrapper">
            <div style="
                background-color: rgba(188, 212, 180, 0.3);
                color: #28a745;
                padding: 0.75rem 1rem;
                border-radius: 0.25rem;
                font-size: 18px;
                margin-bottom: 40px;
                box-shadow:0 2px 8px rgba(0,0,0,.15);
            ">
                Add at least one reactant and one product to compute E factor.
            </div>
            <div class="tooltip-bubble">
                <b>E factor (Environmental factor)</b> is the ratio of the total mass of all waste generated
                (side-products, solvents, auxiliaries, etc.) to the mass of desired product obtained.
                Lower values indicate a greener, more sustainable process, as less waste is produced
                per kilogram of product.
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

# GCP GREEN CHEMISTRY EVALUATION 

st.header("🧑‍🔬 GCP Green Chemistry Evaluation", anchor="GCP_evaluation")

box_style = """
    <div style="
        border: 2px solid #28a745;
        border-radius: 10px;
        padding: 20px;
        margin: 10px;
        background-color: #f9f9f9;
        text-align: center;
    ">
        <h3 style="color: #28a745;">{title}</h3>
        <p style="color: #000000;">{content}</p>
    </div>
"""

box_style_2 = """
<div style="
    border: 2px solid #000000;
    border-radius: 10px;
    padding: 20px;
    margin: 10px;
    background-color: {color};
    text-align: center;
">
    <h3 style="color: #000000;">{title}</h3>
    <p>{content}</p>
</div>
"""

if st.session_state.reactants and st.session_state.products:
    col1, col2 = st.columns(2)

    with col1:
        if st.session_state.extras:
            solv_text, solv_color = get_solvent_info(st.session_state.extras)
            st.markdown(box_style_2.format(title="🌱 Solvent Quality Check", content=solv_text, color=solv_color), unsafe_allow_html=True)
        else:
            st.markdown(box_style.format(title="🌱 Solvent Quality Check", content="Seems like no solvent was used. Nice Job!"), unsafe_allow_html=True)

        Log_P = Descriptors.MolLogP(input_reaction.main_product)
        log_text, log_color = logP_assessment_molecule(Log_P)
        st.markdown(box_style_2.format(title="🌍 Main Product LogP", content=log_text, color=log_color), unsafe_allow_html=True)

        a_ass_text, a_ass_color = atoms_assessment(input_reaction)
        st.markdown(box_style_2.format(title="🔎 Atom Assessment", content=a_ass_text, color=a_ass_color), unsafe_allow_html=True)

        struct_text, struct_color = structural_assessment(input_reaction)
        st.markdown(box_style_2.format(title="🧬 Structural Attributes Analysis", content=struct_text, color=struct_color), unsafe_allow_html=True)

    with col2:
        score_text, color = waste_efficiency(E_result)
        st.markdown(box_style_2.format(title="🚮 E-Factor: Waste Efficiency", content=score_text, color=color), unsafe_allow_html=True)

        pmi_text, pmi_color = PMI_assesment(PMI_result)
        st.markdown(box_style_2.format(title="🅿️ PMI", content=pmi_text, color=pmi_color), unsafe_allow_html=True)

        try:
            atom_economy_m_result = input_reaction.Atom_Economy_M()
            atom_economy_a_result = input_reaction.Atom_Economy_A()
            ae_m_text, ae_m_color = Atom_ec_assesment(atom_economy_m_result)
            ae_a_text, ae_a_color = Atom_ec_assesment(atom_economy_a_result)

            html = f"""<div style="
display: flex;
border: 2px solid #000000;
border-radius: 10px;
overflow: hidden;
margin: 10px;
">
    <div style="
        flex: 1;
        background-color: {ae_m_color};
        padding: 20px;
        text-align: center;
        border-right: 1px solid #000;
    ">
        <h4 style="margin:0; color:#000000;">⚖️ Atom Econ. (Molar Mass)</h4>
        <p style="margin:5px 0 0;">{ae_m_text}</p>
    </div>
    <div style="
        flex: 1;
        background-color: {ae_a_color};
        padding: 20px;
        text-align: center;
    ">
        <h4 style="margin:0; color:#000000;">⚛️ Atom Econ. (Atom no.)</h4>
        <p style="margin:5px 0 0;">{ae_a_text}</p>
    </div>
</div>"""
            st.markdown(html, unsafe_allow_html=True)
        except ValueError as e:
            st.error(f"{e}", icon="🚨")

else:
    col1, col2 = st.columns(2)

    with col1:
        if st.session_state.extras:
            solv_text, solv_color = get_solvent_info(st.session_state.extras)
            st.markdown(box_style_2.format(title="🌱 Solvent Quality Check", content=solv_text, color=solv_color), unsafe_allow_html=True)
        else:
            st.markdown(box_style.format(title="🌱 Solvent Quality Check", content="Seems like no solvent was used. Nice Job!"), unsafe_allow_html=True)
            
        st.markdown(box_style.format(title="🌍 Main Product LogP", content="Add a Reaction!"), unsafe_allow_html=True)
        st.markdown(box_style.format(title="🔎 Atom Assessment", content="Add a Reaction!"), unsafe_allow_html=True)
        st.markdown(box_style.format(title="🧬 Structural Attributes Analysis", content="Add a Reaction!"), unsafe_allow_html=True)

    with col2:
        st.markdown(box_style.format(title="🚮 E-Factor: Waste Efficiency", content="Add a Reaction!"), unsafe_allow_html=True)
        st.markdown(box_style.format(title="🅿️ PMI", content="Add a Reaction!"), unsafe_allow_html=True)
        st.markdown(box_style.format(title="⚖️ Atom Econ. (Molar Mass)", content="Add a Reaction!"), unsafe_allow_html=True)
        st.markdown(box_style.format(title="⚛️ Atom Econ. (Atom no.)", content="Add a Reaction!"), unsafe_allow_html=True)