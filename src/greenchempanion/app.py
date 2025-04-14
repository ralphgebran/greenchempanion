import streamlit as st
import numpy
from rdkit import Chem
from rdkit.Chem import Draw

st.set_page_config(page_title="GreenChemPanion", page_icon="🍃")
st.title("GCP: GreenChemPanion") #TITLE
st.write("Hello Chicho! This is a simple Streamlit app.") #TEXT
name = st.text_input("Enter your name:") #TEXT INPUT
age = st.slider("Select your age", 1, 100, 25, 1) #SLIDER (LABEL, MIN, MAX, VALUE, STEP)

# SMILES TO MOLECULE IMAGE DISPLAY
st.subheader("SMILES TO MOLECULE CONVERTER !⚗️😮")
input_smiles = st.text_input("Enter a SMILES:")

#Developers Easter Egg
if input_smiles == "VRMT":
    st.write("Easter Egg :o")
    st.image("../../assets/banner.png")

elif Chem.MolFromSmiles(input_smiles) is None:
    st.write("⚠️ Enter a valid SMILES format [e.g.: CCO]")

else:
    input_mol = Chem.MolFromSmiles(input_smiles)
    img = Draw.MolToImage(input_mol)
    st.image(img)