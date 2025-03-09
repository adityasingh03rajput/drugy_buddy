import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Streamlit app
st.title("Drug Discovery Assistant")
st.write("Generate and visualize molecules using RDKit.")

# Sidebar for navigation
st.sidebar.title("Navigation")
option = st.sidebar.selectbox(
    "Choose a task:",
    ["Molecule Generation", "About"]
)

# Molecule Generation
if option == "Molecule Generation":
    st.subheader("Molecule Generation")
    smiles_input = st.text_input("Enter a SMILES string (e.g., CCO):")
    
    if smiles_input:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smiles_input)
        if mol:
            # Display molecule image
            st.write("Molecule Visualization:")
            img = Draw.MolToImage(mol)
            st.image(img, caption=f"SMILES: {smiles_input}", use_column_width=True)
            
            # Display basic properties
            st.write("Basic Properties:")
            mol_weight = Chem.Descriptors.MolWt(mol)
            logp = Chem.Descriptors.MolLogP(mol)
            st.write(f"Molecular Weight: {mol_weight:.2f}")
            st.write(f"LogP: {logp:.2f}")
        else:
            st.error("Invalid SMILES string. Please try again.")

# About section
elif option == "About":
    st.subheader("About")
    st.write("This is a basic Drug Discovery Assistant app built using Streamlit and RDKit.")
    st.write("It allows you to generate and visualize molecules using SMILES strings.")
