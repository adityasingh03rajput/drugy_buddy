import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from deepchem.models import GraphConvModel
import pubchempy as pcp
from PIL import Image
import os

# Streamlit app
st.title("Drug Discovery Assistant")
st.write("Generate molecules, predict properties, and analyze drug-like compounds.")

# Sidebar for navigation
st.sidebar.title("Navigation")
option = st.sidebar.selectbox(
    "Choose a task:",
    ["Molecule Generation", "Virtual Screening", "ADMET Prediction", "Drug Search"],
)

# Function to generate molecules using RDKit
def generate_molecules():
    st.subheader("Molecule Generation")
    num_molecules = st.slider("Number of molecules to generate", 1, 10, 5)
    if st.button("Generate Molecules"):
        smiles_list = ["CCO", "CCN", "C1=CC=CC=C1", "C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]  # Example SMILES
        molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list[:num_molecules]]
        for i, mol in enumerate(molecules):
            st.write(f"Molecule {i+1}:")
            st.image(Draw.MolToImage(mol), caption=f"SMILES: {smiles_list[i]}", use_column_width=True)

# Function for virtual screening using DeepChem
def virtual_screening():
    st.subheader("Virtual Screening")
    uploaded_file = st.file_uploader("Upload a CSV file with SMILES strings", type=["csv"])
    if uploaded_file is not None:
        data = pd.read_csv(uploaded_file)
        st.write("Uploaded Data:")
        st.write(data.head())

        if st.button("Run Virtual Screening"):
            # Load a pre-trained DeepChem model
            model = GraphConvModel(n_tasks=1, mode="regression")
            # Example: Predict solubility (replace with your model)
            predictions = model.predict(data["SMILES"].tolist())
            data["Predicted_Solubility"] = predictions
            st.write("Screening Results:")
            st.write(data)

# Function for ADMET prediction
def admet_prediction():
    st.subheader("ADMET Prediction")
    smiles_input = st.text_input("Enter a SMILES string (e.g., CCO):")
    if smiles_input:
        mol = Chem.MolFromSmiles(smiles_input)
        if mol:
            st.image(Draw.MolToImage(mol), caption=f"SMILES: {smiles_input}", use_column_width=True)
            # Example: Calculate molecular weight and logP
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            st.write(f"Molecular Weight: {mw:.2f}")
            st.write(f"LogP: {logp:.2f}")
        else:
            st.error("Invalid SMILES string.")

# Function to search for drugs using PubChem
def drug_search():
    st.subheader("Drug Search")
    drug_name = st.text_input("Enter the name of the drug (e.g., Aspirin):")
    if drug_name:
        compounds = pcp.get_compounds(drug_name, "name")
        if compounds:
            compound = compounds[0]
            st.success(f"Found: {compound.iupac_name}")
            st.write(f"Molecular Formula: {compound.molecular_formula}")
            st.write(f"Molecular Weight: {compound.molecular_weight}")
            st.write(f"Canonical SMILES: {compound.canonical_smiles}")
        else:
            st.error("No matching drug found.")

# Main app logic
if option == "Molecule Generation":
    generate_molecules()
elif option == "Virtual Screening":
    virtual_screening()
elif option == "ADMET Prediction":
    admet_prediction()
elif option == "Drug Search":
    drug_search()
