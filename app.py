import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import requests
import os
import json

# Streamlit app
st.title("Drug Discovery Assistant")
st.write("Generate optimized molecules, visualize them, and get AI explanations.")

# Kaggle API setup
def download_kaggle_dataset(dataset_name, path="data"):
    os.makedirs(path, exist_ok=True)
    os.system(f"kaggle datasets download -d {dataset_name} -p {path} --unzip")

# Load Kaggle dataset
if st.button("Load Kaggle Dataset"):
    dataset_name = "adityasingh03rajput/drug-discovery-dataset"  # Replace with your Kaggle dataset
    download_kaggle_dataset(dataset_name)
    st.success("Dataset downloaded successfully!")

# Molecule generation and optimization
def generate_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Optimize the molecule
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol)
        return img, mol
    else:
        return None, None

# Mini GPT API for explanations
def get_gpt_explanation(text):
    api_key = "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"  # Replace with your Mini GPT API key
    url = "https://api.mini-gpt.com/v1/explain"  # Replace with the actual API endpoint
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }
    data = {
        "text": text
    }
    response = requests.post(url, headers=headers, json=data)
    if response.status_code == 200:
        return response.json().get("explanation", "No explanation available.")
    else:
        return "Failed to get explanation."

# Input SMILES string
smiles_input = st.text_input("Enter a SMILES string (e.g., CCO for ethanol):")

if smiles_input:
    # Generate and visualize molecule
    img, mol = generate_molecule(smiles_input)
    if img:
        st.write("Optimized Molecule Visualization:")
        st.image(img, caption=f"SMILES: {smiles_input}", use_column_width=True)
        
        # Get AI explanation
        explanation = get_gpt_explanation(f"Explain the molecule with SMILES {smiles_input}.")
        st.write("AI Explanation:")
        st.write(explanation)
    else:
        st.error("Invalid SMILES string. Please try again.")
