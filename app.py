import streamlit as st
import tensorflow as tf
import numpy as np
import pubchempy as pcp
from PIL import Image
import os
import zipfile
import kaggle

# Streamlit app
st.title("Drug Detection Machine")
st.write("Upload a microscopic image to analyze or search for a drug by name.")

# Download Kaggle dataset
def download_kaggle_dataset():
    if not os.path.exists("data/drugs"):
        os.makedirs("data/drugs", exist_ok=True)
        kaggle.api.authenticate()
        kaggle.api.dataset_download_files("username/drug-classification-dataset", path="data", unzip=True)
        st.success("Kaggle dataset downloaded successfully!")

# Load a pretrained model (ResNet50)
def load_model():
    model = tf.keras.applications.ResNet50(
        weights="imagenet", include_top=False, pooling="avg"
    )
    return model

# Extract features from an image
def extract_features(image, model):
    image = tf.keras.applications.resnet50.preprocess_input(image)
    image = np.expand_dims(image, axis=0)  # Add batch dimension
    features = model.predict(image)
    return features.flatten()

# Query PubChem for molecular data
def query_pubchem(features):
    try:
        query = " ".join(map(str, features))  # Simplified for example
        compounds = pcp.get_compounds(query, namespace="smiles", limit=1)
        return compounds[0] if compounds else None
    except Exception as e:
        st.error(f"Error querying PubChem: {e}")
        return None

# Preprocess the uploaded image
def preprocess_image(file):
    image = Image.open(file).convert("RGB")
    image = image.resize((224, 224))  # Resize for model input
    return np.array(image)

# Display results
def display_results(compound):
    if compound:
        st.success(f"✅ Identified as {compound.iupac_name}.")
        st.write(f"Molecular Formula: {compound.molecular_formula}")
    else:
        st.error("❌ Not a drug. Likely Water or another substance.")

# Sidebar for drug search
st.sidebar.title("Search Drug by Name")
drug_name = st.sidebar.text_input("Enter the name of the drug (e.g., Aspirin):")
if drug_name:
    st.sidebar.write(f"Searching for {drug_name}...")
    compounds = pcp.get_compounds(drug_name, "name")
    if compounds:
        compound = compounds[0]
        st.sidebar.success(f"Found: {compound.iupac_name}")
        st.sidebar.write(f"Molecular Formula: {compound.molecular_formula}")
    else:
        st.sidebar.error("No matching drug found.")

# Download Kaggle dataset
download_kaggle_dataset()

# File uploader for image analysis
uploaded_file = st.file_uploader("Choose an image...", type=["jpg", "jpeg", "png"])
if uploaded_file is not None:
    # Preprocess the image
    image = preprocess_image(uploaded_file)
    st.image(image, caption="Uploaded Image", use_column_width=True)

    # Load the pretrained model
    model = load_model()

    # Extract features using the pretrained model
    features = extract_features(image, model)

    # Query PubChem for molecular data
    compound = query_pubchem(features)

    # Display results
    display_results(compound)
