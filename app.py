import streamlit as st
import requests
import torch
import torchvision.transforms as transforms
from PIL import Image
import timm
from io import BytesIO
import pubchempy as pcp  # For fetching molecular data
import numpy as np

# Load pre-trained model (EfficientNet) with caching
@st.cache_resource
def load_model():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = timm.create_model("efficientnet_b3", pretrained=True, num_classes=1000)
    model.eval().to(device)
    return model, device

model, device = load_model()

# Image transformation for model input
transform = transforms.Compose([
    transforms.Resize((224, 224)),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])

# Function to classify uploaded image
def classify_image(image):
    img = transform(image).unsqueeze(0).to(device)
    with torch.no_grad():
        output = model(img)
    return output.argmax().item()

# Function to fetch molecular data from PubChem
def fetch_molecular_data(drug_name):
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        if compounds:
            compound = compounds[0]
            return {
                "name": compound.iupac_name,
                "formula": compound.molecular_formula,
                "weight": compound.molecular_weight,
                "smiles": compound.canonical_smiles,
                "image_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={compound.cid}&t=l"
            }
    except Exception as e:
        st.error(f"Error fetching molecular data: {e}")
    return None

# Function to fetch drug images from PubChem
@st.cache_data
def fetch_pubchem_image(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={drug_name}&t=l"
    response = requests.get(url)
    if response.status_code == 200:
        return Image.open(BytesIO(response.content))
    return None

# Streamlit App UI
st.title("🔬 Advanced Drug Detection AI using Microscopic Images")

# Test Form for Scientist Input
with st.form("scientist_input_form"):
    st.write("🔬 Please fill out the following details to help with drug detection:")
    drug_name = st.text_input("Drug Name (e.g., Aspirin, Water)")
    drug_type = st.selectbox("Drug Type", ["Analgesic", "Antibiotic", "Antiviral", "Antifungal", "Other"])
    color = st.color_picker("Color of the Drug")
    solubility = st.selectbox("Solubility", ["Water Soluble", "Fat Soluble", "Insoluble"])
    submitted = st.form_submit_button("Submit Information")

    if submitted:
        st.write(f"Drug Name: {drug_name}")
        st.write(f"Drug Type: {drug_type}")
        st.write(f"Color: {color}")
        st.write(f"Solubility: {solubility}")

# Upload Image
uploaded_file = st.file_uploader("📤 Upload a microscopic image", type=["jpg", "png", "jpeg"])

# Placeholder for real-time output
output_placeholder = st.empty()

if uploaded_file is not None:
    # Display Uploaded Image
    user_image = Image.open(uploaded_file).convert("RGB")
    st.image(user_image, caption="Uploaded Microscopic Image", use_column_width=True)

    # Process Uploaded Image in real-time
    with st.spinner("⏳ Processing Image..."):
        user_prediction = classify_image(user_image)
        output_placeholder.write(f"🔍 Predicted Class: {user_prediction}")

    # Fetch molecular data based on the test form input
    if submitted and drug_name:
        with st.spinner("🔍 Fetching molecular data..."):
            molecular_data = fetch_molecular_data(drug_name)
            if molecular_data:
                st.write("📊 Molecular Data:")
                st.write(f"Name: {molecular_data['name']}")
                st.write(f"Formula: {molecular_data['formula']}")
                st.write(f"Weight: {molecular_data['weight']}")
                st.write(f"SMILES: {molecular_data['smiles']}")

                # Display molecular image
                molecular_image = fetch_pubchem_image(molecular_data['name'])
                if molecular_image:
                    st.image(molecular_image, caption=f"Molecular Structure: {molecular_data['name']}", use_column_width=True)

    # Try matching with known drug images
    with st.spinner("🔍 Searching online databases for a match..."):
        drug_names = ["Aspirin", "Ibuprofen", "Paracetamol", "Morphine", "Water"]  # Example drug list
        found_match = False

        for drug in drug_names:
            online_image = fetch_pubchem_image(drug)
            if online_image:
                online_prediction = classify_image(online_image)
                if user_prediction == online_prediction:
                    output_placeholder.success(f"✅ Drug Detected: **{drug}**")
                    st.image(online_image, caption=f"Matched with: {drug}", use_column_width=True)
                    found_match = True
                    break

        if not found_match:
            output_placeholder.warning("❌ No known drug detected in the image.")
