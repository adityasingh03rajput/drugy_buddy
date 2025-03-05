import streamlit as st
import requests
import torch
import torchvision.transforms as transforms
from PIL import Image
import timm
from io import BytesIO

# Load pre-trained model (EfficientNet)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = timm.create_model("efficientnet_b3", pretrained=True, num_classes=1000)
model.eval().to(device)

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

# Function to fetch drug images from PubChem
def fetch_pubchem_image(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={drug_name}&t=l"
    response = requests.get(url)
    if response.status_code == 200:
        return Image.open(BytesIO(response.content))
    return None

# Streamlit App UI
st.title("🔬 Drug Detection AI using Microscopic Images")

# Upload Image
uploaded_file = st.file_uploader("📤 Upload a microscopic image", type=["jpg", "png", "jpeg"])

if uploaded_file is not None:
    # Display Uploaded Image
    user_image = Image.open(uploaded_file).convert("RGB")
    st.image(user_image, caption="Uploaded Microscopic Image", use_column_width=True)

    # Process Uploaded Image
    st.write("⏳ Processing Image...")
    user_prediction = classify_image(user_image)

    # Try matching with known drug images
    st.write("🔍 Searching online databases for a match...")
    drug_names = ["Aspirin", "Ibuprofen", "Paracetamol", "Morphine"]  # Example drug list

    found_match = False
    for drug in drug_names:
        online_image = fetch_pubchem_image(drug)
        if online_image:
            online_prediction = classify_image(online_image)
            if user_prediction == online_prediction:
                st.success(f"✅ Drug Detected: **{drug}**")
                st.image(online_image, caption=f"Matched with: {drug}", use_column_width=True)
                found_match = True
                break

    if not found_match:
        st.warning("❌ No known drug detected in the image.")
