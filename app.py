import streamlit as st
import os
import requests
import pandas as pd
import turtle

# ✅ Set Kaggle API credentials
os.environ["KAGGLE_USERNAME"] = "adityasingh03rajput"
os.environ["KAGGLE_KEY"] = "d1423178b21e3f6e2ffa3b782bfd1684"

# ✅ Set OpenAI API key
OPENAI_API_KEY = "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"

# ✅ Function to analyze drug requirement
def analyze_requirement(requirement):
    # Step 1: Fetch data from Kaggle
    data = fetch_kaggle_data()
    st.write("Fetched data:", data.head())

    # Step 2: Use OpenAI API to analyze requirement
    response = requests.post(
        "https://api.openai.com/v1/chat/completions",
        headers={"Authorization": f"Bearer {OPENAI_API_KEY}"},
        json={
            "model": "gpt-3.5-turbo",
            "messages": [{"role": "user", "content": f"Analyze this drug requirement: {requirement}"}]
        }
    )
    if response.status_code != 200:
        raise Exception("Failed to analyze requirement.")

    analysis = response.json()["choices"][0]["message"]["content"]
    st.write("Analysis:", analysis)

    # Step 3: Draw molecular structure using Turtle
    draw_molecular_structure()

    return f"Analysis complete. Suggested structure for: {requirement}"

# ✅ Function to fetch data from Kaggle
def fetch_kaggle_data():
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi()
    api.authenticate()
    
    # 🔹 Replace with a drug-related dataset from Kaggle
    dataset_name = "prathamtripathi/drug-classification"  # Example drug dataset
    api.dataset_download_files(dataset_name, path="./data", unzip=True)
    
    # 🔹 Load the dataset (replace "drug200.csv" with the actual file name)
    return pd.read_csv("./data/drug200.csv")

# ✅ Function to draw a molecular structure using Turtle
def draw_molecular_structure():
    screen = turtle.Screen()
    screen.title("Molecular Structure")
    screen.bgcolor("white")

    pen = turtle.Turtle()
    pen.speed(3)

    # Draw a simple molecular representation
    pen.penup()
    pen.goto(-50, 0)
    pen.pendown()
    pen.circle(20)  # First atom

    pen.penup()
    pen.goto(0, 0)
    pen.pendown()
    pen.forward(50)  # Bond

    pen.penup()
    pen.goto(50, 0)
    pen.pendown()
    pen.circle(20)  # Second atom

    pen.penup()
    pen.goto(50, 30)
    pen.pendown()
    pen.circle(10)  # Small branch

    # Close window after 5 seconds
    screen.ontimer(lambda: screen.bye(), 5000)
    screen.mainloop()

# ✅ Streamlit UI for input
st.title("Drug Discovery Lite")

requirement = st.text_input("Enter Drug Requirement:")
if st.button("Submit"):
    if not requirement:
        st.error("Please enter a requirement.")
    else:
        try:
            result = analyze_requirement(requirement)
            st.success(result)
        except Exception as e:
            st.error(str(e))
