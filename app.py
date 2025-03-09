import streamlit as at
import sys
import turtle
import requests
import pandas as pd
import matplotlib.pyplot as plt
import os
from config import KAGGLE_USERNAME, KAGGLE_KEY, OPENAI_API_KEY

# ✅ Correctly setting environment variables
os.environ["KAGGLE_USERNAME"] = KAGGLE_USERNAME
os.environ["KAGGLE_KEY"] = KAGGLE_KEY

# ✅ Function to analyze drug requirement
def analyze_requirement(requirement):
    # Step 1: Fetch data from Kaggle
    data = fetch_kaggle_data()
    print("Fetched data:", data.head())

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
    print("Analysis:", analysis)

    # Step 3: Draw molecular structure using Turtle
    draw_molecular_structure()

    return f"Analysis complete. Suggested structure for: {requirement}"

# ✅ Function to fetch data from Kaggle
def fetch_kaggle_data():
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi()
    api.authenticate()
    
    # 🔹 Replace this with a real dataset name from Kaggle
    dataset_name = "zillow/zecon"  # Example dataset
    api.dataset_download_files(dataset_name, path="./data", unzip=True)
    
    return pd.read_csv("./data/example_dataset.csv")

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
at.title("Drug Discovery Lite")

requirement = at.text_input("Enter Drug Requirement:")
if at.button("Submit"):
    if not requirement:
        at.error("Please enter a requirement.")
    else:
        try:
            result = analyze_requirement(requirement)
            at.success(result)
        except Exception as e:
            at.error(str(e))
