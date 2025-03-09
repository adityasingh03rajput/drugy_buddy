import tkinter as tk
from tkinter import messagebox
import requests
import pandas as pd
import matplotlib.pyplot as plt
from config import KAGGLE_USERNAME, KAGGLE_KEY, OPENAI_API_KEY

# Kaggle API setup
import os
os.environ['KAGGLE_USERNAME'] = KAGGLE_USERNAME
os.environ['KAGGLE_KEY'] = KAGGLE_KEY

# GUI Setup
def create_gui():
    root = tk.Tk()
    root.title("Drug Discovery Lite")

    # Input field for drug requirement
    tk.Label(root, text="Enter Drug Requirement:").pack()
    entry = tk.Entry(root, width=50)
    entry.pack()

    # Button to trigger analysis
    def on_submit():
        requirement = entry.get()
        if not requirement:
            messagebox.showerror("Error", "Please enter a requirement.")
            return

        # Fetch data and generate structure
        try:
            result = analyze_requirement(requirement)
            messagebox.showinfo("Result", result)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    tk.Button(root, text="Submit", command=on_submit).pack()

    root.mainloop()

# Analyze requirement and generate structure
def analyze_requirement(requirement):
    # Step 1: Fetch data from Kaggle (example)
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

    # Step 3: Generate and draw molecular structure (example)
    draw_molecular_structure()

    return f"Analysis complete. Suggested structure for: {requirement}"

# Fetch data from Kaggle (example)
def fetch_kaggle_data():
    # Use Kaggle API to fetch dataset (example)
    # Replace with actual Kaggle API integration
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi()
    api.authenticate()
    api.dataset_download_files("username/dataset-name", path="./data", unzip=True)
    return pd.read_csv("./data/example_dataset.csv")

# Draw molecular structure (example)
def draw_molecular_structure():
    # Example: Draw a simple structure using Matplotlib
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, "C9H8O4", fontsize=12, ha="center")  # Example: Aspirin formula
    ax.axis("off")
    plt.show()

if __name__ == "__main__":
    create_gui()
