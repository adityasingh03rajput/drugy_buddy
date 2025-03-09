import os
import pandas as pd
import requests
import tkinter as tk
from tkinter import messagebox, scrolledtext

# Kaggle API setup
def download_kaggle_dataset(dataset_name, path="data"):
    os.makedirs(path, exist_ok=True)
    os.system(f"kaggle datasets download -d {dataset_name} -p {path} --unzip")

# Load Kaggle dataset
def load_dataset():
    dataset_name = "adityasingh03rajput/drug-discovery-dataset"  # Replace with your Kaggle dataset
    download_kaggle_dataset(dataset_name)
    return pd.read_csv("data/drugs.csv")  # Replace with your dataset file name

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

# Tkinter GUI
class DrugDiscoveryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Drug Discovery Assistant")
        self.root.geometry("600x400")

        # Load dataset
        self.dataset = load_dataset()

        # GUI Components
        self.label = tk.Label(root, text="Enter a drug name:")
        self.label.pack(pady=10)

        self.entry = tk.Entry(root, width=50)
        self.entry.pack(pady=10)

        self.button = tk.Button(root, text="Get Explanation", command=self.show_explanation)
        self.button.pack(pady=10)

        self.text_area = scrolledtext.ScrolledText(root, width=70, height=15)
        self.text_area.pack(pady=10)

    def show_explanation(self):
        drug_name = self.entry.get()
        if drug_name:
            # Find drug in dataset
            drug_info = self.dataset[self.dataset["Drug_Name"] == drug_name]
            if not drug_info.empty:
                # Get GPT explanation
                explanation = get_gpt_explanation(f"Explain the drug {drug_name}.")
                self.text_area.insert(tk.END, f"Drug: {drug_name}\n")
                self.text_area.insert(tk.END, f"Explanation: {explanation}\n\n")
            else:
                messagebox.showerror("Error", "Drug not found in dataset.")
        else:
            messagebox.showerror("Error", "Please enter a drug name.")

# Run the app
if __name__ == "__main__":
    root = tk.Tk()
    app = DrugDiscoveryApp(root)
    root.mainloop()
