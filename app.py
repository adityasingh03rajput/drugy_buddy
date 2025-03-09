import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QMessageBox
from PyQt5.QtCore import Qt
import requests
import pandas as pd
import matplotlib.pyplot as plt
from config import KAGGLE_USERNAME, KAGGLE_KEY, OPENAI_API_KEY

# Kaggle API setup
import os
os.environ['KAGGLE_USERNAME'] = KAGGLE_USERNAME
os.environ['KAGGLE_KEY'] = KAGGLE_KEY

# Main Application Window
class DrugDiscoveryApp(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        # Set up the layout
        layout = QVBoxLayout()

        # Input field for drug requirement
        self.label = QLabel("Enter Drug Requirement:")
        layout.addWidget(self.label)

        self.entry = QLineEdit()
        layout.addWidget(self.entry)

        # Button to trigger analysis
        self.button = QPushButton("Submit")
        self.button.clicked.connect(self.on_submit)
        layout.addWidget(self.button)

        # Set the layout to the window
        self.setLayout(layout)
        self.setWindowTitle("Drug Discovery Lite")
        self.setGeometry(300, 300, 400, 200)

    def on_submit(self):
        requirement = self.entry.text()
        if not requirement:
            QMessageBox.critical(self, "Error", "Please enter a requirement.")
            return

        # Fetch data and generate structure
        try:
            result = self.analyze_requirement(requirement)
            QMessageBox.information(self, "Result", result)
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def analyze_requirement(self, requirement):
        # Step 1: Fetch data from Kaggle (example)
        data = self.fetch_kaggle_data()
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
        self.draw_molecular_structure()

        return f"Analysis complete. Suggested structure for: {requirement}"

    def fetch_kaggle_data(self):
        # Use Kaggle API to fetch dataset (example)
        # Replace with actual Kaggle API integration
        from kaggle.api.kaggle_api_extended import KaggleApi
        api = KaggleApi()
        api.authenticate()
        api.dataset_download_files("username/dataset-name", path="./data", unzip=True)
        return pd.read_csv("./data/example_dataset.csv")

    def draw_molecular_structure(self):
        # Example: Draw a simple structure using Matplotlib
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "C9H8O4", fontsize=12, ha="center")  # Example: Aspirin formula
        ax.axis("off")
        plt.show()

# Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DrugDiscoveryApp()
    window.show()
    sys.exit(app.exec_())
