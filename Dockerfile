FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libsm6 \
    libxext6 \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the app code
COPY . .

# Run the app
CMD ["streamlit", "run", "app.py"]
