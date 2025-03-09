#!/bin/bash

# Install system dependencies
sudo apt-get update
sudo apt-get install -y \
    libxrender1 \
    libsm6 \
    libxext6 \
    libgl1 \
    python3-dev \
    python3-pip \
    cmake \
    git

# Clone RDKit repository
git clone https://github.com/rdkit/rdkit.git
cd rdkit

# Build and install RDKit
mkdir build
cd build
cmake ..
make
sudo make install

# Install Python wrapper
cd ../Code/Python
python3 setup.py install
