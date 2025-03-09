#!/bin/bash

# Function to install libGL on Ubuntu/Debian
install_ubuntu() {
    echo "Detected Ubuntu/Debian-based system."
    sudo apt update
    sudo apt install -y libgl1 mesa-utils
}

# Function to install libGL on Fedora/RHEL
install_fedora() {
    echo "Detected Fedora/RHEL-based system."
    sudo dnf install -y mesa-libGL
}

# Function to install libGL on Arch Linux
install_arch() {
    echo "Detected Arch Linux-based system."
    sudo pacman -S --noconfirm libglvnd
}

# Detect the operating system
if command -v apt &> /dev/null; then
    install_ubuntu
elif command -v dnf &> /dev/null; then
    install_fedora
elif command -v pacman &> /dev/null; then
    install_arch
else
    echo "Unsupported operating system. Please install libGL manually."
    exit 1
fi

# Verify installation
if ldconfig -p | grep libGL.so.1 &> /dev/null; then
    echo "libGL.so.1 installed successfully."
else
    echo "Failed to install libGL.so.1. Please check your system."
    exit 1
fi
