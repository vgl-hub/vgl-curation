#!/bin/bash

# Exit on error
set -e

echo "Creating conda environment 'curationenv'..."
conda create -y -n curationenv python=3.10

eval "$(conda shell.bash hook)"

conda activate curationenv

echo "Installing Python dependencies..."
conda install -y biopython pandas natsort gfastats mashmap

echo "Installation complete. Activate the environment with: conda activate curationenv"
