#!/bin/bash

# Check if we're in a conda environment
if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: No conda environment detected. Please activate your conda environment first."
    echo "Example: conda activate methy3"
    exit 1
fi

echo "Using conda environment: $CONDA_PREFIX"

# Set up paths using current conda environment
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig/"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# Check if required libraries are available
if ! pkg-config --exists htslib; then
    echo "Error: htslib not found. Please install it with:"
    echo "conda install -c bioconda htslib"
    exit 1
fi

if ! pkg-config --exists pbbam; then
    echo "Error: pbbam not found. Please install it with:"
    echo "conda install -c bioconda pbbam"
    exit 1
fi

# Compile the C++ code
echo "Compiling get_control_IPD..."
g++ -o get_control_IPD get_control_IPD.cpp $(pkg-config --cflags --libs htslib pbbam)

if [ $? -eq 0 ]; then
    echo "✅ Compilation successful!"
    echo "Executable created: get_control_IPD"
else
    echo "❌ Compilation failed!"
    exit 1
fi
