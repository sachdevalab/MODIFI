#!/bin/bash

# Compile the C++ code
echo "Compiling get_control_IPD..."
g++ -O2 -o get_control_IPD get_control_IPD.cpp -pthread

if [ $? -eq 0 ]; then
    echo "✅ Compilation successful!"
    echo "Executable created: get_control_IPD"
else
    echo "❌ Compilation failed!"
    exit 1
fi
