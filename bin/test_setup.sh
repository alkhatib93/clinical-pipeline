#!/bin/bash

# Test script for clinical pipeline setup
# This script checks if all required tools are available

echo "=== Clinical Pipeline Setup Test ==="
echo

# Check Nextflow
echo "Checking Nextflow..."
if command -v nextflow &> /dev/null; then
    echo "✓ Nextflow is installed: $(nextflow -version | head -n1)"
else
    echo "✗ Nextflow is not installed or not in PATH"
    exit 1
fi

echo

# Check required tools
echo "Checking required tools..."

tools=("bcftools" "bedtools" "samtools" "tabix")

for tool in "${tools[@]}"; do
    if command -v $tool &> /dev/null; then
        version=$($tool --version 2>/dev/null | head -n1 || echo "version unknown")
        echo "✓ $tool is installed: $version"
    else
        echo "✗ $tool is not installed or not in PATH"
    fi
done

echo

# Check input files
echo "Checking input files..."

if [ -f "data/samplesheet.csv" ]; then
    echo "✓ samplesheet.csv found"
    echo "  Sample count: $(tail -n +2 data/samplesheet.csv | wc -l)"
else
    echo "✗ data/samplesheet.csv not found"
fi

if [ -f "../BED_files/merged_output.bed" ]; then
    echo "✓ merged_output.bed found"
    echo "  Region count: $(wc -l < ../BED_files/merged_output.bed)"
else
    echo "✗ ../BED_files/merged_output.bed not found"
fi

echo

# Check pipeline files
echo "Checking pipeline files..."

if [ -f "main.nf" ]; then
    echo "✓ main.nf found"
else
    echo "✗ main.nf not found"
fi

if [ -f "nextflow.config" ]; then
    echo "✓ nextflow.config found"
else
    echo "✗ nextflow.config not found"
fi

echo

echo "=== Setup Test Complete ==="
echo
echo "To run the pipeline:"
echo "  nextflow run main.nf"
echo
echo "To run with custom parameters:"
echo "  nextflow run main.nf --samplesheet data/my_samples.csv --bed data/my_regions.bed" 