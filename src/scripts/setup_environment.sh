#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
 
if [[ ! -f "$SCRIPT_DIR/../R/03_longitudinal_pipeline.R" ]]; then
    echo "Error: Cannot find R scripts. Make sure you're running from the correct project directory."
    echo "Expected: $SCRIPT_DIR/../R/03_longitudinal_pipeline.R"
    exit 1
fi

echo "=== Longitudinal Neuroanalysis Pipeline Setup ==="

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_status() {
    local status=$1
    local message=$2
    case $status in
        "OK")
            echo -e "${GREEN}✓${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}⚠${NC} $message"
            ;;
        "ERROR")
            echo -e "${RED}✗${NC} $message"
            ;;
    esac
}

echo "Checking system dependencies..."

if command_exists dcm2niix; then
    print_status "OK" "dcm2niix found: $(dcm2niix -h 2>&1 | head -1)"
else
    print_status "ERROR" "dcm2niix not found. Please install dcm2niix."
    exit 1
fi

if [ -n "$FSLDIR" ] && [ -d "$FSLDIR" ]; then
    print_status "OK" "FSL found at: $FSLDIR"
    if command_exists bet; then
        print_status "OK" "FSL tools accessible in PATH"
    else
        print_status "WARN" "FSL directory set but tools not in PATH"
        echo "Add $FSLDIR/bin to your PATH"
    fi
else
    print_status "ERROR" "FSL not found. Please install FSL and set FSLDIR."
    exit 1
fi

if command_exists R; then
    R_VERSION=$(R --version | head -1)
    print_status "OK" "R found: $R_VERSION"
else
    print_status "ERROR" "R not found. Please install R."
    exit 1
fi

echo "Checking R packages..."
REQUIRED_PACKAGES=("oro.dicom" "oro.nifti" "neurobase" "fslr" "scales")

for package in "${REQUIRED_PACKAGES[@]}"; do
    if R -e "library($package)" >/dev/null 2>&1; then
        print_status "OK" "R package '$package' installed"
    else
        print_status "ERROR" "R package '$package' not found"
        echo "Install with: R -e \"install.packages('$package')\""
        if [ "$package" == "neurobase" ] || [ "$package" == "fslr" ]; then
            echo "Or: R -e \"remotes::install_github('muschellij2/$package')\""
        fi
    fi
done

echo "Checking directory structure..."

REQUIRED_DIRS=("data" "results" "src" "config" "logs")

for dir in "${REQUIRED_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        print_status "OK" "Directory '$dir' exists"
    else
        print_status "WARN" "Directory '$dir' not found. Creating..."
        mkdir -p "$dir"
        print_status "OK" "Created directory '$dir'"
    fi
done

echo "Checking script permissions..."

SCRIPTS=("$SCRIPT_DIR/convert_single_patient.sh" "$SCRIPT_DIR/convert_all.sh")

for script in "${SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        if [ -x "$script" ]; then
            print_status "OK" "Script '$script' is executable"
        else
            print_status "WARN" "Script '$script' is not executable. Making executable..."
            chmod +x "$script"
            print_status "OK" "Made '$script' executable"
        fi
    else
        print_status "ERROR" "Script '$script' not found"
    fi
done

echo "Testing basic functionality..."

if dcm2niix -h >/dev/null 2>&1; then
    print_status "OK" "dcm2niix working correctly"
else
    print_status "WARN" "dcm2niix may have issues"
fi

if bet >/dev/null 2>&1; then
    print_status "OK" "FSL BET working correctly"
else
    print_status "WARN" "FSL BET may have issues"
fi

if [ ! -d "data/sample_patient" ]; then
    echo "Creating sample directory structure..."
    mkdir -p data/sample_patient/Visit1/DICOMs
    mkdir -p data/sample_patient/Visit2/DICOMs
    print_status "OK" "Created sample patient directory structure"
fi

echo ""
echo "=== Setup Complete ==="
print_status "OK" "Environment setup finished successfully!"
echo ""
echo "Next steps:"
echo "1. Place your DICOM data in data/PatientX/VisitY/DICOMs/"
echo "2. Run: ./src/scripts/convert_single_patient.sh data/PatientX"
echo "3. Run: Rscript src/R/03_longitudinal_pipeline.R data/PatientX"
echo ""
echo "For batch processing:"
echo "./src/scripts/batch_process.sh -d data/"
echo ""
