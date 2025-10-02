#!/bin/bash
#Convert DICOMs -> NIfTI for a single patient

set -euo pipefail

PATIENT_DIR="${1:-}"

if [ -z "$PATIENT_DIR" ]; then
  echo "Usage: $0 /path/to/PatientFolder"
  exit 1
fi

if ! command -v dcm2niix >/dev/null 2>&1; then
  echo "Error: dcm2niix not found in PATH"
  exit 1
fi

for VISIT_DIR in "$PATIENT_DIR"/*; do
  if [ -d "$VISIT_DIR/DICOMs" ]; then
    echo "Converting: $VISIT_DIR/DICOMs"
    dcm2niix -z y -o "$VISIT_DIR" "$VISIT_DIR/DICOMs"
  fi
done

echo "Completed."
