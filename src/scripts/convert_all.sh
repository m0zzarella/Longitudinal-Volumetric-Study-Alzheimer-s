#!/bin/bash
#Convert DICOMs -> NIfTI for all patients under a root folder

set -euo pipefail

ROOT_DIR="${1:-}"

if [ -z "$ROOT_DIR" ]; then
  echo "Usage: $0 /path/to/root"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

for PATIENT in "$ROOT_DIR"/*; do
  if [ -d "$PATIENT" ]; then
    echo "==> Patient: $PATIENT"
    "$SCRIPT_DIR/convert_single_patient.sh" "$PATIENT"
  fi
done

echo "Completed."
