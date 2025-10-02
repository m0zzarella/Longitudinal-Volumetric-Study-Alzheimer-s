#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if [[ ! -f "$SCRIPT_DIR/../R/03_longitudinal_pipeline.R" ]]; then
    echo "Error: Cannot find R scripts. Make sure you're running from the correct project directory."
    echo "Expected: $SCRIPT_DIR/../R/03_longitudinal_pipeline.R"
    exit 1
fi

ROOT_DIR=""
MAX_PARALLEL=4
SKIP_CONVERSION=false
SKIP_PIPELINE=false
CONFIG_FILE="config/pipeline_config.yaml"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

usage() {
    echo "Usage: $0 -d DATA_DIR [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -d, --data-dir DIR        Root directory containing patient folders"
    echo "  -j, --jobs N             Maximum number of parallel jobs (default: 4)"
    echo "  -c, --config FILE        Configuration file (default: config/pipeline_config.yaml)"
    echo "  --skip-conversion        Skip DICOM to NIfTI conversion"
    echo "  --skip-pipeline          Skip longitudinal pipeline processing"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -d data/ -j 8"
    echo "  $0 -d /path/to/patients --skip-conversion"
}

log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${BLUE}[INFO]${NC} $timestamp - $message" | tee -a logs/batch_process.log
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} $timestamp - $message" | tee -a logs/batch_process.log
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} $timestamp - $message" | tee -a logs/batch_process.log
            ;;
        "SUCCESS")
            echo -e "${GREEN}[SUCCESS]${NC} $timestamp - $message" | tee -a logs/batch_process.log
            ;;
    esac
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--data-dir)
            ROOT_DIR="$2"
            shift 2
            ;;
        -j|--jobs)
            MAX_PARALLEL="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --skip-conversion)
            SKIP_CONVERSION=true
            shift
            ;;
        --skip-pipeline)
            SKIP_PIPELINE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [ -z "$ROOT_DIR" ]; then
    echo "Error: Data directory is required"
    usage
    exit 1
fi

if [ ! -d "$ROOT_DIR" ]; then
    log_message "ERROR" "Data directory does not exist: $ROOT_DIR"
    exit 1
fi

mkdir -p logs

echo "=== Batch Processing Started: $(date) ===" > logs/batch_process.log

log_message "INFO" "Starting batch processing"
log_message "INFO" "Data directory: $ROOT_DIR"
log_message "INFO" "Max parallel jobs: $MAX_PARALLEL"
log_message "INFO" "Configuration file: $CONFIG_FILE"

PATIENT_DIRS=()
for patient_path in "$ROOT_DIR"/*; do
    if [ -d "$patient_path" ]; then
        patient_name=$(basename "$patient_path")
        if ls "$patient_path"/Visit* >/dev/null 2>&1; then
            PATIENT_DIRS+=("$patient_path")
        else
            log_message "WARN" "Skipping $patient_name - no Visit folders found"
        fi
    fi
done

if [ ${#PATIENT_DIRS[@]} -eq 0 ]; then
    log_message "ERROR" "No valid patient directories found in $ROOT_DIR"
    exit 1
fi

log_message "INFO" "Found ${#PATIENT_DIRS[@]} patient directories to process"

process_patient() {
    local patient_dir=$1
    local patient_name=$(basename "$patient_dir")
    
    log_message "INFO" "Processing patient: $patient_name"
    
    if [ "$SKIP_CONVERSION" = false ]; then
        log_message "INFO" "Converting DICOMs for $patient_name"
        if "$SCRIPT_DIR/convert_single_patient.sh" "$patient_dir" 2>&1 | tee -a "logs/${patient_name}_conversion.log"; then
            log_message "SUCCESS" "DICOM conversion completed for $patient_name"
        else
            log_message "ERROR" "DICOM conversion failed for $patient_name"
            return 1
        fi
    fi
    
    if [ "$SKIP_PIPELINE" = false ]; then
        log_message "INFO" "Running longitudinal pipeline for $patient_name"
        if Rscript "$SCRIPT_DIR/../R/03_longitudinal_pipeline.R" "$patient_dir" 2>&1 | tee -a "logs/${patient_name}_pipeline.log"; then
            log_message "SUCCESS" "Pipeline processing completed for $patient_name"
        else
            log_message "ERROR" "Pipeline processing failed for $patient_name"
            return 1
        fi
    fi
    
    log_message "SUCCESS" "Patient $patient_name processed successfully"
    return 0
}

export -f process_patient
export -f log_message
export SKIP_CONVERSION
export SKIP_PIPELINE
export RED GREEN YELLOW BLUE NC

log_message "INFO" "Starting parallel processing with $MAX_PARALLEL jobs"

if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${PATIENT_DIRS[@]}" | parallel -j "$MAX_PARALLEL" process_patient
else
    printf '%s\n' "${PATIENT_DIRS[@]}" | xargs -n 1 -P "$MAX_PARALLEL" -I {} bash -c 'process_patient "$@"' _ {}
fi

log_message "INFO" "Generating processing summary"

TOTAL_PATIENTS=${#PATIENT_DIRS[@]}
SUCCESSFUL=0
FAILED=0

for patient_dir in "${PATIENT_DIRS[@]}"; do
    patient_name=$(basename "$patient_dir")
    if [ -d "results/$patient_name" ]; then
        ((SUCCESSFUL++))
    else
        ((FAILED++))
    fi
done

log_message "INFO" "Processing Summary:"
log_message "INFO" "  Total patients: $TOTAL_PATIENTS"
log_message "INFO" "  Successful: $SUCCESSFUL"
log_message "INFO" "  Failed: $FAILED"

{
    echo "Batch Processing Summary - $(date)"
    echo "=================================="
    echo "Total patients processed: $TOTAL_PATIENTS"
    echo "Successful: $SUCCESSFUL"
    echo "Failed: $FAILED"
    echo ""
    echo "Individual logs can be found in the logs/ directory"
} > logs/batch_summary.txt

log_message "SUCCESS" "Batch processing completed. Summary saved to logs/batch_summary.txt"

if [ $FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
