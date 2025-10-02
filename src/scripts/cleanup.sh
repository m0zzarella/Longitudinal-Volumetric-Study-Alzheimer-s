#!/bin/bash

set -euo pipefail
shopt -s nullglob

RESULTS_DIR="results"
TEMP_DIR="temp"
DRY_RUN=false
KEEP_INTERMEDIATE=false
KEEP_LOGS=true

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -r, --results-dir DIR    Results directory (default: results)"
    echo "  -t, --temp-dir DIR       Temporary directory (default: temp)"
    echo "  --dry-run               Show what would be deleted without actually deleting"
    echo "  --keep-intermediate     Keep intermediate processing files"
    echo "  --remove-logs           Also remove log files"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --dry-run            # Preview cleanup actions"
    echo "  $0 --remove-logs        # Clean including log files"
}

log_message() {
    local level=$1
    local message=$2
    
    case $level in
        "INFO")
            echo -e "${BLUE}[CLEANUP]${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} $message"
            ;;
        "SUCCESS")
            echo -e "${GREEN}[SUCCESS]${NC} $message"
            ;;
        "DRY_RUN")
            echo -e "${YELLOW}[DRY RUN]${NC} Would delete: $message"
            ;;
    esac
}

safe_remove() {
    local target=$1
    local description=$2
    
    if [ -e "$target" ]; then
        if [ "$DRY_RUN" = true ]; then
            log_message "DRY_RUN" "$description: $target"
        else
            rm -rf "$target"
            log_message "SUCCESS" "Removed $description: $target"
        fi
        return 0
    else
        return 1
    fi
}

clean_patient_dir() {
    local patient_dir=$1
    local patient_id=$(basename "$patient_dir")
    
    log_message "INFO" "Cleaning patient: $patient_id"
    
    for visit_dir in "$patient_dir"/Visit*; do
        if [ -d "$visit_dir" ]; then
            local visit_name=$(basename "$visit_dir")
            
            if [ "$KEEP_INTERMEDIATE" = false ]; then
                safe_remove "$visit_dir"/*.mat "transformation matrices"
                safe_remove "$visit_dir"/*_tmp* "temporary files"
                safe_remove "$visit_dir"/*.log "processing logs"
                
                safe_remove "$visit_dir"/*_flirt* "FLIRT intermediate files"
                safe_remove "$visit_dir"/*_fnirt* "FNIRT intermediate files"
                
                safe_remove "$visit_dir"/*.bak "backup files"
                safe_remove "$visit_dir"/*~ "editor backup files"
            fi
            
            safe_remove "$visit_dir"/core.* "core dump files"
            safe_remove "$visit_dir"/.nfs* "NFS temporary files"
        fi
    done
    
    safe_remove "$patient_dir"/*.tmp "patient temporary files"
    safe_remove "$patient_dir"/core.* "core dump files"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--results-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        -t|--temp-dir)
            TEMP_DIR="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --keep-intermediate)
            KEEP_INTERMEDIATE=true
            shift
            ;;
        --remove-logs)
            KEEP_LOGS=false
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

if [ "$DRY_RUN" = true ]; then
    log_message "INFO" "DRY RUN MODE - No files will actually be deleted"
fi

log_message "INFO" "Starting cleanup process"

if [ -d "$TEMP_DIR" ]; then
    safe_remove "$TEMP_DIR"/* "temporary directory contents"
else
    log_message "INFO" "Temporary directory not found: $TEMP_DIR"
fi

if [ -d "$RESULTS_DIR" ]; then
    log_message "INFO" "Cleaning results directory: $RESULTS_DIR"
    
    for patient_dir in "$RESULTS_DIR"/*; do
        if [ -d "$patient_dir" ]; then
            clean_patient_dir "$patient_dir"
        fi
    done
else
    log_message "WARN" "Results directory not found: $RESULTS_DIR"
fi

if [ "$KEEP_LOGS" = false ]; then
    log_message "INFO" "Cleaning log files"
    for f in logs/*.log logs/*_conversion.log logs/*_pipeline.log; do
        safe_remove "$f" "log file"
    done
fi

for f in /tmp/fsl_* /tmp/R_*; do
    safe_remove "$f" "temporary file in /tmp"
done

for f in ./*.tmp ./core.* ./.nfs*; do
    safe_remove "$f" "temporary file in current directory"
done

for f in ./.Rhistory ./.RData; do
    safe_remove "$f" "R workspace/history file"
done

if [ "$DRY_RUN" = true ]; then
    log_message "INFO" "Dry run completed. Use without --dry-run to actually delete files."
else
    log_message "SUCCESS" "Cleanup completed successfully!"
fi

log_message "INFO" "Cleanup settings:"
log_message "INFO" "  Results directory: $RESULTS_DIR"
log_message "INFO" "  Temp directory: $TEMP_DIR"
log_message "INFO" "  Keep intermediate files: $KEEP_INTERMEDIATE"
log_message "INFO" "  Keep log files: $KEEP_LOGS"
