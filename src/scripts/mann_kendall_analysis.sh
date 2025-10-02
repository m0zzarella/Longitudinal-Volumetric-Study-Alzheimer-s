#!/bin/bash

# Mann-Kendall Trend Analysis Script
# Wrapper for the R Mann-Kendall analysis script
#
# Usage: 
#   ./mann_kendall_analysis.sh [results_directory] [patient_id]
#   
# Examples:
#   ./mann_kendall_analysis.sh results/                    #Analyze all patients
#   ./mann_kendall_analysis.sh results/ Patient123        #Analyze specific patient
#   ./mann_kendall_analysis.sh                            #Use default results/ directory

set -euo pipefail

RESULTS_DIR="results"
PATIENT_ID=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT_DIR="$(dirname "$SCRIPT_DIR")/R"
R_SCRIPT="$R_SCRIPT_DIR/06_mann_kendall_analysis.R"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${GREEN}[INFO]${NC} ${timestamp} - $message"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} ${timestamp} - $message"
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} ${timestamp} - $message"
            ;;
        "DEBUG")
            echo -e "${BLUE}[DEBUG]${NC} ${timestamp} - $message"
            ;;
    esac
}

show_help() {
    cat << EOF
Mann-Kendall Trend Analysis for Longitudinal Volume Data

USAGE:
    $0 [OPTIONS] [RESULTS_DIR] [PATIENT_ID]

ARGUMENTS:
    RESULTS_DIR     Directory containing patient results (default: results/)
    PATIENT_ID      Specific patient to analyze (optional, analyzes all if not specified)

OPTIONS:
    -h, --help      Show this help message
    -v, --verbose   Enable verbose output
    --check-deps    Check R dependencies only
    --dry-run       Show what would be analyzed without running

EXAMPLES:
    $0                                    # Analyze all patients in results/
    $0 results/                          # Analyze all patients in results/
    $0 results/ Patient123               # Analyze only Patient123
    $0 --check-deps                      # Check R dependencies
    $0 --dry-run results/                # Preview analysis

DESCRIPTION:
    This script performs Mann-Kendall trend analysis on longitudinal brain volume data.
    It calculates statistical significance of volume changes over time using a modified
    Mann-Kendall test with a 1mL threshold to ignore minor fluctuations.
    
    The analysis produces:
    - Individual patient statistics (S, Z-score, p-value, trend interpretation)
    - Summary statistics across all patients
    - CSV output files for further analysis

REQUIREMENTS:
    - R with data.table package
    - Longitudinal volume data from the processing pipeline
    - At least 2 visits per patient for trend analysis

EOF
}

check_r_dependencies() {
    log_message "INFO" "Checking R dependencies..."
    
    if ! command -v Rscript &> /dev/null; then
        log_message "ERROR" "R/Rscript not found. Please install R."
        return 1
    fi
    
    if ! Rscript -e "library(data.table)" &> /dev/null; then
        log_message "WARN" "data.table package not found. It will be installed automatically."
    fi
    
    if [ ! -f "$R_SCRIPT" ]; then
        log_message "ERROR" "R script not found: $R_SCRIPT"
        return 1
    fi
    
    log_message "INFO" "R dependencies check completed successfully"
    return 0
}

validate_inputs() {
    if [ ! -d "$RESULTS_DIR" ]; then
        log_message "ERROR" "Results directory does not exist: $RESULTS_DIR"
        return 1
    fi
    
    if [ -n "$PATIENT_ID" ] && [ ! -d "$RESULTS_DIR/$PATIENT_ID" ]; then
        log_message "ERROR" "Patient directory does not exist: $RESULTS_DIR/$PATIENT_ID"
        return 1
    fi
    
    local volume_files_found=0
    if [ -n "$PATIENT_ID" ]; then
        if [ -f "$RESULTS_DIR/$PATIENT_ID/${PATIENT_ID}_longitudinal_volumes.csv" ]; then
            volume_files_found=1
        fi
    else
        volume_files_found=$(find "$RESULTS_DIR" -name "*_longitudinal_volumes.csv" | wc -l)
    fi
    
    if [ $volume_files_found -eq 0 ]; then
        log_message "ERROR" "No longitudinal volume files found. Please run the longitudinal pipeline first."
        log_message "INFO" "Expected files: {Patient_ID}_longitudinal_volumes.csv"
        return 1
    fi
    
    log_message "INFO" "Found $volume_files_found longitudinal volume file(s)"
    return 0
}

preview_analysis() {
    log_message "INFO" "=== DRY RUN: Preview of Mann-Kendall Analysis ==="
    
    if [ -n "$PATIENT_ID" ]; then
        log_message "INFO" "Would analyze patient: $PATIENT_ID"
        local vol_file="$RESULTS_DIR/$PATIENT_ID/${PATIENT_ID}_longitudinal_volumes.csv"
        if [ -f "$vol_file" ]; then
            local n_visits=$(tail -n +2 "$vol_file" | wc -l)
            log_message "INFO" "  - Number of visits: $n_visits"
            log_message "INFO" "  - Volume file: $vol_file"
        fi
    else
        log_message "INFO" "Would analyze all patients in: $RESULTS_DIR"
        local patient_dirs=($(find "$RESULTS_DIR" -name "*_longitudinal_volumes.csv" -exec dirname {} \;))
        log_message "INFO" "  - Number of patients: ${#patient_dirs[@]}"
        
        for patient_dir in "${patient_dirs[@]}"; do
            local patient=$(basename "$patient_dir")
            local vol_file="$patient_dir/${patient}_longitudinal_volumes.csv"
            local n_visits=$(tail -n +2 "$vol_file" | wc -l)
            log_message "INFO" "    $patient: $n_visits visits"
        done
    fi
    
    log_message "INFO" "Output files would be created in: $RESULTS_DIR"
    log_message "INFO" "  - mann_kendall_analysis_results.csv"
    log_message "INFO" "  - mann_kendall_summary.txt"
}

run_analysis() {
    log_message "INFO" "=== Starting Mann-Kendall Trend Analysis ==="
    log_message "INFO" "Results directory: $RESULTS_DIR"
    
    if [ -n "$PATIENT_ID" ]; then
        log_message "INFO" "Analyzing specific patient: $PATIENT_ID"
        Rscript "$R_SCRIPT" "$RESULTS_DIR" "$PATIENT_ID"
    else
        log_message "INFO" "Analyzing all patients"
        Rscript "$R_SCRIPT" "$RESULTS_DIR"
    fi
    
    local results_file="$RESULTS_DIR/mann_kendall_analysis_results.csv"
    local summary_file="$RESULTS_DIR/mann_kendall_summary.txt"
    
    if [ -f "$results_file" ]; then
        log_message "INFO" "Analysis completed successfully!"
        log_message "INFO" "Results saved to: $results_file"
        log_message "INFO" "Summary saved to: $summary_file"
        
        local total_analyses=$(tail -n +2 "$results_file" | wc -l)
        local significant_trends=$(tail -n +2 "$results_file" | grep -c "Significant" || echo "0")
        
        log_message "INFO" "Brief summary:"
        log_message "INFO" "  - Total analyses: $total_analyses"
        log_message "INFO" "  - Significant trends: $significant_trends"
        
        return 0
    else
        log_message "ERROR" "Analysis failed - no results file created"
        return 1
    fi
}

VERBOSE=false
CHECK_DEPS=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        --check-deps)
            CHECK_DEPS=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -*)
            log_message "ERROR" "Unknown option: $1"
            show_help
            exit 1
            ;;
        *)
            if [ -z "${RESULTS_DIR_SET:-}" ]; then
                RESULTS_DIR="$1"
                RESULTS_DIR_SET=true
            elif [ -z "$PATIENT_ID" ]; then
                PATIENT_ID="$1"
            else
                log_message "ERROR" "Too many arguments"
                show_help
                exit 1
            fi
            shift
            ;;
    esac
done

main() {
    if [ "$CHECK_DEPS" = true ]; then
        check_r_dependencies
        exit $?
    fi
    
    if ! check_r_dependencies; then
        exit 1
    fi
    
    if ! validate_inputs; then
        exit 1
    fi
    
    if [ "$DRY_RUN" = true ]; then
        preview_analysis
        exit 0
    fi
    
    run_analysis
}

main "$@"
