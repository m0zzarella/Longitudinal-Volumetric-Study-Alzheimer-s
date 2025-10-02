#!/bin/bash

set -euo pipefail

RESULTS_DIR="results"
PATIENT=""
GENERATE_REPORT=true
CHECK_VOLUMES=true
CHECK_REGISTRATION=true

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

usage() {
    echo "Usage: $0 [OPTIONS] [PATIENT_ID]"
    echo ""
    echo "Options:"
    echo "  -r, --results-dir DIR    Results directory (default: results)"
    echo "  --no-report             Don't generate HTML report"
    echo "  --no-volumes            Skip volume validation"
    echo "  --no-registration       Skip registration quality checks"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                      # Check all patients"
    echo "  $0 Patient123           # Check specific patient"
    echo "  $0 -r /path/results     # Use custom results directory"
}

log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${BLUE}[QC]${NC} $message"
            ;;
        "PASS")
            echo -e "${GREEN}[PASS]${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} $message"
            ;;
        "FAIL")
            echo -e "${RED}[FAIL]${NC} $message"
            ;;
    esac
}

check_file() {
    local file=$1
    local description=$2
    
    if [ -f "$file" ] && [ -s "$file" ]; then
        log_message "PASS" "$description exists and is not empty"
        return 0
    else
        log_message "FAIL" "$description missing or empty: $file"
        return 1
    fi
}

check_nifti_integrity() {
    local nifti_file=$1
    local description=$2
    
    if [ ! -f "$nifti_file" ]; then
        log_message "FAIL" "$description file not found: $nifti_file"
        return 1
    fi
    
    if fslinfo "$nifti_file" >/dev/null 2>&1; then
        local dims=$(fslinfo "$nifti_file" | grep "^dim[1-3]" | awk '{print $2}' | tr '\n' 'x' | sed 's/x$//')
        log_message "PASS" "$description valid NIfTI (dimensions: $dims)"
        return 0
    else
        log_message "FAIL" "$description corrupted or invalid NIfTI: $nifti_file"
        return 1
    fi
}

validate_volumes() {
    local volume_file=$1
    local visit_name=$2
    
    if [ ! -f "$volume_file" ]; then
        log_message "FAIL" "Volume file not found: $volume_file"
        return 1
    fi
    
    if head -1 "$volume_file" | grep -qi "CSF\|GM\|WM"; then
        log_message "PASS" "$visit_name volume measurements available"
        read -r csf_vol gm_vol wm_vol < <(awk -F',' 'NR==1{for(i=1;i<=NF;i++){h=tolower($i); if(h~"^csf") c=i; if(h~"^gm") g=i; if(h~"^wm") w=i; } next} {csf=$c; gm=$g; wm=$w} END{print csf, gm, wm}' "$volume_file")
        
        if [ "$(echo "$csf_vol > 0" | bc 2>/dev/null || echo "0")" = "1" ] && \
           [ "$(echo "$gm_vol > 0" | bc 2>/dev/null || echo "0")" = "1" ] && \
           [ "$(echo "$wm_vol > 0" | bc 2>/dev/null || echo "0")" = "1" ]; then
            log_message "PASS" "$visit_name volumes within reasonable range"
        else
            log_message "WARN" "$visit_name volumes may be unrealistic"
        fi
    else
        log_message "FAIL" "$visit_name volume file format invalid"
        return 1
    fi
}

check_registration_quality() {
    local reg_file=$1
    local visit_name=$2
    
    if [ ! -f "$reg_file" ]; then
        log_message "WARN" "Registration quality file not found: $reg_file"
        return 1
    fi
    
    if grep -qi "correlation\|nmi" "$reg_file" 2>/dev/null; then
        local correlation
        correlation=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++){if(tolower($i)~/(^| )correlation( |$)/){c=i;break}}} NR==2{if(c) print $c}' "$reg_file" 2>/dev/null || echo "N/A")
        log_message "PASS" "$visit_name registration metrics available (correlation: $correlation)"
    else
        log_message "WARN" "$visit_name registration quality metrics incomplete"
    fi
}

check_patient() {
    local patient_dir=$1
    local patient_id=$(basename "$patient_dir")
    
    log_message "INFO" "Checking patient: $patient_id"
    
    local issues=0
    
    if [ ! -d "$patient_dir" ]; then
        log_message "FAIL" "Patient directory not found: $patient_dir"
        return 1
    fi
    
    local visit_dirs=($(find "$patient_dir" -maxdepth 1 -type d -name "Visit*" | sort))
    
    if [ ${#visit_dirs[@]} -eq 0 ]; then
        log_message "FAIL" "No visit directories found for $patient_id"
        return 1
    fi
    
    log_message "INFO" "Found ${#visit_dirs[@]} visits for $patient_id"
    
    for visit_dir in "${visit_dirs[@]}"; do
        local visit_name=$(basename "$visit_dir")
        log_message "INFO" "  Checking $visit_name"
        
        local n4_file="$visit_dir/${visit_name}_N4.nii.gz"
        local bet_file="$visit_dir/${visit_name}_BET.nii.gz"
        local volumes_file="$visit_dir/${visit_name}_volumes.csv"
        
        check_nifti_integrity "$n4_file" "$visit_name N4 bias-corrected" || ((issues++))
        check_nifti_integrity "$bet_file" "$visit_name skull-stripped" || ((issues++))
        
        for tissue in 0 1 2; do
            local pve_file="$visit_dir/${visit_name}_BET_pve_${tissue}.nii.gz"
            check_nifti_integrity "$pve_file" "$visit_name tissue $tissue segmentation" || ((issues++))
        done
        
        if [ "$CHECK_VOLUMES" = true ]; then
            validate_volumes "$volumes_file" "$visit_name" || ((issues++))
        fi
        
        if [[ "$visit_name" != "Visit1" ]] && [ "$CHECK_REGISTRATION" = true ]; then
            local reg_file="$visit_dir/${visit_name}_BET_reg.nii.gz"
            local reg_quality="$visit_dir/${visit_name}_registration_quality.csv"
            
            check_nifti_integrity "$reg_file" "$visit_name registered" || ((issues++))
            check_registration_quality "$reg_quality" "$visit_name" || ((issues++))
        fi
    done
    
    local long_volumes="$patient_dir/${patient_id}_longitudinal_volumes.csv"
    local long_deltas="$patient_dir/${patient_id}_longitudinal_deltas.csv"
    
    check_file "$long_volumes" "Longitudinal volumes summary" || ((issues++))
    check_file "$long_deltas" "Longitudinal deltas summary" || ((issues++))
    
    if [ $issues -eq 0 ]; then
        log_message "PASS" "Patient $patient_id passed all quality checks"
        return 0
    else
        log_message "FAIL" "Patient $patient_id has $issues issues"
        return 1
    fi
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--results-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        --no-report)
            GENERATE_REPORT=false
            shift
            ;;
        --no-volumes)
            CHECK_VOLUMES=false
            shift
            ;;
        --no-registration)
            CHECK_REGISTRATION=false
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            PATIENT="$1"
            shift
            ;;
    esac
done

if [ ! -d "$RESULTS_DIR" ]; then
    log_message "FAIL" "Results directory not found: $RESULTS_DIR"
    exit 1
fi

log_message "INFO" "Starting quality control checks"
log_message "INFO" "Results directory: $RESULTS_DIR"

total_patients=0
passed_patients=0
failed_patients=0

if [ -n "$PATIENT" ]; then
    patient_dir="$RESULTS_DIR/$PATIENT"
    total_patients=1
    if check_patient "$patient_dir"; then
        ((passed_patients++))
    else
        ((failed_patients++))
    fi
else
    for patient_dir in "$RESULTS_DIR"/*; do
        if [ -d "$patient_dir" ]; then
            ((total_patients++))
            if check_patient "$patient_dir"; then
                ((passed_patients++))
            else
                ((failed_patients++))
            fi
        fi
    done
fi

log_message "INFO" "Quality Control Summary:"
log_message "INFO" "  Total patients: $total_patients"
log_message "INFO" "  Passed: $passed_patients"
log_message "INFO" "  Failed: $failed_patients"

if [ $failed_patients -eq 0 ]; then
    log_message "PASS" "All quality checks passed!"
    exit 0
else
    log_message "FAIL" "$failed_patients patient(s) failed quality checks"
    exit 1
fi
