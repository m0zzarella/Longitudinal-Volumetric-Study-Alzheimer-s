DATA_DIR ?= data
RESULTS_DIR ?= results
PATIENT ?= 
MAX_JOBS ?= 4

BLUE = \033[0;34m
GREEN = \033[0;32m
YELLOW = \033[1;33m
RED = \033[0;31m
NC = \033[0m 

.DEFAULT_GOAL := help

help:
	@echo "$(BLUE)Longitudinal Neuroanalysis Pipeline$(NC)"
	@echo "=================================="
	@echo ""
	@echo "Available targets:"
	@echo ""
	@grep -E '^##' $(MAKEFILE_LIST) | sed 's/##//' | column -t -s ':'
	@echo ""
	@echo "Examples:"
	@echo "  make setup                    #Initial environment setup"
	@echo "  make convert PATIENT=Patient1 #Convert single patient"
	@echo "  make process PATIENT=Patient1 #Process single patient"
	@echo "  make full PATIENT=Patient1    #Convert + process single patient"
	@echo "  make batch DATA_DIR=data/     #Batch process all patients"
	@echo "  make qc                       #Quality control check"
	@echo "  make mann-kendall             #Mann-Kendall trend analysis"


PATIENT_DIR := $(DATA_DIR)/$(PATIENT)
RESULT_PATIENT_DIR := $(RESULTS_DIR)/$(PATIENT)
VOLUME_CSV := $(RESULTS_DIR)/$(PATIENT)/$(PATIENT)_longitudinal_volumes.csv

.PHONY: check-patient
check-patient:
	@if [ -z "$(PATIENT)" ]; then \
		echo "$(RED)Error: PATIENT variable is required$(NC)"; \
		echo "Usage: make <target> PATIENT=Patient123"; \
		exit 1; \
	fi

setup:
	@echo "$(BLUE)Setting up environment...$(NC)"
	@mkdir -p data results logs temp config
	@chmod +x src/scripts/*.sh 2>/dev/null || echo "Note: chmod not available on Windows"
	@bash src/scripts/setup_environment.sh

convert: check-patient $(DATA_DIR)/$(PATIENT)/.conversion.done
	@echo "$(GREEN)Conversion up-to-date for $(PATIENT)$(NC)"

process: check-patient $(VOLUME_CSV)
	@echo "$(GREEN)Processing up-to-date for $(PATIENT)$(NC)"

full: check-patient $(DATA_DIR)/$(PATIENT)/.conversion.done $(VOLUME_CSV)
	@echo "$(GREEN)Full pipeline up-to-date for $(PATIENT)$(NC)"

batch:
	@echo "$(BLUE)Batch processing all patients in $(DATA_DIR)...$(NC)"
	@bash src/scripts/batch_process.sh -d "$(DATA_DIR)" -j $(MAX_JOBS)

convert-all:
	@echo "$(BLUE)Converting all patients in $(DATA_DIR)...$(NC)"
	@bash src/scripts/convert_all.sh "$(DATA_DIR)"

qc:
	@echo "$(BLUE)Running quality control checks...$(NC)"
	@bash src/scripts/quality_check.sh -r "$(RESULTS_DIR)"

clean:
	@echo "$(BLUE)Cleaning temporary files...$(NC)"
	@bash src/scripts/cleanup.sh --dry-run

clean-force:
	@echo "$(YELLOW)Cleaning temporary files (force)...$(NC)"
	@bash src/scripts/cleanup.sh

clean-all:
	@echo "$(RED)Warning: This will remove all results!$(NC)"
	@read -p "Are you sure? (y/N): " confirm && [ "$$confirm" = "y" ]
	@rm -rf results/* logs/* temp/*
	@echo "$(GREEN)All results cleaned$(NC)"



install-r:
	@echo "$(BLUE)Installing R dependencies...$(NC)"
	@R -e "install.packages(c('oro.dicom', 'oro.nifti', 'scales', 'devtools', 'remotes'), repos='https://cran.r-project.org')"
	@R -e "remotes::install_github('muschellij2/neurobase')"
	@R -e "remotes::install_github('muschellij2/fslr')"

validate:
	@echo "$(BLUE)Validating environment...$(NC)"
	@bash src/scripts/setup_environment.sh

status:
	@echo "$(BLUE)Pipeline Status$(NC)"
	@echo "==============="
	@echo "Data directory: $(DATA_DIR)"
	@echo "Results directory: $(RESULTS_DIR)"
	@echo ""
	@if [ -d "$(DATA_DIR)" ]; then \
		echo "Patients in data directory: $$(find $(DATA_DIR) -maxdepth 1 -type d -name 'Patient*' | wc -l)"; \
	else \
		echo "$(YELLOW)Data directory not found$(NC)"; \
	fi
	@if [ -d "$(RESULTS_DIR)" ]; then \
		echo "Processed patients: $$(find $(RESULTS_DIR) -maxdepth 1 -type d | wc -l)"; \
	else \
		echo "$(YELLOW)Results directory not found$(NC)"; \
	fi

report:
	@echo "$(BLUE)Generating processing report...$(NC)"
	@bash src/scripts/quality_check.sh -r "$(RESULTS_DIR)" > logs/quality_report.txt
	@echo "$(GREEN)Report saved to logs/quality_report.txt$(NC)"

mann-kendall: $(RESULTS_DIR)/mann_kendall_analysis_results.csv
	@echo "$(GREEN)Mann-Kendall analysis up-to-date: $(RESULTS_DIR)/mann_kendall_analysis_results.csv$(NC)"

mann-kendall-patient: check-patient $(VOLUME_CSV)
	@echo "$(BLUE)Running Mann-Kendall analysis for $(PATIENT)...$(NC)"
	@bash src/scripts/mann_kendall_analysis.sh "$(RESULTS_DIR)" "$(PATIENT)"
	@echo "$(GREEN)Analysis complete for $(PATIENT). Check $(RESULTS_DIR)/mann_kendall_analysis_results.csv$(NC)"

mann-kendall-check:
	@echo "$(BLUE)Checking Mann-Kendall analysis dependencies...$(NC)"
	@bash src/scripts/mann_kendall_analysis.sh --check-deps

mann-kendall-preview:
	@echo "$(BLUE)Previewing Mann-Kendall analysis...$(NC)"
	@bash src/scripts/mann_kendall_analysis.sh --dry-run "$(RESULTS_DIR)"

sample:
	@echo "$(BLUE)Creating sample data structure...$(NC)"
	@mkdir -p data/sample_patient/Visit1/DICOMs
	@mkdir -p data/sample_patient/Visit2/DICOMs
	@echo "$(GREEN)Sample structure created in data/sample_patient/$(NC)"


.SECONDEXPANSION:

$(RESULTS_DIR)/%:
	@mkdir -p "$@"

$(DATA_DIR)/%/.conversion.done: $$(wildcard $$(@D)/DICOMs/*)
	@echo "$(BLUE)Converting DICOMs -> NIfTI for $*...$(NC)"
	@bash src/scripts/convert_single_patient.sh "$(DATA_DIR)/$*"
	@printf "" > "$@"

$(RESULTS_DIR)/%/%_longitudinal_volumes.csv: $(DATA_DIR)/%/.conversion.done | $(RESULTS_DIR)/%
	@echo "$(BLUE)Running longitudinal pipeline for $*...$(NC)"
	@Rscript src/R/03_longitudinal_pipeline.R "$(DATA_DIR)/$*"
	@echo "$(GREEN)Wrote: $(RESULTS_DIR)/$*/$*_longitudinal_volumes.csv$(NC)"

$(RESULTS_DIR)/mann_kendall_analysis_results.csv: $(wildcard $(RESULTS_DIR)/*/*_longitudinal_volumes.csv)
	@echo "$(BLUE)Running Mann-Kendall trend analysis...$(NC)"
	@bash src/scripts/mann_kendall_analysis.sh "$(RESULTS_DIR)"
	@echo "$(GREEN)Analysis complete. Check $(RESULTS_DIR)/mann_kendall_analysis_results.csv$(NC)"

.PHONY: help setup check-patient convert process full batch convert-all qc clean clean-force clean-all \
        install-r validate status report sample mann-kendall mann-kendall-patient \
        mann-kendall-check mann-kendall-preview
