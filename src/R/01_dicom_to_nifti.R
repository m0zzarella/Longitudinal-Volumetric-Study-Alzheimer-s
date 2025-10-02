#!/usr/bin/env Rscript

suppressMessages({
  library(oro.dicom)
  library(oro.nifti)
})

args <- commandArgs(trailingOnly = TRUE)
patient_dir <- args[1]
if (is.na(patient_dir)) stop("Usage: Rscript 01_dicom_to_nifti.R /path/to/Patient1")
if (!dir.exists(patient_dir)) stop("Patient directory does not exist: ", patient_dir)

visit_dirs <- list.dirs(patient_dir, recursive = FALSE, full.names = TRUE)
visit_dirs <- visit_dirs[dir.exists(file.path(visit_dirs, "DICOMs"))]
if (length(visit_dirs) == 0) stop("No visit directories with 'DICOMs' folders found in: ", patient_dir)

for (visit in visit_dirs) {
  dicom_dir <- file.path(visit, "DICOMs")
  message("Reading DICOMs from: ", dicom_dir)
  
  tryCatch({
    dcm <- readDICOM(dicom_dir)
    nii <- dicom2nifti(dcm)
    
    out_base <- file.path(visit, paste0(basename(visit), "_scan"))
    writeNIfTI(nii, filename = out_base)
    message("Wrote: ", out_base, ".nii.gz")
  }, error = function(e) {
    warning("Failed to process ", visit, ": ", e$message)
  })
}

message("Done.")
