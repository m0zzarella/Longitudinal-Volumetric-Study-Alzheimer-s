#!/usr/bin/env Rscript

suppressMessages({
  library(oro.nifti)
  library(neurobase)
})

args <- commandArgs(trailingOnly = TRUE)
patient_dir <- args[1]
if (is.na(patient_dir)) stop("Usage: Rscript 02_visualize_qc.R /path/to/Patient1")
if (!dir.exists(patient_dir)) stop("Patient directory does not exist: ", patient_dir)

qc_dir <- file.path("results", basename(patient_dir), "qc")
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

niftis <- list.files(patient_dir, pattern = "\\.nii(\\.gz)?$", recursive = TRUE, full.names = TRUE)
if (length(niftis) == 0) stop("No NIfTI files found under: ", patient_dir)

for (nii_path in niftis) {
  message("QC view for: ", nii_path)
  img <- readNIfTI(nii_path, reorient = FALSE)
  basefile <- sub("\\.nii(\\.gz)?$", "", basename(nii_path))
  qc_name <- paste0(basefile, "_qc.png")
  png_path <- file.path(qc_dir, qc_name)

  png(png_path, width = 1000, height = 800)
  try({
    orthographic(img, xyz = c(128,128,170))
  }, silent = TRUE)
  dev.off()
}

message("PNGs written to: ", qc_dir)
