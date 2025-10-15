#!/usr/bin/env Rscript

# Expects visit folders:
#   Patient#/Visit1/
#   Patient#/Visit2/
#   Patient#/Visit3
#   ...
# Each visit should contain a NIfTI (e.g., scan.nii.gz or any .nii/.nii.gz)
# If starting from DICOMs, run convert_patient.sh or 01_dicom_to_nifti.R first.

suppressMessages({
  library(oro.nifti)
  library(neurobase)
  library(fslr)
  library(scales)
  library(yaml)

  get_script_dir <- function() {
    script_path <- tryCatch({
      sys.frame(1)$ofile
    }, error = function(e) NULL)
    
    if (is.null(script_path)) {
      script_path <- tryCatch({
        args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("--file=", args, value = TRUE)
        if (length(file_arg) > 0) {
          sub("--file=", "", file_arg[1])
        } else {
          NULL
        }
      }, error = function(e) NULL)
    }
    
    if (is.null(script_path)) {
      file.path(getwd(), "src", "R")
    } else {
      dirname(script_path)
    }
  }
  
  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))
  source(file.path(script_dir, "04_intensity_normalization.R"))
  source(file.path(script_dir, "05_registration.R"))

  config <- load_config()
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
patient_dir <- args[1]
if (is.na(patient_dir)) stop("Usage: Rscript 03_longitudinal_pipeline.R /path/to/Patient1")
if (!dir.exists(patient_dir)) stop("Patient directory does not exist: ", patient_dir)

visit_dirs <- list.dirs(patient_dir, recursive = FALSE, full.names = TRUE)
visit_dirs <- visit_dirs[grepl("Visit", basename(visit_dirs), ignore.case = TRUE)]
if (length(visit_dirs) == 0) stop("No Visit* folders found under ", patient_dir)

patient_id <- basename(normalizePath(patient_dir))
out_root <- file.path("results", patient_id)
ensure_dir(out_root)


volumes <- list()
visit_names <- c()
baseline_img <- NULL  
baseline_bet <- NULL  
failed_visits <- list()  

for (visit in sort(visit_dirs)) {
  vname <- basename(visit)
  visit_names <- c(visit_names, vname)
  message("\n== Processing ", vname, " ==")

  visit_out <- file.path(out_root, vname)
  ensure_dir(visit_out)

  nifti_path <- find_visit_nifti(visit)
  if (is.na(nifti_path)) {
    warning("No NIfTI found in ", visit, " — skipping")
    next
  }

  img <- readNIfTI(nifti_path, reorient = FALSE)

  #N4 bias field correction
  message("Bias correction (N4 ITK): ")
  bc_img <- fslr::fsl_biascorrect(img, retimg = TRUE)
  writeNIfTI(bc_img, file.path(visit_out, paste0(vname, "_N4")))

  message("Intensity normalization: ")
  if (needs_intensity_normalization(bc_img)) {
    norm_img <- normalize_intensity(bc_img, method = "zscore")
    writeNIfTI(norm_img, file.path(visit_out, paste0(vname, "_N4_norm")))
    message("Applied z-score normalization")
  } else {
    norm_img <- bc_img
    message("Skipped normalization (CV within acceptable range)")
  }

  #Skull stripping using FSL BET
  message("Skull stripping: ")
  bet_img <- run_bet(norm_img, cog_init = TRUE)
  writeNIfTI(bet_img, file.path(visit_out, paste0(vname, "_BET")))
  
  if (vname == "Visit1" || (is.null(baseline_img) && length(visit_names) == 1)) {
    baseline_img <- norm_img
    baseline_bet <- bet_img
    message("Set as baseline for registration")
  }
  
  if (!is.null(baseline_img) && vname != "Visit1" && !(is.null(baseline_img) && length(visit_names) == 1)) {
    message("Registration to baseline: ")
    
    brain_mask <- create_brain_mask(bet_img)
    baseline_mask <- create_brain_mask(baseline_bet)
    
    reg_result <- register_with_mask(
      moving_img = bet_img,
      fixed_img = baseline_bet,
      moving_mask = brain_mask,
      fixed_mask = baseline_mask,
      registration_type = "linear",
      cost_function = "corratio"
    )

    reg_img <- reg_result$final_registered
    writeNIfTI(reg_img, file.path(visit_out, paste0(vname, "_BET_reg")))

    reg_quality <- validate_registration(reg_img, baseline_bet, baseline_mask)

    if (reg_quality$is_valid) {
      message("Registration quality: valid")
      bet_img <- reg_img
    } else {
      message("Linear registration quality: using original image")
      message("Linear registration - Correlation: ", round(reg_quality$metrics$correlation, 3))

      if (config$quality_control$enable_registration_retry) {
        message("Linear registration failed, attempting ", config$quality_control$retry_registration_type, "...")

        reg_result_retry <- register_with_mask(
          moving_img = bet_img,
          fixed_img = baseline_bet,
          moving_mask = brain_mask,
          fixed_mask = baseline_mask,
          registration_type = config$quality_control$retry_registration_type,
          cost_function = "corratio"
        )

        reg_img_retry <- reg_result_retry$final_registered

        retry_suffix <- paste0("_BET_reg_", config$quality_control$retry_registration_type)
        writeNIfTI(reg_img_retry, file.path(visit_out, paste0(vname, retry_suffix)))

        reg_quality_retry <- validate_registration(reg_img_retry, baseline_bet, baseline_mask)

        if (reg_quality_retry$is_valid) {
          message(config$quality_control$retry_registration_type, " registration succeeded")
          reg_img <- reg_img_retry
          reg_quality <- reg_quality_retry
          bet_img <- reg_img
        } else {
          warning("Both linear and ", config$quality_control$retry_registration_type, " registration failed - skipping visit")
          message(config$quality_control$retry_registration_type, " registration - Correlation: ",
                  round(reg_quality_retry$metrics$correlation, 3))

          failed_visits[[vname]] <- list(
            reason = paste("Both linear and", config$quality_control$retry_registration_type, "registration failed"),
            correlation = reg_quality_retry$metrics$correlation,
            mse = reg_quality_retry$metrics$mse,
            registration_type = config$quality_control$retry_registration_type
          )

          next
        }
      } else {
        message("Registration retry disabled, using original image")
      }
    }
    
    reg_quality_csv <- file.path(visit_out, paste0(vname, "_registration_quality.csv"))
    write.csv(data.frame(
      Visit = vname,
      Correlation = reg_quality$metrics$correlation,
      NMI = reg_quality$metrics$normalized_mutual_information,
      MSE = reg_quality$metrics$mse,
      MAE = reg_quality$metrics$mae,
      SSIM = reg_quality$metrics$ssim,
      Quality_Check = if(reg_quality$is_valid) "PASSED" else "WARNING"
    ), reg_quality_csv, row.names = FALSE)
    
    reg_qc_png <- file.path(visit_out, paste0(vname, "_registration_qc.png"))
    create_registration_qc(bet_img, baseline_bet, reg_img, reg_qc_png)
  }

  message("FAST segmentation ...")
  fast_out <- fslr::fast(
    file = bet_img,
    outfile = file.path(visit_out, paste0(vname, "_BET"))
  )

  pve_csf_path <- file.path(visit_out, paste0(vname, "_BET_pve_0.nii.gz"))
  pve_gm_path  <- file.path(visit_out, paste0(vname, "_BET_pve_1.nii.gz"))
  pve_wm_path  <- file.path(visit_out, paste0(vname, "_BET_pve_2.nii.gz"))

  if (!file.exists(pve_csf_path) || !file.exists(pve_gm_path) || !file.exists(pve_wm_path)) {
    warning("FAST PVE outputs not found for ", vname)
    next
  }

  pve_csf <- readNIfTI(pve_csf_path, reorient = FALSE)
  pve_gm  <- readNIfTI(pve_gm_path,  reorient = FALSE)
  pve_wm  <- readNIfTI(pve_wm_path,  reorient = FALSE)

  vol_csf <- compute_volume_ml(pve_csf, threshold = config$processing$pve_threshold)
  vol_gm  <- compute_volume_ml(pve_gm,  threshold = config$processing$pve_threshold)
  vol_wm  <- compute_volume_ml(pve_wm,  threshold = config$processing$pve_threshold)


  visit_csv <- file.path(visit_out, paste0(vname, "_volumes.csv"))
  write.csv(
    data.frame(Visit = vname, CSF_ml = vol_csf, GM_ml = vol_gm, WM_ml = vol_wm),
    visit_csv, row.names = FALSE
  )

  volumes[[vname]] <- c(CSF_ml = vol_csf, GM_ml = vol_gm, WM_ml = vol_wm)

  #Create QC overlay images showing tissue segmentation

  qc_png <- file.path(visit_out, paste0(vname, "_qc_overlay.png"))
  png(qc_png, width = 1400, height = 900)
  par(mfrow = c(1,3), mar = c(2,2,2,2))
  mtext(paste("QC Overlay:", patient_id, "-", vname), side = 3, line = -2, outer = TRUE, cex = 1.5, font = 2)

  try(ortho2(bet_img, pve_csf > config$processing$pve_threshold, col.y = alpha("red", 0.5),  text = "CSF", xyz = c(128,128,170)), silent = TRUE)
  try(ortho2(bet_img, pve_gm  > config$processing$pve_threshold, col.y = alpha("blue", 0.5), text = "GM",  xyz = c(128,128,170)), silent = TRUE)
  try(ortho2(bet_img, pve_wm  > config$processing$pve_threshold, col.y = alpha("green", 0.5), text = "WM",  xyz = c(128,128,170)), silent = TRUE)
  dev.off()
}

if (length(volumes) >= 1) {
  long_csv <- file.path(out_root, paste0(patient_id, "_longitudinal_volumes.csv"))
  vol_mat <- do.call(rbind, volumes)
  vol_df <- data.frame(Visit = names(volumes), vol_mat, row.names = NULL)
  write.csv(vol_df, long_csv, row.names = FALSE)
  message("\nWrote longitudinal table: ", long_csv)

  if (nrow(vol_df) >= 2) {
    baseline <- vol_df[1, ]
    deltas <- vol_df
    deltas$CSF_delta_ml <- deltas$CSF_ml - baseline$CSF_ml
    deltas$GM_delta_ml  <- deltas$GM_ml  - baseline$GM_ml
    deltas$WM_delta_ml  <- deltas$WM_ml  - baseline$WM_ml
    delta_csv <- file.path(out_root, paste0(patient_id, "_longitudinal_deltas.csv"))
    write.csv(deltas[, c("Visit", "CSF_delta_ml", "GM_delta_ml", "WM_delta_ml")], delta_csv, row.names = FALSE)
    message("Wrote deltas table: ", delta_csv)
  }
}

if (length(failed_visits) > 0) {
  failed_csv <- file.path(out_root, paste0(patient_id, "_failed_visits.csv"))
  failed_df <- do.call(rbind, lapply(names(failed_visits), function(v) {
    data.frame(
      Visit = v,
      Reason = failed_visits[[v]]$reason,
      Correlation = failed_visits[[v]]$correlation,
      MSE = failed_visits[[v]]$mse,
      Registration_Type = failed_visits[[v]]$registration_type,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(failed_df, failed_csv, row.names = FALSE)
  message("\nFailed visits logged to: ", failed_csv)
  message("Total failed visits: ", nrow(failed_df))
}

message("\nDone.")
