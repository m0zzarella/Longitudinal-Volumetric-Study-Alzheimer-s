#!/usr/bin/env Rscript

suppressMessages({
  library(oro.nifti)
  library(neurobase)
  library(fslr)

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
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

register_to_baseline <- function(moving_img, fixed_img, registration_type = "linear", 
                                cost_function = "corratio", search_range = 90,
                                dof = 12, interpolation = "trilinear") {
  # Register followup visits to baseline scan using FSL FLIRT/FNIRT
  #Args:
  #  moving_img: NIfTI image object
  #  fixed_img: NIfTI image object
  #  registration_type: registration type
  #  cost_function: cost function
  #  search_range: search range
  #  dof: degrees of freedom
  #  interpolation: interpolation method
  #Returns:
  #  list of registration results
  result <- list()
  
  if (registration_type %in% c("linear", "both")) {
    message("Performing linear registration...")
    
    linear_reg <- fslr::flirt(
      infile = moving_img,
      reffile = fixed_img,
      dof = dof,
      cost = cost_function,
      searchrx = c(-search_range, search_range),
      searchry = c(-search_range, search_range),
      searchrz = c(-search_range, search_range),
      interp = interpolation,
      retimg = TRUE
    )
    
    result$linear_registered <- linear_reg
    result$linear_matrix <- attr(linear_reg, "transformation_matrix")
    
    if (registration_type == "linear") {
      result$final_registered <- linear_reg
      return(result)
    }
  }
  
  if (registration_type %in% c("nonlinear", "both")) {
    message("Performing nonlinear registration...")
    
    if (exists("linear_reg")) {
      input_img <- linear_reg
    } else {
      input_img <- moving_img
    }
    
    tryCatch({
      nonlinear_reg <- fslr::fnirt(
        infile = input_img,
        reffile = fixed_img,
        retimg = TRUE
      )
      
      result$nonlinear_registered <- nonlinear_reg
      result$final_registered <- nonlinear_reg
    }, error = function(e) {
      warning("Nonlinear registration failed: ", e$message)
      if (exists("linear_reg")) {
        message("Using linear registration result instead")
        result$final_registered <- linear_reg
      } else {
        stop("Both linear and nonlinear registration failed")
      }
    })
  }
  
  return(result)
}

register_with_mask <- function(moving_img, fixed_img, moving_mask = NULL, fixed_mask = NULL,
                              registration_type = "linear", cost_function = "corratio") {
  # Register images with optional brain masks for improved registration
  
  if (!is.null(moving_mask) && !is.null(fixed_mask)) {
    message("Using brain masks for registration...")
    
    masked_moving <- moving_img
    masked_fixed <- fixed_img
    masked_moving[moving_mask == 0] <- 0
    masked_fixed[fixed_mask == 0] <- 0
    
    result <- register_to_baseline(
      masked_moving, masked_fixed,
      registration_type = registration_type,
      cost_function = cost_function
    )
  } else {
    result <- register_to_baseline(
      moving_img, fixed_img,
      registration_type = registration_type,
      cost_function = cost_function
    )
  }
  
  return(result)
}

create_registration_qc <- function(moving_img, fixed_img, registered_img, output_path) {
  #Create quality control images for registration assessment
  
  png(output_path, width = 1800, height = 600)
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
  
  try(ortho2(moving_img, text = "Followup Visit", xyz = c(128, 128, 170)), silent = TRUE)
  
  try(ortho2(fixed_img, text = "Baseline", xyz = c(128, 128, 170)), silent = TRUE)
  
  try(ortho2(registered_img, text = "Registered", xyz = c(128, 128, 170)), silent = TRUE)
  
  dev.off()
  message("Registration QC image saved: ", output_path)
}

process_visit_registration <- function(visit_path, baseline_path, output_dir, 
                                     registration_type = "linear", create_qc = TRUE) {
  #Process registration for a single visit
  
  visit_name <- basename(dirname(visit_path))
  
  message("Loading images...")
  moving_img <- readNIfTI(visit_path, reorient = FALSE)
  fixed_img <- readNIfTI(baseline_path, reorient = FALSE)
  
  message("Registering ", visit_name, " to baseline...")
  reg_result <- register_to_baseline(
    moving_img, fixed_img,
    registration_type = registration_type
  )

  output_file <- file.path(output_dir, paste0(visit_name, "_registered.nii.gz"))
  writeNIfTI(reg_result$final_registered, output_file)
  message("Registered image saved: ", output_file)
  
  if (!is.null(reg_result$linear_matrix)) {
    matrix_file <- file.path(output_dir, paste0(visit_name, "_transform.mat"))
    write.table(reg_result$linear_matrix, matrix_file, row.names = FALSE, col.names = FALSE)
  }
  
  quality_metrics <- calculate_image_similarity(reg_result$final_registered, fixed_img)
  quality_file <- file.path(output_dir, paste0(visit_name, "_registration_quality.csv"))
  write.csv(data.frame(
    Visit = visit_name,
    Correlation = quality_metrics$correlation,
    NMI = quality_metrics$normalized_mutual_information,
    MSE = quality_metrics$mse,
    SSIM = quality_metrics$ssim
  ), quality_file, row.names = FALSE)
  
  if (create_qc) {
    qc_file <- file.path(output_dir, paste0(visit_name, "_registration_qc.png"))
    create_registration_qc(moving_img, fixed_img, reg_result$final_registered, qc_file)
  }
  
  return(list(
    registered_image = reg_result$final_registered,
    quality_metrics = quality_metrics,
    output_file = output_file
  ))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 3) {
    cat("Usage: Rscript 05_registration.R <moving_image> <fixed_image> <output_file> [registration_type] [cost_function]\n")
    cat("registration_type: linear, nonlinear, both (default: linear)\n")
    cat("cost_function: corratio, mutualinfo, normmi, normcorr, leastsq (default: corratio)\n")
    quit(status = 1)
  }
  
  moving_file <- args[1]
  fixed_file <- args[2]
  output_file <- args[3]
  registration_type <- if (length(args) >= 4) args[4] else "linear"
  cost_function <- if (length(args) >= 5) args[5] else "corratio"
  
  if (!file.exists(moving_file)) {
    stop("Moving image file does not exist: ", moving_file)
  }
  
  if (!file.exists(fixed_file)) {
    stop("Fixed image file does not exist: ", fixed_file)
  }

  moving_img <- readNIfTI(moving_file, reorient = FALSE)
  fixed_img <- readNIfTI(fixed_file, reorient = FALSE)
  
  message("Performing registration...")
  reg_result <- register_to_baseline(
    moving_img, fixed_img,
    registration_type = registration_type,
    cost_function = cost_function
  )

  writeNIfTI(reg_result$final_registered, output_file)
  
  quality_metrics <- calculate_image_similarity(reg_result$final_registered, fixed_img)
  quality_file <- paste0(tools::file_path_sans_ext(output_file), "_quality.csv")
  write.csv(data.frame(
    Correlation = quality_metrics$correlation,
    NMI = quality_metrics$normalized_mutual_information,
    MSE = quality_metrics$mse,
    SSIM = quality_metrics$ssim
  ), quality_file, row.names = FALSE)
  
  message("Quality metrics saved: ", quality_file)
}
