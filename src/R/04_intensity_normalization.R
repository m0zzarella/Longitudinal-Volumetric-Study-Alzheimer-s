#!/usr/bin/env Rscript

# Intensity Normalization Script
# This script performs intensity normalization on NIfTI files using various methods
# It can be run standalone or integrated into the longitudinal pipeline

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

normalize_intensity <- function(img, method = "zscore", mask = NULL, reference_img = NULL) {
  #Normalize image intensities using various methods
  #Args:
  #  img: NIfTI image object
  #  method: normalization method
  #  mask: NIfTI image object
  #  reference_img: NIfTI image object
  #Returns:
  #  normalized NIfTI image object
  if (!is.null(mask)) {
    # Apply mask to get brain-only voxels
    brain_voxels <- img[mask > 0]
  } else {
    # Use all non-zero voxels
    brain_voxels <- img[img > 0]
  }
  
  if (length(brain_voxels) == 0) {
    warning("No brain voxels found for normalization")
    return(img)
  }
  
  normalized_img <- img
  
  switch(method,
    "zscore" = {
      #Z-score normalization (mean=0, std=1)
      mean_val <- mean(brain_voxels, na.rm = TRUE)
      std_val <- sd(brain_voxels, na.rm = TRUE)
      
      if (std_val > 0) {
        if (!is.null(mask)) {
          normalized_img[mask > 0] <- (img[mask > 0] - mean_val) / std_val
        } else {
          normalized_img[img > 0] <- (img[img > 0] - mean_val) / std_val
        }
      }
    },
    
    "minmax" = {
      #Min-max normalization (0-1 range)
      min_val <- min(brain_voxels, na.rm = TRUE)
      max_val <- max(brain_voxels, na.rm = TRUE)
      
      if (max_val > min_val) {
        if (!is.null(mask)) {
          normalized_img[mask > 0] <- (img[mask > 0] - min_val) / (max_val - min_val)
        } else {
          normalized_img[img > 0] <- (img[img > 0] - min_val) / (max_val - min_val)
        }
      }
    },
    
    "histogram_match" = {
      # Histogram matching to reference image
      if (is.null(reference_img)) {
        stop("Reference image required for histogram matching")
      }
      
      if (!is.null(mask)) {
        ref_voxels <- reference_img[mask > 0]
      } else {
        ref_voxels <- reference_img[reference_img > 0]
      }
      
      normalized_img <- histogram_match_image(img, brain_voxels, ref_voxels, mask)
    },
    
    "nyul" = {
      # Nyul and Udupa intensity normalization
      normalized_img <- nyul_normalize(img, mask)
    },
    
    {
      warning("Unknown normalization method: ", method, ". Using zscore.")
      # Default to z-score
      mean_val <- mean(brain_voxels, na.rm = TRUE)
      std_val <- sd(brain_voxels, na.rm = TRUE)
      
      if (std_val > 0) {
        if (!is.null(mask)) {
          normalized_img[mask > 0] <- (img[mask > 0] - mean_val) / std_val
        } else {
          normalized_img[img > 0] <- (img[img > 0] - mean_val) / std_val
        }
      }
    }
  )
  
  return(normalized_img)
}

# Histogram matching implementation
histogram_match_image <- function(source_img, source_voxels, reference_voxels, mask = NULL) {
  #Match histogram of source image to reference image
  #Args:
  #  source_img: NIfTI image object
  #  source_voxels: vector of voxel values
  #  reference_voxels: vector of voxel values
  #  mask: NIfTI image object
  #Returns:
  #  matched NIfTI image object
  
  # Get CDFs
  source_hist <- hist(source_voxels, breaks = 256, plot = FALSE)
  ref_hist <- hist(reference_voxels, breaks = 256, plot = FALSE)
  
  source_cdf <- cumsum(source_hist$counts) / sum(source_hist$counts)
  ref_cdf <- cumsum(ref_hist$counts) / sum(ref_hist$counts)
  
  # Create mapping
  matched_img <- source_img
  
  for (i in seq_along(source_hist$mids)) {
    source_val <- source_hist$mids[i]
    source_prob <- source_cdf[i]
    
    # Find closest CDF value in reference
    ref_idx <- which.min(abs(ref_cdf - source_prob))
    ref_val <- ref_hist$mids[ref_idx]
    
    # Map values
    if (!is.null(mask)) {
      matched_img[mask > 0 & abs(source_img - source_val) < diff(source_hist$breaks)[1]/2] <- ref_val
    } else {
      matched_img[abs(source_img - source_val) < diff(source_hist$breaks)[1]/2] <- ref_val
    }
  }
  
  return(matched_img)
}

# Nyul normalization (simplified implementation)
nyul_normalize <- function(img, mask = NULL) {
  # Nyul and Udupa intensity normalization using percentile-based landmarks
  #Args:
  #  img: NIfTI image object
  #  mask: NIfTI image object
  #Returns:
  #  normalized NIfTI image object
  if (!is.null(mask)) {
    brain_voxels <- img[mask > 0]
  } else {
    brain_voxels <- img[img > 0]
  }
  
  if (length(brain_voxels) == 0) {
    return(img)
  }

  percentiles <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
  landmarks <- quantile(brain_voxels, percentiles, na.rm = TRUE)

  standard_landmarks <- c(10, 20, 40, 80, 120, 160, 200, 220, 240)

  normalized_img <- img
  
  if (!is.null(mask)) {
    brain_mask <- mask > 0
    brain_intensities <- img[brain_mask]
  } else {
    brain_mask <- img > 0
    brain_intensities <- img[brain_mask]
  }
  
  # Apply piecewise linear transformation
  normalized_intensities <- approx(
    x = landmarks,
    y = standard_landmarks,
    xout = brain_intensities,
    method = "linear",
    rule = 2
  )$y
  
  normalized_img[brain_mask] <- normalized_intensities
  
  return(normalized_img)
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript 04_intensity_normalization.R <input_file> <output_file> [method] [mask_file] [reference_file]\n")
    cat("Methods: zscore, minmax, histogram_match, nyul (default: zscore)\n")
    quit(status = 1)
  }
  
  input_file <- args[1]
  output_file <- args[2]
  method <- if (length(args) >= 3) args[3] else "zscore"
  mask_file <- if (length(args) >= 4) args[4] else NULL
  reference_file <- if (length(args) >= 5) args[5] else NULL
  
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  message("Loading input image: ", input_file)
  img <- readNIfTI(input_file, reorient = FALSE)
  
  mask <- NULL
  if (!is.null(mask_file) && file.exists(mask_file)) {
    message("Loading mask: ", mask_file)
    mask <- readNIfTI(mask_file, reorient = FALSE)
  }
  
  reference_img <- NULL
  if (!is.null(reference_file) && file.exists(reference_file)) {
    message("Loading reference image: ", reference_file)
    reference_img <- readNIfTI(reference_file, reorient = FALSE)
  }
  
  message("Normalizing intensities using method: ", method)
  normalized_img <- normalize_intensity(img, method = method, mask = mask, reference_img = reference_img)
  
  message("Saving normalized image: ", output_file)
  writeNIfTI(normalized_img, output_file)
  
  message("Intensity normalization completed successfully.")
}
