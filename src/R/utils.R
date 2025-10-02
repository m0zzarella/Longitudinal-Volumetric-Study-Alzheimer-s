
suppressMessages({
  library(oro.nifti)
  library(neurobase)
  library(fslr)
  library(scales)
})

#compute volume from a PVE map with thresholding
compute_volume_ml <- function(pve_img, threshold = 0.33) {
  vdim <- prod(voxdim(pve_img))         #mm^3 per voxel
  nvox <- sum(pve_img > threshold)      #voxel count
  vol_ml <- (vdim * nvox) / 1000        #convert mm^3 to mL
  return(vol_ml)
}

#robust skull stripping with optional COG initialization
run_bet <- function(nifti_img, cog_init = TRUE) {
  if (cog_init) {
    cog_coords <- cog(nifti_img, ceil = TRUE)
    opts <- paste("-c", paste(cog_coords, collapse = " "))
  } else {
    opts <- ""
  }
  fslr::fslbet(infile = nifti_img, retimg = TRUE, opts = opts)
}


find_visit_nifti <- function(visit_dir) {
  preferred <- file.path(visit_dir, "scan.nii.gz")
  if (file.exists(preferred)) return(preferred)

  candidates <- list.files(visit_dir, pattern = "\\.nii(\\.gz)?$", full.names = TRUE)
  if (length(candidates) == 0) return(NA_character_)
  return(candidates[1])
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

needs_intensity_normalization <- function(img, threshold_cv = 0.5) {
  #Check if image needs intensity normalization based on coefficient of variation
  #Args:
  #  img: NIfTI image object
  #  threshold_cv: coefficient of variation threshold
  
  #Returns:
  #  TRUE if normalization is needed, FALSE otherwise
  brain_voxels <- img[img > 0]
  if (length(brain_voxels) == 0) return(FALSE)
  
  cv <- sd(brain_voxels, na.rm = TRUE) / mean(brain_voxels, na.rm = TRUE)
  return(cv > threshold_cv)
}

get_image_stats <- function(img, mask = NULL) {
  #Basic image statistics
  if (!is.null(mask)) {
    voxels <- img[mask > 0]
  } else {
    voxels <- img[img > 0]
  }
  
  if (length(voxels) == 0) {
    return(list(
      mean = NA, median = NA, sd = NA, min = NA, max = NA,
      q25 = NA, q75 = NA, cv = NA, n_voxels = 0
    ))
  }
  
  stats <- list(
    mean = mean(voxels, na.rm = TRUE),
    median = median(voxels, na.rm = TRUE),
    sd = sd(voxels, na.rm = TRUE),
    min = min(voxels, na.rm = TRUE),
    max = max(voxels, na.rm = TRUE),
    q25 = quantile(voxels, 0.25, na.rm = TRUE),
    q75 = quantile(voxels, 0.75, na.rm = TRUE),
    n_voxels = length(voxels)
  )
  
  stats$cv <- stats$sd / stats$mean
  
  return(stats)
}

create_brain_mask <- function(img, threshold = 0.1) {
  mask <- img
  mask[mask > threshold] <- 1
  mask[mask <= threshold] <- 0
  return(mask)
}

apply_transformation <- function(coords, transform_matrix) {
  #Apply transformation matrix to image coordinates (Nx3 via 4x4 homogeneous transform)
  #Args:
  #  coords: Nx3 matrix of coordinates
  #  transform_matrix: 4x4 transformation matrix
  #Returns:
  #  Transformed coordinates
  
  homo_coords <- cbind(coords, 1)
  transformed <- t(transform_matrix %*% t(homo_coords))
  return(transformed[, 1:3])
}

calculate_normalized_mutual_information <- function(x, y, bins = 50) {
  #Calculate normalized mutual information between two vectors
  x_range <- range(x, na.rm = TRUE)
  y_range <- range(y, na.rm = TRUE)
  
  x_bins <- seq(x_range[1], x_range[2], length.out = bins + 1)
  y_bins <- seq(y_range[1], y_range[2], length.out = bins + 1)
  
  x_idx <- cut(x, x_bins, include.lowest = TRUE, labels = FALSE)
  y_idx <- cut(y, y_bins, include.lowest = TRUE, labels = FALSE)
  
  joint_hist <- table(x_idx, y_idx)
  joint_prob <- joint_hist / sum(joint_hist)
  
  px <- rowSums(joint_prob)
  py <- colSums(joint_prob)
  
  hx <- -sum(px * log(px + 1e-10))
  hy <- -sum(py * log(py + 1e-10))
  hxy <- -sum(joint_prob * log(joint_prob + 1e-10))
  
  nmi <- 2 * (hx + hy - hxy) / (hx + hy)
  return(nmi)
}

calculate_image_similarity <- function(img1, img2, mask = NULL) {
  #Calculate similarity metrics between two images  
  #Args:
  #  img1: NIfTI image object
  #  img2: NIfTI image object
  #  mask: NIfTI image object
  #Returns:
  #  list of similarity metrics
  if (!is.null(mask)) {
    vox1 <- img1[mask > 0]
    vox2 <- img2[mask > 0]
  } else {
    non_zero <- img1 > 0 & img2 > 0
    vox1 <- img1[non_zero]
    vox2 <- img2[non_zero]
  }
  
  if (length(vox1) == 0 || length(vox2) == 0) {
    return(list(correlation = NA, mse = NA, mae = NA, ssim = NA, normalized_mutual_information = NA))
  }
  
  correlation <- cor(vox1, vox2, use = "complete.obs")
  mse <- mean((vox1 - vox2)^2, na.rm = TRUE)
  mae <- mean(abs(vox1 - vox2), na.rm = TRUE)
  nmi <- calculate_normalized_mutual_information(vox1, vox2)
  
  mu1 <- mean(vox1, na.rm = TRUE)
  mu2 <- mean(vox2, na.rm = TRUE)
  sigma1 <- sd(vox1, na.rm = TRUE)
  sigma2 <- sd(vox2, na.rm = TRUE)
  sigma12 <- cov(vox1, vox2, use = "complete.obs")
  
  c1 <- 0.01^2
  c2 <- 0.03^2
  
  ssim <- ((2 * mu1 * mu2 + c1) * (2 * sigma12 + c2)) / 
          ((mu1^2 + mu2^2 + c1) * (sigma1^2 + sigma2^2 + c2))
  
  return(list(
    correlation = correlation,
    mse = mse,
    mae = mae,
    ssim = ssim,
    normalized_mutual_information = nmi
  ))
}

validate_registration <- function(registered_img, reference_img, mask = NULL, 
                                min_correlation = 0.7, max_mse = 1000) {
  #Validate registration quality based on similarity metrics
  #Args:
  #  registered_img: NIfTI image object
  #  reference_img: NIfTI image object
  #  mask: NIfTI image object
  #  min_correlation: minimum correlation threshold
  #  max_mse: maximum MSE threshold
  #Returns:
  #  list of validation results
  similarity <- calculate_image_similarity(registered_img, reference_img, mask)
  
  is_valid <- !is.na(similarity$correlation) && 
              similarity$correlation >= min_correlation &&
              !is.na(similarity$mse) &&
              similarity$mse <= max_mse
  
  return(list(
    is_valid = is_valid,
    metrics = similarity,
    thresholds = list(min_correlation = min_correlation, max_mse = max_mse)
  ))
}

create_log_entry <- function(step_name, input_file, output_file, parameters = NULL, 
                           quality_metrics = NULL, timestamp = Sys.time()) {
  #Create standardized log entry for processing steps
  
  log_entry <- list(
    timestamp = as.character(timestamp),
    step = step_name,
    input_file = input_file,
    output_file = output_file,
    parameters = parameters,
    quality_metrics = quality_metrics
  )
  
  return(log_entry)
}


# TODO: Segmentation quality (try if possible)
#
# maybe provide intrinsic metrics to assess the FAST segmentations
# (CSF/GM/WM) without manual annotations reliability
#
# some rough basic ideas for helper functions to implement in later
#
# 1) PosteriorEntropy(pve_csf, pve_gm, pve_wm)
#      Inputs: three NIfTI probability maps aligned voxelwise.
#      For each voxel v: H(v) = -sum_c p_c(v) * log(p_c(v) + 1e-12).
#      Return: list(mean_entropy, median_entropy, p95_entropy, entropy_map [maybe return the original NIfTI as well]).
#      Pseudocode:
#        p_stack <- abind::abind(pve_csf[], pve_gm[], pve_wm[], along=4)
#        H <- -rowSums(p_stack * log(p_stack + 1e-12))
#        stats <- c(mean=mean(H[brain_mask]), median=median(...), p95=quantile(...))
#        return(list(stats=stats, map=entropy_img))
#
# 2) BhattacharyyaDistanceFromHist(x, y, bins=128)
#    - bbuild normalized histograms for x and y; compute DB = -ln(sum_sqrt(p_i * q_i)).
#    - use brain-masked intensities or tissue-weighted samples.
#
# 3) ClassSeparability(intensity_img, pve_csf, pve_gm, pve_wm, bins=128)
#    - compute pairwise for bhattacharyya distance, js divergence etc
#    - return minimal separability across pairs and per-pair values
#
# 4) BoundaryGradientScore(intensity_img, pve_csf, pve_gm, pve_wm, iso=0.33)
#    - Derive binary masks at iso thresholds; compute boundary voxels via morphological edge.
#    - Compute image gradient magnitude (3D Sobel or finite diff - idk which better) at boundary voxels.
#    - Score = mean gradient along all tissue boundaries (higher is better).
#
# 5) ICVSanity(pve_csf, pve_gm, pve_wm, threshold=0.33)
#    - Estimate ICV by thresholding each PVE and add the 3 volumes
#    - Return the total intracranial volume and plausibility flags for outliers etc
#
# 6) SoftDiceProb(p_c_native, p_c_ref)
#    -softdice for probabilities: 2 * sum(p*q) / (sum(p^2) + sum(q^2) + 1e-12).
#    - return per-tissue and average vectorized across all voxels
#
# 7) can compute the surface hausdorff distance too
#
