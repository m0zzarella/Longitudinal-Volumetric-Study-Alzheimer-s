#!/usr/bin/env Rscript

# Mann-Kendall Test for Longitudinal Volume Analysis
# 
# This script calculates the Mann-Kendall test statistics for each patient's 
# longitudinal volume data to determine if there are statistically significant 
# increasing or decreasing trends over time.
#
# The test uses a modified sign function with a 1mL threshold to ignore 
# minute changes and focus on clinically meaningful volume differences.
#
# Usage: Rscript 06_mann_kendall_analysis.R /path/to/results/directory [patient_id]
#        If patient_id is provided, analyzes only that patient
#        If not provided, analyzes all patients in the results directory

suppressMessages({
  library(data.table)
  if (!require("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cran.r-project.org/")
    library(data.table)
  }
  library(yaml)
})

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

config <- load_config()

modified_sign <- function(x_j, x_k, threshold = NULL) {
  if (is.null(threshold)) {
    threshold <- config$processing$mann_kendall_threshold
  }

  diff <- x_j - x_k
  ifelse(diff > threshold, 1, ifelse(diff < -threshold, -1, 0))
}

calculate_mk_statistic <- function(volumes) {
  n <- length(volumes)
  if (n < 2) return(list(S = NA, n = n))
  
  S <- 0
  for (k in 1:(n-1)) {
    for (j in (k+1):n) {
      S <- S + modified_sign(volumes[j], volumes[k])
    }
  }
  
  return(list(S = S, n = n))
}

calculate_mk_variance <- function(n) {
  if (n < 2) return(NA)
  var_S <- (n * (n - 1) * (2 * n + 5)) / 18
  return(var_S)
}

calculate_z_score <- function(S, var_S) {
  if (is.na(S) || is.na(var_S) || var_S <= 0) return(NA)
  
  sqrt_var <- sqrt(var_S)
  
  if (S > 0) {
    Z <- (S - 1) / sqrt_var
  } else if (S < 0) {
    Z <- (S + 1) / sqrt_var
  } else {
    Z <- 0
  }
  
  return(Z)
}

interpret_trend <- function(Z_vol, alpha = 0.05) {
  if (is.na(Z_vol)) return("Insufficient data")
  
  # Critical value for two-tailed test at alpha = 0.05 is ±1.96
  critical_value <- qnorm(1 - alpha/2)
  
  if (abs(Z_vol) >= critical_value) {
    if (Z_vol > 0) {
      return("Significant increase")
    } else {
      return("Significant decrease")
    }
  } else {
    if (Z_vol > 0) {
      return("Non-significant increase")
    } else if (Z_vol < 0) {
      return("Non-significant decrease")
    } else {
      return("No trend")
    }
  }
}

calculate_p_value <- function(Z_vol) {
  if (is.na(Z_vol)) return(NA)
  p_value <- 2 * (1 - pnorm(abs(Z_vol)))
  return(p_value)
}

calculate_confidence_percentage <- function(Z_vol) {
  if (is.na(Z_vol)) return(NA)
  
  #For one-tailed test (directional hypothesis)
  #Confidence = P(Z <= |Z_vol|) * 100
  confidence_pct <- pnorm(abs(Z_vol)) * 100
  
  return(confidence_pct)
}

analyze_patient_tissue <- function(volumes, tissue_name, patient_id, visit_names) {
  mk_result <- calculate_mk_statistic(volumes)
  S <- mk_result$S
  n <- mk_result$n
  
  var_S <- calculate_mk_variance(n)
  Z_vol <- calculate_z_score(S, var_S)
  p_value <- calculate_p_value(Z_vol)
  confidence_pct <- calculate_confidence_percentage(Z_vol)
  trend <- interpret_trend(Z_vol)
  
  tau <- ifelse(n >= 2, S / (n * (n - 1) / 2), NA)
  
  result <- data.frame(
    Patient_ID = patient_id,
    Tissue_Type = tissue_name,
    N_Visits = n,
    S_Statistic = S,
    Variance_S = var_S,
    Z_Score = Z_vol,
    P_Value = p_value,
    Confidence_Percentage = confidence_pct,
    Kendall_Tau = tau,
    Trend_Interpretation = trend,
    Visit_Sequence = paste(visit_names, collapse = ", "),
    Volume_Sequence = paste(round(volumes, 2), collapse = ", "),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

analyze_patient_mk <- function(patient_dir) {
  patient_id <- basename(patient_dir)
  message("\n=== Analyzing patient: ", patient_id, " ===")
  
  long_volumes_file <- file.path(patient_dir, paste0(patient_id, "_longitudinal_volumes.csv"))
  
  if (!file.exists(long_volumes_file)) {
    warning("Longitudinal volumes file not found for patient ", patient_id, ": ", long_volumes_file)
    return(NULL)
  }
  
  vol_data <- read.csv(long_volumes_file, stringsAsFactors = FALSE)
  
  if (nrow(vol_data) < 2) {
    warning("Patient ", patient_id, " has fewer than 2 visits. Cannot perform Mann-Kendall test.")
    return(NULL)
  }
  
  visit_order <- order(as.numeric(gsub("Visit", "", vol_data$Visit)))
  vol_data <- vol_data[visit_order, ]
  
  message("Found ", nrow(vol_data), " visits: ", paste(vol_data$Visit, collapse = ", "))
  message("Volume data preview:")
  print(vol_data)
  
  results_list <- list()
  
  if ("CSF_ml" %in% colnames(vol_data)) {
    csf_result <- analyze_patient_tissue(
      vol_data$CSF_ml, "CSF", patient_id, vol_data$Visit
    )
    results_list[["CSF"]] <- csf_result
    message("CSF: S=", csf_result$S_Statistic, ", Z=", round(csf_result$Z_Score, 3), 
            ", p=", round(csf_result$P_Value, 4), ", confidence=", round(csf_result$Confidence_Percentage, 2), "%")
  }
  
  if ("GM_ml" %in% colnames(vol_data)) {
    gm_result <- analyze_patient_tissue(
      vol_data$GM_ml, "GM", patient_id, vol_data$Visit
    )
    results_list[["GM"]] <- gm_result
    message("GM:  S=", gm_result$S_Statistic, ", Z=", round(gm_result$Z_Score, 3), 
            ", p=", round(gm_result$P_Value, 4), ", confidence=", round(gm_result$Confidence_Percentage, 2), "%")
  }
  
  if ("WM_ml" %in% colnames(vol_data)) {
    wm_result <- analyze_patient_tissue(
      vol_data$WM_ml, "WM", patient_id, vol_data$Visit
    )
    results_list[["WM"]] <- wm_result
    message("WM:  S=", wm_result$S_Statistic, ", Z=", round(wm_result$Z_Score, 3), 
            ", p=", round(wm_result$P_Value, 4), ", confidence=", round(wm_result$Confidence_Percentage, 2), "%")
  }
  
  if (length(results_list) > 0) {
    combined_results <- do.call(rbind, results_list)
    return(combined_results)
  } else {
    warning("No valid tissue volume columns found for patient ", patient_id)
    return(NULL)
  }
}

generate_summary <- function(all_results) {
  if (is.null(all_results) || nrow(all_results) == 0) {
    message("No results to summarize.")
    return(NULL)
  }
  
  total_analyses <- nrow(all_results)
  significant_results <- sum(grepl("Significant", all_results$Trend_Interpretation), na.rm = TRUE)
  
  message("Total analyses performed: ", total_analyses)
  message("Significant trends found: ", significant_results, " (", 
          round(100 * significant_results / total_analyses, 1), "%)")
  
  message("\nBy tissue type:")
  tissue_summary <- aggregate(cbind(Z_Score, P_Value) ~ Tissue_Type, 
                             data = all_results, 
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                               median = median(x, na.rm = TRUE)))
  
  for (tissue in unique(all_results$Tissue_Type)) {
    tissue_data <- all_results[all_results$Tissue_Type == tissue, ]
    n_patients <- nrow(tissue_data)
    n_sig <- sum(grepl("Significant", tissue_data$Trend_Interpretation), na.rm = TRUE)
    n_inc <- sum(grepl("increase", tissue_data$Trend_Interpretation), na.rm = TRUE)
    n_dec <- sum(grepl("decrease", tissue_data$Trend_Interpretation), na.rm = TRUE)
    
    message(sprintf("  %s: %d patients, %d significant (%d increasing, %d decreasing)", 
                   tissue, n_patients, n_sig, n_inc, n_dec))
  }
  
  patient_sig_counts <- aggregate(grepl("Significant", all_results$Trend_Interpretation) ~ Patient_ID, 
                                 data = all_results, FUN = sum)
  multi_sig_patients <- patient_sig_counts[patient_sig_counts[,2] > 1, ]
  
  if (nrow(multi_sig_patients) > 0) {
    message("\nPatients with multiple significant trends:")
    for (i in 1:nrow(multi_sig_patients)) {
      message("  ", multi_sig_patients[i,1], ": ", multi_sig_patients[i,2], " significant trends")
    }
  }
  
  return(tissue_summary)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    stop("Usage: Rscript 06_mann_kendall_analysis.R /path/to/results/directory [patient_id]")
  }
  
  results_dir <- args[1]
  specific_patient <- if (length(args) >= 2) args[2] else NULL
  
  if (!dir.exists(results_dir)) {
    stop("Results directory does not exist: ", results_dir)
  }
  
  if (!is.null(specific_patient)) {
    patient_dirs <- file.path(results_dir, specific_patient)
    if (!dir.exists(patient_dirs)) {
      stop("Patient directory not found: ", patient_dirs)
    }
    patient_dirs <- list(patient_dirs)
  } else {
    all_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
    patient_dirs <- all_dirs[file.exists(file.path(all_dirs, paste0(basename(all_dirs), "_longitudinal_volumes.csv")))]
    
    if (length(patient_dirs) == 0) {
      stop("No patient directories with longitudinal volume data found in: ", results_dir)
    }
  }
  
  message("Found ", length(patient_dirs), " patient(s) to analyze")
  
  all_results <- list()
  
  for (patient_dir in patient_dirs) {
    result <- analyze_patient_mk(patient_dir)
    if (!is.null(result)) {
      all_results[[basename(patient_dir)]] <- result
    }
  }
  
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    rownames(combined_results) <- NULL
    
    output_file <- file.path(results_dir, "mann_kendall_analysis_results.csv")
    write.csv(combined_results, output_file, row.names = FALSE)
    
    summary_stats <- generate_summary(combined_results)
    
    summary_file <- file.path(results_dir, "mann_kendall_summary.txt")
    sink(summary_file)
    cat("Mann-Kendall Trend Analysis Summary\n")
    cat("Analysis performed on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Results directory:", results_dir, "\n")
    cat("Number of patients analyzed:", length(unique(combined_results$Patient_ID)), "\n\n")
    
    print(combined_results)
    sink()
    
    message("Summary written to: ", summary_file)
    
    return(combined_results)
  } else {
    message("No valid results obtained.")
    return(NULL)
  }
}

if (!interactive()) {
  main()
}
