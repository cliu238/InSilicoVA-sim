# Load packages
library(openVA)
library(MCMCpack)
source("functions.R")

# Create output directory if it does not exist
if (!dir.exists("rda")) {
  dir.create("rda")
}

# Read PHMRC data from local CSV file
PHMRC <- read.csv("../PHMRC/IHME_PHMRC_VA_DATA_ADULT_Y2013M09D11_0.csv")
causes <- as.character(unique(PHMRC$gs_text34))
C <- length(causes)  # number of causes

# New splitting logic based on 'site' column:
#   - 50% of records with site "UP" are used for test set.
#   - 10% of records with site "UP" are used for validation set.
#   - The remaining 40% of records with site "UP" are used for unlabelled training.
#   - Records with sites in c("AP", "Bohol", "Dar", "Mexico", "Pemba") are used for labelled training.
#   The final training set is the union of the labelled training and unlabelled "UP" records.

for(rep in 1:2) {
  
  # For the first replication, use HCE; for later replications, do not.
  if(rep > 1) {
    rep_index <- rep - 1
    has.HCE <- FALSE
  } else {
    rep_index <- rep
    has.HCE <- TRUE
  }
  
  set.seed(123 + rep_index)
  
  ## Splitting based on 'site'
  up_indices <- which(PHMRC$site == "UP")
  n_up <- length(up_indices)
  
  # 50% of "UP" records for test
  test_up <- sample(up_indices, size = round(0.5 * n_up))
  
  # 10% of "UP" records for validation (sampled from remaining "UP" records)
  remaining_up <- setdiff(up_indices, test_up)
  validation_up <- sample(remaining_up, size = round(0.1 * n_up))
  
  # Remaining "UP" records for unlabelled training
  unlabelled_up <- setdiff(up_indices, c(test_up, validation_up))
  
  # Labelled training: records with sites in the specified list
  labelled_idx <- which(PHMRC$site %in% c("AP", "Bohol", "Dar", "Mexico", "Pemba"))
  
  # Final training indices: union of labelled and unlabelled "UP" records
  train_idx <- c(labelled_idx, unlabelled_up)
  
  # Create final datasets
  train <- PHMRC[train_idx, ]
  test <- PHMRC[test_up, ]
  validation <- PHMRC[validation_up, ]
  
  ## Resample within test data to obtain a new CSMF:
  csmf.resample <- rdirichlet(1, rep(1, C))
  names(csmf.resample) <- causes
  
  # Remove causes not present in the test data
  missing_causes <- which(!(causes %in% test$gs_text34))
  if(length(missing_causes) > 0) {
    csmf.resample[missing_causes] <- 0
    csmf.resample <- csmf.resample / sum(csmf.resample)
  }
  
  # Create new cause assignment for test data based on the resampled CSMF
  causes.test <- sample(size = nrow(test), x = 1:length(causes), replace = TRUE, prob = csmf.resample)
  
  # Build a list of indices for each cause in test data
  testlist <- vector("list", C)
  for(i in 1:C) {
    testlist[[i]] <- which(test$gs_text34 == causes[i])
  }
  
  # For each test observation, sample one index from the appropriate cause pool
  test.index <- rep(NA, nrow(test))
  for(i in 1:nrow(test)) {
    test.index[i] <- sample(testlist[[causes.test[i]]], 1)
  }
  
  # Subset test data and compute the true CSMF from test data
  test <- test[test.index, ]
  csmf.true <- table(test$gs_text34)
  csmf.true <- csmf.true[causes]
  csmf.true[is.na(csmf.true)] <- 0
  csmf.true <- csmf.true / sum(csmf.true)
  
  ## Run InSilicoVA depending on HCE usage
  if(!has.HCE) {
    binary <- ConvertData.phmrc(train, test, phmrc.type = "adult", cause = "gs_text34")
    reduced <- paste0("a1_01_", 1:14)
    binary$output <- binary$output[, -which(colnames(binary$output) %in% reduced)]
    binary$output.test <- binary$output.test[, -which(colnames(binary$output.test) %in% reduced)]
    
    insilico <- codeVA(data = binary$output.test, 
                       data.type = "customize", 
                       model = "InSilicoVA",
                       data.train = binary$output, 
                       causes.train = "Cause", 
                       jump.scale = 0.05, 
                       convert.type = "fixed",
                       Nsim = 10000, auto.length = FALSE)
  } else {
    insilico <- codeVA(data = test,
                       data.type = "PHMRC",
                       model = "InSilicoVA",
                       data.train = train,
                       causes.train = "gs_text34",
                       phmrc.type = "adult",
                       jump.scale = 0.05,
                       convert.type = "fixed",
                       Nsim = 10000, auto.length = FALSE)
  }
  
  ## Retrieve estimated CSMF and top cause assignments
  csmf <- getCSMF(insilico)[, "Mean"]
  csmf <- csmf[causes]  # ensure the order follows the original cause list
  topcause <- getTopCOD(insilico)
  
  ## Align estimated CSMF and true CSMF to have identical cause names:
  all_causes <- union(names(csmf), names(csmf.true))
  
  missing_in_est <- setdiff(all_causes, names(csmf))
  if(length(missing_in_est) > 0) {
    csmf[missing_in_est] <- 0
  }
  
  missing_in_true <- setdiff(all_causes, names(csmf.true))
  if(length(missing_in_true) > 0) {
    csmf.true[missing_in_true] <- 0
  }
  
  # Order both vectors by the same cause names
  csmf <- csmf[all_causes]
  csmf.true <- csmf.true[all_causes]
  
  ## Calculate CSMF accuracy and normalized CCCSMF accuracy
  csmfacc <- getCSMF_accuracy(csmf = csmf, truth = csmf.true)
  cccsmfacc <- (csmfacc - (1 - exp(-1))) / (1 - (1 - exp(-1)))
  
  ## Compute Cause-specific Concordance (CCC) for individual assignments
  ccc <- getCCC(est = topcause[,2], truth = as.character(test$gs_text34), causes = causes)
  names(ccc) <- causes
  
  # Compute individual-level accuracy: the proportion of individuals with correct top cause assignment
  indiv_accuracy <- mean(topcause[,2] == as.character(test$gs_text34))
  
  metrics <- list(ccc = ccc, csmfacc = csmfacc, cccsmfacc = cccsmfacc, indiv_accuracy = indiv_accuracy)
  
  # --- Output accuracy metrics to the console ---
  cat("CSMF Accuracy:", csmfacc, "\n")
  cat("Normalized CCCSMF Accuracy:", cccsmfacc, "\n")
  cat("Cause-specific Concordance (CCC):\n")
  print(ccc)
  cat("Individual Level Accuracy:", indiv_accuracy, "\n")
  # --- End output ---
  
  # Optionally, to check convergence issues, you can run:
  # diag_info <- csmf.diag(insilico)
  # print(diag_info)
  
  # Save the results to an .rda file based on whether HCE was used
  if(has.HCE) {
    save(metrics, file = paste0("rda/full-HCE-rep", rep_index, ".rda"))
  } else {
    save(metrics, file = paste0("rda/full-noHCE-rep", rep_index, ".rda"))
  }
}