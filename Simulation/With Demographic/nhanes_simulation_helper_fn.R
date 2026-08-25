library(tidyverse)

# -----------------------------------------------------------------------------
# Helper Functions for Non-Parametric Simulation
# 
# Description:
# Provides utility functions to calculate dissimilarity scores and probability 
# matrices for individuals based on demographic, ordinal, and categorical variables.
# Also includes functions to resample dietary data according to these probabilities.
# 
# Requirements:
# - Required packages: tidyverse
# -----------------------------------------------------------------------------



# Compute (dis)similarity score.
# Takes as input: 
#     @data (the dataframe to build from)
#     @base_row (the row index the get similarities to),
#     @compare_rows (the row indices to compare @base_row to),
#     @numeric_var_names (the names of continuous vars for comparison)
#     @ordinal_var_names (the names of ordinal vars for comparison)
#     @categorical_var_names (the names of categorical vars for comparison)
# 
# Returns
#     A numeric vector with a length equal to the number of compare_rows which
#     computes a (standardized) dissimilarity score -- higher values are less
#     similar.
row_dissimilarity_score <- function(data,
                                    base_row, 
                                    compare_rows, 
                                    numeric_var_names, 
                                    ordinal_var_names,
                                    categorical_var_names) {
  
  m <- length(compare_rows)
  
  # Loop over each of the numeric variables that are under consideration.
  # Find the Scaling and Imputation Factors, and then compute the difference
  # Between @base_row and @compare_rows
  numeric_comparisons <- vector(mode = "list", length(numeric_var_names))
  names(numeric_comparisons) <- numeric_var_names
  
  for(numvar in numeric_var_names) {
    varvals <- as.numeric(data |> pull(numvar)) # Ensure it is read as numeric
    
    # Get Scaling and Imputation Factors
    max_v <- max(varvals, na.rm = TRUE)
    min_v <- min(varvals, na.rm = TRUE)
    imp_val <- (max_v + min_v) / 2
    
    varvals <- replace_na(varvals, imp_val)
    
    # Calculate Difference then Scale.
    numeric_comparisons[[numvar]] <- abs(varvals[base_row] - varvals[compare_rows])/(max_v - min_v)
  }
  
  
  # Find the Scaling Factors for the Ordinal Variables
  factor_level_counts <- unlist(lapply(data[,ordinal_var_names], FUN = function(x) { nlevels(as.factor(x)) } ))
  factor_counts <- matrix(factor_level_counts,
                          nrow = m,
                          ncol = length(factor_level_counts), 
                          byrow = TRUE)

  # Compute the Ordinal Differences and then scale by the matrix that was 
  # computed above. First form a matrix of the base_row repeated often enough, 
  # and then consider the subtraction of the compare_rows, scaling as needed.
  r1_ordinal <- matrix(replace_na(as.numeric(data[base_row,ordinal_var_names]), 0), 
                       ncol = length(ordinal_var_names),
                       nrow = m,
                       byrow = TRUE)
  
  r2_ordinal <- data.matrix(data[compare_rows, ordinal_var_names])
  r2_ordinal[is.na(r2_ordinal)] <- 0
  
  # Ordinal Comparisons
  ordinal_compare <- abs(r1_ordinal - r2_ordinal)/factor_counts

  # Categorical Variable Comparisons.
  # These are compared at equality versus inequality. 
  r1_categorical <- replace_na(as.numeric(data[base_row, categorical_var_names]), -1)
  r1_categorical <- matrix(r1_categorical, 
                           nrow = m, 
                           ncol = length(r1_categorical), 
                           byrow = TRUE)
  
  r2_categorical <- data.matrix(data[compare_rows, categorical_var_names])
  r2_categorical[is.na(r2_categorical)] <- -1
  
  categorical <- r1_categorical != r2_categorical

  # Return the Row Sums
  as.numeric(
    rowSums(cbind(do.call(cbind, numeric_comparisons), 
                  ordinal_compare, 
                  categorical))
  )
}

# Computes a (dis)similarity matrix from data
# Takes as input: 
#     @data (the data frame to build from)
#     @numeric_var_names (the names of continuous vars for comparison)
#     @ordinal_var_names (the names of ordinal vars for comparison)
#     @categorical_var_names (the names of categorical vars for comparison)
#     @verbose(default:TRUE) (whether to echo the progress or not)
# 
# Returns
#     A matrix which is n x n with each value is how dissimilar (i, j) are. 
#     This will be a symmetric matrix, with zeroes on the diagonal (for perfect
#     similarity)

compute_dissimilarity_matrix <- function(data,
                                         numeric_var_names,
                                         ordinal_var_names,
                                         categorical_var_names,
                                         verbose=TRUE) {
  n <- nrow(data)
  
  # Compute Similarity Matrix
  sim_matrix <- matrix(0, nrow = n, ncol = n)
  for( ii in 1:(n - 1) ) {
    if(verbose) { cat(paste0("\r Working on ", ii, "/", n-1)) }
    
    row_diffs <- row_dissimilarity_score(data = data,
                                         base_row = ii, 
                                         compare_rows = (ii+1):n, 
                                         numeric_var_names = numeric_var_names,
                                         ordinal_var_names = ordinal_var_names,
                                         categorical_var_names = categorical_var_names)

    sim_matrix[ii,(ii+1):n] <- row_diffs
    sim_matrix[(ii+1):n, ii] <- row_diffs
  }
  
  # Return the Matrix
  sim_matrix
}

# Compute a probability matrix from the data
# Takes as input:
#     @dissimilarity_matrix (the matrix which gives dissimilarity scores)
#     @transformation (a function which converts a dissimilarity to similarity)
#
# Returns
#   A matrix with the same dimensions as the dissimilarity matrix, with 
#   each entry being a normalized probability (such that rows sum to 1)
compute_probability_matrix <- function(dissimilarity_matrix,
                                      transformation) {
  
  sim_matrix <- transformation(dissimilarity_matrix)
  
  prob_matrix <- sim_matrix/(rowSums(sim_matrix)-1)
  diag(prob_matrix) <- 0
  
  prob_matrix
}


# Select the usual intake and recorded options for a particular row of data
# Takes as input
#     @dietary_data (a data frame with SEQN + identifiers, then all dietary 
#                    variables; note variables from day1 should end _D1 and 
#                     variables for day2 should end _D2).
#     @row_idx (the index for the row you want to select from)
#     @probability_matrix (the matrix of probabilities for re-sampling)
#     @K(default=1000) (the total number of individuals whose data should be 
#                       selected to compute usual intake)
#     @L(default=2) (the number of days of recalls to keep)
# 
# Returns
#   A data frame with 1 row which has SEQN followed by L measurements of each
#   variable in the dietary data, as well as a "usual intake" for each of these
#   measurements. 
row_resample_diet <- function(dietary_data,
                              row_idx,
                              probability_matrix,
                              K = 1000,
                              L = 2) {

  # First we sample which items will actually be borrowed
  loaner_idx <- sample(1:nrow(probability_matrix),
                       K,
                       prob = probability_matrix[row_idx,],
                       replace = TRUE)
  
  # Form the usual intake matrix
  usual_intake_mat <- dietary_data[loaner_idx, ]
  
  # Distort the values up or down by 10%
  d1_cols_filter <- which(endsWith(colnames(dietary_data), "_D1"))
  d2_cols_filter <- which(endsWith(colnames(dietary_data), "_D2"))
  diet_cols <- colnames(dietary_data)[c(d1_cols_filter, d2_cols_filter)]
  
  distortion <- matrix(runif(K*length(diet_cols), 0.9, 1.1),
                       nrow = K,
                       ncol = length(diet_cols))
  
  usual_intake_mat[,diet_cols] <- usual_intake_mat[,diet_cols] * distortion
  
  # Compute Usual Intakes by Combining Day 1 + Day 2 Matrices
  ui_d1 <- usual_intake_mat[,d1_cols_filter]
  colnames(ui_d1) <- str_replace_all(colnames(ui_d1), "_D1", "")
  
  ui_d2 <- usual_intake_mat[,d2_cols_filter]
  colnames(ui_d2) <- str_replace_all(colnames(ui_d2), "_D2", "")
  
  ui_only <- colMeans(rbind(ui_d1, ui_d2))
  names(ui_only) <- paste0(names(ui_only), "_usual_intake")
  
  # Construct Row to Return
  return_row <- dietary_data[row_idx, ]
  return_row[,diet_cols] <- usual_intake_mat[sample(1:nrow(usual_intake_mat),
                                                    1),
                                             diet_cols]
  cbind(return_row, t(ui_only))
}

# Calculates the dietary matrix with usual intake and days measurements 
# Takes as input
#     @dietary_data (a data frame with SEQN + identifiers, then all dietary 
#                    variables; note variables from day1 should end _D1 and 
#                     variables for day2 should end _D2).
#     @n(default=nrow(dietary_data)) (The number of individuals to create the
#                     dataset for. By default this will be resampled based on
#                     the complete nhanes. If n < nrow(dietary_data) sampling
#                     is done without replacement; otherwise, sampling is done
#                     with replacement.)
#     @probability_matrix (the matrix of probabilities for re-sampling)
#     @K(default=1000) (the total number of individuals whose data should be 
#                       selected to compute usual intake)
#     @L(default=2) (the number of days of recalls to keep)
# 
# Returns
#   A data frame with n rows which has SEQN followed by L measurements of each
#   variable in the dietary data, as well as a "usual intake" for each of these
#   measurements. 
#
# Important: The row_idx should math both the probability matrix and the dietary
#            data. As a quick check on this, the function enforces the same number
#            of rows between them.
#

build_resampled_diet_matrix <- function(dietary_data,
                                        probability_matrix,
                                        n = nrow(dietary_data),
                                        K = 1000,
                                        L = 2) {
  
  if(nrow(dietary_data) != nrow(probability_matrix)) {
    stop(paste0("Please ensure that the probability matrix corresponds 1:1", 
                " with the dietary data."))
  }
  
  if(n == nrow(dietary_data)) {
    sample_idx <- 1:n
  } else if (n < nrow(dietary_data)) {
    sample_idx <- sample(1:nrow(dietary_data), n)
  } else {
    sample_idx <- sample(1:nrow(dietary_data), n, replace = TRUE)
  }

  
  do.call(rbind,
          apply(matrix(sample_idx), MARGIN = 1, FUN = function(idx){
              row_resample_diet(dietary_data,
                                idx,
                                probability_matrix,
                                K = 1000,
                                L = 2)
            }))
  
}
