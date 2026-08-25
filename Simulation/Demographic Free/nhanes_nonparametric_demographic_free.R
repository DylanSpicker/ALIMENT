library(tidyverse)


# -----------------------------------------------------------------------------
# Demographic Free Non-Parametric Simulation
# 
# Description:
# This script creates a synthetic dataset resembling NHANES dietary data without 
# imposing parametric modelling assumptions. It samples individuals from the 
# NHANES dataset and then resamples their 24-hour recalls to generate a "usual 
# intake" along with two days of recorded intake, introducing some random distortion.
# 
# Requirements:
# - NHANES and FPED data files placed in respective directories:
#   - NHANES/DR1TOT_J.xpt, NHANES/DR2TOT_J.xpt
#   - FPED/fped_dr1tot_1718.sas7bdat, FPED/fped_dr2tot_1718.sas7bdat
# - Required packages: tidyverse, haven, labelled
# 
# Output:
# - simulated_nhanes_complete_no_demographic.csv: Synthetic dataset
# - TrueDistCSV/: Directory containing true usual intake population subsets
# - TrueDISTCSV/NHANES_Simulation_True_Population.csv: Combined true population data
# -----------------------------------------------------------------------------


################################################################################

set.seed(31415)

diet1_df <- haven::read_xpt("NHANES/DR1TOT_J.xpt", NULL) |> labelled::remove_labels() |> 
              select(SEQN, "DR1T_KCAL"=DR1TKCAL, "DR1T_PFAT"=DR1TPFAT, "DR1T_MFAT"=DR1TMFAT, "DR1T_SFAT"=DR1TSFAT, "DR1T_SODI"=DR1TSODI)
diet2_df <- haven::read_xpt("NHANES/DR2TOT_J.xpt", NULL) |> labelled::remove_labels() |> 
              select(SEQN, "DR2T_KCAL"=DR2TKCAL, "DR2T_PFAT"=DR2TPFAT, "DR2T_MFAT"=DR2TMFAT, "DR2T_SFAT"=DR2TSFAT, "DR2T_SODI"=DR2TSODI)

FPED_D1 <- haven::read_sas("FPED/fped_dr1tot_1718.sas7bdat", NULL) %>% 
            filter(DR1DRSTZ == 1,     
                   RIDAGEYR >= 2) %>% 
            select(SEQN, 
                   WTDRD1, 
                   WTDR2D, 
                   starts_with("DR1T"),
                   -DR1TNUMF) |>
              left_join(diet1_df, by = "SEQN")

FPED_D2 <- haven::read_sas("FPED/fped_dr2tot_1718.sas7bdat", NULL) %>% 
            filter(DR2DRSTZ == 1,     
                   RIDAGEYR >= 2) %>% 
            select(SEQN, 
                   starts_with("DR2T"),
                   -DR2TNUMF) |>
              left_join(diet2_df, by = "SEQN")

colnames(FPED_D1) <- c("SEQN", "WTDRD1", "WTDR2D", paste0(str_replace_all(colnames(FPED_D1)[4:ncol(FPED_D1)], "DR1T_", ""), "_D1"))
colnames(FPED_D2) <- c("SEQN", paste0(str_replace_all(colnames(FPED_D2)[2:ncol(FPED_D1)], "DR2T_", ""), "_D2"))

complete_df <- FPED_D1 %>% inner_join(FPED_D2, by = "SEQN")
difference_dist_mat <- complete_df %>% select(ends_with("_D1")) %>% as.matrix() - complete_df %>% select(ends_with("_D2")) %>% as.matrix()

d1_col_names <- colnames(complete_df %>% select(ends_with("_D1")))
d2_col_names <- colnames(complete_df %>% select(ends_with("_D2")))
usual_intake_names <- str_replace_all(d2_col_names, "_D2", "_truth")

n <- nrow(complete_df)
K <- 1000
L <- 2

# 1. Resample rows from complete_df_D1
complete_df_D1 <- complete_df %>% select(ends_with("_D1"))

simulated_df <- complete_df_D1[sample(1:nrow(complete_df), n, TRUE), ]

simulated_df[,c(d2_col_names, usual_intake_names)] <- 0

# 2. Loop over each observation; generate usual intake and 2 days

for(idx in 1:nrow(simulated_df)) {
    individual_obs <- simulated_df[rep(idx, K), d1_col_names] + difference_dist_mat[sample(1:nrow(difference_dist_mat), K, TRUE), ]
    distortion <- matrix(runif(K*ncol(individual_obs), 0.9, 1.1), nrow = K)
    individual_obs <- individual_obs * distortion

    individual_obs[individual_obs < 0] <- 0

    usual_intake <- colMeans(rbind(simulated_df[idx, d1_col_names], individual_obs))


    day2 <- individual_obs[sample(1:K, 1), ]

    simulated_df[idx, c(d2_col_names, usual_intake_names)] <- cbind(day2, t(usual_intake))
}

simulated_df %>% write.csv(file = "simulated_nhanes_complete_no_demographic.csv", row.names = FALSE)

set.seed(31415)

n <- 10000
K <- 1000

POP <- 10000000

for(idx in 1:(POP/n)) {
    cat(paste0("\r Progress: ", round(((idx-1) * n) / POP * 100,3), "%"))

    # 1. Resample rows from complete_df_D1
    full_truth <- complete_df_D1[sample(1:nrow(complete_df), n, TRUE), ]
    full_truth$rowID <- 1:n

    expanded_truth <- full_truth[rep(1:nrow(full_truth), each = K), ]
    additive_error <- difference_dist_mat[sample(1:nrow(difference_dist_mat), K * n, TRUE), ]
    distortion <- matrix(runif(K*n*ncol(additive_error), 0.9, 1.1), nrow = K * n)
    expanded_truth[,d1_col_names] <- (expanded_truth[,d1_col_names] + additive_error) * distortion
    expanded_truth[expanded_truth < 0] <- 0

    usual_intake_df <- rbind(full_truth, expanded_truth) %>%
                            group_by(rowID) %>%
                            summarize(across(everything(), mean)) %>%
                            select(-rowID)

    colnames(usual_intake_df) <- usual_intake_names

    usual_intake_df %>% write.csv(file = paste0("TrueDistCSV/true_pop_group_", idx, ".csv"), row.names = FALSE)
}

library(dplyr)
library(readr)

list.files(path="TrueDistCSV", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows() %>% 
  write.csv(file = "TrueDISTCSV/NHANES_Simulation_True_Population.csv")
