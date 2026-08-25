library(tidyverse)

source("./nhanes_simulation_helper_fn.R")

# -----------------------------------------------------------------------------
# Non-Parametric Simulation with Demographic Factors
# 
# Description:
# This script creates a synthetic NHANES-like dataset for dietary data without 
# specific parametric modelling assumptions. It incorporates demographic variables 
# to build a dissimilarity/probability matrix. Based on this matrix, it borrows 
# dietary values from similar individuals to construct a "usual intake" and 
# simulated recall days.
# 
# Requirements:
# - NHANES and FPED data files placed in respective directories:
#   - NHANES/DEMO_J.xpt, NHANES/DR1TOT_J.xpt, NHANES/DR2TOT_J.xpt
#   - FPED/fped_dr1tot_1718.sas7bdat, FPED/fped_dr2tot_1718.sas7bdat
# - Required packages: tidyverse, haven, labelled
# 
# Output:
# - probability_matrix.RDa: Cached probability matrix to speed up future runs
# - true_nhanes_complete.csv: The joined real NHANES/FPED dataset used as a base
# - simulated_nhanes_complete.csv: The synthetic output dataset
# -----------------------------------------------------------------------------

# 
# This file will simulate the creation of an NHANES-like data set for dietary
# data, without imposing specific parametric modelling assumptions on the nature
# of the variates or errors. To do so we first create a similarity matrix based
# on a set of variates between all individuals in the data. This matrix will 
# mark how close individuals are to one another in terms of any demographic, or
# health factors to one another. 
#
# Then, a synthetic dataset is created by first randomly sampling (with replace-
# ment) a set number (n) of individuals from the individuals within NHANES. Each 
# individual has their relevant demographic, health, and dietary factors 
# included. Then, based on the similarity matrix, a large number (K) of other 
# individuals' dietary data are selected. These intake values are multiplied by
# a randomly selected number between 0.9 and 1.1.
# 
# Then 2K + 2 total dietary values are averaged to give the "true" usual intake
# for each category, and a certain number of days (L) are selected from these 
# intakes to be the "recorded" intakes in the data file.
# 
# This process repeated across the full data produces the synthetic data file. 
#
################################################################################

set.seed(31415)
find_prob_matrix <- FALSE

# Convert Ages to all be measured in months for easier comparison
demo_df <- haven::read_xpt("NHANES/DEMO_J.xpt", NULL) |> labelled::remove_labels()

diet1_df <- haven::read_xpt("NHANES/DR1TOT_J.xpt", NULL) |> labelled::remove_labels() |> 
              select(SEQN, "DR1T_KCAL"=DR1TKCAL, "DR1T_PFAT"=DR1TPFAT, "DR1T_MFAT"=DR1TMFAT, "DR1T_SFAT"=DR1TSFAT, "DR1T_SODI"=DR1TSODI)
diet2_df <- haven::read_xpt("NHANES/DR2TOT_J.xpt", NULL) |> labelled::remove_labels() |> 
              select(SEQN, "DR2T_KCAL"=DR2TKCAL, "DR2T_PFAT"=DR2TPFAT, "DR2T_MFAT"=DR2TMFAT, "DR2T_SFAT"=DR2TSFAT, "DR2T_SODI"=DR2TSODI)

# Use FPED: note these files need to be downloaded from
# https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/
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

# Complete the Joining together of the Data Frames
complete_df <- demo_df %>% 
                inner_join(FPED_D1, by = "SEQN") %>% 
                inner_join(FPED_D2, by = "SEQN") %>%
                mutate(RIDAGEMN_MOD = 12*RIDAGEYR + ifelse(is.na(RIDAGEMN), 0, RIDAGEMN)) |> 
                select(-c("RIDAGEYR","RIDAGEMN", "RIDSTATR", "RIDEXAGM", 
                          "SDDSRVYR", "WTINT2YR", "WTMEC2YR", "SDMVPSU", 
                          "SDMVSTRA")) 

if(find_prob_matrix) {
       # Decide which columns count for which variable types 
       numeric_cols <- c( "RIDAGEMN_MOD", "INDFMPIR")
       ordinal_cols <- c("DMDYRSUS", "DMDEDUC3", "DMDEDUC2", "DMDHHSIZ", "DMDFMSIZ", 
                     "DMDHHSZA", "DMDHHSZB", "DMDHHSZE", "DMDHRAGZ", "DMDHSEDZ", 
                     "INDHHIN2", "INDFMIN2")

       categorical_cols <- intersect(setdiff(colnames(demo_df), 
                                          c("SEQN", numeric_cols, ordinal_cols)),
                                   colnames(complete_df))

       # Compute the Dissimilarity Matrix
       dissimilarity_matrix <- compute_dissimilarity_matrix(data = complete_df,
                                                        numeric_var_names = numeric_cols,
                                                        ordinal_var_names = ordinal_cols,
                                                        categorical_var_names = categorical_cols)

       probability_matrix <- compute_probability_matrix(dissimilarity_matrix,
                                                        transformation = function(x) { exp(-x) } )

       save(probability_matrix, file = "probability_matrix.RDa")

       remove(dissimilarity_matrix)
} else {
       load("probability_matrix.RDa")
}


# Create Resampled Matrix
# ~ 2.4 minutes to run on 7122 observations (2024-06-19)
resampled_matrix <- build_resampled_diet_matrix(complete_df,
                                                probability_matrix,
                                                n = nrow(complete_df),
                                                K = 1000,
                                                L = 2)

# Save as .csv
complete_df %>% write.csv(file = "true_nhanes_complete.csv", row.names = FALSE)
resampled_matrix %>% write.csv(file = "simulated_nhanes_complete.csv", row.names = FALSE)
