library(tidyverse)

# -----------------------------------------------------------------------------
# Parametric Dietary Data Simulation
# 
# Description:
# This script simulates dietary data sets using specific parametric distributions 
# based on the NCI method settings. It provides two separate simulation scenarios:
# 1. Gamma Distribution Setting: Simulates data akin to Vitamin A intake using a 
#    gamma distribution with a log link.
# 2. Normal Distribution Setting: Simulates data akin to Calcium intake using a 
#    normal distribution, which is then transformed using an inverse Box-Cox transformation.
# 
# Requirements:
# - Required packages: tidyverse, stringr
# - ParametricCSV/ directory must exist prior to running
# 
# Output:
# - ParametricCSV/VitaminA_*.csv: Simulated Vitamin A datasets
# - ParametricCSV/Calcium_*.csv: Simulated Calcium datasets
# -----------------------------------------------------------------------------


# Gamma Setting
# A mixed-effects model approach for estimating the distribution of usual intake of nutrients: The NCI method
#   Section 5.2: Data were simulated from the vitamin A data for women in EATS, using a gamma distribution to simulate data for the model with a log link; the intercept was 8.7, random effect variance 0.37, and scale parameter 1.26. Five hundred datasets of size 500 were simulated. Truth was estimated by 365-d means for all 250 000 simulated individuals. The ISU method and the NCI method were ﬁt to 2 days of data.

set.seed(31415)

total_days <- 365
n <- 500 * 500

intercept <- 8.7
re_var <- 0.37
scale <- 1.26
re_s <- rnorm(n, 0, sqrt(re_var))

all_data <- as.data.frame(do.call(cbind, lapply(1:total_days, function(x){
    obs <- rgamma(n, shape = 1/scale, rate = 1/(scale * exp(intercept + re_s)))
})))

all_data$Truth <- rowMeans(all_data)

# Return a Single Set of 500
all_data[,"Truth"] %>% write.csv(file = "ParametricCSV/VitaminA_000_all.csv", row.names = FALSE)

for(ii in 1:500) {
    cat(paste0("\r Working on ", ii, " / ", 500))
    my_df <- all_data[((ii - 1)*500 + 1):(ii*500), c(paste0("V", sample(1:(ncol(all_data)-1), 2)))] 
    colnames(my_df) <- c("D1", "D2")
    my_df %>% write.csv(file = paste0("ParametricCSV/VitaminA_", stringr::str_pad(ii, 3, pad="0"), ".csv"), row.names = FALSE)
}

# Normal Setting
# Three hundred data sets were simulated with 500 individuals each. Data were simulated from a normal distribution basedon the distribution of transformed calcium in the EATS data, and then transformed using the inverse of the Box–Coxtransformation with a lambda of 0.3. Simple means were calculated at the individual level for 2 and 1000 days ofintake. Truth was obtained from the percentiles of the 1000-d means for all 300 data sets combined. The ISU methodand the NCI method were ﬁt to the two days of data, and the mean of the 5th, 10th, 25th, 50th, 75th, 90th, and 95thpercentiles were estimated and compared with truth and with 2-d means

set.seed(31415)

total_days <- 1000
n <- 300 * 500

intercept <- 20.2
re_sd <- 2.85
error_sd <- 0.9
lambda <- 0.3
re_s <- rnorm(n, 0, re_sd)  

all_data <- as.data.frame(do.call(cbind, lapply(1:total_days, function(x){
    obs <- (lambda * (intercept + re_s + rnorm(n, 0, error_sd)) + 1)^(1/lambda)
})))


all_data$Truth <- rowMeans(all_data)


all_data$Truth <- rowMeans(all_data)

# Return a Single Set of 500
all_data[,"Truth"] %>% write.csv(file = "ParametricCSV/Calcium_000_all.csv", row.names = FALSE)

for(ii in 1:500) {
    cat(paste0("\r Working on ", ii, " / ", 500))
    my_df <- all_data[((ii - 1)*300 + 1):(ii*300), c(paste0("V", sample(1:(ncol(all_data)-1), 2)))] 
    colnames(my_df) <- c("D1", "D2")
    my_df %>% write.csv(file = paste0("ParametricCSV/Calcium_", stringr::str_pad(ii, 3, pad="0"), ".csv"), row.names = FALSE)
}
