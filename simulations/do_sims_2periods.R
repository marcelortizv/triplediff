# ------------------------------------------------------------------------
# Purpose: Simulate finite sample performance of estimators based on DGPs
# Date: 06/05/2024
# ------------------------------------------------------------------------


library(foreach)
library(parallel)
library(doParallel)
library(dplyr)
library(data.table)
library(DRDID)
PAT <- Sys.getenv("GITHUB_PAT")
install_github("marcelortizv/triplediff", auth_token = PAT, ref = "monte-carlo-simulations")
library(triplediff)

# Clear memory
rm(list=ls())

# Setting parameters
nrep <- 1000    # Monte Carlo replications

# Set the Working Directory
# address <- "set/here/your/working/directory/for/panel"
address <- getwd()
setwd(address)

# Load the functions for simulations
source(paste0(address, "/simulations/", "dgps_sims.R", sep=""))


do_one_sim <- function(seed, n, dgp_type){
  set.seed(seed)
  # get efficiency bound
  dgp <- gen_dgp(n, dgp_type)
  data <- dgp$data
  att.true <- dgp$att
  att.unf <- dgp$att.unf
  eff <- dgp$eff

  # -----------------------------------------
  # DRDDD Estimator
  # -----------------------------------------
  drddd <- triplediff::ddd(yname = "y", tname = "time", idname = "id", dname = "state",
               gname = NULL, partition_name = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
               data = data, est_method = "trad", skip_data_checks = TRUE)

  # Whether the CI covers the true ATT (coverage probability)
  cp_drddd <- as.numeric((drddd$lci <= att.true) * (drddd$uci >= att.true))
  #Length of confidence interval
  len_drddd <- drddd$uci - drddd$lci

  # -----------------------------------------
  # TWFE Estimator
  # -----------------------------------------
  twfe <- triplediff::twfe_ddd(yname = "y", tname = "time", dname = "state",
                   partition_name = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
                   data = data)
  # Whether the CI covers the true ATT (coverage probability)
  cp_twfe <- as.numeric((twfe$lci <= att.true) * (twfe$uci >= att.true))
  #Length of confidence interval
  len_twfe <- twfe$uci - twfe$lci

  # -----------------------------------------
  # Diff of 2 DRDID
  # -----------------------------------------
  #subset data based on partition
  drdid_partition1 <- data[data$partition == 1,]
  drdid_partition0 <- data[data$partition == 0,]

  # run drdid among obs with partition value equal to 1
  drdid_2 <- DRDID::drdid(yname = "y", tname="time", idname="id", dname="state",
                          xformla = ~cov1 + cov2 + cov3 + cov4,
                          data = drdid_partition1, panel = TRUE)
  # run drdid among obs with partition value equal to 0
  drdid_1 <-  DRDID::drdid(yname = "y", tname="time", idname="id", dname="state",
                           xformla = ~cov1 + cov2 + cov3 + cov4,
                           data = drdid_partition0, panel = TRUE)
  drdid <- drdid_2$ATT - drdid_1$ATT


  summary <- cbind(
            # true ATT
            att.true = att.true,
            # unfeasible ATT
            att.unf = att.unf,

            # DRDDD
            drddd_att = drddd$ATT,
            var_drddd = (drddd$se * sqrt(n))^2,
            cp_drddd = cp_drddd,
            len_drddd = len_drddd,

            # TWFE
            twfe_att = twfe$ATT,
            var_twfe = (twfe$se * sqrt(n))^2,
            cp_twfe = cp_twfe,
            len_twfe = len_twfe,

            # Diff of 2 DRDID
            drdid = drdid,

            # Semiparametric Efficiency Bound
            eff_bound = eff)

  return(summary)
}

# Function to calculate the performance metrics
calculate_metrics <- function(estimator, true_att, var_estimator, cp, len) {
  avg_estimator <- mean(estimator, na.rm = TRUE)
  bias <- avg_estimator - mean(true_att, na.rm = TRUE)
  rmse <- sqrt(mean((estimator - true_att)^2, na.rm = TRUE))
  mean_var <- mean(var_estimator, na.rm = TRUE)
  mean_cp <- mean(cp, na.rm = TRUE)
  mean_len <- mean(len, na.rm = TRUE)
  c(avg_estimator, bias, rmse, mean_var, mean_cp, mean_len)
}


# make a parameter grid
param_grid <- expand.grid(
  seed = 1:nrep,
  #set sample size
  n = c(1000, 5000, 10000),
  dgp_type = c(1, 2, 3, 4)
)


# DOING SIMULATION
# Detect the number of cores
numCores <- detectCores()
# Register doParallel as the backend for foreach
cl <- makeCluster(numCores - 1) # leave one core free for other tasks
registerDoParallel(cl)

results_sims <- foreach(i = 1:nrow(param_grid), .combine = rbind, .init = data.frame(),
                        .packages = c('dplyr', 'data.table', 'DRDID', 'triplediff')) %dopar% {
  row <- param_grid[i, ]
  est <- do_one_sim(row$seed + 711*row$dgp_type, row$n, row$dgp_type)
  # Create a data frame to save results for each iteration
  data.frame(
    seed = row$seed + 711*row$dgp_type,
    n = row$n,
    dgp_type = row$dgp_type,
    att_true = est[[1]],
    att_unf = est[[2]],
    drddd_att = est[[3]],
    var_drddd = est[[4]],
    cp_drddd = est[[5]],
    len_drddd = est[[6]],
    twfe_att = est[[7]],
    var_twfe = est[[8]],
    cp_twfe = est[[9]],
    len_twfe = est[[10]],
    drdid = est[[11]],
    eff_bound = est[[12]]
  )
}

# Stop the parallel cluster
stopCluster(cl)


# ------------------------------------------------------------------------
# Creating summary table
# ------------------------------------------------------------------------

# Initialize a list to store the summary tables
summary_tables <- list()

# Iterate through each dgp_type
for (dgp in unique(results_sims$dgp_type)) {
  # Filter data for the current dgp_type
  dgp_data <- results_sims[results_sims$dgp_type == dgp,]

  # Initialize a data frame to store summary for current dgp
  dgp_summary <- data.frame()

  # Iterate through each sample size
  for (n in unique(dgp_data$n)) {
    # Filter data for the current sample size
    n_data <- dgp_data[dgp_data$n == n, ]

    # Compute the metrics for each estimator
    drddd_metrics <- calculate_metrics(n_data$drddd_att, n_data$att_true, n_data$var_drddd, n_data$cp_drddd, n_data$len_drddd)
    twfe_metrics <- calculate_metrics(n_data$twfe_att, n_data$att_true, n_data$var_twfe, n_data$cp_twfe, n_data$len_twfe)
    # Note: DRDID does not have a variance column, using twfe var as a placeholder
    drdid_metrics <- calculate_metrics(n_data$drdid, n_data$att_true, n_data$var_twfe, n_data$cp_twfe, n_data$len_twfe)
    # Combine the metrics into a summary table
    n_summary <- data.frame(
      estimator = c("DRDDD", "TWFE", "DRDID"),
      sample_size = n,
      avg_estimator = c(drddd_metrics[1], twfe_metrics[1], drdid_metrics[1]),
      bias = c(drddd_metrics[2], twfe_metrics[2], drdid_metrics[2]),
      rmse = c(drddd_metrics[3], twfe_metrics[3], drdid_metrics[3]),
      mean_variance = c(drddd_metrics[4], twfe_metrics[4], drdid_metrics[4]),
      mean_coverage_prob = c(drddd_metrics[5], twfe_metrics[5], drdid_metrics[5]),
      mean_length = c(drddd_metrics[6], twfe_metrics[6], drdid_metrics[6])
    )

    # Append the summary for the current sample size to the dgp summary
    dgp_summary <- rbind(dgp_summary, n_summary)
  }

  # Add the dgp summary to the list
  summary_tables[[paste0("DGP_Type_", dgp)]] <- dgp_summary
}

# If you want to save the summary tables as CSV files
for (dgp in names(summary_tables)) {
  write.csv(summary_tables[[dgp]], paste0(address, "/simulations/", dgp, "_summary.csv"), row.names = FALSE)
}





