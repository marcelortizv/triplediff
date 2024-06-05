# ------------------------------------------------------------------------
# Purpose: Simulate finite sample performance of estimators based on DGPs
# Date: 05/29/2024
# ------------------------------------------------------------------------


library(foreach)
library(parallel)
library(doParallel)
library(dplyr)
library(data.table)
library(DRDID)

# Clear memory
rm(list=ls())

# Setting parameters
nrep <- 10    # Monte Carlo replications

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
  drddd <- ddd(yname = "y", tname = "time", idname = "id", dname = "state",
               gname = NULL, partition_name = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
               data = data, est_method = "trad", skip_data_checks = TRUE)

  # Whether the CI covers the true ATT (coverage probability)
  cp_drddd <- as.numeric((drddd$lci <= att.true) * (drddd$uci >= att.true))
  #Length of confidence interval
  len_drddd <- drddd$uci - drddd$lci

  # -----------------------------------------
  # TWFE Estimator
  # -----------------------------------------
  twfe <- twfe_ddd(yname = "y", tname = "time", dname = "state",
                   partition_name = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
                   data = data)
  # Whether the CI covers the true ATT (coverage probability)
  cp_twfe <- as.numeric((twfe$lci <= att.true) * (twfe$uci >= att.true))
  #Length of confidence interval
  len_twfe <- twfe$uci - twfe$lci

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

            # Semiparametric Efficiency Bound
            eff_bound = eff)

  return(summary)
}


# make a parameter grid
param_grid <- expand.grid(
  seed = 1:nrep,
  #set sample size
  #n = c(500, 1000, 100000),
  n = 500,
  dgp_type = c(1, 2, 3, 4)
)


# DOING SIMULATION
# Detect the number of cores
numCores <- detectCores()
# Register doParallel as the backend for foreach
cl <- makeCluster(numCores - 1) # leave one core free for other tasks
registerDoParallel(cl)

results_sims <- foreach(i = 1:nrow(param_grid), .combine = rbind, .init = data.frame(),
                        .packages = c('dplyr', 'data.table', 'DRDID')) %dopar% {
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
    eff_bound = est[[11]]
  )
}

# Stop the parallel cluster
stopCluster(cl)


# ------------------------------------------------------------------------
# Creating summary table
# ------------------------------------------------------------------------




# # save data
# write.csv(results_efficiency, file = paste0(address, "/simulations/", "efficiency_bound_n-1m.csv"), row.names = FALSE)
#
# mean_efficiency <- results_efficiency %>%
#   group_by(dgp_type) %>%
#   summarise(mean_efficiency = mean(eff_bound, na.rm = TRUE))
#
# # Save the result to a CSV file
# write.csv(mean_efficiency, file = paste0(address, "/simulations/", "eff_bounds_by_dgp_n-1m.csv"), row.names = FALSE)


