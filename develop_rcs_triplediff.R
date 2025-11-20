# Script to develop and test Repeated Cross-Section (RCS) Triple Diff
# This script mimics the workflow for RCS data using the newly created att_dr_rc function.

rm(list=ls())
devtools::load_all() # Load the package

set.seed(123)

# 1. Generate Data (Balanced Panel initially)
# We use the 2-period DGP
# dgp_type=1: Homogeneous effects (ATT=0)
dgp <- gen_dgp_2periods(size = 2000, dgp_type = 1)
dt_panel <- dgp$data

# 2. Convert to Repeated Cross-Section
# To simulate RCS, we will sample different individuals for each time period.
# Time periods are usually labeled, let's assume min(time) is pre and max(time) is post.
times <- sort(unique(dt_panel$time))
t1 <- times[1]
t2 <- times[2]

# Split IDs into two disjoint sets
all_ids <- unique(dt_panel$id)
ids_t1 <- sample(all_ids, length(all_ids)/2)
ids_t2 <- setdiff(all_ids, ids_t1)

# Create RCS dataset
dt_rc <- rbind(
  dt_panel[time == t1 & id %in% ids_t1],
  dt_panel[time == t2 & id %in% ids_t2]
)

# Verify it's RCS (no ID should appear twice)
if(any(table(dt_rc$id) > 1)) stop("Data is not Repeated Cross-Section!")

message("RCS Data generated. N = ", nrow(dt_rc))

# 3. Manual Preprocessing
# Since we haven't updated the main ddd() function or preprocess.R to handle RC yet,
# we manually construct the object required by att_dr_rc.

# Define variables
yname <- "y"
tname <- "time"
idname <- "id"
treatname <- "state" # In this DGP, state indicates treatment group assignment
partitionname <- "partition"
xformula <- ~ x1 + x2 + x3 + x4

# Create 'post' variable
dt_rc[, post := as.numeric(time == t2)]

# Create 'treat' variable (already exists as state, but let's standardize)
dt_rc[, treat := get(treatname)]

# Create 'weights'
dt_rc[, weights := 1]

# Create 'subgroup'
# 4: Treated (treat=1) & Eligible (partition=1)
# 3: Treated (treat=1) & Ineligible (partition=0)
# 2: Control (treat=0) & Eligible (partition=1)
# 1: Control (treat=0) & Ineligible (partition=0)
dt_rc[, subgroup := fifelse(partition == 1 & treat == 1, 4,
                    fifelse(partition == 0 & treat == 1, 3,
                    fifelse(partition == 1 & treat == 0, 2, 1)))]

# Check subgroup counts
subgroup_counts <- dt_rc[, .N, by = subgroup][order(-subgroup)]
print(subgroup_counts)

# Prepare Covariates
# We need to cbind the model matrix to the data.table (excluding intercept as per convention in this package?)
# att_dr_rc uses xformula to compute model matrices internally in compute_pscore/reg.
# However, att_dr seems to expect `did_preprocessed$preprocessed_data` to contain the columns?
# Let's check att_dr.R -> compute_pscore -> it uses `xformula` on `data`.
# So we don't strictly need to cbind columns if they are in the dataframe, 
# BUT preprocess.R does cbind them.
# Let's rely on columns x1, x2, x3, x4 being present in dt_rc.

# Construct the input list
did_preprocessed <- list(
  preprocessed_data = dt_rc,
  est_method = "dr", # Doubly Robust
  xformula = xformula,
  boot = FALSE,      # Analytical SEs (Influence Function)
  nboot = NULL,
  alpha = 0.05,
  use_parallel = FALSE,
  cores = 1,
  cband = FALSE,
  inffunc = TRUE,
  subgroup_counts = subgroup_counts
)

# 4. Run Estimation using the internal function
message("Running att_dr_rc...")
# We need to use triplediff::: because it's not exported yet
result <- triplediff:::att_dr_rc(did_preprocessed)

# 5. Display Results
print(result)

message("ATT Estimate: ", round(result$ATT, 4))
message("Standard Error: ", round(result$se, 4))
message("95% CI: [", round(result$lci, 4), ", ", round(result$uci, 4), "]")
message("True ATT is 0")

# 6. Compare with simple means (for sanity check)
# Simple DDD = (dT_post - dT_pre) - (dC_post - dC_pre)
# where d = (y_elig - y_inelig)
means <- dt_rc[, .(mu = mean(y)), by = .(treat, partition, post)]

get_mu <- function(tr, pa, po) means[treat==tr & partition==pa & post==po, mu]

# Treated Eligible
y_111 <- get_mu(1, 1, 1)
y_110 <- get_mu(1, 1, 0)
# Treated Ineligible
y_101 <- get_mu(1, 0, 1)
y_100 <- get_mu(1, 0, 0)
# Control Eligible
y_011 <- get_mu(0, 1, 1)
y_010 <- get_mu(0, 1, 0)
# Control Ineligible
y_001 <- get_mu(0, 0, 1)
y_000 <- get_mu(0, 0, 0)

simple_ddd <- ((y_111 - y_110) - (y_101 - y_100)) - ((y_011 - y_010) - (y_001 - y_000))
message("Simple DDD (Naive): ", round(simple_ddd, 4))

# ---------------------------------------------------------
# 7. VALIDATION AGAINST DRDID
# ---------------------------------------------------------
message("\n--- VALIDATION AGAINST DRDID ---")

if (!requireNamespace("DRDID", quietly = TRUE)) {
  message("DRDID package not installed. Skipping validation.")
} else {
  library(DRDID)
  
  # Extract data for Subgroups 4 (Treated) and 3 (Control)
  # Subgroup 4: Treated=1, Partition=1
  # Subgroup 3: Treated=1, Partition=0 (Acts as control in this comparison)
  
  # In DRDID context:
  # "D" indicates treatment group. 
  # For comparison 4 vs 3: 
  #   - Unit is "Treated" if subgroup == 4
  #   - Unit is "Control" if subgroup == 3
  
  # Filter data
  data_check <- dt_rc[subgroup %in% c(3, 4)]
  
  # Setup variables for DRDID
  # y: outcome
  # post: post indicator
  # D: 1 if subgroup 4, 0 if subgroup 3
  y_val <- data_check$y
  post_val <- data_check$post
  D_val <- ifelse(data_check$subgroup == 4, 1, 0)
  
  # Covariates (including intercept is handled by DRDID usually, but let's pass matrix)
  # DRDID usually wants a matrix of covariates
  # We need to define the matrix explicitly using xformula and the data_check
  cov_val <- model.matrix(xformula, data = data_check)
  
  # Run DRDID
  message("Running DRDID::drdid_rc...")
  # Note: if cov_val has intercept, we should check if drdid_rc adds another one.
  # The function documentation says "Please add a vector of constants if you want to include an intercept".
  # So model.matrix output is perfect.
  res_drdid <- DRDID::drdid_rc(y = y_val, post = post_val, D = D_val, covariates = cov_val)
  
  # Run internal compute_did_rc
  # We need to prep the arguments exactly as att_dr_rc does
  message("Running triplediff:::compute_did_rc...")
  
  # Pre-compute nuisances using internal functions
  # 1. Propensity Score
  ps_list <- list()
  # In att_dr_rc logic:
  # If comparison is subgroup 3 (control), we put it in pscores[[1]]?
  # compute_did_rc(..., condition_subgroup = 3, ...)
  # inside compute_did_rc for condition_subgroup=3:
  # idx <- 1 -> pscores[[1]]
  
  # compute_pscore_rc expects 'data' with 'subgroup' column
  # It computes pscore for D=1 (subgroup 4) vs D=0 (condition_subgroup)
  ps_res <- triplediff:::compute_pscore_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
  # We need to wrap it in a list structure that compute_did_rc expects
  # pscores arg is a list of results. compute_did_rc selects pscores[[idx]]
  # So we create a list where the first element is our result
  pscores_arg <- list(ps_res, NULL, NULL) 
  
  # 2. Outcome Regression
  mu_res <- triplediff:::compute_outcome_regression_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
  reg_adjust_arg <- list(mu_res, NULL, NULL)
  
  # Run compute_did_rc
  # Note: data argument for compute_did_rc is the full data? 
  # Yes, it subsets internally: condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  res_triplediff <- triplediff:::compute_did_rc(
    data = dt_rc, 
    condition_subgroup = 3, 
    pscores = pscores_arg, 
    reg_adjustment = reg_adjust_arg, 
    xformula = xformula, 
    est_method = "dr"
  )
  
  message("\n--- COMPARISON ---")
  message(sprintf("DRDID ATT:        %.10f", res_drdid$ATT))
  message(sprintf("triplediff ATT:   %.10f", res_triplediff$dr_att))
  message(sprintf("Difference:       %.10f", res_drdid$ATT - res_triplediff$dr_att))
  
  if(abs(res_drdid$ATT - res_triplediff$dr_att) < 1e-8) {
    message("SUCCESS: Results match!")
  } else {
    message("WARNING: Results differ.")
  }
}
