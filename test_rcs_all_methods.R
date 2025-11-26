# Comprehensive Test Script for RCS Triple Diff
# Tests all three estimation methods (DR, IPW, REG) against DRDID
# Tests both with and without covariates

rm(list=ls())
devtools::load_all()

set.seed(123)

# Generate Data
dgp <- gen_dgp_2periods(size = 2000, dgp_type = 1)
dt_panel <- dgp$data

# Convert to Repeated Cross-Section
times <- sort(unique(dt_panel$time))
t1 <- times[1]
t2 <- times[2]

all_ids <- unique(dt_panel$id)
ids_t1 <- sample(all_ids, length(all_ids)/2)
ids_t2 <- setdiff(all_ids, ids_t1)

dt_rc <- rbind(
  dt_panel[time == t1 & id %in% ids_t1],
  dt_panel[time == t2 & id %in% ids_t2]
)

# Verify it's RCS
if(any(table(dt_rc$id) > 1)) stop("Data is not Repeated Cross-Section!")

message("RCS Data generated. N = ", nrow(dt_rc))

# Prepare data
dt_rc[, post := as.numeric(time == t2)]
dt_rc[, treat := state]
dt_rc[, weights := 1]
dt_rc[, subgroup := fifelse(partition == 1 & treat == 1, 4,
                    fifelse(partition == 0 & treat == 1, 3,
                    fifelse(partition == 1 & treat == 0, 2, 1)))]

# Extract data for validation (Subgroups 4 vs 3)
data_check <- dt_rc[subgroup %in% c(3, 4)]
y_val <- data_check$y
post_val <- data_check$post
D_val <- ifelse(data_check$subgroup == 4, 1, 0)

#---------------------------------------------------------
# TEST 1: DR Method with intercept only
#---------------------------------------------------------
message("\n=== TEST 1: DR Method (Intercept Only) ===")
xformula <- ~ 1
cov_val <- model.matrix(xformula, data = data_check)

# DRDID
res_drdid <- DRDID::drdid_rc(y = y_val, post = post_val, D = D_val, covariates = cov_val, inffunc = TRUE)

# triplediff
ps_res <- triplediff:::compute_pscore_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
mu_res <- triplediff:::compute_outcome_regression_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
res_triplediff <- triplediff:::compute_did_rc(
  data = dt_rc,
  condition_subgroup = 3,
  pscores = list(ps_res, NULL, NULL),
  reg_adjustment = list(mu_res, NULL, NULL),
  xformula = xformula,
  est_method = "dr"
)

message(sprintf("DRDID ATT:        %.10f", res_drdid$ATT))
message(sprintf("triplediff ATT:   %.10f", res_triplediff$dr_att))
message(sprintf("Difference:       %.10f", abs(res_drdid$ATT - res_triplediff$dr_att)))
if(abs(res_drdid$ATT - res_triplediff$dr_att) < 1e-8) {
  message("✓ SUCCESS: Results match!")
} else {
  message("✗ FAILED: Results differ.")
}

#---------------------------------------------------------
# TEST 2: IPW Method with intercept only
#---------------------------------------------------------
message("\n=== TEST 2: IPW Method (Intercept Only) ===")
xformula <- ~ 1
cov_val <- model.matrix(xformula, data = data_check)

# DRDID
res_drdid_ipw <- DRDID::std_ipw_did_rc(y = y_val, post = post_val, D = D_val, covariates = cov_val, inffunc = TRUE)

# triplediff
ps_res <- triplediff:::compute_pscore_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
mu_res <- triplediff:::compute_outcome_regression_null_rc(dt_rc, condition_subgroup = 3)
res_triplediff_ipw <- triplediff:::compute_did_rc(
  data = dt_rc,
  condition_subgroup = 3,
  pscores = list(ps_res, NULL, NULL),
  reg_adjustment = list(mu_res, NULL, NULL),
  xformula = xformula,
  est_method = "ipw"
)

message(sprintf("DRDID ATT:        %.10f", res_drdid_ipw$ATT))
message(sprintf("triplediff ATT:   %.10f", res_triplediff_ipw$dr_att))
message(sprintf("Difference:       %.10f", abs(res_drdid_ipw$ATT - res_triplediff_ipw$dr_att)))
if(abs(res_drdid_ipw$ATT - res_triplediff_ipw$dr_att) < 1e-8) {
  message("✓ SUCCESS: Results match!")
} else {
  message("✗ FAILED: Results differ.")
}

#---------------------------------------------------------
# TEST 3: REG Method with intercept only
#---------------------------------------------------------
message("\n=== TEST 3: REG Method (Intercept Only) ===")
xformula <- ~ 1
cov_val <- model.matrix(xformula, data = data_check)

# DRDID
res_drdid_reg <- DRDID::reg_did_rc(y = y_val, post = post_val, D = D_val, covariates = cov_val, inffunc = TRUE)

# triplediff
ps_res <- triplediff:::compute_pscore_null_rc(dt_rc, condition_subgroup = 3)
mu_res <- triplediff:::compute_outcome_regression_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
res_triplediff_reg <- triplediff:::compute_did_rc(
  data = dt_rc,
  condition_subgroup = 3,
  pscores = list(ps_res, NULL, NULL),
  reg_adjustment = list(mu_res, NULL, NULL),
  xformula = xformula,
  est_method = "reg"
)

message(sprintf("DRDID ATT:        %.10f", res_drdid_reg$ATT))
message(sprintf("triplediff ATT:   %.10f", res_triplediff_reg$dr_att))
message(sprintf("Difference:       %.10f", abs(res_drdid_reg$ATT - res_triplediff_reg$dr_att)))
if(abs(res_drdid_reg$ATT - res_triplediff_reg$dr_att) < 1e-8) {
  message("✓ SUCCESS: Results match!")
} else {
  message("✗ FAILED: Results differ.")
}

#---------------------------------------------------------
# TEST 4: DR Method with covariates
#---------------------------------------------------------
message("\n=== TEST 4: DR Method (With Covariates) ===")
xformula <- ~ cov1 + cov2 + cov3 + cov4
cov_val <- model.matrix(xformula, data = data_check)

# DRDID
res_drdid_cov <- DRDID::drdid_rc(y = y_val, post = post_val, D = D_val, covariates = cov_val, inffunc = TRUE)

# triplediff
ps_res <- triplediff:::compute_pscore_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
mu_res <- triplediff:::compute_outcome_regression_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
res_triplediff_cov <- triplediff:::compute_did_rc(
  data = dt_rc,
  condition_subgroup = 3,
  pscores = list(ps_res, NULL, NULL),
  reg_adjustment = list(mu_res, NULL, NULL),
  xformula = xformula,
  est_method = "dr"
)

message(sprintf("DRDID ATT:        %.10f", res_drdid_cov$ATT))
message(sprintf("triplediff ATT:   %.10f", res_triplediff_cov$dr_att))
message(sprintf("Difference:       %.10f", abs(res_drdid_cov$ATT - res_triplediff_cov$dr_att)))
if(abs(res_drdid_cov$ATT - res_triplediff_cov$dr_att) < 1e-8) {
  message("✓ SUCCESS: Results match!")
} else {
  message("✗ FAILED: Results differ.")
}

#---------------------------------------------------------
# TEST 5: IPW Method with covariates
#---------------------------------------------------------
message("\n=== TEST 5: IPW Method (With Covariates) ===")
xformula <- ~ cov1 + cov2 + cov3 + cov4
cov_val <- model.matrix(xformula, data = data_check)

# DRDID
res_drdid_ipw_cov <- DRDID::std_ipw_did_rc(y = y_val, post = post_val, D = D_val, covariates = cov_val, inffunc = TRUE)

# triplediff
ps_res <- triplediff:::compute_pscore_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
mu_res <- triplediff:::compute_outcome_regression_null_rc(dt_rc, condition_subgroup = 3)
res_triplediff_ipw_cov <- triplediff:::compute_did_rc(
  data = dt_rc,
  condition_subgroup = 3,
  pscores = list(ps_res, NULL, NULL),
  reg_adjustment = list(mu_res, NULL, NULL),
  xformula = xformula,
  est_method = "ipw"
)

message(sprintf("DRDID ATT:        %.10f", res_drdid_ipw_cov$ATT))
message(sprintf("triplediff ATT:   %.10f", res_triplediff_ipw_cov$dr_att))
message(sprintf("Difference:       %.10f", abs(res_drdid_ipw_cov$ATT - res_triplediff_ipw_cov$dr_att)))
if(abs(res_drdid_ipw_cov$ATT - res_triplediff_ipw_cov$dr_att) < 1e-8) {
  message("✓ SUCCESS: Results match!")
} else {
  message("✗ FAILED: Results differ.")
}

#---------------------------------------------------------
# TEST 6: REG Method with covariates
#---------------------------------------------------------
message("\n=== TEST 6: REG Method (With Covariates) ===")
xformula <- ~ cov1 + cov2 + cov3 + cov4
cov_val <- model.matrix(xformula, data = data_check)

# DRDID
res_drdid_reg_cov <- DRDID::reg_did_rc(y = y_val, post = post_val, D = D_val, covariates = cov_val, inffunc = TRUE)

# triplediff
ps_res <- triplediff:::compute_pscore_null_rc(dt_rc, condition_subgroup = 3)
mu_res <- triplediff:::compute_outcome_regression_rc(dt_rc, condition_subgroup = 3, xformula = xformula)
res_triplediff_reg_cov <- triplediff:::compute_did_rc(
  data = dt_rc,
  condition_subgroup = 3,
  pscores = list(ps_res, NULL, NULL),
  reg_adjustment = list(mu_res, NULL, NULL),
  xformula = xformula,
  est_method = "reg"
)

message(sprintf("DRDID ATT:        %.10f", res_drdid_reg_cov$ATT))
message(sprintf("triplediff ATT:   %.10f", res_triplediff_reg_cov$dr_att))
message(sprintf("Difference:       %.10f", abs(res_drdid_reg_cov$ATT - res_triplediff_reg_cov$dr_att)))
if(abs(res_drdid_reg_cov$ATT - res_triplediff_reg_cov$dr_att) < 1e-8) {
  message("✓ SUCCESS: Results match!")
} else {
  message("✗ FAILED: Results differ.")
}

message("\n=== ALL TESTS COMPLETE ===")
