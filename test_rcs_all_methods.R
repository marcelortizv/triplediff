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
# TEST 1: DR Method with intercept only (local DID with group 3 vs 4)
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
# TEST 2: IPW Method with intercept only (local DID with group 3 vs 4)
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
# TEST 3: REG Method with intercept only (local DID with group 3 vs 4)
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
# TEST 3.5: OLS Regression Triple Interaction Check (No Covariates)
#---------------------------------------------------------
message("\n=== TEST 3.5: OLS Triple Interaction (No Covariates) ===")
message("Testing that ATT = coefficient on Si*Qi*f2t in regression")

# Create regression variables
dt_rc[, Si := as.numeric(treat == 1)]  # Treatment indicator
dt_rc[, Qi := as.numeric(partition == 1)]  # Eligibility indicator
dt_rc[, f2t := as.numeric(post == 1)]  # Time period 2 indicator

# Run OLS regression: Y ~ 1 + Si + Qi + f2t + Si*Qi + Si*f2t + Qi*f2t + Si*Qi*f2t
ols_model <- lm(y ~ Si + Qi + f2t + Si:Qi + Si:f2t + Qi:f2t + Si:Qi:f2t, data = dt_rc)
ols_coef <- coef(ols_model)["Si:Qi:f2t"]
ols_se <- summary(ols_model)$coefficients["Si:Qi:f2t", "Std. Error"]

message(sprintf("OLS Si*Qi*f2t coefficient: %.10f", ols_coef))
message(sprintf("OLS Si*Qi*f2t SE:          %.10f", ols_se))

# Compare with triplediff DR estimator (no covariates)
# Now test att_dr_rc() with proper preprocessing structure
# Create subgroup counts (total observations per subgroup)
subgroup_counts <- dt_rc[, .(count = .N), by = subgroup][order(-subgroup)]

# Create did_preprocessed structure
did_preprocessed <- list(
  preprocessed_data = dt_rc,
  est_method = "dr",
  xformula = xformula,
  boot = FALSE,
  nboot = NULL,
  alpha = 0.05,
  use_parallel = FALSE,
  cores = 1,
  cband = FALSE,
  inffunc = TRUE,
  subgroup_counts = subgroup_counts
)

# Call att_dr_rc()
res_att_dr_rc_nocov <- triplediff:::att_dr_rc(did_preprocessed)

message(sprintf("triplediff DDD ATT:        %.10f", res_att_dr_rc_nocov$ATT))
message(sprintf("Difference from OLS:       %.10f", abs(ols_coef - res_att_dr_rc_nocov$ATT)))

if(abs(ols_coef - res_att_dr_rc_nocov$ATT) < 1e-8) {
  message("✓ SUCCESS: triplediff ATT matches OLS triple interaction coefficient!")
} else {
  message("✗ FAILED: triplediff ATT differs from OLS.")
}

# Also compare with DRDID manual calculation (Diff of two DRDIDs)
data_partition1 <- dt_rc[partition == 1]
data_partition0 <- dt_rc[partition == 0]

y_nocov1 <- data_partition1$y
post_nocov1 <- data_partition1$post
D_nocov1 <- ifelse(data_partition1$subgroup == 4, 1, 0)
cov_nocov1 <- model.matrix(~ 1, data = data_partition1)
res_drdid_nocov1 <- DRDID::drdid_rc(y = y_nocov1, post = post_nocov1, D = D_nocov1,
                                   covariates = cov_nocov1, inffunc = TRUE)
y_nocov0 <- data_partition0$y
post_nocov0 <- data_partition0$post
D_nocov0 <- ifelse(data_partition0$subgroup == 3, 1, 0)
cov_nocov0 <- model.matrix(~ 1, data = data_partition0)
res_drdid_nocov0 <- DRDID::drdid_rc(y = y_nocov0, post = post_nocov0, D = D_nocov0,
                                   covariates = cov_nocov0, inffunc = TRUE)
# Compute DDD manually: ATT_3 = ATT_4vs3 - ATT_4vs2
manual_drdid_nocov <- res_drdid_nocov1$ATT - res_drdid_nocov0$ATT

message(sprintf("Diff of two DRDIDs:   %.10f", manual_drdid_nocov))

message(sprintf("Difference from OLS:       %.10f", abs(ols_coef - manual_drdid_nocov)))

if(abs(ols_coef - manual_drdid_nocov) < 1e-8) {
  message("✓ SUCCESS: Diff of 2 DRDIDs ATT matches OLS triple interaction coefficient!")
} else {
  message("✗ FAILED: Diff of 2 DRDIDs ATT differs from OLS.")
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

#---------------------------------------------------------
# TEST 7: Full att_dr_rc() function - DR Method
#---------------------------------------------------------
message("\n=== TEST 7: Full att_dr_rc() Function (DR Method) ===")

# Manually compute the DDD using DRDID for all three comparisons
xformula <- ~ cov1 + cov2 + cov3 + cov4

# Comparison 1: Subgroup 4 vs Subgroup 1 (Treated+Eligible vs Untreated+Ineligible)
data_comp1 <- dt_rc[subgroup %in% c(1, 4)]
y_comp1 <- data_comp1$y
post_comp1 <- data_comp1$post
D_comp1 <- ifelse(data_comp1$subgroup == 4, 1, 0)
cov_comp1 <- model.matrix(xformula, data = data_comp1)
res_comp1 <- DRDID::drdid_rc(y = y_comp1, post = post_comp1, D = D_comp1,
                             covariates = cov_comp1, inffunc = TRUE)

# Comparison 2: Subgroup 4 vs Subgroup 2 (Treated+Eligible vs Eligible+Untreated)
data_comp2 <- dt_rc[subgroup %in% c(2, 4)]
y_comp2 <- data_comp2$y
post_comp2 <- data_comp2$post
D_comp2 <- ifelse(data_comp2$subgroup == 4, 1, 0)
cov_comp2 <- model.matrix(xformula, data = data_comp2)
res_comp2 <- DRDID::drdid_rc(y = y_comp2, post = post_comp2, D = D_comp2,
                             covariates = cov_comp2, inffunc = TRUE)

# Comparison 3: Subgroup 4 vs Subgroup 3 (Treated+Eligible vs Treated+Ineligible)
data_comp3 <- dt_rc[subgroup %in% c(3, 4)]
y_comp3 <- data_comp3$y
post_comp3 <- data_comp3$post
D_comp3 <- ifelse(data_comp3$subgroup == 4, 1, 0)
cov_comp3 <- model.matrix(xformula, data = data_comp3)
res_comp3 <- DRDID::drdid_rc(y = y_comp3, post = post_comp3, D = D_comp3,
                             covariates = cov_comp3, inffunc = TRUE)

# Compute manual DDD: ATT_3 + ATT_2 - ATT_1
manual_ddd <- res_comp3$ATT + res_comp2$ATT - res_comp1$ATT

message(sprintf("Comparison 1 ATT (4 vs 1): %.10f", res_comp1$ATT))
message(sprintf("Comparison 2 ATT (4 vs 2): %.10f", res_comp2$ATT))
message(sprintf("Comparison 3 ATT (4 vs 3): %.10f", res_comp3$ATT))
message(sprintf("Manual DDD = ATT_3 + ATT_2 - ATT_1: %.10f", manual_ddd))

# Now test att_dr_rc() with proper preprocessing structure
# Create subgroup counts (total observations per subgroup)
subgroup_counts <- dt_rc[, .(count = .N), by = subgroup][order(-subgroup)]

# Create did_preprocessed structure
did_preprocessed <- list(
  preprocessed_data = dt_rc,
  est_method = "dr",
  xformula = xformula,
  boot = FALSE,
  nboot = NULL,
  alpha = 0.05,
  use_parallel = FALSE,
  cores = 1,
  cband = FALSE,
  inffunc = TRUE,
  subgroup_counts = subgroup_counts
)

# Call att_dr_rc()
res_att_dr_rc <- triplediff:::att_dr_rc(did_preprocessed)

message(sprintf("att_dr_rc() ATT:  %.10f", res_att_dr_rc$ATT))
message(sprintf("Difference:       %.10f", abs(manual_ddd - res_att_dr_rc$ATT)))
message(sprintf("DEBUG: inf_func class: %s, length: %d",
                class(res_att_dr_rc$inf_func), length(res_att_dr_rc$inf_func)))

if(abs(manual_ddd - res_att_dr_rc$ATT) < 1e-8) {
  message("✓ SUCCESS: att_dr_rc() matches manual DDD calculation!")
} else {
  message("✗ FAILED: att_dr_rc() differs from manual DDD.")
}

# Verify influence function properties
if (!is.null(res_att_dr_rc$inf_func)) {
  message(sprintf("Influence function length: %d (should equal nrow(dt_rc) = %d)",
                  length(res_att_dr_rc$inf_func), nrow(dt_rc)))
  message(sprintf("Mean of inf_func: %.10f (should be ≈ 0)", mean(res_att_dr_rc$inf_func)))
  message(sprintf("SE from influence function: %.10f", res_att_dr_rc$se))
  if(length(res_att_dr_rc$inf_func) == nrow(dt_rc)) {
    message("✓ Influence function has correct length!")
  } else {
    message("✗ WARNING: Influence function length mismatch!")
  }
} else {
  message("✗ WARNING: Influence function is NULL despite inffunc = TRUE")
}

#---------------------------------------------------------
# TEST 8: Full att_dr_rc() function - IPW Method
#---------------------------------------------------------
message("\n=== TEST 8: Full att_dr_rc() Function (IPW Method) ===")

# Compute manual DDD using IPW
res_comp1_ipw <- DRDID::std_ipw_did_rc(y = y_comp1, post = post_comp1, D = D_comp1,
                                       covariates = cov_comp1, inffunc = TRUE)
res_comp2_ipw <- DRDID::std_ipw_did_rc(y = y_comp2, post = post_comp2, D = D_comp2,
                                       covariates = cov_comp2, inffunc = TRUE)
res_comp3_ipw <- DRDID::std_ipw_did_rc(y = y_comp3, post = post_comp3, D = D_comp3,
                                       covariates = cov_comp3, inffunc = TRUE)

manual_ddd_ipw <- res_comp3_ipw$ATT + res_comp2_ipw$ATT - res_comp1_ipw$ATT

did_preprocessed$est_method <- "ipw"
res_att_dr_rc_ipw <- triplediff:::att_dr_rc(did_preprocessed)

message(sprintf("Manual DDD (IPW): %.10f", manual_ddd_ipw))
message(sprintf("att_dr_rc() ATT:  %.10f", res_att_dr_rc_ipw$ATT))
message(sprintf("Difference:       %.10f", abs(manual_ddd_ipw - res_att_dr_rc_ipw$ATT)))

if(abs(manual_ddd_ipw - res_att_dr_rc_ipw$ATT) < 1e-8) {
  message("✓ SUCCESS: att_dr_rc() IPW matches manual DDD calculation!")
} else {
  message("✗ FAILED: att_dr_rc() IPW differs from manual DDD.")
}

#---------------------------------------------------------
# TEST 9: Full att_dr_rc() function - REG Method
#---------------------------------------------------------
message("\n=== TEST 9: Full att_dr_rc() Function (REG Method) ===")

# Compute manual DDD using REG
res_comp1_reg <- DRDID::reg_did_rc(y = y_comp1, post = post_comp1, D = D_comp1,
                                   covariates = cov_comp1, inffunc = TRUE)
res_comp2_reg <- DRDID::reg_did_rc(y = y_comp2, post = post_comp2, D = D_comp2,
                                   covariates = cov_comp2, inffunc = TRUE)
res_comp3_reg <- DRDID::reg_did_rc(y = y_comp3, post = post_comp3, D = D_comp3,
                                   covariates = cov_comp3, inffunc = TRUE)

manual_ddd_reg <- res_comp3_reg$ATT + res_comp2_reg$ATT - res_comp1_reg$ATT

did_preprocessed$est_method <- "reg"
res_att_dr_rc_reg <- triplediff:::att_dr_rc(did_preprocessed)

message(sprintf("Manual DDD (REG): %.10f", manual_ddd_reg))
message(sprintf("att_dr_rc() ATT:  %.10f", res_att_dr_rc_reg$ATT))
message(sprintf("Difference:       %.10f", abs(manual_ddd_reg - res_att_dr_rc_reg$ATT)))

if(abs(manual_ddd_reg - res_att_dr_rc_reg$ATT) < 1e-8) {
  message("✓ SUCCESS: att_dr_rc() REG matches manual DDD calculation!")
} else {
  message("✗ FAILED: att_dr_rc() REG differs from manual DDD.")
}

message("\n=== ALL TESTS COMPLETE ===")
