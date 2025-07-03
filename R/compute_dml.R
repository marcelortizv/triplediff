#' Utility functions to compute nuisances parameters based on DML.
#' @importFrom stats update
#' @import parglm
#' @import speedglm
#' @import data.table
#' @import mlr3
#' @import mlr3learners
#' @import mlr3tuning
#' @import DoubleML
#' @import caret
#' @noRd
NULL
# -------------------------------------------------
# FUNCTIONS FOR DML
# -------------------------------------------------

# Utility function to reshape the dataset from long to wide format
get_wide_data <- function(dt, xformula = ~ 1) {

  if (!inherits(dt, "data.table"))
    dt <- data.table::as.data.table(dt)

  if ("period" %in% names(dt)) dt[, period := NULL]

  x_cols <- attr(terms(xformula), "term.labels")
  inv_cols <- setdiff(names(dt), c("id", "y", "post", x_cols))

  ## wide Y
  cast_formula <- as.formula(
    paste("id +", paste(inv_cols, collapse = " +"), "~ post")
  )
  wide_y <- data.table::dcast(dt, cast_formula, value.var = "y")
  data.table::setnames(wide_y, c("0", "1"), c("y0", "y1"))

  ## covariates at post == 0
  if (length(x_cols)) {
    base_x <- dt[post == 0, .SD[1L], by = id, .SDcols = x_cols]
    data.table::setkey(base_x, id)     # keyed join is fastest
    wide_y <- base_x[wide_y]           # left join keeps order
  }

  wide_y[]
}


#' Stratified K-fold splits for cross-fitting
#' Quickly builds *train* and *test* index lists that preserve the class
#' proportions of a categorical variable in every fold
#' @param data  A `data.frame` / `data.table` with at least the column named
#'              in `strat_col`. This comes from preprocessing steps.
#' @param strat_col Character scalar. Column used for stratification (e.g.
#'                  `"subgroup"`).
#' @param n_folds   Integer ≥ 2. Number of folds. Default is 2.
#' @param seed      Optional integer. If supplied, the function calls
#'                  `set.seed()` for reproducible splits.
#'
#' @return A list with two components — `test_ids` and `train_ids` — each a
#'         length-K list of integer vectors.  This is the exact structure
#'         expected by `DoubleML`:
#'         \preformatted{
#'         splits$train_ids[[j]]  # rows used to fit learners   (fold j)
#'         splits$test_ids [[j]]  # rows used to compute scores (fold j)
#'         }
#' @keywords internal
make_stratified_folds <- function(data,
                                  strat_col,
                                  n_folds = 2,
                                  seed    = NULL) {
  if (!is.null(seed)) set.seed(seed)

  ## ----- stratified test indices ----------------------------------------
  fold_test <- unname(caret::createFolds(data[[strat_col]], k = n_folds))

  ## ----- stratified train indices  -------------------------------------
  fold_train <- lapply(fold_test,
                       function(test) setdiff(seq_len(nrow(data)), test))

  list(fold_test  = fold_test,
       fold_train = fold_train)
}


#' Function to get local sample indices for "local" DiDs comparisons.
#' @param units_selected Vector of units selected for the local sample
#' @param n_folds Number of folds
#' @param global_fold_test List of test indices for global folds
#' @param global_fold_train List of train indices for global folds
#' @return A list containing train and test indices for each fold at the local level.
#' @keywords internal
get_local_folds <- function(units_selected,
                            n_folds,
                            global_fold_test,
                            global_fold_train) {

  ## ---------- helper to remap one fold -------------------------------------
  remap <- function(test_vec, train_vec) {
    # keep only indices that occur in the subset, then translate
    list(
      test  = match(intersect(test_vec,  units_selected), units_selected),
      train = match(intersect(train_vec, units_selected), units_selected)
    )
  }

  folds_local <- lapply(
    seq_len(n_folds),
    function(j) remap(global_fold_test[[j]], global_fold_train[[j]])
  )

  list(list(
    train_ids = lapply(folds_local, `[[`, "train"),
    test_ids  = lapply(folds_local, `[[`, "test")
  ))
}


#' Function to compute the DML ATT + score function
#' @param data A data.table containing the data processed
#' @param xformula The formula for the propensity score model
#' @param ml_pa `mlr3` learner for propensity score estimation
#' @param ml_md `mlr3` learner for regression adjustment estimation
#' @param fold_test List of test indices for global folds
#' @param fold_train List of train indices for global folds
#' @param n_folds The number of folds for cross-fitting
#' @return A list containing ids, estimator for each k fold, and scores
#' @keywords internal
compute_dml_nuisances <- function(data,
                                  xformula,
                                  ml_pa,
                                  ml_md,
                                  fold_test,
                                  fold_train,
                                  n_folds){
  # Initialize lists to store results
  # index of comparison groups;
  # 3:= S=g, Q=0
  # 2:= S=0, Q=1
  # 1:= S=0, Q=0
  compare_values = c(3, 2, 1)

  att_vec  <- numeric(length(compare_values))
  names(att_vec) <- as.character(compare_values)

  psi_mat  <- matrix(0, nrow = nrow(data), ncol = length(compare_values))
  colnames(psi_mat) <- paste0("psi", compare_values)

  # get covariates list names
  x_cols <- attr(terms(xformula), "term.labels")

  ## ---- loop over comparison groups ------------------------------------
  for (idx in seq_along(compare_values)){
    comp <- compare_values[idx]

    # Note: processed data already has the column "subgroup". "4" is the treated group
    units_keep <- which(data[["subgroup"]] %in% c(4, comp))

    # build local folds
    smpls <- get_local_folds(units_keep, n_folds, fold_test, fold_train)

    # prepare DoubleML data: build a *pairwise* DiD inside the subset
    dtmp  <- data[units_keep, c("deltaY", "D", x_cols), with = FALSE]

    dml_data <- DoubleML::DoubleMLData$new(
      data  = dtmp,
      y_col  = "deltaY",
      d_cols = "D",
      x_cols = x_cols
    )

    obj <- DoubleML::DoubleMLIRM$new(
      data  = dml_data,
      ml_g  = ml_md, # outcome model
      ml_m  = ml_pa, # propensity score model
      score = "ATTE",
      draw_sample_splitting = FALSE
    )
    obj$set_sample_splitting(smpls) # we provide our own sampling scheme adapted to DDD
    obj$fit()

    ## ---- store results ------------------------------------------------
    att_vec[idx]            <- obj$coef[1]
    psi_mat[units_keep, idx] <- obj$psi[, 1, 1] # AIPW score
  }

  # assemble influence matrix
  influence_dt <- data.table::data.table(
    id  = data$id,
    trt = data$D,
    psi_mat
  )

  # return results
  return(list(att = att_vec,
              influence_matrix = influence_dt
              )
         )
}

# Function to compute se for DML ATT
compute_se_dml <- function(influence_matrix) {
  n <- nrow(influence_matrix) # data is already wide data
  influence_matrix[, inf := psi3 + psi2 - psi1]

  se_dml_att <- stats::sd(influence_matrix[, inf])/sqrt(n)
  inf <- influence_matrix[, inf]
  return(list(
    se = se_dml_att,
    inf_func = inf
  ))
}
