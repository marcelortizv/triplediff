#' @title Multiplier Bootstrap
#'
#' @description This function take an influence function and use the
#'  multiplier bootstrap to compute standard errors and critical values for
#'  uniform confidence bands.
#'
#' @param inf_func an influence function
#' @param did_preprocessed A `dp` object obtained after preprocess
# #' @param multiple_periods Boolean of whether or not to use multiple periods. This is useful to get the right tlist. Default is \code{multiple_periods=TRUE}.
#' @param use_parallel Boolean of whether or not to use parallel processing in the multiplier
#'  bootstrap, default is \code{use_parallel=FALSE}
#' @param cores the number of cores to use with parallel processing, default is \code{cores=1}
#'
#' @return A list with following elements
#' \item{bres}{results from each bootstrap iteration}
#' \item{V}{variance matrix}
#' \item{se}{standard errors}
#' \item{crit_val}{a critical value for computing uniform confidence bands}
#'
#' @export
mboot <- function(inf_func, did_preprocessed, use_parallel = FALSE, cores = 1) {

  # setup needed variables
  data <- did_preprocessed$preprocessed_data # we only need data for first period
  cluster <- did_preprocessed$cluster
  biters <- did_preprocessed$nboot
  alpha <- did_preprocessed$alpha
  tlist <- sort(unique(data$period))


  # just get n observations. Only for panel data
  dta <- data[period == tlist[1]]

  # Make sure inf_func is matrix because we need this for computing n below
  inf_func <- as.matrix(inf_func)

  # set correct number of units
  n <- nrow(inf_func)

  # multiplier bootstrap
  n_clusters <- n
  if (length(cluster)==0) {
    bres <- sqrt(n) * run_multiplier_bootstrap(inf_func, biters, use_parallel, cores)
  } else {
    # Compute multiplier bootstrap for clustered standard errors

    # Extract the unique clusters along with their IDs
    unique_clusters <- dta[, c("id", "cluster")]
    # Count the number of unique clusters
    n_clusters <- length(unique(unique_clusters$cluster))
    # Count the number of observations in each cluster
    cluster_counts <- as.vector(table(unique_clusters$cluster))
    # Compute the mean influence function per cluster
    cluster_mean_if <- rowsum(inf_func, unique_clusters$cluster, reorder = TRUE) / cluster_counts

    # Run the bootstrap procedure
    bres <- sqrt(n_clusters) * run_multiplier_bootstrap(cluster_mean_if, biters, use_parallel, cores)

  }


  if (isTRUE(class(bres) == "numeric")) bres <- as.matrix(bres)

  # Non-degenerate dimensions
  ndg.dim <- (!is.na(colSums(bres))) & (base::colSums(bres^2) > sqrt(.Machine$double.eps)*10)
  bres <- as.matrix(bres[ , ndg.dim])

  # bootstrap variance matrix (this matrix can be defective because of degenerate cases)
  V <- cov(bres)
  # bootstrap standard error
  bSigma <- apply(bres, 2,
                  function(b) (quantile(b, .75, type=1, na.rm = T) -
                                 quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))

  bSigma[bSigma <= sqrt(.Machine$double.eps)*10] <- NA

  # critical value for uniform confidence band
  bT <- base::suppressWarnings(apply(bres, 1, function(b) max(abs(b / bSigma), na.rm = T)))
  bT <- bT[is.finite(bT)]
  crit_val <- quantile(bT, 1-alpha, type=1, na.rm = T) # uniform critical value
  se <- rep(NA, length(ndg.dim))
  se[ndg.dim] <- as.numeric(bSigma) / sqrt(n_clusters)

  return(list(bres = bres, V = V, se = se, bT= bT, unif_crit_val = crit_val))
}

run_multiplier_bootstrap <- function(inf_func, biters, use_parallel = FALSE, cores = 1) {
  ngroups = ceiling(biters/cores)
  chunks = rep(ngroups, cores)
  # Round down so you end up with the right number of biters
  chunks[1] = chunks[1] + biters - sum(chunks)

  n <- nrow(inf_func)
  parallel.function <- function(biters) {
    BMisc::multiplier_bootstrap(inf_func, biters)
  }
  # From tests, this is about where it becomes worth it to parallelize
  if(n > 2500 & use_parallel == TRUE & cores > 1) {
    results = parallel::mclapply(
      chunks,
      FUN = parallel.function,
      mc.cores = cores
    )
    results = do.call(rbind, results)
  } else {
    results = BMisc::multiplier_bootstrap(inf_func, biters)
  }
  return(results)
}
