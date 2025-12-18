#' @useDynLib triplediff
#' @importFrom Rcpp sourceCpp
NULL
utils::globalVariables(c("first_treat", "treat", "tail", "period", "id", "subgroup", "partition",
                         "post", "deltaY", "y1", "y0", "y", "..cols_to_keep","asif_never_treated","treated_first_period",
                         "V1", "..tname", "..pname", "..idname", "..gname", "psi1", "psi2", "psi3", ".", ".rcs_id", "N", "control"))
