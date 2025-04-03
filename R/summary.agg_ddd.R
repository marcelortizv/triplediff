#' @title Summary
#' @description Summary of an aggddd object
#' @param object A aggddd object
#' @param ... Other params (required as generic function, but not used)
#' @export
#' @noRd
# Define new summary function
summary.agg_ddd <- function(object, ...){
  agg_ddd.obj <- object
  print(agg_ddd.obj)
}
