#' @title Summary
#' @description Summary of an aggddd object
#' @param object A aggddd object
#' @param ... Other params (required as generic function, but not used)
#' @export
#' @noRd
# Define new summary function
summary.aggddd <- function(object, ...){
  aggddd.obj <- object
  print(aggddd.obj)
}
