#' @title Summary
#' @description Summary of a ddd object
#' @param object A ddd object
#' @param ... Other params (required as generic function, but not used)
#' @export
#' @noRd
# Define new summary function
summary.ddd <- function(object, ...){
  ddd.obj <- object
  print(ddd.obj)
}
