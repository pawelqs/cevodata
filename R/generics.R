#' Split object
#' @param object object to split
#' @param ... other arguments
#' @export
split_by <- function(object, ...) {
  UseMethod("split_by")
}
