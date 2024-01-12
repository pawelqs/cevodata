
#' @rdname cevodata_components
#' @export
SNVs <- function(object, ...) {
  UseMethod("SNVs")
}


#' @rdname cevodata_components
#' @export
CNAs <- function(object, ...) {
  UseMethod("CNAs")
}


#' Split object
#' @param object object to split
#' @param ... other arguments
#' @export
split_by <- function(object, ...) {
  UseMethod("split_by")
}
