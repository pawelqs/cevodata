
# --------------------------- cevodata functions ------------------------------

#' @describeIn assays Add new CNAs to cevodata
#' @export
add_CNA_data <- function(object, cnas, name = NULL) {
  if (is.null(name)) {
    n <- length(object$CNAs)
    name <- if (n == 0) "cnas" else str_c("cnas", n)
  }
  object$CNAs[[name]] <- validate_CNAs(cnas)
  default_CNAs(object) <- name
  meta <- tibble(sample_id = unique(cnas$sample_id))
  object <- add_metadata(object, meta)
  object
}


#' @describeIn active_assays Get default CNAs assay of cevodata
#' @export
default_CNAs <- function(object, ...) {
  object$settings$active_CNAs
}


#' @describeIn active_assays Set default CNAs assay of cevodata
#' @export
`default_CNAs<-` <- function(object, ..., value) {
  if (value %not in% names(object$CNAs)) {
    stop("Chosen CNA assay must exist in object$CNAs")
  }
  object$settings$active_CNAs <- value
  object
}


# --------------------------- Other functions --------------------------------

validate_CNAs <- function(cnas) {
  require_columns(cnas, c("sample_id", "chrom", "start", "end"))
  cnas
}


#' Annotate chromosome ploidies in CNA data
#'
#' Adds the normal_cn column to the data. This column is required for e.g.
#' by Dentro CCF calculation method. Requires 'sex' column in the metadata.
#' Males should be encoded by 'M' or "male'.
#'
#' @param object <cevodata> object
#' @param which_cnas Name of the CNAs slot
#'
#' @return <cevodata> object
#' @export
annotate_normal_cn <- function(object, which_cnas = default_CNAs(object)) {
  msg("Assuming human genome")
  sex <- get_metadata(object) |>
    select("sample_id", "sex")
  cnas <- CNAs(object, which_cnas) |>
    left_join(sex, by = "sample_id") |>
    mutate(
      normal_cn = if_else(
        .data$sex %in% c("male", "M") & .data$chrom %in% c("chrX", "chrY"), 1, 2
      )
    )
  object |>
    add_CNA_data(cnas, name = which_cnas)
}
