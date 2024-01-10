
validate_CNAs <- function(cnas) {
  required_cols <- c(
    "sample_id", "chrom", "start", "end"
    # "log_ratio", "BAF", "total_cn", "major_cn", "minor_cn"
  )
  missing_cols <- setdiff(required_cols, names(cnas))
  if (length(missing_cols)) {
    stop(str_c("cnas object is missing the following columns:", str_c(missing_cols, collapse = ", ")))
  }
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
  cnas <- CNAs(object, which_cnas) |>
    left_join(get_patient_sex(object), by = "sample_id") |>
    mutate(
      normal_cn = if_else(
        .data$sex %in% c("male", "M") & .data$chrom %in% c("chrX", "chrY"), 1, 2
      )
    )
  object |>
    add_CNA_data(cnas, name = which_cnas)
}
