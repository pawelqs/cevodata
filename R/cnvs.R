
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
