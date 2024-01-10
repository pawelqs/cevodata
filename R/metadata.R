
# --------------------------- cevodata functions ------------------------------

#' Add metadata to the cevodata object
#' @param object object
#' @param data name of new default assay
add_metadata <- function(object, data) {
  id_cols <- c("patient_id", "sample_id")
  if (!any(id_cols %in% colnames(data))) {
    stop("Metadata must have patient_id or sample_id column!")
  }

  if (is.null(object$metadata)) {
    object$metadata <- data
  } else {
    meta <- get_metadata(object)
    keys <- intersect(colnames(meta), colnames(data))
    meta <- full_join(meta, data, by = keys)
  }

  object$metadata <- meta |>
    select(any_of(id_cols), everything())

  object
}


#' Get sample metadata
#' @param cd cevodata object
#' @export
get_metadata <- function(cd) {
  cd$metadata
}


#' Choose purity measure
#'
#' <cevodata> metadata can contain purity measures in columns other than 'purity'.
#' T his function can be used to set 'purity' values using values from requested
#' column
#'
#' @param cd <cevodata> object
#' @param name Name of the metadata column with chosen purity values
#' @param verbose Verbose?
#' @export
use_purity <- function(cd, name, verbose = verbose::verbose("cevoverse")) {
  if (name %not in% names(cd$metadata)) {
    stop(
      "`name` should be a name of the column in the metadata tibble, ",
      "which should be used as purity measure"
    )
  } else {
    msg("Using '", name, "' as default purity measure", verbose = verbose)
    if (!is.null(cd$metadata[["purity"]])) {
      cd$metadata$prev_purity <- cd$metadata$purity
    }
    cd$metadata$purity <- cd$metadata[[name]]
    cd
  }
}


# ---------------------------------- Other -----------------------------------

get_purities <- function(cd) {
  cd$metadata |>
    select("sample_id", "purity")
}


get_patients_data <- function(metadata) {
  patient_data_cols <- metadata |>
    group_by(.data$patient_id) |>
    summarise_all(n_distinct) |>
    map(~all(.x == 1)) |>
    keep(~.x) |>
    names()
  metadata |>
    select("patient_id", all_of(patient_data_cols))
}


get_sample_ids <- function(cd) {
  cd$metadata$sample_id
}


# get_patient_sex <- function(cd) {
#   cd$metadata |>
#     select("sample_id", "sex")
# }
