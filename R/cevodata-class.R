
# ---------------------------- cevodata constructors ---------------------------

new_cevodata <- function(name) {
  cd <- list(
    name = name,
    metadata = NULL,
    SNVs = list(),
    CNAs = list(),
    models = list(),
    misc = list(),
    settings = list(
      active_SNVs = NULL,
      active_CNAs = NULL,
      active_models = NULL
    )
  )
  structure(cd, class = "cevodata")
}


#' Create new cevomod dataset object
#'
#' @param name dataset name
#' @param snvs tibble with SNVs
#' @param snvs_name name for SNVs assay
#' @param cnas tibble with CNAs
#' @param cnas_name name for CNAs assay
#' @return `cevodata` object
#'
#' @export
init_cevodata <- function(name = "Unnamed dataset",
                          snvs = NULL, snvs_name = NULL,
                          cnas = NULL, cnas_name = NULL) {
  cd <- new_cevodata(name)
  if (!is.null(snvs)) {
    cd <- add_SNV_data(cd, snvs, snvs_name)
  }
  if (!is.null(cnas)) {
    cd <- add_CNA_data(cd, cnas, cnas_name)
  }
  cd
}


# --------------------------------- Add data -----------------------------------

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


# --------------------------------- Misc ---------------------------------------

#' Update to cevodata to v3
#' @param cd <cevodata> object
#' @export
update_cevodata_v2_to_v3 <- function(cd) {
  new <- init_cevodata(cd$name)
  new$metadata <- cd$metadata
  new$SNVs <- cd$SNVs
  new$CNAs <- cd$CNVs
  new$settings <- list(
    active_SNVs = cd$active_SNVs,
    active_CNAs = cd$active_CNAs,
    active_models = cd$active_models
  )
  new
}


is_cevodata_singlepatient <- function(object) {
  n_patients <- count_patients(object)
  if (is.na(n_patients)) {
    return(FALSE)
  } else if (n_patients > 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


count_patients <- function(cd) {
  if (!is.null(cd$metadata[["patient_id"]])) {
    n_distinct(cd$metadata$patient_id)
  } else {
    NA_integer_
  }
}
