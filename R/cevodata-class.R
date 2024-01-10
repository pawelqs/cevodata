
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

#' Get/Add SNV/CNA data from the cevodata dataset
#' @param object object
#' @param snvs tibble with SNVs
#' @param cnas tibble with CNAs
#' @param name name for SNVs/CNAs assay
#' @param which assay to use - uses active_SNVs if none
#' @name assays
NULL


#' @describeIn assays Add new SNVs to cevodata
#' @export
add_SNV_data <- function(object, snvs, name = NULL) {
  if (is.null(name)) {
    n <- length(object$SNVs)
    name <- if (n == 0) "snvs" else str_c("snvs", n)
  }
  object$SNVs[[name]] <- as_cevo_snvs(snvs)
  default_SNVs(object) <- name
  meta <- tibble(sample_id = unique(snvs$sample_id))
  object <- add_sample_data(object, meta)
  object
}


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
  object <- add_sample_data(object, meta)
  object
}


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
    keys <- intersect(id_cols, colnames(data))
    object$metadata <- full_join(object$metadata, data, by = keys)
  }

  object$metadata <- object$metadata |>
    select(all_of(id_cols), everything())

  object
}


# -------------------------------- Defaults ------------------------------------

#' Get/Set active assays of the cevodata object
#' @param object object
#' @param value name of new default assay
#' @param ... other arguments
#' @name active_assays
NULL


#' @describeIn active_assays Get default SNVs assay of cevodata
#' @export
default_SNVs <- function(object, ...) {
  object$settings$active_SNVs
}


#' @describeIn active_assays Set default SNVs assay of cevodata
#' @export
`default_SNVs<-` <- function(object, ..., value) {
  if (value %not in% names(object$SNVs)) {
    stop("Chosen SNV assay must exist in object$SNVs")
  }
  object$settings$active_SNVs <- value
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
use_purity <- function(cd, name, verbose = get_cevomod_verbosity()) {
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
