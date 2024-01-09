
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
init_cevodata <- function(name = "cevodata object",
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


#' Get/Add SNV/CNA data from the cevodata dataset
#' @param object object
#' @param snvs tibble with SNVs
#' @param cnas tibble with CNAs
#' @param name name for SNVs/CNAs assay
#' @param which assay to use - uses active_SNVs if nonei
#' @param ... other arguments
#' @name assays
NULL


#' @describeIn assays Add new SNVs to cevodata
#' @export
add_SNV_data <- function(object, snvs, name = NULL, ...) {
  if(is.null(name)) {
    n <- length(object$SNVs)
    name <- if (n == 0) "snvs" else str_c("snvs", n)
  }
  snvs <- as_cevo_snvs(snvs)
  object$SNVs[[name]] <- snvs
  default_SNVs(object) <- name
  meta <- snvs |>
    select("sample_id") |>
    unique() |>
    as_tibble()
  object <- add_sample_data(object, meta)
  object
}


#' @describeIn assays Add new CNAs to cevodata
#' @export
add_CNA_data <- function(object, cnas, name = NULL, ...) {
  if(is.null(name)) {
    n <- length(object$CNAs)
    name <- if (n == 0) "cnas" else str_c("cnas", n)
  }
  validate_CNAs(cnas)
  object$CNAs[[name]] <- cnas
  default_CNAs(object) <- name
  meta <- cnas |>
    select("sample_id") |>
    unique() |>
    as_tibble()
  object <- add_sample_data(object, meta)
  object
}


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


#' Add metadata to the cevodata object
#' @param object object
#' @param data name of new default assay
#' @param ... other arguments
#' @name cevo_metadata
NULL


#' @describeIn cevo_metadata Add patient data to cevodata object
#' @export
add_patient_data <- function(object, data, ...) {
  if ("patient_id" %not in% colnames(data)) {
    stop("Data must have patient_id column!")
  }
  if (is.null(object$metadata)) {
    object$metadata <- data
  } else {
    keys <- intersect(c("patient_id", "sample_id"), colnames(data))
    object$metadata <- full_join(object$metadata, data, by = keys)
  }
  object
}


#' @describeIn cevo_metadata Add samples' data to cevodata object
#' @export
add_sample_data <- function(object, data, ...) {
  if ("sample_id" %not in% colnames(data)) {
    stop("Data must have sample_id column!")
  }
  if (is.null(object$metadata)) {
    object$metadata <- data
  } else {
    keys <- c("patient_id", "sample_id", "sample") |>
      intersect(colnames(data)) |>
      intersect(colnames(object$metadata))
    object$metadata <- full_join(object$metadata, data, by = keys)
  }
  if (all(c("patient_id", "sample_id", "sample") %in% colnames(object$metadata))) {
    object$metadata <- object$metadata |>
      select("patient_id", "sample_id", "sample", everything())
  }
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
