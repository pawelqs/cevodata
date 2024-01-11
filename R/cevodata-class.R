
# ---------------------------- cevodata constructors ---------------------------

new_cevodata <- function(name) {
  cd <- list(
    name = name,
    metadata = NULL,
    SNVs = list(),
    CNAs = list(),
    stats = list(),
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


# ------------------------------ Adding data -----------------------------------

#' Add metadata to the cevodata object
#' @param object object
#' @param data name of new default assay
add_metadata <- function(object, data) {
  id_cols <- c("patient_id", "sample_id")
  if (!any(id_cols %in% colnames(data))) {
    stop("Metadata must have patient_id or sample_id column!")
  }

  if (is.null(object$metadata)) {
    meta <- data
  } else {
    meta <- get_metadata(object)
    keys <- intersect(colnames(meta), colnames(data))
    meta <- full_join(meta, data, by = keys)
  }

  object$metadata <- meta |>
    select(any_of(id_cols), everything())

  object
}


#' Get/Add SNV/CNA data from the cevodata dataset
#' @param object cevodata object
#' @param snvs Tibble with SNVs
#' @param cnas Tibble with CNAs
#' @param models cv_subitem object with models
#' @param name Name for SNVs/CNAs/models
#' @param ... Other arguments
#' @name cevodata_components
NULL


#' @describeIn cevodata_components Add new SNVs to cevodata
#' @export
add_SNV_data <- function(object, snvs, name = NULL) {
  if (is.null(name)) {
    n <- length(object$SNVs)
    name <- if (n == 0) "snvs" else str_c("snvs", n)
  }
  object$SNVs[[name]] <- as_cevo_snvs(snvs)
  default_SNVs(object) <- name
  meta <- tibble(sample_id = unique(snvs$sample_id))
  object <- add_metadata(object, meta)
  object
}


#' @describeIn cevodata_components Add new CNAs to cevodata
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


add_stats <- function(object, stats, name) {
  object$stats[[name]] <- stats
  object
}


#' @describeIn cevodata_components Add new models to cevodata
#' @export
add_models <- function(object, models, name) {
  object$models[[name]] <- as_cv_subitem(models)
  active_models(object) <- name
  object
}


# ------------------------------ Settings --------------------------------------

#' Get/Set active assays of the cevodata object
#' @param object object
#' @param value name of new default assay
#' @param ... other arguments
#' @name settings
NULL


#' @describeIn settings Get default SNVs assay of cevodata
#' @export
default_SNVs <- function(object, ...) {
  object$settings$active_SNVs
}


#' @describeIn settings Set default SNVs assay of cevodata
#' @export
`default_SNVs<-` <- function(object, ..., value) {
  if (value %not in% names(object$SNVs)) {
    stop("Chosen SNV assay must exist in object$SNVs")
  }
  object$settings$active_SNVs <- value
  object
}


#' @describeIn settings Get default CNAs assay of cevodata
#' @export
default_CNAs <- function(object, ...) {
  object$settings$active_CNAs
}


#' @describeIn settings Set default CNAs assay of cevodata
#' @export
`default_CNAs<-` <- function(object, ..., value) {
  if (value %not in% names(object$CNAs)) {
    stop("Chosen CNA assay must exist in object$CNAs")
  }
  object$settings$active_CNAs <- value
  object
}


#' @describeIn settings Get active models from cevodata
#' @export
active_models <- function(object, ...) {
  object$settings$active_models
}


#' @describeIn settings Set active models in cevodata
#' @export
`active_models<-` <- function(object, ..., value) {
  if (value %not in% names(object$models)) {
    stop("Chosen models must exist in object$models")
  }
  object$settings$active_models <- value
  object
}


# ------------------------------ Basic getters ---------------------------------
# more getters in other files

#' @return tibble
#' @describeIn cevodata_components Get SNVs from cevodata dataset
#' @export
SNVs.cevodata <- function(object, name = default_SNVs(object), ...) {
  if (name %not in% names(object$SNVs)) {
    stop(str_c(name, " does not exist in object$SNVs"))
  }
  object$SNVs[[name]]
}


#' @describeIn cevodata_components Get CNAs from cevodata dataset
#' @export
CNAs.cevodata <- function(object, name = default_CNAs(object), ...) {
  if (name %not in% names(object$CNAs)) {
    stop(str_c(name, " does not exist in object$CNAs"))
  }
  object$CNAs[[name]]
}


#' Get sample metadata
#' @param cd cevodata object
#' @export
get_metadata <- function(cd) {
  cd$metadata
}


#' @describeIn cevodata_components Get models cevodata dataset
#' @export
get_models <- function(object, name = active_models(object)) {
  if (name %not in% names(object$models)) {
    stop(str_c(name, " does not exist in object$models"))
  }
  object$models[[name]]
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
  new$transformations <- list()
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
