#' Get extra sampe/patient stats
#' @param object cevodata object
#' @param snvs Which SNVs to use
#' @name misc_stats
NULL


## ------------------------------- Sample Stats --------------------------------

#' @describeIn misc_stats Get misc sample level stats
#' @export
calc_all_sample_stats <- function(object, snvs = default_SNVs(object)) {
  stats <- list(
    calc_sample_mutation_burden(object, snvs = snvs)
  )

  stats_merged <- if (length(stats) == 1) {
    stats[[1]]
  } else {
    reduce(stats, full_join, by = "sample_id")
  }

  add_stats(object, stats_merged, name = "sample_stats")
}


calc_sample_mutation_burden <- function(object, snvs = default_SNVs(object), ...) {
  snvs <- SNVs(object, which = which)
  if ("mutation_id" %not in% colnames(snvs)) {
    snvs <- unite_mutation_id(snvs)
  }

  SNVs(object, which = snvs) |>
    filter(.data$VAF > 0, !is.na(.data$VAF)) |>
    select("sample_id", "mutation_id") |>
    unique() |>
    group_by(.data$sample_id) |>
    summarise(mutation_burden = n(), .groups = "drop")
}


# ------------------------------- Patient Stats --------------------------------

#' @describeIn misc_stats Get misc patient level stats
#' @export
calc_all_patient_stats <- function(object, snvs = default_SNVs(object)) {
  stats <- list(
    calc_patient_mutation_burden(object, snvs = snvs)
  )

  stats_merged <- if (length(stats) == 1) {
    stats[[1]]
  } else {
    reduce(stats, full_join, by = "patient_id")
  }

  add_stats(object, stats_merged, name = "patient_stats")
}


calc_patient_mutation_burden <- function(object, snvs = default_SNVs(object)) {
  snvs <- SNVs(object, which = which)
  if ("mutation_id" %not in% colnames(snvs)) {
    snvs <- unite_mutation_id(snvs)
  }

  snvs |>
    filter(.data$VAF > 0) |>
    join_metadata(object) |>
    group_by(.data$patient_id) |>
    select("patient_id", "mutation_id") |>
    unique() |>
    summarise(mutation_burden = n(), .groups = "drop")
}
