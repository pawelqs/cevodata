#' Calc extra sample/patient stats
#'
#' Calculated stats are saved to appropriate cevodata$stats slots and can be
#' obtained with get_stats() function.
#'
#' @param object cevodata object
#' @param snvs_name Which SNVs to use
#' @name misc_stats
NULL


## ------------------------------- Sample Stats --------------------------------

#' @describeIn misc_stats Calc misc sample level stats
#'
#' Currently calculates only the sample mutation burden
#'
#' @export
calc_all_sample_stats <- function(object, snvs_name = default_SNVs(object)) {
  stats <- list(
    calc_sample_mutation_burden(object, snvs_name = snvs_name)
  )

  stats_merged <- if (length(stats) == 1) {
    stats[[1]]
  } else {
    reduce(stats, full_join, by = "sample_id")
  }

  add_stats(object, stats_merged, name = "sample_stats")
}


calc_sample_mutation_burden <- function(object, snvs_name = default_SNVs(object), ...) {
  snvs <- SNVs(object, name = snvs_name)
  if ("mutation_id" %not in% colnames(snvs)) {
    snvs <- unite_mutation_id(snvs)
  }

  snvs |>
    filter(.data$VAF > 0, !is.na(.data$VAF)) |>
    select("sample_id", "mutation_id") |>
    unique() |>
    group_by(.data$sample_id) |>
    summarise(mutation_burden = n(), .groups = "drop")
}


# ------------------------------- Patient Stats --------------------------------

#' @describeIn misc_stats Calc misc patient level stats
#'
#' Currently calculates only the patient mutation burden
#'
#' @export
calc_all_patient_stats <- function(object, snvs_name = default_SNVs(object)) {
  stats <- list(
    calc_patient_mutation_burden(object, snvs_name = snvs_name)
  )

  stats_merged <- if (length(stats) == 1) {
    stats[[1]]
  } else {
    reduce(stats, full_join, by = "patient_id")
  }

  add_stats(object, stats_merged, name = "patient_stats")
}


calc_patient_mutation_burden <- function(object, snvs_name = default_SNVs(object)) {
  snvs <- SNVs(object, name = snvs_name)
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


# ---------------------------- Sample-sample stats -----------------------------

#' @describeIn misc_stats Calc sample-sample stats
#'
#' Currently calculates only the Jaccard index
#'
#' @export
calc_sample_sample_stats <- function(object, snvs_name = default_SNVs(object)) {
  stats <- list(
    if ("sample" %in% names(get_metadata(object))) calc_Jaccard_indexes(object, snvs_name = snvs_name)
  )
  stats <- stats[!map_lgl(stats, is.null)]

  if (length(stats) == 0) {
    warning("No sample-sample stats calculated")
    return(object)
  } else if (length(stats) == 1) {
    stats_merged <- stats[[1]]
  } else {
    stats_merged <- reduce(stats, full_join, by = "patient_id")
  }

  add_stats(object, stats_merged, name = "sample_sample_stats")
}


calc_Jaccard_indexes <- function(object, snvs_name = default_SNVs(object)) {
  require_columns(get_metadata(object), c("patient_id", "sample"))
  patients_to_keep <- get_metadata(object) |>
    count(.data$patient_id) |>
    filter(n > 1) |>
    pull("patient_id")

  snvs <- SNVs(object, name = snvs_name) |>
    filter(.data$VAF > 0) |>
    join_metadata(object) |>
    filter(.data$patient_id %in% patients_to_keep)

  if ("mutation_id" %not in% colnames(snvs)) {
    snvs <- unite_mutation_id(snvs)
  }

  dt <- snvs |>
    select("patient_id", "sample", "mutation_id") |>
    unique() |>
    mutate(sample = as.character(.data$sample)) |>
    nest_by(.data$patient_id)

  Jaccard <- dt |>
    reframe(.calc_Jaccard_index(.data$data)) |>
    rename(sample1 = "group1", sample2 = "group2")

  if (any(duplicated(Jaccard$patient_id))) {
    Jaccard |>
      nest_by(.data$patient_id) |>
      rename(Jaccard_indexes = "data")
  } else {
    Jaccard
  }
}


.calc_Jaccard_index <- function(tbl) {
  groups_vec <- tbl[[1]]
  items_vec <- tbl[[2]]
  groups <- unique(tbl[[1]])
  res <- get_combinations_tbl(groups) |>
    set_names(c("group1", "group2")) |>
    mutate(Jaccard_index = NA_real_)
  for(i in 1:nrow(res)) {
    A <- items_vec[groups_vec == res$group1[[i]]]
    B <- items_vec[groups_vec == res$group2[[i]]]
    res$Jaccard_index[[i]] <- length(intersect(A, B)) / length(union(A, B))
  }
  res
}


get_combinations_tbl <- function(items) {
  utils::combn(items, m = 2) |>
    t() |>
    `colnames<-`(c("item1", "item2")) |>
    as_tibble()
}
