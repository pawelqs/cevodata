
#' Intervalize the mutation frequencies
#'
#' Intervalize the mutation frequencies for the subsequent analyses and plots.
#' Adds f_interval column to the SNV tibble
#'
#' @param object <cevodata> object
#' @param which_snvs Which SNVs to use
#' @param column Which frequency measure column to intervalize? By default, uses
#'   the CCF/2  if it is found in the SNV tibble, and VAF otherwise
#' @param bins Number of interval bins
#' @param ... Other args
#' @param verbose Verbose?
#'
#' @return <cevodata> object
#' @export
intervalize_mutation_frequencies <- function(object, ...) {
  UseMethod("intervalize_mutation_frequencies")
}


#' @rdname intervalize_mutation_frequencies
#' @export
intervalize_mutation_frequencies.cevodata <- function(object,
                                                      which_snvs = default_SNVs(object),
                                                      column = get_snvs_frequency_measure(object, which_snvs),
                                                      bins = NULL,
                                                      verbose = get_verbosity(),
                                                      ...) {
  snvs <- SNVs(object, which_snvs) |>
    intervalize_mutation_frequencies(column, bins, verbose)
  object |>
    add_SNV_data(snvs, name = which_snvs)
}


#' @rdname intervalize_mutation_frequencies
#' @export
intervalize_mutation_frequencies.cevo_snvs <- function(object,
                                                       column = get_snvs_frequency_measure(object),
                                                       bins = NULL,
                                                       verbose = get_verbosity(),
                                                       ...) {
  msg("Calculating f intervals, using ", column, " column", verbose = verbose)
  object |>
    cut_f_intervals(bins = bins, column = column) |>
    as_cevo_snvs()
}


# Get mutation frequency intervals
#
# Cuts requested column e.g. CCF or VAF into intervals. Adds f_interval and f
# columns to the input SNVs tibble
cut_f_intervals <- function(snvs, column, bins = NULL) {
  breaks <- get_interval_breaks(snvs, bins = bins) |>
    enframe(name = "sample_id", value = "breaks")
  res <- snvs |>
    nest_by(.data$sample_id) |>
    left_join(breaks, by = "sample_id") |>
    mutate(data = list(cut_f(.data$data, .data$breaks, column = column))) |>
    ungroup()
  interval_levels <- res |>
    select("sample_id", "data") |>
    deframe() |>
    map("f_interval") |>
    map(levels)
  res$data <- res$data |>
    map(~mutate(.x, across(f_interval, as.character)))
  res <- res |>
    select("sample_id", "data") |>
    unnest("data") |>
    mutate(f = get_interval_centers(.data$f_interval), .after = "f_interval")
  attr(res, "intervals") <- interval_levels
  attr(res, "f_column") <- column
  as_cevo_snvs(res)
}


cut_f <- function(tbl, breaks, column) {
  tbl |>
    mutate(
      f_interval = cut(.data[[column]], breaks = breaks),
      .after = all_of(column)
    )
}


get_interval_breaks <- function(object, bins = NULL, sample_id = NULL) {
  if (is.null(bins)) {
    bins_by_sample <- get_sample_sequencing_depths(object) |>
      transmute(
        .data$sample_id,
        bins = round(.data$median_DP),
        bins = if_else(.data$bins > 100, 100, bins)
      )
  } else {
    bins_by_sample <- tibble(
      sample_id = unique(object$sample_id),
      bins = bins
    )
  }
  breaks <- bins_by_sample |>
    deframe() |>
    map(~c(-1/.x, seq(0, 1, length.out = .x + 1)))

  if (is.null(sample_id)) {
    breaks
  } else {
    breaks[[sample_id]]
  }
}


get_sample_sequencing_depths <- function(snvs) {
  sequencing_depths <- snvs |>
    group_by(.data$sample_id) |>
    summarise(
      mean_DP = mean(.data$DP, na.rm = TRUE),
      median_DP = stats::median(.data$DP, na.rm = TRUE),
      sd_DP = sd(.data$DP, na.rm = TRUE),
      .groups = "drop"
    )
  sequencing_depths
}


complete_missing_f_intervals <- function(tbl, intervals) {
  intervals <- intervals |>
    enframe(name = "sample_id", value = "f_interval") |>
    unnest("f_interval")
  missing_intervals <- intervals |>
    anti_join(tbl, by = c("sample_id", "f_interval"))
  tbl |>
    bind_rows(missing_intervals) |>
    arrange(.data$sample_id, .data$f_interval)
}


get_interval_centers <- function(intervals) {
  res <- intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    as_tibble_col("from_and_to") |>
    separate_wider_delim("from_and_to", names = c("from", "to"), delim = ",") |>
    map_df(parse_double) |>
    mutate(centers = .data$from + (.data$to - .data$from) / 2)
  res$centers
}


get_interval_width <- function(intervals) {
  res <- intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    as_tibble_col("from_and_to") |>
    separate_wider_delim("from_and_to", names = c("from", "to"), delim = ",") |>
    map_df(parse_double) |>
    mutate(width = .data$to - .data$from)
  stats::median(res$width)
}

