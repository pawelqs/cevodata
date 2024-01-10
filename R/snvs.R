
# --------------------------- cevo_snvs class --------------------------------

new_cevo_snvs <- function(tbl) {
  structure(tbl, class = c("cevo_snvs", class(tibble::tibble())))
}


#' Create cevo_snvs tibble
#' @param snvs tibble
#' @export
as_cevo_snvs <- function(snvs) {
  validate_SNVs(snvs)
  columns_to_order <- c(
    "sample_id", "chrom", "pos", "gene_symbol",
    "ref", "alt", "ref_reads", "alt_reads", "impact", "VAF"
  )
  snvs |>
    select(any_of(columns_to_order), everything()) |>
    new_cevo_snvs()
}


validate_SNVs <- function(snvs) {
  # snvs is not empty
  if (length(snvs) == 0) {
    stop("snvs cannot be empty tibble")
  }

  # snvs contains sample_id and required mutation columns
  required_cols <- list(
    option1 = c("sample_id", "chrom", "pos", "alt", "ref"),
    option2 = c("sample_id", "mutation_id")
  )
  pass <- map_lgl(required_cols, ~all(.x %in% names(snvs)))
  if (!any(pass)) {
    stop(
      "snvs must contain either: \n",
      required_cols |>
        map_chr(str_c, collapse = ", ") |>
        str_c(collapse = "\nor\n")
    )
  }

  # snvs contains VAF or read counts
  required_cols_2 <- list(
    option1 = c("alt_reads", "ref_reads"),
    option2 = c("VAF")
  )
  pass <- map_lgl(required_cols_2, ~all(.x %in% names(snvs)))
  if (!any(pass)) {
    stop(
      "snvs must contain either: \n",
      required_cols_2 |>
        map_chr(str_c, collapse = ", ") |>
        str_c(collapse = "\nor\n")
    )
  }
}


# ------------------------------- More getters ---------------------------------

#' @describeIn assays Get SNVs in the wide table form
#' @param object cevodata object
#' @param fill_na fill missing with this value
#' @export
get_SNVs_wider <- function(object, fill_na = NULL) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

  snvs <- SNVs(object) |>
    select("sample_id", "chrom":"alt", "VAF") |>
    left_join(patients_to_samples, by = "sample_id") |>
    unite_mutation_id() |>
    select(-"sample_id") |>
    select("patient_id", everything()) |>
    pivot_wider(names_from = "sample", values_from = "VAF")

  if (!is.null(fill_na)) {
    snvs[is.na(snvs)] <- fill_na
  }
  snvs
}


#' @describeIn assays Get SNVs matrix (for single-patient cevodata only)
#'
#' Requires 'sample' column in metadata.
#'
#' @param object cevodata object
#' @param rows_sample Sample to use as rows
#' @param cols_sample Sample to use as cols
#' @param bins Number of bins to use for VAF intervals
#' @param verbose Verbose?
#' @export
get_SNVs_2d_matrix <- function(object,
                               rows_sample = NULL, cols_sample = NULL,
                               bins = NULL,
                               verbose = verbose::verbose("cevoverse")) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

  if (is.null(rows_sample) || is.null(cols_sample)) {
    rows_sample <- patients_to_samples$sample[[1]]
    cols_sample <- patients_to_samples$sample[[2]]
    msg("Using '", rows_sample, "' as rows and '", cols_sample, "' as cols", verbose = verbose)
  }
  if (n_distinct(patients_to_samples$patient_id) > 1) {
    stop("This function works only for single sample objects")
  }
  if (rows_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for rows is not present for this patient")
  }
  if (cols_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for cols is not present for this patient")
  }

  breaks <- object |>
    SNVs() |>
    get_interval_breaks(bins = bins)
  sample_ids <- patients_to_samples$sample_id |>
    set_names(patients_to_samples$sample)
  rowsample_breaks <- breaks[[sample_ids[rows_sample]]]
  colsample_breaks <- breaks[[sample_ids[cols_sample]]]

  # Prepare SNVs wider tibble
  mutations <-get_SNVs_wider(object, fill_na = 0)
  mutations <- mutations[c("mutation_id", rows_sample, cols_sample)]
  colnames(mutations) <- c("mutation_id", "rows_sample", "cols_sample")

  # Cut intervals
  mutations <- mutations |>
    mutate(
      rows_sample = cut(rows_sample, breaks = rowsample_breaks),
      cols_sample = cut(cols_sample, breaks = colsample_breaks)
    )
  row_intervals <- levels(mutations$rows_sample)
  col_intervals <- levels(mutations$cols_sample)

  # Prepare square matrix
  incomplete_mat <- mutations |>
    mutate(across(c("rows_sample", "cols_sample"), as.character)) |>
    group_by(.data$rows_sample, .data$cols_sample) |>
    count() |>
    ungroup() |>
    pivot_wider(names_from = "cols_sample", values_from = "n") |>
    column_to_rownames("rows_sample")
  incomplete_mat[is.na(incomplete_mat)] <- 0

  mat <- matrix(0, nrow = length(row_intervals), ncol = length(col_intervals))
  rownames(mat) <- row_intervals
  colnames(mat) <- col_intervals
  mat[rownames(incomplete_mat), colnames(incomplete_mat)] <- as.matrix(incomplete_mat)

  attr(mat, "rows_sample") <- rows_sample
  attr(mat, "cols_sample") <- cols_sample
  attr(mat, "rows_sample_id") <- sample_ids[[rows_sample]]
  attr(mat, "cols_sample_id") <- sample_ids[[cols_sample]]
  mat
}


get_SNVs_wider_intervals <- function(object, fill_na = NULL, bins = NULL) {
  metadata <- object$metadata
  if (n_distinct(metadata$patient_id) > 1) {
    stop("This function works only for single sample objects")
  }

  breaks <- object |>
    SNVs() |>
    get_interval_breaks(bins = bins)
  breaks <- breaks[metadata$sample_id]
  names(breaks) <- metadata$sample

  mutations <-get_SNVs_wider(object, fill_na = 0)
  sample_intervals <- list()
  for (sample in metadata$sample) {
    VAF_intervals <- cut(mutations[[sample]], breaks = breaks[[sample]])
    sample_intervals[[sample]] <- levels(VAF_intervals)
    mutations[[sample]] <- as.character(VAF_intervals)
  }

  attr(mutations, "sample_intervals") <- sample_intervals
  mutations
}


#' Get SNVs with merged CNAs
#' @param object cevodata object with SNVs and CNAs
#' @export
SNVs_CNAs <- function(object) {
  SNVs(object) |>
    join_CNAs(CNAs(object))
}


join_CNAs <- function(snvs, cnas) {
  left_join(
    snvs, cnas,
    by = join_by("sample_id", "chrom", "pos" >= "start", "pos" <= "end"),
    relationship = "many-to-one"
  )
}


# ---------------------------- Misc functions --------------------------------

#' Unite many columns to create mutation_id column
#' @param snvs SNVs
#' @param sep Separator
#' @param remove Remove united columns?
#' @export
unite_mutation_id <- function(snvs, sep = "-", include_gene_symbol = FALSE, remove = TRUE) {
  cols <- c("chrom", "pos", if (include_gene_symbol) "gene_symbol" else NULL, "ref", "alt")
  unite(snvs, "mutation_id", all_of(cols), sep = sep, remove = remove)
}


#' Filter SNVs by position: using regions tbl or bed file
#' @param snvs snvs tbl with columns: sample_id, chrom, pos
#' @param regions regions tbl with columns chrom, start, end
#' @param bed_file bed file
#' @export
filter_SNVs_by_regions <- function(snvs, regions = NULL, bed_file = NULL) {
  if (is.null(regions) && is.null(bed_file)) {
    stop("Provide one of: regions, bed_file")
  }

  if (!is.null(bed_file)) {
    # bed_file = "tests/testdata/regions.tsv" # for tests only
    regions <- bed_file |>
      read_tsv(col_types = "cii", col_names = c("chrom", "start", "end"))

    # Unlike the coordinate system used by other standards such as GFF, the system
    # used by the BED format is zero-based for the coordinate start and one-based
    # for the coordinate end.
    regions <- regions |>
      mutate(start = .data$start + 1)
  }
  snv_classes <- class(snvs)

  regions_gr <- regions |>
    rename(seqnames = "chrom") |>
    plyranges::as_granges()
  snvs_gr <- snvs |>
    mutate(
      seqnames = .data$chrom,
      start = .data$pos,
      end = .data$pos,
      .before = "sample_id"
    ) |>
    plyranges::as_granges()

  filtered_snvs <- plyranges::filter_by_overlaps(snvs_gr, regions_gr) |>
    as_tibble() |>
    select(-("seqnames":"strand"))
  class(filtered_snvs) <- snv_classes

  filtered_snvs
}
