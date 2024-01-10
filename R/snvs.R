
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
#' @param include_gene_symbol Include gene symbol in mutation_id?
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
