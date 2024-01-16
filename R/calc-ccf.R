
#' Calc mutation frequencies
#'
#' Calculates the CNV-corrected mutation frequencies. Implemented methods:
#' - Dentro - calculates the Cancer Cell Fraction (CCF) using the formulas from
#'   [Dentro et al. *Principles of Reconstructing the Subclonal Architecture of Cancers* (2015)](https://doi.org/10.1101/cshperspect.a026625)
#'
#' @param object <cevodata> object
#' @param method Available methods: Dentro
#' @param which_snvs Which SNVs to use
#' @param which_cnas Which CNAs to use
#' @param rm_intermediate_cols Should the columns used to get CCF be removed?
#' @param verbose Verbose?
#'
#' @return <cevodata> object
#' @export
calc_mutation_frequencies <- function(object,
                                      method = "Dentro",
                                      which_snvs = default_SNVs(object),
                                      which_cnas = default_CNAs(object),
                                      rm_intermediate_cols = TRUE,
                                      verbose = get_verbosity()) {
  if (method == "Dentro") {
    cnas <- CNAs(object, which_cnas) |>
      select("sample_id", "chrom", "start", "end", "total_cn", "normal_cn")
    purities <- get_purities(object)
    snvs <- SNVs(object, which_snvs) |>
      join_CNAs(cnas) |>
      left_join(purities, by = "sample_id") |>
      dentro_2015_correction() |>
      set_snvs_frequency_measure("CCF/2")
    intermediate_cols <- c("start", "end", "total_cn", "normal_cn", "purity", "u", "m")
  } else {
    stop("Currently supported methods: 'Dentro'")
  }

  if (rm_intermediate_cols) {
    snvs <- snvs |>
      select(-any_of(intermediate_cols))
  }

  nas <- snvs |>
    filter(is.na(.data$CCF)) |>
    nrow()
  nas_pct <- round(nas*100/nrow(snvs), digits = 2)
  msg(nas, " variants (", nas_pct, " %), have NA CCF value", verbose = verbose)

  object |>
    add_SNV_data(snvs, name = which_snvs)
}


#' @describeIn calc_mutation_frequencies Implements the CNV-based frequency
#' correction method described in [Dentro et al. 'Principles of Reconstructing the Subclonal Architecture of Cancers' (2015)](https://doi.org/10.1101/cshperspect.a026625)
#'
#' @param tbl tibble that contains columns: VAF, total_cn, normal_cn, purity
#'
#' @return The same tibble with new columns: u, m, CCF (Cancer Cell Fraction)
#' @export
dentro_2015_correction <- function(tbl) {
  tbl |>
    mutate(
      u = .data$VAF * (1 / .data$purity) * (.data$purity * .data$total_cn + (1 - .data$purity) * .data$normal_cn),
      m = if_else(.data$u < 1, 1, round(.data$u)),
      CCF = .data$u / .data$m,
      `CCF/2` = .data$CCF / 2
    ) |>
    relocate("CCF", "CCF/2", .after = "VAF")
}
