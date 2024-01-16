data("tcga_brca_test")

cd <- tcga_brca_test |>
  intervalize_mutation_frequencies(column = "VAF")


# ---------------------------------- Calc -------------------------------------


test_that("Calculation of cumulative tails works", {
  cd <- calc_cumulative_tails(cd)
  path <- test_path("testdata", "tcga_brca_cumulative_tails.tsv")
  expected <- read_tsv(path, col_types = "ccdiid")
  # write_tsv(cd$models$cumulative_tails, path)
  class(expected) <- c("cevo_cumulative_tails_tbl", class(expected))
  attr(expected, "f_column") <- "VAF"
  expect_equal(get_cumulative_tails(cd), expected)
})


# ---------------------------------- Get -------------------------------------


test_that("get_SFS() runs calc_SFS() if SFS is missing", {
  object <- test_data
  object$stats[["cumulative_tails"]] <- NULL
  tails <- get_cumulative_tails(object)
  expect_s3_class(tails, "cevo_cumulative_tails_tbl")
  expect_equal(attr(tails, "f_column"), "VAF")

  tails2 <- object |>
    calc_cumulative_tails() |>
    get_cumulative_tails()
  expect_equal(tails, tails2)
})


# ---------------------------------- Plot -------------------------------------


test_that("Plotting of cumulative tails works", {
  p <- cd |>
    calc_cumulative_tails() |>
    plot_cumulative_tails()
  expect_s3_class(p, c("gg", "ggplot"))

  p2 <- SNVs(cd) %>%
    calc_cumulative_tails() %>%
    plot()
  expect_s3_class(p2, c("gg", "ggplot"))
})
