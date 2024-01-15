verbose::verbose(cevoverse = 0)
data("tcga_brca_test")


# ------------------------------ Calc ------------------------------------------


test_that("Calculation of Mf_1f works", {
  cd <- calc_Mf_1f(tcga_brca_test, verbose = FALSE)
  path <- test_path("testdata", "tcga_brca_Mf_1f.tsv")
  expected <- read_tsv(path, col_types = "ccdiid")
  class(expected) <- c("cevo_Mf_1f_tbl", class(expected))
  attr(expected, "f_column") <- "VAF"
  # write_tsv(get_Mf_1f(cd), path)
  expect_equal(get_Mf_1f(cd), expected)
})


# ------------------------------- Get ------------------------------------------


test_that("get_Mf_1f() runs calc_Mf_1f() if Mf_1f is missing", {
  object <- test_data
  Mf_1f <- get_Mf_1f(object)
  expect_s3_class(Mf_1f, "cevo_Mf_1f_tbl")
  expect_equal(attr(Mf_1f, "f_column"), "VAF")

  Mf_1f2 <- object |>
    calc_Mf_1f() |>
    get_Mf_1f()
  expect_equal(Mf_1f, Mf_1f2)
})


# ------------------------------ Plot ------------------------------------------


test_that("plot_Mf_1f() works", {
  p <- tcga_brca_test |>
    calc_Mf_1f() |>
    plot_Mf_1f()
  vdiffr::expect_doppelganger("plot_Mf_1f", p)
})


test_that("plot(calc_Mf_1f()) works", {
  p <- SNVs(tcga_brca_test) |>
    calc_Mf_1f(verbose = FALSE) |>
    plot()
  vdiffr::expect_doppelganger("plot(calc_Mf_1f())", p)
})
