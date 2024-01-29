data("test_data")


# ---------------------------- Calc SFS ----------------------------------------

test_that("Calculation of SFS works", {
  object <- test_data
  cd <- calc_SFS(object)
  expected <- test_path("testdata", "SFS.tsv") |>
    read_tsv(col_types = "ccdid")
  class(expected) <- c("cevo_SFS_tbl", class(expected))
  attr(expected, "f_column") <- "VAF"
  # write_tsv(get_SFS(cd), test_path("testdata", "SFS.tsv"))
  expect_equal(get_SFS(cd), expected)
})


# ---------------------------- Get SFS ----------------------------------------


test_that("get_SFS() runs calc_SFS() if SFS is missing", {
  object <- test_data
  sfs <- get_SFS(object)
  expect_s3_class(sfs, "cevo_SFS_tbl")
  expect_equal(attr(sfs, "f_column"), "VAF")

  sfs2 <- object |>
    calc_SFS() |>
    get_SFS()
  expect_equal(sfs, sfs2)
})


# ---------------------------- Plot SFS ----------------------------------------

test_that("plot.cevo_snvs() works", {
  dt <- SNVs(test_data)
  expect_warning({
    p <- dt %>%
      calc_SFS() %>%
      plot()
  })
  expect_s3_class(p, c("gg", "ggplot"))
  vdiffr::expect_doppelganger("test_data_sfs_1", p)
})


test_that("plot_SFS() works", {
  expect_warning({
    p <- test_data |>
      calc_SFS() |>
      plot_SFS()
  })
  vdiffr::expect_doppelganger("test_data_sfs_2", p)
})
