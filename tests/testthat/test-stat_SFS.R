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


# test_that("Calculation of SFS with resampling works", {
#   object <- test_data
#   sfs_resamples <- calc_SFS_resamples(object, times = 2)
#
#   expect_type(sfs_resamples, "list")
#   expect_equal(length(sfs_resamples), 4)
#
#   sfs_resamples |>
#     walk(expect_s3_class, "cevo_SFS_bootstraps")
#
#   attribs <- attributes(sfs_resamples$`Sample 1`$sfs[[1]])
#   expect_equal(attribs$f_column, "VAF")
# })
#
#
# test_that("Calculation of SFS with resampling works with CCF/2", {
#   object <- test_data |>
#     calc_mutation_frequencies()
#   sfs_resamples <- calc_SFS_resamples(object, times = 2)
#
#   expect_type(sfs_resamples, "list")
#   expect_equal(length(sfs_resamples), 4)
#
#   sfs_resamples |>
#     walk(expect_s3_class, "cevo_SFS_bootstraps")
#
#   attribs <- attributes(sfs_resamples$`Sample 1`$sfs[[1]])
#   expect_equal(attribs$f_column, "CCF/2")
# })


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


# test_that("get_non_zero_SFS_range works", {
#   SFS <- tibble(
#     sample_id = "S1",
#     f = 1:100,
#     y = 10
#   ) |>
#     mutate(
#       y = case_when(
#         f < 12 ~ 0,
#         f == 40 ~ 0,
#         f %in% c(54, 55) ~ 1,
#         f > 65 ~ 0,
#         TRUE ~ y
#       ),
#       f = f / 100
#     )
#   allowed_zero_bins <- 1
#   y_treshold <- 1
#
#   expected <- tibble(
#     sample_id = "S1",
#     from = 0.12,
#     to = 0.65
#   )
#   res <- get_non_zero_SFS_range(SFS)
#
#   expect_identical(res, expected)
# })
