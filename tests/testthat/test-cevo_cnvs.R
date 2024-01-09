verbose(cevo = 0)


# test_that("annotate_normal_cn() works", {
#   cnas <- tibble(
#     sample_id = c("A", "A", "B", "B", "C", "C"),
#     chrom = c("ch1", "chrX", "chr1", "chrY", "chr1", "chrX"),
#     start = 1,
#     end = 100
#   )
#   cd <- init_cevodata("test", cnas = cnas)
#   cd$metadata$sex <- c("M", "male", "female")
#
#   res <- cd |>
#     annotate_normal_cn() |>
#     CNAs()
#   expected <- c(2, 1, 2, 1, 2, 2)
#   expect_true("normal_cn" %in% names(res))
#   expect_equal(res$normal_cn, expected)
# })
