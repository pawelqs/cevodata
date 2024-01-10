
test_that("as_cevo_snvs() works", {
  snvs <- tibble(sample_id = "A", mutation_id = str_c("chr", 1:3), VAF = 0.5)
  res <- as_cevo_snvs(snvs)
  expect_s3_class(res, c("cevo_snvs", "tbl_df", "tbl", "data.frame"))
})
