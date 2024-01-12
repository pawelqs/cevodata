data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnas <- CNAs(tcga_brca_test)
meta <- tcga_brca_test$metadata


test_that("print.cevodata runs without error works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_metadata(meta)
  output <- print(cd) |>
    capture.output()
  expect_true(length(cd) > 1)
})


test_that("print.cevodata runs without error foe empty object", {
  cd <- init_cevodata("TCGA BRCA small")
  output <- print(cd) |>
    capture.output()
  expect_true(length(cd) > 1)
})
