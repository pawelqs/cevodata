data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnas <- CNAs(tcga_brca_test)
meta <- tcga_brca_test$metadata


test_that("init_cevodata works", {
  cd <- init_cevodata("TCGA BRCA small")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 0)
  expect_equal(length(cd$cNVs), 0)
  cd <- init_cevodata("TCGA BRCA small")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 0)
  expect_equal(length(cd$cNVs), 0)
})


test_that("init_cevodata with SNV assay works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs)
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 1)
  expect_s3_class(cd$SNVs[[1]], "tbl_df")
  expect_equal(cd$active_SNVs, "snvs")
})


test_that("Adding SNVs to cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = snvs, "head")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 2)
  expect_s3_class(cd$SNVs[[1]], "tbl_df")
  expect_s3_class(cd$SNVs[[2]], "tbl_df")
  expect_equal(cd$active_SNVs, "head")
})


test_that("Adding SNVs to cevodata extends metadata", {
  small_snvs <- snvs |>
    filter(sample_id %in% c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- init_cevodata("TCGA BRCA small", snvs = small_snvs)
  expect_equal(cd$metadata$sample_id, c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- add_SNV_data(cd, snvs = snvs, "head")
  expect_equal(nrow(cd$metadata), 4)
})


test_that("Adding incorrect SNV to cevodata throws error", {
  expect_error(init_cevodata("TCGA BRCA small", snvs = cnas))
})


test_that("Getting SNV from cevodata works", {
  tcga2 <- head(snvs)
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = tcga2, "head")
  expect_identical(SNVs(cd), tcga2)
  expect_equal(SNVs(cd) |> nrow(), 6)
  expect_identical(SNVs(cd, "head"), tcga2)
  expect_equal(SNVs(cd, "snvs") |> nrow(), 21570)
  expect_error(SNVs(cd, "xxx"))
})


test_that("Getting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = snvs, "head")
  expect_equal(default_SNVs(cd), "head")
})


test_that("Setting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = snvs, "head")
  default_SNVs(cd) <- "snvs"
  expect_equal(default_SNVs(cd), "snvs")
  expect_error(default_SNVs(cd) <- "xxx")
})


test_that("Adding CNAs to cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_CNA_data(cnas, "tcga")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$CNAs), 1)
  expect_s3_class(cd$CNAs[[1]], "tbl_df")
  expect_equal(cd$active_CNAs, "tcga")
})


test_that("Adding CNAs to cevodata extends metadata", {
  small_cnas <- cnas |>
    filter(sample_id %in% c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- init_cevodata("TCGA BRCA small", cnas = small_cnas)
  expect_equal(cd$metadata$sample_id, c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- add_CNA_data(cd, cnas = cnas, "head")
  expect_equal(nrow(cd$metadata), 4)
})


test_that("Adding incorrect CNA to cevodata throws error", {
  expect_error(init_cevodata("TCGA BRCA small", cnas = snvs))
})


test_that("Getting CNA from cevodata works", {
  cnas2 <- head(cnas)
  cd <- init_cevodata("TCGA BRCA small", cnas = cnas) |>
    add_CNA_data(cnas2, "tcga2")
  expect_identical(CNAs(cd), cnas2)
  expect_equal(CNAs(cd) |> nrow(), 6)
  expect_identical(CNAs(cd, "tcga2"), cnas2)
  expect_identical(CNAs(cd, "cnas"), cnas)
  expect_equal(CNAs(cd, "cnas") |> nrow(), 714)
  expect_error(CNAs(cd, "xxx"))
})


test_that("Getting active CNA assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnas = cnas) |>
    add_CNA_data(cnas = cnas, "tcga")
  expect_equal(default_CNAs(cd), "tcga")
})


test_that("Setting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnas = cnas) |>
    add_CNA_data(cnas = cnas, "tcga2")
  default_CNAs(cd) <- "cnas"
  expect_equal(default_CNAs(cd), "cnas")
  expect_error(default_CNAs(cd) <- "xxx")
})


test_that("Setting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnas = cnas) |>
    add_CNA_data(cnas = cnas, "tcga2")
  default_CNAs(cd) <- "cnas"
  expect_equal(default_CNAs(cd), "cnas")
  expect_error(default_CNAs(cd) <- "xxx")
})


test_that("Adding sample data to cevodata works", {
  meta <- tibble(
    sample_id =  c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"),
    meta1 = 1:2
  )
  cd <- init_cevodata("TCGA BRCA small") |>
    add_sample_data(meta)
  expect_equal(cd$metadata$sample_id, meta$sample_id)
  expect_equal(cd$metadata$meta1, meta$meta1)

  meta2 <- tibble(
    sample_id =  c("TCGA-D8-A27V-01", "TCGA-AN-A046-01"),
    meta2 = c("a", "b")
  )
  cd <- add_sample_data(cd, meta2)
  expect_equal(cd$metadata$sample_id, c("TCGA-AC-A23H-01", "TCGA-AN-A046-01", "TCGA-D8-A27V-01"))
  expect_equal(cd$metadata$meta1, c(1, 2, NA_real_))
  expect_equal(cd$metadata$meta2, c(NA_character_, "b", "a"))

  cd <- add_SNV_data(cd, snvs)
  expect_equal(nrow(cd$metadata), 4)
  expect_named(cd$metadata, c("sample_id", "meta1", "meta2"))
  expect_true(all(unique(snvs$sample_id) %in% cd$metadata$sample_id))
  expect_equal(cd$metadata$meta1, c(1, 2, NA_real_, NA_real_))
})
