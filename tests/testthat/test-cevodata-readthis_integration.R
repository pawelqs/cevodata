test_that("adding ASCAT data works", {
  ascat_dir <- system.file("extdata", "ASCAT", package = "readthis")
  data <- readthis::read_ascat_files(ascat_dir, sample_id_pattern = "(?<=ASCAT\\/)[:alnum:]*(?=\\.)")
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(CNAs(cd), "tbl")
  expect_equal(default_CNAs(cd), "ASCAT")
  expect_equal(dim(CNAs(cd)), c(20, 8))
  expect_equal(get_metadata(cd)$purity, c(0.99322, 0.99322))
  expect_equal(get_metadata(cd)$purity, get_metadata(cd)$ascat_purity)
})



test_that("adding FACETS data works", {
  facets_dir <- system.file("extdata", "FACETS", package = "readthis")
  data <- readthis::read_facets_cnas(facets_dir)
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(CNAs(cd), "tbl")
  expect_equal(default_CNAs(cd), "FACETS")
  expect_equal(dim(CNAs(cd)), c(128, 16))
  expect_equal(get_metadata(cd)$purity, c(0.3, 0.3))
  expect_equal(get_metadata(cd)$purity, get_metadata(cd)$facets_purity)
})



test_that("adding Mutect2 data works", {
  path <- system.file("extdata", "Mutect", package = "readthis")
  data <- readthis::read_mutect_snvs(
    path,
    patient_id_pattern = "(?<=Mutect\\/)[:alnum:]*(?=\\.)",
    verbose = FALSE
  )
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(SNVs(cd), "tbl")
  expect_equal(default_SNVs(cd), "Mutect")
  expect_equal(dim(SNVs(cd)), c(16, 12))
  expect_equal(get_metadata(cd)$sample_id, c("S1_L1", "S1_P1", "S2_L1", "S2_P1"))
  expect_equal(get_metadata(cd)$patient_id, c("S1", "S1", "S2", "S2"))
})



test_that("adding Strelka data works", {
  path <- system.file("extdata", "Strelka", package = "readthis")
  data <- readthis::read_strelka_somatic_snvs(
    path,
    patient_id_pattern = "(?<=Strelka\\/)[:alnum:]*(?=\\.)",
    verbose = FALSE
  ) |>
    mutate(sample_id = str_c(patient_id, sample_id, sep = "_"))
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(SNVs(cd), "tbl")
  expect_equal(default_SNVs(cd), "Strelka")
  expect_equal(dim(SNVs(cd)), c(18, 9))
  expect_equal(get_metadata(cd)$sample_id, c("S1_TUMOR", "S2_TUMOR"))
})
