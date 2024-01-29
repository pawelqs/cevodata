# -------------------------------- Sample Stats --------------------------------


test_that("calc_all_sample_stats() works", {
  snvs <- tibble(
    sample_id = c("sample1", "sample1", "sample2", "sample2", "sample3"),
    mutation_id = str_c("mutation", 0:4),
    VAF = 0.5
  )

  res <- init_cevodata(snvs = snvs) |>
    calc_all_sample_stats() |>
    get_stats("sample_stats")

  expect_equal(nrow(res), 3)
  expect_true("mutation_burden" %in% colnames(res))
})


test_that("calc_sample_mutation_burden() works", {
  snvs <- tibble(
    sample_id = c("sample1", "sample1", "sample2", "sample2", "sample3"),
    mutation_id = str_c("mutation", 0:4),
    VAF = 0.5
  )

  res <- init_cevodata(snvs = snvs) |>
    calc_sample_mutation_burden()
  expected <- tibble(
    sample_id = c("sample1", "sample2", "sample3"),
    mutation_burden = c(2, 2, 1)
  )

  expect_equal(res, expected)
})


# -------------------------------- Patient Stats --------------------------------


test_that("calc_all_patient_stats() works", {
  snvs <- tibble(
    sample_id = c("sampleA1", "sampleA1", "sampleA2", "sampleA2", "sampleB1"),
    mutation_id = str_c("mutation", 0:4),
    VAF = 0.5
  )
  meta <- tibble(
    sample_id = c("sampleA1", "sampleA2", "sampleB1"),
    patient_id = c("patientA", "patientA", "patientB")
  )

  object <- init_cevodata(snvs = snvs) |>
    add_metadata(meta)
  snvs <- default_SNVs(object)

  res <- object |>
    calc_all_patient_stats() |>
    get_stats("patient_stats")

  expect_equal(nrow(res), 2)
  expect_true("mutation_burden" %in% colnames(res))
})


test_that("calc_patient_mutation_burden() works", {
  snvs <- tibble(
    sample_id = c("sampleA1", "sampleA1", "sampleA2", "sampleA2", "sampleB1"),
    mutation_id = str_c("mutation", 0:4),
    VAF = 0.5
  )
  meta <- tibble(
    sample_id = c("sampleA1", "sampleA2", "sampleB1"),
    patient_id = c("patientA", "patientA", "patientB")
  )

  res <- init_cevodata(snvs = snvs) |>
    add_metadata(meta) |>
    calc_patient_mutation_burden()
  expected <- tibble(
    patient_id = c("patientA", "patientB"),
    mutation_burden = c(4, 1)
  )

  expect_equal(res, expected)
})


# -------------------------- Sample-sample stats -------------------------------

test_that("calc_sample_sample_stats() works", {
  snvs <- bind_rows(
    tibble(sample_id = "sampleA1", mutation_id = str_c("m", 1:4), VAF = 0.5),
    tibble(sample_id = "sampleA2", mutation_id = str_c("m", 3:4), VAF = 0.5),
    tibble(sample_id = "sampleA3", mutation_id = str_c("m", 1:4), VAF = 0.5),
    tibble(sample_id = "sampleB1", mutation_id = str_c("m", 1:4), VAF = 0.5)
  )
  meta <- tibble(
    patient_id = c("patientA", "patientA", "patientA", "patientB"),
    sample_id = c("sampleA1", "sampleA2", "sampleA3", "sampleB1"),
    sample = c("primary", "relapse", "metastasis", "primary")
  )

  object <- init_cevodata(snvs = snvs) |>
    add_metadata(meta)
  res <- object |>
    calc_sample_sample_stats() |>
    get_stats("sample_sample_stats")

  expected <- tibble(
    sample1 = c("primary", "primary", "relapse"),
    sample2 = c("relapse", "metastasis", "metastasis"),
    Jaccard_index = c(0.5, 1, 0.5)
  )

  expect_equal(nrow(res), 1)
  expect_equal(res$patient_id, "patientA")
  expect_equal(res$Jaccard_indexes[[1]], expected)
})


test_that("calc_Jaccard_indexes() works with paired samples", {
  snvs <- bind_rows(
    tibble(sample_id = "sampleA1", mutation_id = str_c("m", 1:4), VAF = 0.5),
    tibble(sample_id = "sampleA2", mutation_id = str_c("m", 3:4), VAF = 0.5),
    tibble(sample_id = "sampleB1", mutation_id = str_c("m", 1:4), VAF = 0.5),
    tibble(sample_id = "sampleB2", mutation_id = str_c("m", 1:4), VAF = 0.5)
  )
  meta <- tibble(
    patient_id = c("patientA", "patientA", "patientB", "patientB"),
    sample_id = c("sampleA1", "sampleA2", "sampleB1", "sampleB2"),
    sample = c("primary", "relapse", "primary", "relapse")
  )

  res <- init_cevodata(snvs = snvs) |>
    add_metadata(meta) |>
    calc_Jaccard_indexes()
  expected <- tibble(
    patient_id = c("patientA", "patientB"),
    sample1 = c("primary", "primary"),
    sample2 = c("relapse", "relapse"),
    Jaccard_index = c(0.5, 1)
  )

  expect_equal(nrow(res), 2)
})


test_that("calc_Jaccard_indexes() works with unpaired samples", {
  snvs <- bind_rows(
    tibble(sample_id = "sampleA1", mutation_id = str_c("m", 1:4), VAF = 0.5),
    tibble(sample_id = "sampleA2", mutation_id = str_c("m", 3:4), VAF = 0.5),
    tibble(sample_id = "sampleA3", mutation_id = str_c("m", 1:4), VAF = 0.5),
    tibble(sample_id = "sampleB1", mutation_id = str_c("m", 1:4), VAF = 0.5)
  )
  meta <- tibble(
    patient_id = c("patientA", "patientA", "patientA", "patientB"),
    sample_id = c("sampleA1", "sampleA2", "sampleA3", "sampleB1"),
    sample = c("primary", "relapse", "metastasis", "primary")
  )

  res <- init_cevodata(snvs = snvs) |>
    add_metadata(meta) |>
    calc_Jaccard_indexes()
  expected <- tibble(
    sample1 = c("primary", "primary", "relapse"),
    sample2 = c("relapse", "metastasis", "metastasis"),
    Jaccard_index = c(0.5, 1, 0.5)
  )

  expect_equal(nrow(res), 1)
  expect_equal(res$patient_id, "patientA")
  expect_equal(res$Jaccard_indexes[[1]], expected)
})
