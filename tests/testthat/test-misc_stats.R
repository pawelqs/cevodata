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

  res <- init_cevodata(snvs = snvs) |>
    add_metadata(meta) |>
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


