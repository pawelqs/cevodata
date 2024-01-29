test_that("use_purity() works", {
  object <- test_data
  object$metadata$purity_2 <- 1
  res <- use_purity(object, "purity_2")
  expect_equal(get_purities(res)$purity, rep(1, 4))
})

test_that("use_purity() throws error if nonexisting purity column is requested", {
  object <- test_data
  expect_error(use_purity(object, "purity_2"))
})


test_that("get_patients_data() works", {
  sample_data <- tibble(
    sample_id = c("sampleA1", "sampleA2", "sampleB1", "sampleB2"),
    purity = c(0.5, 0.6, 0.5, 0.5)
  )
  samples_to_patients <- tibble(
    sample_id = c("sampleA1", "sampleA2", "sampleB1", "sampleB2"),
    patient_id = c("patientA", "patientA", "patientB", "patientB")
  )
  patients_data <- tibble(
    patient_id = c("patientA", "patientB"),
    sex = c("male", "female")
  )

  object <- init_cevodata() |>
    add_metadata(sample_data) |>
    add_metadata(samples_to_patients) |>
    add_metadata(patients_data)

  res <- get_patients_data(object)
  expect_equal(res, patients_data)
})
