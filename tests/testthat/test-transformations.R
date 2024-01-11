verbose(cevoverse = 0)

data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnas <- CNAs(tcga_brca_test)
meta <- get_metadata(tcga_brca_test) |>
  mutate(
    patient_id = c("A", "A", "B", "B"),
    sample = c("S1", "S2", "S1", "S2")
  )

sample_ids <- unique(meta$sample_id)
test_models <- list(
  coefs = tibble(sample_id = sample_ids, coef = 1:4),
  residuals = tibble(sample_id = sample_ids, resid = 11:14),
  info = list(source_stat = "xxx")
)
simple_stats <- tibble(
  sample_id = sample_ids,
  a = 1:4
)


# ------------------------- Filter --------------------------------------------

test_that("Filtering cevodata works", {
  patients <- unique(meta$patient_id)
  cd <- init_cevodata("TCGA BRCA small", cnas = cnas) |>
    add_CNA_data(cnas = cnas, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga") |>
    add_metadata(meta) |>
    add_stats(simple_stats, name = "simple_stats") |>
    add_models(test_models, name = "simple_models")

  res <- filter(cd, sample_id == "TCGA-AC-A23H-01")

  # SNVs
  expect_equal(nrow(res$SNVs$tcga), 6419)
  expect_equal(unique(res$SNVs$tcga$sample_id), "TCGA-AC-A23H-01")

  # CNAs
  expect_equal(nrow(res$CNAs$cnas), 285)
  expect_equal(nrow(res$CNAs$tcga2), 285)
  expect_equal(unique(res$CNAs$cnas$sample_id), "TCGA-AC-A23H-01")
  expect_equal(unique(res$CNAs$tcga2$sample_id), "TCGA-AC-A23H-01")

  # stats
  expect_equal(nrow(res$stats$simple_stats), 1)

  # Models
  res_models <- get_models(res)
  expect_equal(nrow(res_models$coefs), 1)
  expect_equal(nrow(res_models$residuals), 1)
  expect_equal(res_models$info, list(source_stat = "xxx"))
})


# ------------------------------- Split by -------------------------------------

test_that("Splitting cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnas = cnas) |>
    add_CNA_data(cnas = cnas, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga") |>
    add_stats(simple_stats, name = "simple_stats") |>
    add_models(test_models, name = "simple_models")
  cd$metadata$sex <- c("M", "M", "M", "F")
  splits <- split_by(cd, sex, verbose = FALSE)
  expect_named(splits, c("M", "F"))
  expect_equal(get_metadata(splits$M)$sex, c("M", "M", "M"))
  expect_equal(get_metadata(splits$F)$sex, c("F"))

  # SNVs are splitted correctly
  expect_equal(nrow(splits$M$SNVs$tcga), 17338)
  expect_equal(nrow(splits$`F`$SNVs$tcga), 4232)

  # CNAs are splitted correctly
  expect_equal(nrow(splits$M$CNAs$cnas), 545)
  expect_equal(nrow(splits$M$CNAs$tcga2), 545)
  expect_equal(nrow(splits$`F`$CNAs$tcga2), 169)
  expect_equal(nrow(splits$`F`$CNAs$tcga2), 169)

  # Stats are splitted correctly
  expect_equal(nrow(splits$M$stats$simple_stat), 3)
  expect_equal(nrow(splits$`F`$stats$simple_stat), 1)

  # Models are splitted correctly
  res_models <- map(splits, get_models)
  expect_equal(nrow(res_models$M$coefs), 3)
  expect_equal(nrow(res_models$`F`$coefs), 1)
  expect_equal(nrow(res_models$M$residuals), 3)
  expect_equal(nrow(res_models$`F`$residuals), 1)
})


# ------------------------------- Merge -------------------------------------

test_that("Merging cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnas = cnas) |>
    add_CNA_data(cnas = cnas, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga") |>
    add_models(test_models, name = "simple_models")
  cd1 <- filter(cd, sample_id == "TCGA-AC-A23H-01")
  cd2 <- filter(cd, sample_id != "TCGA-AC-A23H-01")
  cd_merged <- merge(cd1, cd2, name = "TCGA BRCA small")
  output <- print(cd_merged) |>
    capture.output()
  expect_equal(object.size(cd_merged), object.size(cd))
})

