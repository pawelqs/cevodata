cd <- test_data
new_cnas <- CNAs(cd) |>
  filter(sample_id != "Sample 3")
cd <- cd |>
  add_CNA_data(new_cnas, "new_cnas")
res <- export_to_clip(cd)


test_that("export_to_clip() creates correct tibbles for each sample", {
  expect_named(res, str_c("Sample ", 1:4))
  res |>
    map(\(x) names(x) == c("snvs", "cnas", "purities")) |>
    map_lgl(all) |>
    all() |>
    expect_true()
})


test_that("export_to_clip() creates output files with correct column names", {
  expect_named(res$`Sample 1`$snvs, c("chromosome_index", "position", "alt_count", "ref_count"))
  expect_named(
    res$`Sample 1`$cnas,
    c("chromosome_index", "start_position", "end_position", "major_cn", "minor_cn", "total_cn")
  )
  expect_type(res$`Sample 1`$purities, "double")
})
