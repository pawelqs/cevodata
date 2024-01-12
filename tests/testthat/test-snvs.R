
test_that("as_cevo_snvs() works", {
  snvs <- tibble(sample_id = "A", mutation_id = str_c("chr", 1:3), VAF = 0.5)
  res <- as_cevo_snvs(snvs)
  expect_s3_class(res, c("cevo_snvs", "tbl_df", "tbl", "data.frame"))
})


test_that("SNVs_CNAs() works", {
  snvs <- tribble(
    ~sample_id, ~chrom, ~pos, ~VAF,
    "S1",       "chr1", 100,  0.1,
    "S1",       "chr1", 200,  0.1,
    "S1",       "chr1", 300,  0.1,
    "S1",       "chr1", 400,  0.1,
    "S1",       "chr2", 100,  0.1,
    "S2",       "chr1", 100,  0.1,
    "S2",       "chr1", 300,  0.1
  )
  snvs$ref <- "C"
  snvs$alt <- "T"
  cnas <- tribble(
    ~sample_id, ~chrom, ~start, ~end, ~total_cn, ~minor_cn, ~normal_cn,
    "S1",       "chr1", 1,      250,  3,         1,         2,
    "S1",       "chr1", 350,    450,  4,         1,         2,
    "S1",       "chr2", 1,      250,  2,         1,         2,
    "S2",       "chr1", 1,      250,  3,         1,         2
  )
  cd <- init_cevodata("test", snvs = snvs, cnas = cnas)
  expected <- tribble(
    ~sample_id, ~chrom, ~pos, ~VAF, ~total_cn, ~minor_cn, ~normal_cn,
    "S1",       "chr1", 100,  0.1,  3,         1,         2,
    "S1",       "chr1", 200,  0.1,  3,         1,         2,
    "S1",       "chr1", 300,  0.1,  NA_real_,  NA_real_,  NA_real_,
    "S1",       "chr1", 400,  0.1,  4,         1,         2,
    "S1",       "chr2", 100,  0.1,  2,         1,         2,
    "S2",       "chr1", 100,  0.1,  3,         1,         2,
    "S2",       "chr1", 300,  0.1,  NA_real_,  NA_real_,  NA_real_
  )
  res <- SNVs_CNAs(cd)
  expect_identical(res$total_cn, expected$total_cn)
  expect_identical(res$minor_cn, expected$minor_cn)
})


test_that("filter_SNVs_by_regions works", {
  snvs <- tibble(
    sample_id = "S1",
    chrom = str_c("chr", c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)),
    pos = c(100:104, 100:104)
  )
  regions <- tibble(
    chrom = c("chr1", "chr2"),
    start = c(103, 100),
    end = c(200, 102)
  )
  bed_file <- test_path("testdata", "regions.tsv")
  expected <- tibble(
    sample_id = "S1",
    chrom = str_c("chr", c(1, 1, 2, 2, 2)),
    pos = c(103L, 104L, 100L, 101L, 102L)
  )

  expect_identical(filter_SNVs_by_regions(snvs, regions = regions), expected)
  expect_identical(filter_SNVs_by_regions(snvs, bed_file = bed_file), expected)
})
