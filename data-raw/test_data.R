
## --------------------------------- Data ------------------------------------

data("test_data", package = "cevodatasets")
test_data <- update_cevodata_v2_to_v3(test_data)

sample_data <- test_data$metadata |>
  mutate(purity = c(1, 0.7, 1, 1))

snvs <- SNVs(test_data)
snvs <- snvs |>
  mutate(
    chrom = if_else(VAF > 0.75, "chr2", "chr1"),
    alt_reads = if_else(sample_id == "Sample 2", floor(alt_reads * 0.7), alt_reads),
    ref_reads = if_else(sample_id == "Sample 2", DP - alt_reads, ref_reads),
    VAF = if_else(sample_id == "Sample 2", alt_reads / (ref_reads + alt_reads), VAF)
  )

cnas <- tribble(
  ~sample_id,  ~chrom, ~start, ~end,  ~total_cn, ~major_cn, ~minor_cn, ~normal_cn,
   "Sample 1",  "chr1", 1,      4000,  2,         1,         0,         2,
   "Sample 1",  "chr2", 1,      4000,  1,         1,         0,         2,
   "Sample 2",  "chr1", 1,      4000,  2,         1,         0,         2,
   "Sample 2",  "chr2", 1,      4000,  1,         1,         0,         2,
   "Sample 3",  "chr1", 1,      4000,  2,         1,         0,         2,
   "Sample 3",  "chr2", 1,      4000,  1,         1,         0,         2,
   "Sample 4",  "chr1", 1,      4000,  2,         1,         0,         2,
   "Sample 4",  "chr2", 1,      4000,  1,         1,         0,         2
)



test_data <- init_cevodata(name = "test_data") |>
  add_SNV_data(snvs) |>
  add_CNA_data(cnas) |>
  add_metadata(sample_data)

use_data(test_data, overwrite = TRUE)
