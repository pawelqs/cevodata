## code to prepare `sysdata` goes here
library(tidyverse)


driver_genes <- openxlsx::read.xlsx(
    "data-raw/driver_genes.xlsx",
    sheet = "Table S1", startRow = 4
  ) |>
  as_tibble() |>
  rename(
    Type = `Tumor.suppressor.or.oncogene.prediction.(by.20/20+)`
  )


usethis::use_data(
  driver_genes,
  overwrite = TRUE, internal = TRUE
)
usethis::use_data(driver_genes, overwrite = TRUE)
