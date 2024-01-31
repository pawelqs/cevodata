
# cevodata

<!-- badges: start -->
[![R-CMD-check](https://github.com/pawelqs/cevodata/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pawelqs/cevodata/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`cevodata` R package provides a `cevodata` class - a data structure designed to store the mutational data (SNVs, CNAs), metadata, and the analytical results (e.g. mathematical models) on a cohort of tumor samples. It also facilitates the data manipulation (e.g. filtering or splitting the data object) and visualization. The package also provides a set of functions to calculate the variant's Cancer Cell Fractions (CCF), and a number of sample-level statistics, such as the VAF/CCF spectra, cumulative tails, or [Williams's](https://doi.org/10.1038/ng.3489) $M(f) \sim 1/f$ statistics.


## Installation

You can install the development version of cevodata from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pawelqs/cevodata")
```


## Example

This is a basic example which shows how to create a `cevodata` object and to plot the VAF spectra:

``` r
library(cevodata)

cd <- init_cevodata() |> 
  add_snv_data(snv_tbl) |>
  add_cna_data(cna_tbl) |>
  calc_mutation_frequencies()

plot_SFS(cd)
```

For more examples go to the vignettes.

## Last changes

* **v1.0.0** - cevodata package extracted from the original cevomod R package


## Help and support

[GitHub Issues](https://github.com/pawelqs/cevomod/issues)


