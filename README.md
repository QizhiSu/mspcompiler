
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mspcompiler

<!-- badges: start -->
<!-- badges: end -->

The goal of mspcompiler is to offer ways to compile either EI or tandem
mass spectral libraries from various sources, such as NIST (if you have
it installed), MoNA, and GPNS, and organize them into a neat and
up-to-date msp file that can be used in MS-DIAL.

## Installation

Install by trying the following code.

``` r
# install.packages("devtools")
devtools::install_github("QizhiSu/mspcompiler", build_vignettes = TRUE)
```

If and fail, please try:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("ChemmineR")

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("ChemmineOB")
```

Please read the vignettes to learn how to compile EI and tandem mass
spectral libraries.

``` r
browseVignettes("mspcompiler")
```
