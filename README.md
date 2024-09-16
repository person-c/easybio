
<!-- README.md is generated from README.Rmd. Please edit that file -->

# easybio

<!-- badges: start -->

[![R-CMD-check](https://github.com/person-c/easybio/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/person-c/easybio/actions/workflows/check-standard.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/easybio)](https://CRAN.R-project.org/package=easybio)
<!-- badges: end -->

`easybio` provides a comprehensive toolkit for single-cell RNA-seq
annotation using the CellMarker2.0 database. It streamlines the process
of assigning biological labels in scRNA-seq data, integrating seamlessly
with tools like Seurat. While the package includes additional
bioinformatics workflows, such as handling TCGA and GEO datasets,
differential expression analysis, and enrichment analysis visualization,
for details specifically on the single-cell annotation functionality,
please refer to the [bioRxiv
preprint](https://doi.org/10.1101/2024.09.14.609619).

## Download and Usage

You can install the development version of `easybio` from
[GitHub](https://github.com/) with:

``` r
devtools::install("person-c/easybio", build_vignettes = TRUE)
```

To know how to use this package, please see the
[wiki](https://github.com/person-c/easybio/wiki) or run:

``` r
vignette(topic = "example-bulk-rna-seq-workflow", package = 'easybio')
vignette(topic = "example-single-cell-annotation", package = "easybio")
```

To learn the difference between development version and CRAN version,
see [NEWS](./NEWS.md)

## Citation

If you use the single-cell annotation functionality from `easybio`,
consider cite:

Wei, Cui. (2024). *easybio: an R Package for Single-Cell Annotation with
CellMarker2.0*. bioRxiv. <https://doi.org/10.1101/2024.09.14.609619>

C. Hu, T. Li, Y. Xu, X. Zhang, F. Li, J. Bai, J. Chen, W. Jiang, K.
Yang, Q. Ou, X. Li, P. Wang, Y. Zhang, CellMarker 2.0: an updated
database of manually curated cell markers in human/mouse and web tools
based on scRNA-seq data, Nucleic Acids Research 51 (D1) (2022)
D870–D876. <doi:10.1093/nar/gkac947>
