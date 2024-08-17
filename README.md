
<!-- README.md is generated from README.Rmd. Please edit that file -->

# easybio

<!-- badges: start -->

[![R-CMD-check](https://github.com/person-c/easybio/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/person-c/easybio/actions/workflows/check-standard.yaml)
<!-- badges: end -->

You can install the development version of easybio from
[GitHub](https://github.com/) with:

``` r
# install.packages("BiocManger")
BiocManger::install("person-c/easybio")
```

To know how to use this package, please see the
[wiki](https://github.com/person-c/easybio/wiki).

## Coding Principles Reminders

- Use `[[` and `slot` to access child nodes instead of `$` or `@`.
- Ensure all operations are written in scripts, not executed directly in
  the terminal.
- Save all raw data used for your visualizations.
- Detailed script name to describle its target and number it by order.
- Use `lintr`.

An example project looks like this:

``` r
r$> fs::dir_tree('example-project/')
example-project/
+-- .gitignore
+-- main.R
+-- patch1-GSEA4MF.R
+-- Rawdata
+-- README.qmd
\-- Result
    +-- main
    +-- patch1
    \-- README
```
