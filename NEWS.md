# Future

* It seems that there is no need to use `GEOquery` package to prepare GEO series data...maybe I can remove this dependency.
* The `Artist` should have a default field `data` to show a data with different method.
* All `plot*` function should have more user-customized arguments.
* A new s3 class like `dgeList`, but will more customized and faster.

# Developed

* `prepare_geo`: 
  
1. convert non-character ID column in gpl data to character to avoid the ID column is like numeric 1, 2, 3, 4...
2. read potential expression data from supplementary files...

# easybio 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Avoid error of function `uniprot_id_map()` in latest mac in examples.
* update vignettes.