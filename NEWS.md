# Future Developments

1. It appears that the GEOquery package may not be necessary for preparing GEO series data. We might consider removing this dependency.
2. All `plot*` functions should be enhanced to allow for greater user customization.
3. Development of a new S3 class similar to dgeList, which would offer more customization and improved performance.

# Recently Developed Features

enhancements:

1. Now, the `Artist`  class has a default data argument which can be used to explore the data in a new way.
2. The `Artist` class now will save all your command and result, you can recheck those results anytime you want.

`prepare_geo` function: 

1. Converts non-character ID columns in gpl data to character to prevent issues with IDs being read as numeric (e.g., 1, 2, 3, 4…).
2. Reads potential expression data from supplementary files…

`plotORA` function:

1. Fix bug when you map a variable to fill aesthetic parameters, the plot doesn't show legend.

`get_marker` function:

1. add extra messages when you search  undefined cells in database.
  
# easybio 1.0.1 Release Notes

1. Added a NEWS.md file to document changes and updates to the package.
2. Fixed an error in the `uniprot_id_map()` function that occurred on the latest macOS in example usage.
3. Updated vignettes for optimization.