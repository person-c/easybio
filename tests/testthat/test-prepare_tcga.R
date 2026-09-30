# prepare_tcga() takes a SummarizedExperiment, which would be a heavy
# dependency for a test. The four slots it reads are what the stand-ins below
# reproduce: colData, rowRanges (a data.frame that also carries the names of
# its ranges), and assays@data.
setClass("fakeNames", representation(NAMES = "character"))
setClass("fakeRowRanges", contains = "data.frame", representation(ranges = "fakeNames"))
setClass("fakeAssays", representation(data = "list"))
setClass("fakeSE", representation(colData = "data.frame", rowRanges = "fakeRowRanges", assays = "fakeAssays"))

fake_se <- function() {
  dimnames <- list(c("g1", "g2"), c("s1", "s2", "s3"))
  new("fakeSE",
    colData = data.frame(
      sample_type = c("Primary Tumor", "Primary Tumor", "Solid Tissue Normal"),
      days_to_death = c(10, NA, NA),
      days_to_last_follow_up = c(100, 30, 50),
      row.names = dimnames[[2]]
    ),
    rowRanges = new("fakeRowRanges", data.frame(symbol = dimnames[[1]]),
      ranges = new("fakeNames", NAMES = dimnames[[1]])
    ),
    assays = new("fakeAssays", data = list(
      unstranded = matrix(1:6, 2, dimnames = dimnames),
      fpkm_unstrand = matrix(7:12, 2, dimnames = dimnames)
    ))
  )
}

# the warning is once per session, so the tests that expect it have to start
# from a clean slate rather than from whatever ran before them
reset_field_warning <- function() {
  # the namespace is locked, so empty the environment instead of replacing it
  warned <- get(".tcga_renamed_warned", envir = asNamespace("easybio"))
  rm(list = ls(warned), envir = warned)
}

test_that("prepare_tcga() returns snake_case tables", {
  lt <- prepare_tcga(fake_se())

  expect_identical(names(lt$all), c("expr_count", "features_info", "sample_info"))
  expect_identical(names(lt$tumor), c("expr_fpkm", "features_info", "sample_info"))
  expect_identical(rownames(lt$all$expr_count), c("g1", "g2"))
  expect_identical(colnames(lt$all$expr_count), c("s1", "s2", "s3"))
  # only the two tumor samples are kept in the tumor table
  expect_identical(colnames(lt$tumor$expr_fpkm), c("s1", "s2"))
  expect_identical(rownames(lt$all$sample_info), c("s1", "s2", "s3"))
  # OS falls back to the last follow-up when the death date is missing
  expect_identical(lt$all$sample_info[["OS"]], c(10, 30, 50))
})

test_that("the renamed expression fields still answer to their old names", {
  reset_field_warning()
  lt <- prepare_tcga(fake_se())

  expect_warning(old <- lt$all$exprCount, "renamed")
  expect_identical(old, lt$all$expr_count)
  expect_warning(old <- lt$tumor[["exprFpkm"]], "renamed")
  expect_identical(old, lt$tumor$expr_fpkm)

  # and each name warns only once per session, however often it is read
  expect_silent(lt$all$exprCount)
  expect_silent(lt$tumor$exprFpkm)
})

test_that("assigning an old field name writes the new one", {
  reset_field_warning()
  lt <- prepare_tcga(fake_se())

  expect_warning(lt$all$exprCount <- 1:6, "renamed")
  expect_identical(lt$all$expr_count, 1:6)
  expect_false("exprCount" %in% names(lt$all))

  expect_warning(lt$tumor[["exprFpkm"]] <- 7:10, "renamed")
  expect_identical(lt$tumor$expr_fpkm, 7:10)
  expect_false("exprFpkm" %in% names(lt$tumor))
})

test_that("a classed table behaves like the list it is", {
  lt <- prepare_tcga(fake_se())

  # the class only carries the accessors above: printing, positions and
  # missing names should look like a plain list
  expect_false(any(grepl("attr", capture.output(print(lt$all)))))
  expect_identical(lt$all[[1]], lt$all$expr_count)
  expect_error(lt$all[["nope"]], "subscript out of bounds")
  expect_null(lt$all$nope)
})
