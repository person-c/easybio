test_that("available_tissue_type() lists the types recorded under one class", {
  all_types <- available_tissue_type("Human")
  blood <- available_tissue_type("Human", tissue_class = "Blood")

  expect_true(length(blood) > 0)
  expect_true(all(blood %in% all_types))
  expect_true("Peripheral blood" %in% blood)
  # a class is one of two labels on an entry, not its parent, so restricting to
  # one class drops the types recorded under the others
  expect_lt(length(blood), length(all_types))
  # the classes are not disjoint, so a class can add types another one has
  expect_gt(
    length(available_tissue_type("Human", tissue_class = c("Blood", "Bone marrow"))),
    length(blood)
  )
})

test_that("available_tissue_type() defaults to every class of the species", {
  expect_setequal(
    available_tissue_type("Human"),
    available_tissue_type("Human", tissue_class = available_tissue_class("Human"))
  )
})

test_that("available_tissue_type() rejects a class the species does not have", {
  expect_error(available_tissue_type("Human", tissue_class = "NoSuchTissue"), "tissue_class")
  # the labels are matched as they are stored, so a stray space is a typo, not
  # an empty result
  expect_error(available_tissue_type("Human", tissue_class = "Blood "), "tissue_class")
})
