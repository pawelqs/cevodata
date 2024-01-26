test_that("use_purity() works", {
  object <- test_data
  object$metadata$purity_2 <- 1
  res <- use_purity(object, "purity_2")
  expect_equal(get_purities(res)$purity, rep(1, 4))
})

test_that("use_purity() throws error if nonexisting purity column is requested", {
  object <- test_data
  expect_error(use_purity(object, "purity_2"))
})
