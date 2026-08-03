test_that("normalize_samples: rank scales to [0, 10000] per column, preserving NAs", {
  m <- matrix(c(5, 1, 3, NA, 10, 20, 15, 5), ncol = 2)
  out <- normalize_samples(m, "rank")

  ## column 1: values 5,1,3,NA -> ranks 3,1,2,NA -> *10000/4
  expect_equal(out[, 1], c(3, 1, 2, NA) * 10000 / 4)
  expect_true(is.na(out[4, 1]))
  expect_false(anyNA(out[, 2]))
})

test_that("normalize_samples: log.rank applies log(rank-scaled + e) on top of rank", {
  m <- matrix(c(5, 1, 3), ncol = 1)
  out <- normalize_samples(m, "log.rank")
  ranks <- c(3, 1, 2)
  expect_equal(out[, 1], log(10000 * ranks / 3 + exp(1)))
})

test_that("normalize_samples: log clamps values below 1 before log-transforming", {
  m <- matrix(c(0.5, 2, -3), ncol = 1)
  out <- normalize_samples(m, "log")
  expect_equal(out[, 1], log(c(1, 2, 1) + exp(1)))
})

test_that("normalize_samples: none leaves the matrix unchanged", {
  m <- matrix(c(5, 1, 3, 8), ncol = 2)
  expect_equal(normalize_samples(m, "none"), m)
})
