test_that("package scaffold loads and bundled data is intact", {
  expect_true(is.character(hallmark_categories))
  expect_true(is.character(hallmark_category_colors))
  expect_true(all(unique(hallmark_categories) %in% names(hallmark_category_colors)))
})
