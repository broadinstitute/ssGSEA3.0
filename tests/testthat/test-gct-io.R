test_that("read_gct_safe: flags 'deduped' when rids were made unique, not otherwise", {
  wd <- withr::local_tempdir()
  f.ok <- file.path(wd, "ok.gct")
  writeLines(c("#1.3", "2\t1\t0\t0", "id\tS1", "G1\t1", "G2\t2"), f.ok)
  ok <- read_gct_safe(f.ok)
  expect_null(attr(ok, "deduped"))

  f.dup <- file.path(wd, "dup.gct")
  writeLines(c("#1.3", "2\t1\t0\t0", "id\tS1", "G1\t1", "G1\t2"), f.dup)
  expect_warning(dup <- read_gct_safe(f.dup), "not unique")
  expect_true(isTRUE(attr(dup, "deduped")))
  expect_equal(dup@rid, c("G1", "G1_1"))
})

test_that("read_gct_safe: parses a single-data-row, multi-rdesc-column GCT (cmapR parse_gctx bug workaround)", {
  wd <- withr::local_tempdir()
  f <- file.path(wd, "single.gct")
  writeLines(c(
    "#1.3", "1\t2\t2\t0",
    "id\tdescA\tdescB\tS1\tS2",
    "GS1\tvalA\tvalB\t1.5\t2.5"
  ), f)
  g <- read_gct_safe(f)
  expect_equal(dim(g@mat), c(1, 2))
  expect_equal(unname(g@mat["GS1", ]), c(1.5, 2.5))
  expect_equal(g@rdesc$descA, "valA")
  expect_equal(g@rdesc$descB, "valB")
})
