make_gmt_file <- function(lines) {
  f <- tempfile(fileext = ".gmt")
  writeLines(lines, f)
  f
}

test_that("read_gene_sets: filters by size and merges multiple databases", {
  gmt1 <- make_gmt_file(c(
    "SET_SMALL\tdesc1\tA\tB", ## size 2, below thres.min
    "SET_OK\tdesc2\tA\tB\tC\tD\tE"  ## size 5
  ))
  gmt2 <- make_gmt_file(c(
    "SET_OTHER_DB\tdesc3\tX\tY\tZ\tW\tV" ## size 5
  ))

  res <- read_gene_sets(c(gmt1, gmt2), thres.min = 3, thres.max = 2000)

  expect_setequal(res$gs.names, c("SET_OK", "SET_OTHER_DB"))
  expect_equal(res$gs$SET_OK, c("A", "B", "C", "D", "E"))
  expect_null(res$gs.direction)
})

test_that("read_gene_sets: de-duplicates members within a gene set with a warning", {
  gmt <- make_gmt_file(c("SET_DUP\tdesc\tA\tA\tB\tC"))
  expect_warning(res <- read_gene_sets(gmt, thres.min = 1), "Duplicated gene set members")
  expect_equal(sort(res$gs$SET_DUP), c("A", "B", "C"))
  expect_equal(unname(res$size.G["SET_DUP"]), 3)
})

test_that("read_gene_sets: detects PTMsigDB-style directionality and strips suffixes", {
  gmt <- make_gmt_file(c("SIG1\tdesc\tP1;u\tP2;u\tP3;d"))
  res <- read_gene_sets(gmt, thres.min = 1)

  expect_equal(res$gs$SIG1, c("P1", "P2", "P3"))
  expect_equal(res$gs.direction$SIG1, c("u", "u", "d"))
})

test_that("match_gene_sets: drops gene sets below min.overlap and intersects members", {
  gmt <- make_gmt_file(c(
    "SET_A\tdesc\tG1\tG2\tG3\tG4\tG5",
    "SET_B\tdesc\tG1\tGxx\tGyy\tGzz\tGww" ## only 1 real overlap
  ))
  gsdb <- read_gene_sets(gmt, thres.min = 1)
  gene.names <- paste0("G", 1:10)

  matched <- match_gene_sets(gsdb, gene.names, min.overlap = 3)

  expect_equal(unname(matched$gs.names), "SET_A")
  expect_equal(sort(matched$gs$SET_A), c("G1", "G2", "G3", "G4", "G5"))
})

test_that("match_gene_sets: errors clearly when no gene set overlaps at all", {
  gmt <- make_gmt_file(c("SET_A\tdesc\tXX1\tXX2\tXX3"))
  gsdb <- read_gene_sets(gmt, thres.min = 1)
  expect_error(match_gene_sets(gsdb, paste0("G", 1:10), min.overlap = 1), "No overlap")
})

test_that("match_gene_sets: errors clearly when overlap exists but below min.overlap", {
  gmt <- make_gmt_file(c("SET_A\tdesc\tG1\tXX2\tXX3"))
  gsdb <- read_gene_sets(gmt, thres.min = 1)
  expect_error(match_gene_sets(gsdb, paste0("G", 1:10), min.overlap = 2), "minimum-overlap")
})

test_that("match_gene_sets: collapses duplicate gene set names, keeping the first", {
  gmt <- make_gmt_file(c(
    "SET_A\tdesc1\tG1\tG2\tG3",
    "SET_A\tdesc2\tG1\tG2\tG4" ## same name, different members
  ))
  gsdb <- read_gene_sets(gmt, thres.min = 1)
  matched <- match_gene_sets(gsdb, paste0("G", 1:10), min.overlap = 1)

  expect_equal(length(matched$gs.names), 1)
  expect_equal(sort(matched$gs$SET_A), c("G1", "G2", "G3")) ## first occurrence kept
})
