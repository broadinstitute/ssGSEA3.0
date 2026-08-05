## ##########################################################################
##            .align_geneset_to_direction() (the desync fix)
## ##########################################################################

test_that(".align_geneset_to_direction: drops absent members without desyncing direction labels", {
  gene.names <- paste0("G", 1:10)
  gene.set <- c("G1", "NOTPRESENT1", "G2", "NOTPRESENT2", "G9", "G10")
  gs.dir <- c("u", "u", "u", "u", "d", "d")

  aligned <- .align_geneset_to_direction(gene.set, gs.dir, gene.names)

  expect_equal(aligned$gene.set2, c(1, 2, 9, 10))
  expect_equal(aligned$gene.set.direction, c("u", "u", "d", "d"))
  expect_equal(length(aligned$gene.set2), length(aligned$gene.set.direction))
})

test_that(".align_geneset_to_direction: expands direction labels for redundant (multi-probe) gene names", {
  gene.names <- c("G1", "G1", "G2") ## G1 measured twice (e.g. multi-probe data)
  gene.set <- c("G1", "G2")
  gs.dir <- c("u", "d")

  aligned <- .align_geneset_to_direction(gene.set, gs.dir, gene.names)

  expect_equal(aligned$gene.set2, c(1, 2, 3))
  expect_equal(aligned$gene.set.direction, c("u", "u", "d")) ## G1's "u" label repeated for both its rows
})

test_that(".align_geneset_to_direction: non-directional (NULL) case is an unchanged passthrough", {
  gene.names <- paste0("G", 1:5)
  gene.set <- c("G1", "NOTPRESENT", "G3")
  aligned <- .align_geneset_to_direction(gene.set, NULL, gene.names)
  expect_equal(aligned$gene.set2, c(1, 3))
  expect_null(aligned$gene.set.direction)
})

test_that("project_geneset: a directional gene set with a partially-absent member no longer crashes", {
  gene.names <- paste0("G", 1:20)
  data.expr <- 20:1
  mat <- matrix(data.expr, ncol = 1, dimnames = list(gene.names, "S1"))

  ## 3 "u" members but only 2 exist in the data; 2 "d" members, both exist
  gene.set <- c("G1", "NOTPRESENT", "G2", "G19", "G20")
  gs.dir <- c("u", "u", "u", "d", "d")

  res <- project_geneset(data.array = mat, gene.names = gene.names, n.cols = 1, n.rows = 20,
                          weight = 0, statistic = "area.under.RES", gene.set = gene.set,
                          nperm = 0, correl.type = "rank", gene.set.direction = gs.dir,
                          min.overlap = 3, size.G.current = 5)

  expect_false(is.na(res$ES.vector[1]))
  expect_equal(res$ES.vector[1], 20) ## verified independently: ES.u=10, ES.d=-10, ES=ES.u-ES.d=20
})

test_that("project_geneset: bootstrap resampling of a directional gene set is crash-free and reproducible", {
  set.seed(1)
  gene.names <- paste0("G", 1:30)
  mat <- matrix(rnorm(30 * 2), ncol = 2, dimnames = list(gene.names, c("S1", "S2")))
  gene.set <- gene.names[1:10]
  gs.dir <- rep(c("u", "d"), 5)

  set.seed(42)
  res1 <- project_geneset(data.array = mat, gene.names = gene.names, n.cols = 2, n.rows = 30,
                           weight = 0, statistic = "area.under.RES", gene.set = gene.set,
                           nperm = 300, correl.type = "rank", gene.set.direction = gs.dir,
                           min.overlap = 5, size.G.current = 10)
  set.seed(42)
  res2 <- project_geneset(data.array = mat, gene.names = gene.names, n.cols = 2, n.rows = 30,
                           weight = 0, statistic = "area.under.RES", gene.set = gene.set,
                           nperm = 300, correl.type = "rank", gene.set.direction = gs.dir,
                           min.overlap = 5, size.G.current = 10)

  expect_false(anyNA(res1$NES.vector))
  expect_identical(res1$NES.vector, res2$NES.vector)
  expect_identical(res1$p.val.vector, res2$p.val.vector)
})

## ##########################################################################
##                          gsea_score()
## ##########################################################################

test_that("gsea_score (non-directional): a gene set concentrated at the top ranks scores as hand-derived", {
  gene.names <- paste0("G", 1:10)
  data.expr <- 10:1 ## G1 highest, G10 lowest -> already rank-ordered
  gene.list <- order(data.expr, decreasing = TRUE)
  gene.set2 <- c(1, 2, 3) ## top 3 genes

  r_auc <- gsea_score(gene.list, gene.set2, weight = 0, n.rows = 10, correl.type = "rank",
                       gene.set.direction = NULL, data.expr = data.expr,
                       statistic = "area.under.RES", min.overlap = 1)
  r_ks <- gsea_score(gene.list, gene.set2, weight = 0, n.rows = 10, correl.type = "rank",
                      gene.set.direction = NULL, data.expr = data.expr,
                      statistic = "Kolmogorov-Smirnov", min.overlap = 1)

  expect_equal(r_auc$ES, 5)
  expect_equal(r_ks$ES, 1)
})

test_that("gsea_score (non-directional): a gene set at the bottom ranks scores negatively", {
  gene.names <- paste0("G", 1:10)
  data.expr <- 10:1
  gene.list <- order(data.expr, decreasing = TRUE)
  gene.set2 <- c(8, 9, 10) ## bottom 3 genes

  r_auc <- gsea_score(gene.list, gene.set2, weight = 0, n.rows = 10, correl.type = "rank",
                       gene.set.direction = NULL, data.expr = data.expr,
                       statistic = "area.under.RES", min.overlap = 1)

  expect_lt(r_auc$ES, 0)
  expect_equal(r_auc$ES, -5) ## symmetric to the top-3 case by construction
})

test_that("gsea_score (directional): up-at-top + down-at-bottom combine as hand-derived", {
  gene.names <- paste0("G", 1:20)
  data.expr <- 20:1
  gene.list <- order(data.expr, decreasing = TRUE)
  gene.set2 <- c(1, 2, 3, 4, 5, 16, 17, 18, 19, 20)
  gs.dir <- c("u", "u", "u", "u", "u", "d", "d", "d", "d", "d")

  d_auc <- gsea_score(gene.list, gene.set2, weight = 0, n.rows = 20, correl.type = "rank",
                       gene.set.direction = gs.dir, data.expr = data.expr,
                       statistic = "area.under.RES", min.overlap = 1)
  d_ks <- gsea_score(gene.list, gene.set2, weight = 0, n.rows = 20, correl.type = "rank",
                      gene.set.direction = gs.dir, data.expr = data.expr,
                      statistic = "Kolmogorov-Smirnov", min.overlap = 1)

  expect_equal(d_auc$ES.all$u, 10)
  expect_equal(d_auc$ES.all$d, -10)
  expect_equal(d_auc$ES, 20) ## ES.u - ES.d = 10 - (-10)

  expect_equal(d_ks$ES.all$u, 1)
  expect_equal(d_ks$ES.all$d, -1)
  expect_equal(d_ks$ES, 2)
})

test_that("gsea_score (directional): mismatched direction (up-at-bottom, down-at-top) scores negatively", {
  gene.names <- paste0("G", 1:20)
  data.expr <- 20:1
  gene.list <- order(data.expr, decreasing = TRUE)
  gene.set2 <- c(1, 2, 3, 4, 5, 16, 17, 18, 19, 20)
  gs.dir <- c("d", "d", "d", "d", "d", "u", "u", "u", "u", "u") ## swapped

  d_auc <- gsea_score(gene.list, gene.set2, weight = 0, n.rows = 20, correl.type = "rank",
                       gene.set.direction = gs.dir, data.expr = data.expr,
                       statistic = "area.under.RES", min.overlap = 1)

  expect_equal(d_auc$ES, -20)
})

## ##########################################################################
##                        project_geneset()
## ##########################################################################

test_that("project_geneset (nperm=0): matches gsea_score's ES exactly and sets NES=ES, p=1", {
  mat <- matrix(10:1, ncol = 1, dimnames = list(paste0("G", 1:10), "S1"))
  res <- project_geneset(data.array = mat, gene.names = rownames(mat), n.cols = 1, n.rows = 10,
                          weight = 0, statistic = "area.under.RES", gene.set = c("G1", "G2", "G3"),
                          nperm = 0, correl.type = "rank", gene.set.direction = NULL,
                          min.overlap = 1, size.G.current = 3)

  expect_equal(res$ES.vector[1], 5)
  expect_equal(res$NES.vector[1], 5)
  expect_equal(res$p.val.vector[1], 1)
  expect_equal(res$OL.numb.vector[1], 3)
  expect_equal(res$OL.perc.vector[1], 100)
})

test_that("project_geneset: missing values are excluded from ranking, not left as NA scores", {
  expr <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
  expr[4] <- NA ## drop G4
  mat <- matrix(expr, ncol = 1, dimnames = list(paste0("G", 1:10), "S1"))

  res <- project_geneset(data.array = mat, gene.names = rownames(mat), n.cols = 1, n.rows = 10,
                          weight = 0, statistic = "area.under.RES", gene.set = c("G1", "G2", "G3"),
                          nperm = 0, correl.type = "rank", gene.set.direction = NULL,
                          min.overlap = 1, size.G.current = 3)

  expect_false(is.na(res$ES.vector[1]))
  expect_equal(res$OL.numb.vector[1], 3) ## none of G1-G3 were the dropped gene
})

test_that("project_geneset: insufficient overlap reports NA scores and overlap diagnostics, not an error", {
  mat <- matrix(10:1, ncol = 1, dimnames = list(paste0("G", 1:10), "S1"))
  res <- project_geneset(data.array = mat, gene.names = rownames(mat), n.cols = 1, n.rows = 10,
                          weight = 0, statistic = "area.under.RES", gene.set = c("G1", "Gxx", "Gyy"),
                          nperm = 0, correl.type = "rank", gene.set.direction = NULL,
                          min.overlap = 2, size.G.current = 3)

  expect_true(is.na(res$ES.vector[1]))
  expect_true(is.na(res$random.walk[[1]]))
  expect_equal(res$OL.numb.vector[1], 1) ## only G1 overlaps
  expect_true(grepl("G1", res$OL.vector[1]))
})

test_that("project_geneset: directional gene set with an unusable up/down split (lone member on one side, other side < min.overlap) is NA, not a fake zero score", {
  gene.names <- paste0("G", 1:20)
  mat <- matrix(20:1, ncol = 1, dimnames = list(gene.names, "S1"))

  ## 1 "u" member + 4 "d" members = 5 total (clears min.overlap on raw count),
  ## but gsea_score()'s lone-member rule zeroes both sides since the "d" side
  ## (4) is itself below min.overlap (5) -- this used to report ES=0/p=1 as
  ## if a real (null) score had been computed, instead of NA/unscored
  gene.set <- c("G1", "G10", "G11", "G12", "G13")
  gs.dir <- c("u", "d", "d", "d", "d")

  res <- project_geneset(data.array = mat, gene.names = gene.names, n.cols = 1, n.rows = 20,
                          weight = 0, statistic = "area.under.RES", gene.set = gene.set,
                          nperm = 0, correl.type = "rank", gene.set.direction = gs.dir,
                          min.overlap = 5, size.G.current = 5)

  expect_true(is.na(res$ES.vector[1]))
  expect_true(is.na(res$NES.vector[1]))
  expect_true(is.na(res$p.val.vector[1]))
  expect_true(is.na(res$random.walk[[1]]))
  expect_equal(res$OL.numb.vector[1], 5) ## raw overlap still reported for diagnostics
  expect_true(grepl("G1;u", res$OL.vector[1]))
})

test_that("project_geneset: a directional gene set with a *usable* up/down split still scores normally (no regression)", {
  gene.names <- paste0("G", 1:20)
  mat <- matrix(20:1, ncol = 1, dimnames = list(gene.names, "S1"))

  ## same total size (5) and same 1-member "u" side as the unusable case above,
  ## but here the "d" side (4) is reduced to fit within a smaller min.overlap,
  ## so it clears the lone-member rule and should score as usual
  gene.set <- c("G1", "G10", "G11", "G12", "G13")
  gs.dir <- c("u", "d", "d", "d", "d")

  res <- project_geneset(data.array = mat, gene.names = gene.names, n.cols = 1, n.rows = 20,
                          weight = 0, statistic = "area.under.RES", gene.set = gene.set,
                          nperm = 0, correl.type = "rank", gene.set.direction = gs.dir,
                          min.overlap = 4, size.G.current = 5)

  expect_false(is.na(res$ES.vector[1]))
  expect_true(is.list(res$random.walk[[1]])) ## a real GSEA.results list, not the NA sentinel
})

test_that("project_geneset (bootstrap): reproducible under a fixed seed, NES same sign as ES, p in [0,1]", {
  set.seed(42)
  mat <- matrix(c(rnorm(15, mean = 3), rnorm(15, mean = -1)), ncol = 2,
                dimnames = list(paste0("G", 1:15), c("S1", "S2")))
  gene.set <- paste0("G", 1:5)

  set.seed(123)
  res1 <- project_geneset(data.array = mat, gene.names = rownames(mat), n.cols = 2, n.rows = 15,
                           weight = 0, statistic = "area.under.RES", gene.set = gene.set,
                           nperm = 50, correl.type = "rank", gene.set.direction = NULL,
                           min.overlap = 3, size.G.current = 5)
  set.seed(123)
  res2 <- project_geneset(data.array = mat, gene.names = rownames(mat), n.cols = 2, n.rows = 15,
                           weight = 0, statistic = "area.under.RES", gene.set = gene.set,
                           nperm = 50, correl.type = "rank", gene.set.direction = NULL,
                           min.overlap = 3, size.G.current = 5)

  expect_identical(res1$NES.vector, res2$NES.vector)
  expect_identical(res1$p.val.vector, res2$p.val.vector)
  expect_true(all(res1$p.val.vector >= 0 & res1$p.val.vector <= 1))
  expect_equal(sign(res1$NES.vector), sign(res1$ES.vector))
})
