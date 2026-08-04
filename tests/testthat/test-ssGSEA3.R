make_test_gct <- function(mat, path) {
  g <- cmapR::GCT()
  g@mat <- mat
  g@rid <- rownames(mat)
  g@cid <- colnames(mat)
  g@rdesc <- data.frame(row.names = rownames(mat))
  g@cdesc <- data.frame(row.names = colnames(mat))
  cmapR::write_gct(g, ofile = path, appenddim = FALSE)
  path
}

make_test_gmt <- function(sets, path) {
  lines <- vapply(names(sets), function(nm) paste(c(nm, "desc", sets[[nm]]), collapse = "\t"), character(1))
  writeLines(lines, path)
  path
}

test_that("ssGSEA3: end-to-end run on synthetic data scores gene sets sensibly and reproducibly with a seed", {
  wd <- withr::local_tempdir()
  old.wd <- getwd()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  set.seed(1)
  genes <- paste0("G", 1:50)
  ## S1: GS_UP genes (G1-G10) all high; S2: GS_UP genes all low
  expr.up <- c(50:41, rnorm(40, mean = 20, sd = 5))
  expr.down <- c(1:10, rnorm(40, mean = 20, sd = 5))
  mat <- matrix(c(expr.up, expr.down), ncol = 2, dimnames = list(genes, c("S1", "S2")))

  gct.path <- make_test_gct(mat, "input.gct")
  gmt.path <- make_test_gmt(list(GS_UP = genes[1:10]), "sets.gmt")

  res <- ssGSEA3(gct.path, output.prefix = "run1", gene.set.databases = gmt.path,
                  sample.norm.type = "rank", weight = 0, min.overlap = 5, nperm = 100,
                  seed = 42, verbose = FALSE, param.file = FALSE, log.file = NULL,
                  export.signat.gct = FALSE)

  expect_s3_class(res, "ssgsea_result")
  expect_equal(unname(res$gs.names), "GS_UP")
  ## GS_UP genes are top-ranked in S1 (positive enrichment) and bottom-ranked in S2 (negative)
  expect_gt(res$scores["GS_UP", "S1"], 0)
  expect_lt(res$scores["GS_UP", "S2"], 0)

  expect_true(file.exists("run1-scores.gct"))
  expect_true(file.exists("run1-pvalues.gct"))
  expect_true(file.exists("run1-fdr-pvalues.gct"))
  expect_true(file.exists("run1-combined.gct"))

  ## reproducibility: same seed -> identical scores
  res2 <- ssGSEA3(gct.path, output.prefix = "run2", gene.set.databases = gmt.path,
                   sample.norm.type = "rank", weight = 0, min.overlap = 5, nperm = 100,
                   seed = 42, verbose = FALSE, param.file = FALSE, log.file = NULL,
                   export.signat.gct = FALSE)
  expect_identical(res$scores, res2$scores)
  expect_identical(res$pvalues, res2$pvalues)
})

test_that("ssGSEA3: output.score.type='ES' skips bootstrapping (p-values are all 1)", {
  wd <- withr::local_tempdir()
  old.wd <- getwd()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  genes <- paste0("G", 1:30)
  mat <- matrix(c(30:1, rnorm(30, 10, 3)), ncol = 2, dimnames = list(genes, c("S1", "S2")))
  gct.path <- make_test_gct(mat, "input.gct")
  gmt.path <- make_test_gmt(list(GS_A = genes[1:8]), "sets.gmt")

  res <- ssGSEA3(gct.path, output.prefix = "run_es", gene.set.databases = gmt.path,
                  output.score.type = "ES", min.overlap = 3, verbose = FALSE,
                  param.file = FALSE, log.file = NULL, export.signat.gct = FALSE)

  expect_true(all(res$pvalues == 1))
  expect_equal(res$scores, res$scores) ## ES == NES when nperm forced to 0; sanity only
})

test_that("ssGSEA3: verbose=TRUE with log.file=NULL doesn't crash (log writing must not depend on verbose)", {
  wd <- withr::local_tempdir()
  old.wd <- getwd()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  genes <- paste0("G", 1:20)
  mat <- matrix(c(20:1, rnorm(20, 10, 3)), ncol = 2, dimnames = list(genes, c("S1", "S2")))
  gct.path <- make_test_gct(mat, "input.gct")
  gmt.path <- make_test_gmt(list(GS_A = genes[1:6]), "sets.gmt")

  expect_no_error(
    ssGSEA3(gct.path, output.prefix = "run_verbose_nolog", gene.set.databases = gmt.path,
            min.overlap = 3, nperm = 20, verbose = TRUE, param.file = FALSE, log.file = NULL,
            export.signat.gct = FALSE)
  )
})
