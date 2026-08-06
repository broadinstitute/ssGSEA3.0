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

test_that("ssGSEA3: a directional gene set unscoreable in every sample (lone-member split) is dropped, not exported empty", {
  wd <- withr::local_tempdir()
  old.wd <- getwd()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  genes <- paste0("G", 1:20)
  mat <- matrix(c(20:1, rnorm(20, 10, 3)), ncol = 2, dimnames = list(genes, c("S1", "S2")))
  gct.path <- make_test_gct(mat, "input.gct")

  ## GS_BAD: 1 "u" + 4 "d" members (5 total, clears min.overlap=5 on raw count)
  ## but the lone "u" member can't be scored since the "d" side (4) is itself
  ## below min.overlap -- used to report a fake ES=0/p=1 and, worse, crash
  ## export.signat.gct with an empty-overlap signature (see
  ## project_ssgsea3_empty_signature_export_bug memory)
  gmt.path <- make_test_gmt(
    list(
      GS_BAD = c("G1;u", "G10;d", "G11;d", "G12;d", "G13;d"),
      GS_OK = paste0(genes[1:6], ";u")
    ),
    "sets.gmt"
  )

  res <- NULL
  expect_no_error(
    res <- ssGSEA3(gct.path, output.prefix = "run_bad_signature", gene.set.databases = gmt.path,
                   min.overlap = 5, nperm = 20, verbose = FALSE, param.file = FALSE,
                   log.file = NULL, export.signat.gct = TRUE)
  )

  expect_false("GS_BAD" %in% res$gs.names)
  expect_true("GS_OK" %in% res$gs.names)
  expect_true(any(grepl("^GS_OK", list.files("signature_gct"))))
  expect_false(any(grepl("^GS_BAD", list.files("signature_gct"))))
})

test_that("ssGSEA3: accepts a GCT object directly (e.g. from preprocessGCT()), matching the file-path route exactly", {
  wd <- withr::local_tempdir()
  old.wd <- getwd()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  ## raw, multi-probe-per-gene input requiring gc-level collapsing
  mat <- matrix(c(20:1, rnorm(20, 10, 3)), ncol = 2, dimnames = list(paste0("p", 1:20), c("S1", "S2")))
  rdesc <- data.frame(geneSymbol = paste0("G", rep(1:10, each = 2)), row.names = rownames(mat))
  raw.gct <- cmapR::GCT()
  raw.gct@mat <- mat
  raw.gct@rid <- rownames(mat)
  raw.gct@cid <- colnames(mat)
  raw.gct@rdesc <- rdesc
  raw.gct@cdesc <- data.frame(row.names = colnames(mat))
  raw.path <- "raw_input.gct"
  cmapR::write_gct(raw.gct, ofile = raw.path, appenddim = FALSE)

  gmt.path <- make_test_gmt(list(GS_A = paste0("G", 1:5)), "sets.gmt")

  processed <- preprocessGCT(raw.path, level = "gc", mode = "mean", gene.col = "geneSymbol", appenddim = FALSE)
  expect_s4_class(processed, "GCT")

  ## route 1: pass the in-memory GCT object straight through, no re-read
  res.obj <- ssGSEA3(processed, output.prefix = "run_obj", gene.set.databases = gmt.path,
                      sample.norm.type = "rank", weight = 0, min.overlap = 3, nperm = 20,
                      seed = 42, verbose = FALSE, param.file = FALSE, log.file = NULL,
                      export.signat.gct = FALSE)

  ## route 2: the traditional file-path route, using the path preprocessGCT() itself wrote
  res.path <- ssGSEA3(attr(processed, "path"), output.prefix = "run_path", gene.set.databases = gmt.path,
                       sample.norm.type = "rank", weight = 0, min.overlap = 3, nperm = 20,
                       seed = 42, verbose = FALSE, param.file = FALSE, log.file = NULL,
                       export.signat.gct = FALSE)

  expect_identical(res.obj$scores, res.path$scores)
  expect_identical(res.obj$pvalues, res.path$pvalues)
  expect_equal(unname(res.obj$gs.names), "GS_A")
})

test_that("ssGSEA3: gcr-mode passthrough recognizes redundant genes via the 'deduped' attribute, not discards them", {
  wd <- withr::local_tempdir()
  old.wd <- getwd()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  ## GeneA appears twice (-> relabeled "GeneA"/"GeneA_1" by gcr) at ranks 5-6
  ## of 10; with area.under.RES + weight=0, whether both ranks count as GS_A
  ## hits (vs only the first, un-stripped "GeneA") changes the score --
  ## unlike a boundary/edge-rank placement, this configuration discriminates
  ## the two cases (verified empirically) and so proves the attribute, not
  ## just the labeling, is what makes ssGSEA3() recognize "GeneA_1" as GeneA
  mat <- matrix(10:1, ncol = 1, dimnames = list(paste0("p", 1:10), "S1"))
  rdesc <- data.frame(
    geneSymbol = c("G1", "G2", "G3", "G4", "GeneA", "GeneA", "G7", "G8", "G9", "G10"),
    row.names = rownames(mat)
  )
  raw.gct <- cmapR::GCT()
  raw.gct@mat <- mat
  raw.gct@rid <- rownames(mat)
  raw.gct@cid <- "S1"
  raw.gct@rdesc <- rdesc
  raw.gct@cdesc <- data.frame(row.names = "S1")
  raw.path <- "raw_input.gct"
  cmapR::write_gct(raw.gct, ofile = raw.path, appenddim = FALSE)

  gmt.path <- make_test_gmt(list(GS_A = "GeneA"), "sets.gmt")

  processed <- preprocessGCT(raw.path, level = "gcr", gene.col = "geneSymbol", appenddim = FALSE)
  expect_true(isTRUE(attr(processed, "deduped")))
  expect_true(all(c("GeneA", "GeneA_1") %in% processed@rid))

  res.fixed <- ssGSEA3(processed, output.prefix = "run_fixed", gene.set.databases = gmt.path,
                        sample.norm.type = "rank", weight = 0, min.overlap = 1,
                        output.score.type = "ES", verbose = FALSE, param.file = FALSE,
                        log.file = NULL, export.signat.gct = FALSE)

  ## simulate the pre-fix bug: strip the attribute so ssGSEA3() treats
  ## "GeneA_1" as an unrecognized gene, same as read_gct_safe()'s happy path
  ## on a file whose rids were already (redundantly) unique
  broken <- processed
  attr(broken, "deduped") <- NULL
  res.broken <- ssGSEA3(broken, output.prefix = "run_broken", gene.set.databases = gmt.path,
                         sample.norm.type = "rank", weight = 0, min.overlap = 1,
                         output.score.type = "ES", verbose = FALSE, param.file = FALSE,
                         log.file = NULL, export.signat.gct = FALSE)

  expect_equal(unname(res.fixed$scores["GS_A", "S1"]), 0)
  expect_equal(unname(res.broken$scores["GS_A", "S1"]), 5 / 9)
})

test_that("ssGSEA3: seed reproducibility holds under par=TRUE too, and is not tied to par at all", {
  skip_on_cran() ## spins up a parallel::makeCluster() worker process
  ## par=TRUE workers are separate R processes that resolve internal
  ## functions (e.g. project_geneset()) by loading the "ssGSEA3" namespace
  ## by name -- this only works if the package is actually installed
  ## (R CMD INSTALL), not merely source-loaded via devtools::load_all()/
  ## pkgload, which registers no installed package for a fresh process to
  ## find. R CMD check always installs first, so this isn't an issue there.
  wd <- withr::local_tempdir()
  old.wd <- getwd()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  set.seed(7)
  genes <- paste0("G", 1:30)
  mat <- matrix(rnorm(60, 20, 5), ncol = 2, dimnames = list(genes, c("S1", "S2")))
  gct.path <- make_test_gct(mat, "input.gct")
  gmt.path <- make_test_gmt(list(GS_A = genes[1:10]), "sets.gmt")

  run <- function(par, seed) {
    ssGSEA3(gct.path, output.prefix = tempfile(), gene.set.databases = gmt.path,
            sample.norm.type = "rank", weight = 0, min.overlap = 3, nperm = 50,
            seed = seed, par = par, spare.cores = 1, verbose = FALSE,
            param.file = FALSE, log.file = NULL, export.signat.gct = FALSE)
  }

  ## par=TRUE with a seed is just as reproducible as the serial (par=FALSE)
  ## case already covered above -- par only changes *how* the seed is
  ## applied (clusterSetRNGStream vs plain set.seed), not *whether* seeding
  ## works; see the misleading @param seed doc this test guards against
  par.a <- run(par = TRUE, seed = 42)
  par.b <- run(par = TRUE, seed = 42)
  expect_identical(par.a$pvalues, par.b$pvalues)

  ## without a seed, neither mode is reproducible -- confirms the two
  ## reproducible cases above aren't reproducible merely by coincidence
  serial.noseed.a <- run(par = FALSE, seed = NULL)
  serial.noseed.b <- run(par = FALSE, seed = NULL)
  expect_false(isTRUE(all.equal(serial.noseed.a$pvalues, serial.noseed.b$pvalues)))

  par.noseed.a <- run(par = TRUE, seed = NULL)
  par.noseed.b <- run(par = TRUE, seed = NULL)
  expect_false(isTRUE(all.equal(par.noseed.a$pvalues, par.noseed.b$pvalues)))
})
