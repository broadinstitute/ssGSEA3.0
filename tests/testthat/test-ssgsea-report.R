make_tmp_entry <- function(es, nes, pval, ol, ol.numb, ol.perc) {
  list(ES.vector = es, NES.vector = nes, p.val.vector = pval,
       OL.vector = ol, OL.numb.vector = ol.numb, OL.perc.vector = ol.perc,
       gene.set.size = length(ol), random.walk = as.list(rep(NA, length(es))))
}

test_that("build_ssgsea_result: FDR is computed per-sample-column by default, globally when requested", {
  tmp <- list(
    GS1 = make_tmp_entry(c(1, 1), c(1, 1), c(0.01, 0.20), c("G1", "G1"), c(1, 1), c(50, 50)),
    GS2 = make_tmp_entry(c(2, 2), c(2, 2), c(0.02, 0.30), c("G2", "G2"), c(1, 1), c(50, 50)),
    GS3 = make_tmp_entry(c(3, 3), c(3, 3), c(0.03, 0.40), c("G3", "G3"), c(1, 1), c(50, 50))
  )
  gs.names <- c("GS1", "GS2", "GS3")
  gs <- list(GS1 = "G1", GS2 = "G2", GS3 = "G3")
  sample.names <- c("S1", "S2")
  sample.descs <- data.frame(row.names = sample.names)

  res.col <- build_ssgsea_result(tmp, gs.names, c("d1", "d2", "d3"), c(1, 1, 1), NULL, gs,
                                  sample.names, sample.descs, output.score.type = "NES", global.fdr = FALSE)
  expect_equal(unname(res.col$fdr[, "S1"]), stats::p.adjust(c(0.01, 0.02, 0.03), method = "fdr"))
  expect_equal(unname(res.col$fdr[, "S2"]), stats::p.adjust(c(0.20, 0.30, 0.40), method = "fdr"))

  res.glob <- build_ssgsea_result(tmp, gs.names, c("d1", "d2", "d3"), c(1, 1, 1), NULL, gs,
                                   sample.names, sample.descs, output.score.type = "NES", global.fdr = TRUE)
  expect_equal(as.vector(t(as.matrix(res.glob$fdr))),
               stats::p.adjust(c(0.01, 0.20, 0.02, 0.30, 0.03, 0.40), method = "fdr"))
})

test_that("build_ssgsea_result: gene sets that scored in zero samples are dropped, keeping gs/random.walk in sync", {
  tmp <- list(
    GS1 = make_tmp_entry(c(1, 2), c(1.5, 2.5), c(0.01, 0.02), c("G1|G2", "G1|G2"), c(2, 2), c(100, 100)),
    GS2 = make_tmp_entry(c(NA, NA), c(NA, NA), c(NA, NA), c(NA, NA), c(0, 0), c(0, 0)),
    GS3 = make_tmp_entry(c(3, NA), c(3.5, NA), c(0.03, NA), c("G3", NA), c(1, 0), c(50, 0))
  )
  gs.names <- c("GS1", "GS2", "GS3")
  gs <- list(GS1 = c("G1", "G2"), GS2 = c("G4", "G5"), GS3 = c("G3", "G6"))
  sample.names <- c("S1", "S2")
  sample.descs <- data.frame(row.names = sample.names)

  res <- build_ssgsea_result(tmp, gs.names, c("d1", "d2", "d3"), c(2, 2, 2), NULL, gs,
                              sample.names, sample.descs, output.score.type = "NES")

  expect_equal(res$gs.names, c("GS1", "GS3"))
  expect_equal(names(res$gs), c("GS1", "GS3"))
  expect_equal(res$gs$GS3, c("G3", "G6"))
  expect_equal(nrow(res$scores), 2)
  expect_equal(unname(res$scores["GS1", ]), c(1.5, 2.5))
  expect_equal(length(res$random.walk), 2)
})

test_that("build_ssgsea_result: overlap/gs.descs columns are named and populated as expected", {
  tmp <- list(GS1 = make_tmp_entry(c(1), c(1.5), c(0.04), c("G1|G2"), c(2), c(66.7)))
  gs.names <- "GS1"
  gs <- list(GS1 = c("G1", "G2", "G3"))
  sample.names <- "S1"
  sample.descs <- data.frame(row.names = sample.names)

  res <- build_ssgsea_result(tmp, gs.names, "a hallmark set", 3, NULL, gs,
                              sample.names, sample.descs, output.score.type = "NES", extended.output = TRUE)

  expect_true("Signature.set.overlap.S1" %in% colnames(res$gs.descs))
  expect_equal(res$gs.descs$Signature.set.description, "a hallmark set")
  expect_equal(res$gs.descs$Signature.set.size, 3)
  expect_equal(res$overlap["GS1", 1], "G1|G2")
})

test_that("build_ssgsea_result: errors clearly if every gene set scored zero samples", {
  tmp <- list(GS1 = make_tmp_entry(c(NA), c(NA), c(NA), c(NA), c(0), c(0)))
  expect_error(
    build_ssgsea_result(tmp, "GS1", "d", 1, NULL, list(GS1 = "G1"), "S1", data.frame(row.names = "S1")),
    "nothing to report"
  )
})

test_that("write_ssgsea_gcts: writes scores/pvalues/fdr/combined + per-signature GCTs that round-trip", {
  tmp <- list(GS1 = make_tmp_entry(c(1, 2), c(1.5, 2.5), c(0.01, 0.02), c("G1|G2", "G1|G2"), c(2, 2), c(100, 100)))
  gs.names <- "GS1"
  gs <- list(GS1 = c("G1", "G2"))
  sample.names <- c("S1", "S2")
  sample.descs <- data.frame(row.names = sample.names)

  result <- build_ssgsea_result(tmp, gs.names, "desc", 2, NULL, gs, sample.names, sample.descs,
                                 output.score.type = "NES")

  m.org <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, dimnames = list(c("G1", "G2", "G3"), sample.names))
  gene.descs <- data.frame(row.names = rownames(m.org))

  old.wd <- getwd()
  wd <- withr::local_tempdir()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  written <- write_ssgsea_gcts(result, output.prefix = "out", m.org = m.org, gene.names = rownames(m.org),
                                gene.descs = gene.descs, export.signat.gct = TRUE, appenddim = FALSE)

  expect_true(file.exists("out-scores.gct"))
  expect_true(file.exists("out-pvalues.gct"))
  expect_true(file.exists("out-fdr-pvalues.gct"))
  expect_true(file.exists("out-combined.gct"))
  expect_true(length(dir("signature_gct")) >= 1)

  ## NOTE: with only 1 surviving gene set, this GCT has exactly 1 data row and
  ## >1 row-descriptor columns -- exactly the shape that triggers a cmapR
  ## parse_gctx() bug (see read_gct_safe()), so read it back with the safe
  ## wrapper rather than cmapR::parse_gctx() directly.
  scores.gct <- read_gct_safe("out-scores.gct")
  expect_equal(unname(scores.gct@mat["GS1", ]), c(1.5, 2.5))

  sig.file <- dir("signature_gct", full.names = TRUE)[1]
  sig.gct <- read_gct_safe(sig.file)
  expect_setequal(sig.gct@rid, c("G1", "G2")) ## only signature members, not G3
})

test_that("write_ssgsea_gcts: a single surviving gene set (1 data row, >1 rdesc cols) round-trips via read_gct_safe", {
  tmp <- list(GS1 = make_tmp_entry(c(1, 2), c(1.5, 2.5), c(0.01, 0.02), c("G1|G2", "G1|G2"), c(2, 2), c(100, 100)))
  result <- build_ssgsea_result(tmp, "GS1", "desc", 2, NULL, list(GS1 = c("G1", "G2")),
                                 c("S1", "S2"), data.frame(row.names = c("S1", "S2")), output.score.type = "NES")

  m.org <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, dimnames = list(c("G1", "G2", "G3"), c("S1", "S2")))
  gene.descs <- data.frame(row.names = rownames(m.org))

  old.wd <- getwd()
  wd <- withr::local_tempdir()
  setwd(wd)
  on.exit(setwd(old.wd), add = TRUE)

  write_ssgsea_gcts(result, output.prefix = "single", m.org = m.org, gene.names = rownames(m.org),
                     gene.descs = gene.descs, export.signat.gct = FALSE, appenddim = FALSE)

  reparsed <- read_gct_safe("single-combined.gct")
  expect_equal(dim(reparsed@mat), c(1, 2))
  expect_equal(unname(reparsed@mat["GS1", ]), c(1.5, 2.5))
})
