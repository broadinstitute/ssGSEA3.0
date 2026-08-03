## Shared helper: build a minimal, valid GCT object and write it to a temp file
make_gct_file <- function(mat, rid, cid, rdesc = NULL) {
  gct <- cmapR::GCT()
  dimnames(mat) <- list(rid, cid)
  gct@mat <- mat
  gct@rid <- rid
  gct@cid <- cid
  gct@rdesc <- if (is.null(rdesc)) data.frame(row.names = rid) else rdesc
  gct@cdesc <- data.frame(row.names = cid)
  f <- tempfile(fileext = ".gct")
  cmapR::write_gct(gct, ofile = f, appenddim = FALSE)
  f
}

## #####################################################################
##                              gcr
## #####################################################################
test_that("gcr keeps every row and relabels by gene symbol", {
  mat <- matrix(1:6, nrow = 3, ncol = 2)
  rdesc <- data.frame(geneSymbol = c("GeneA", "GeneA", "GeneB"), row.names = c("p1", "p2", "p3"))
  f <- make_gct_file(mat, c("p1", "p2", "p3"), c("S1", "S2"), rdesc)

  out.f <- preprocessGCT(f, level = "gcr", gene.col = "geneSymbol", appenddim = FALSE)
  out <- cmapR::parse_gctx(out.f)

  expect_equal(nrow(out@mat), 3)
  expect_equal(sort(out@rid), sort(c("GeneA", "GeneA_1", "GeneB")))
})

## #####################################################################
##                     gc: mean / median / abs.max
## #####################################################################
test_that("gc mean/median/abs.max collapse multiple rows per gene", {
  mat <- matrix(c(1, 3, 10, 10,   ## S1: probe1(GeneA)=1, probe2(GeneA)=3, probe1(GeneB)=10, probe2(GeneB)=10
                  2, 4, 5, 100),  ## S2
                nrow = 4, ncol = 2,
                dimnames = list(c("pA1", "pA2", "pB1", "pB2"), c("S1", "S2")))
  rdesc <- data.frame(geneSymbol = c("GeneA", "GeneA", "GeneB", "GeneB"), row.names = rownames(mat))
  f <- make_gct_file(mat, rownames(mat), colnames(mat), rdesc)

  out.f <- preprocessGCT(f, level = "gc", mode = "mean", gene.col = "geneSymbol", appenddim = FALSE)
  out <- cmapR::parse_gctx(out.f)
  expect_equal(nrow(out@mat), 2)
  expect_equal(unname(out@mat["GeneA", "S1"]), mean(c(1, 3)))
  expect_equal(unname(out@mat["GeneB", "S2"]), mean(c(5, 100)))

  out.f2 <- preprocessGCT(f, level = "gc", mode = "abs.max", gene.col = "geneSymbol", appenddim = FALSE)
  out2 <- cmapR::parse_gctx(out.f2)
  expect_equal(unname(out2@mat["GeneB", "S2"]), 100) ## |100| > |5|
})

## #####################################################################
##                     gc: mode = "sd" (the fix)
## #####################################################################
test_that("gc mode='sd' keeps the most cross-sample-variable probe per gene (previously crashed)", {
  ## GeneA: probe1 is flat (sd=0), probe2 is variable (sd>0) -> keep probe2
  mat <- matrix(c(
    10, 10, 10,   ## pA1 (flat)
    1, 50, 100,   ## pA2 (variable) -- keep this one
    5, 5, 5,      ## pB1 (flat)     -- keep this one (ties broken by which.max -> first)
    5, 5, 5       ## pB2 (flat, identical to pB1)
  ), nrow = 4, byrow = TRUE, dimnames = list(c("pA1", "pA2", "pB1", "pB2"), c("S1", "S2", "S3")))
  rdesc <- data.frame(geneSymbol = c("GeneA", "GeneA", "GeneB", "GeneB"), row.names = rownames(mat))
  f <- make_gct_file(mat, rownames(mat), colnames(mat), rdesc)

  out.f <- preprocessGCT(f, level = "gc", mode = "sd", gene.col = "geneSymbol", appenddim = FALSE)
  out <- cmapR::parse_gctx(out.f)

  expect_equal(nrow(out@mat), 2)
  expect_equal(unname(out@mat["GeneA", ]), c(1, 50, 100)) ## pA2 kept, not pA1
})

## #####################################################################
##                     gc: mode = "SGT"
## #####################################################################
test_that("gc mode='SGT' keeps the lowest subgroup per gene", {
  mat <- matrix(1:6, nrow = 3, ncol = 2, dimnames = list(c("p1", "p2", "p3"), c("S1", "S2")))
  rdesc <- data.frame(
    geneSymbol = c("GeneA", "GeneA", "GeneB"),
    subgroupNum = c("100.2", "100.1", "200.1"), ## GeneA: subgroup 1 < subgroup 2 -> keep p2
    row.names = rownames(mat)
  )
  f <- make_gct_file(mat, rownames(mat), colnames(mat), rdesc)

  out.f <- preprocessGCT(f, level = "gc", mode = "SGT", gene.col = "geneSymbol", SGT.col = "subgroupNum", appenddim = FALSE)
  out <- cmapR::parse_gctx(out.f)

  expect_equal(unname(out@mat["GeneA", "S1"]), mat["p2", "S1"])
})

## #####################################################################
##                  ssc (seqwin path -- no ID mapping needed)
## #####################################################################
test_that("ssc/seqwin explodes multi-site peptides and averages duplicate site detections", {
  ## PEP1 covers 2 sites (SITEAAAA, SITEBBBB); PEP2 separately/singly covers
  ## SITEAAAA (so only SITEAAAA -- not SITEBBBB -- is excluded from PEP1's
  ## multi-site list, leaving PEP1 to represent SITEBBBB alone). PEP3 and
  ## PEP4 both singly cover SITECCCC (a genuine duplicate site detection,
  ## to be averaged). PEP5 singly covers SITEDDDD (no duplication at all).
  mat <- matrix(c(1, 2, 10, 20, 100, 110, 200, 220, 1000, 1010), nrow = 5, byrow = TRUE,
                 dimnames = list(c("PEP1", "PEP2", "PEP3", "PEP4", "PEP5"), c("S1", "S2")))
  rdesc <- data.frame(
    VMsiteFlanks = c("SITEAAAA|SITEBBBB", "SITEAAAA", "SITECCCC", "SITECCCC", "SITEDDDD"),
    row.names = rownames(mat)
  )
  f <- make_gct_file(mat, rownames(mat), colnames(mat), rdesc)

  out.f <- preprocessGCT(f, level = "ssc", id.type = "sm", id.type.out = "seqwin",
                          seqwin.col = "VMsiteFlanks", mod.type = "p", mode = "mean", appenddim = FALSE)
  out <- cmapR::parse_gctx(out.f)

  expect_setequal(out@rid, c("SITEAAAA-p", "SITEBBBB-p", "SITECCCC-p", "SITEDDDD-p"))
  expect_equal(unname(out@mat["SITEBBBB-p", ]), c(1, 2))         ## PEP1, with SITEAAAA excluded
  expect_equal(unname(out@mat["SITEAAAA-p", ]), c(10, 20))       ## PEP2 alone
  expect_equal(unname(out@mat["SITECCCC-p", ]), c(mean(c(100, 200)), mean(c(110, 220)))) ## PEP3 + PEP4 averaged
  expect_equal(unname(out@mat["SITEDDDD-p", ]), c(1000, 1010))  ## PEP5 alone, no duplication
})

## #####################################################################
##                        error conditions
## #####################################################################
test_that("id.type='wg' and id.type.out='psp' are no longer accepted (previously crashed or silently misbehaved)", {
  expect_error(match.arg("wg", c("sm", "ph")), "should be one of")
  expect_error(match.arg("psp", c("uniprot", "refseq", "seqwin")), "should be one of")
})

test_that("mode='sd'/'SGT' error clearly (not cryptically) for level='ssc'", {
  mat <- matrix(1:2, nrow = 1, dimnames = list("PEP1", c("S1", "S2")))
  rdesc <- data.frame(VMsiteFlanks = "SITEAAAA", row.names = "PEP1")
  f <- make_gct_file(mat, rownames(mat), colnames(mat), rdesc)

  expect_error(
    preprocessGCT(f, level = "ssc", id.type = "sm", id.type.out = "seqwin",
                  seqwin.col = "VMsiteFlanks", mode = "sd", appenddim = FALSE),
    "not supported for level = 'ssc'"
  )
})

test_that("missing organism errors clearly when UniProt mapping is required", {
  mat <- matrix(1:2, nrow = 1, dimnames = list("NP_000001_1_1_1_S100s", c("S1", "S2")))
  rdesc <- data.frame(row.names = rownames(mat))
  f <- make_gct_file(mat, rownames(mat), colnames(mat), rdesc)

  expect_error(
    preprocessGCT(f, level = "ssc", id.type = "sm", id.type.out = "uniprot",
                  acc.type = "refseq", organism = NULL, appenddim = FALSE),
    "`organism` must be specified"
  )
})

## #####################################################################
##                    RefSeq2UniProt (explicit organism)
## #####################################################################
test_that("RefSeq2UniProt maps real human RefSeq accessions with organism='human'", {
  skip_if_not_installed("org.Hs.eg.db")
  res <- RefSeq2UniProt(c("NP_005219.2", "NP_000537.3"), keytype = "REFSEQ", organism = "human")

  expect_equal(res$organism, "human")
  expect_equal(res$id.map$id.query, c("NP_005219", "NP_000537"))
  expect_true(all(grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]", res$id.map$id.mapped)))
})

test_that("RefSeq2UniProt requires an installed annotation package for the requested organism", {
  skip_if_not_installed("org.Hs.eg.db") ## sanity: this test only makes sense if human works
  ## zebrafish annotation package is unlikely to be installed in CI; if it
  ## is, this test is skipped rather than giving a false failure
  skip_if(requireNamespace("org.Dr.eg.db", quietly = TRUE))
  expect_error(
    RefSeq2UniProt("NP_005219.2", keytype = "REFSEQ", organism = "zebrafish"),
    "is required to map organism"
  )
})
