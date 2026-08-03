make_result <- function(scores, fdr, gs.names, gs, gs.direction = NULL, sample.names, data = NULL, gene.names = NULL) {
  structure(
    list(scores = scores, fdr = fdr, gs.names = gs.names, gs = gs, gs.direction = gs.direction,
         sample.names = sample.names, data = data, gene.names = gene.names),
    class = "ssgsea_result"
  )
}

test_that("select_top_terms: picks top-N significant terms per direction", {
  ## A: sig, best (max-|.|) score = +5. B: sig, best score = -8. C: not
  ## significant (min fdr = 0.3 > default 0.05).
  scores <- rbind(A = c(1, 5), B = c(-8, 2), C = c(-0.1, 3))
  colnames(scores) <- c("S1", "S2")
  fdr <- rbind(A = c(0.2, 0.01), B = c(0.01, 0.2), C = c(0.3, 0.3))
  colnames(fdr) <- c("S1", "S2")
  res <- make_result(scores, fdr, c("A", "B", "C"), NULL, sample.names = c("S1", "S2"))

  sel <- select_top_terms(res, top_n = 10, fdr.max = 0.05)
  expect_setequal(sel, c("A", "B"))

  sel.pos <- select_top_terms(res, direction = "pos", fdr.max = 0.05)
  expect_equal(sel.pos, "A")

  sel.neg <- select_top_terms(res, direction = "neg", fdr.max = 0.05)
  expect_equal(sel.neg, "B")
})

test_that("select_top_terms: returns character(0) when nothing is significant", {
  scores <- matrix(1, nrow = 1, dimnames = list("A", "S1"))
  fdr <- matrix(0.9, nrow = 1, dimnames = list("A", "S1"))
  res <- make_result(scores, fdr, "A", NULL, sample.names = "S1")
  expect_equal(select_top_terms(res, fdr.max = 0.05), character(0))
})

test_that("select_top_terms: top_n caps the number of terms per direction", {
  scores <- matrix(c(5, 4, 3, 2, 1), ncol = 1, dimnames = list(paste0("T", 1:5), "S1"))
  fdr <- matrix(0.01, nrow = 5, ncol = 1, dimnames = list(paste0("T", 1:5), "S1"))
  res <- make_result(scores, fdr, paste0("T", 1:5), NULL, sample.names = "S1")

  sel <- select_top_terms(res, top_n = 2, direction = "pos", fdr.max = 0.05)
  expect_equal(sel, c("T1", "T2")) ## highest two scores
})

test_that("plot_rank_diagnostics: errors clearly if result lacks the original data", {
  res <- make_result(matrix(1, dimnames = list("A", "S1")), matrix(0.01, dimnames = list("A", "S1")),
                      "A", list(A = "G1"), sample.names = "S1")
  expect_error(plot_rank_diagnostics(res, terms = "A", file = NULL), "original expression data")
})

test_that("plot_rank_diagnostics: errors clearly for an unknown term name", {
  res <- make_result(matrix(1, dimnames = list("A", "S1")), matrix(0.01, dimnames = list("A", "S1")),
                      "A", list(A = "G1"), sample.names = "S1",
                      data = matrix(1:3, dimnames = list(c("G1", "G2", "G3"), "S1")),
                      gene.names = c("G1", "G2", "G3"))
  expect_error(plot_rank_diagnostics(res, terms = "NOT_A_TERM", file = NULL), "not found")
})

test_that("plot_rank_diagnostics: draws a non-directional gene set to a PDF without error", {
  wd <- withr::local_tempdir()
  set.seed(1)
  genes <- paste0("G", 1:30)
  mat <- matrix(c(30:1, rnorm(30, 15, 5)), ncol = 2, dimnames = list(genes, c("S1", "S2")))
  scores <- matrix(c(4, -1), nrow = 1, dimnames = list("GS_A", c("S1", "S2")))
  fdr <- matrix(c(0.01, 0.2), nrow = 1, dimnames = list("GS_A", c("S1", "S2")))
  res <- make_result(scores, fdr, "GS_A", list(GS_A = genes[1:8]), sample.names = c("S1", "S2"),
                      data = mat, gene.names = genes)

  out <- plot_rank_diagnostics(res, terms = "GS_A", file = file.path(wd, "out.pdf"))
  expect_true(file.exists(file.path(wd, "out.pdf")))
  expect_equal(out, "GS_A")
})

test_that("plot_rank_diagnostics: draws a directional (PTMsigDB-style) gene set without error", {
  wd <- withr::local_tempdir()
  set.seed(2)
  genes <- paste0("G", 1:30)
  mat <- matrix(c(30:1, rnorm(30, 15, 5)), ncol = 2, dimnames = list(genes, c("S1", "S2")))
  scores <- matrix(c(6, -1), nrow = 1, dimnames = list("GS_DIR", c("S1", "S2")))
  fdr <- matrix(c(0.001, 0.2), nrow = 1, dimnames = list("GS_DIR", c("S1", "S2")))
  gs <- list(GS_DIR = c(genes[1:3], genes[28:30]))
  gs.dir <- list(GS_DIR = c("u", "u", "u", "d", "d", "d"))
  res <- make_result(scores, fdr, "GS_DIR", gs, gs.direction = gs.dir, sample.names = c("S1", "S2"),
                      data = mat, gene.names = genes)

  out <- plot_rank_diagnostics(res, terms = "GS_DIR", file = file.path(wd, "out.pdf"))
  expect_true(file.exists(file.path(wd, "out.pdf")))
})
