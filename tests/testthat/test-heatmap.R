make_result <- function(terms, scores, fdr, sample.names = c("S1", "S2")) {
  structure(list(scores = scores, fdr = fdr, gs.names = terms, sample.names = sample.names,
                 output.score.type = "NES"),
            class = "ssgsea_result")
}

test_that("plot_ssgsea_heatmap: HALLMARK terms are split/colored by category with prefix stripped", {
  terms <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_E2F_TARGETS")
  scores <- matrix(c(5, -3, 2, 1, 6, -2), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(c(0.01, 0.2, 0.3, 0.5, 0.02, 0.6), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms)

  expect_s4_class(ht, "Heatmap")
  expect_equal(rownames(ht@matrix), terms) ## underlying matrix keeps full names
  expect_equal(ht@row_names_param$labels, c("TNFA_SIGNALING_VIA_NFKB", "HYPOXIA", "E2F_TARGETS"))
  split.col <- ht@matrix_param$row_split[[1]]
  expect_equal(as.character(split.col), c("signaling", "pathway", "proliferation"))
  expect_false(is.null(ht@left_annotation))
})

test_that("plot_ssgsea_heatmap: HALLMARK terms not in the bundled lookup get an 'other' bucket, not an error", {
  terms <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_NOT_A_REAL_SET")
  scores <- matrix(c(5, -3, 1, 6), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(c(0.01, 0.2, 0.5, 0.02), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms)
  split.col <- as.character(ht@matrix_param$row_split[[1]])
  expect_equal(split.col, c("signaling", ""))
})

test_that("plot_ssgsea_heatmap: PTMSigDB-style multi-category prefixes are split, no HALLMARK coloring", {
  terms <- c("KINASE-PSP_CDK1", "KINASE-PSP_CDK2", "PERT-PSP_EGF", "PATH-WP_Ras_signaling")
  scores <- matrix(rnorm(8), nrow = 4, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(8), nrow = 4, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms)
  split.col <- as.character(ht@matrix_param$row_split[[1]])
  expect_equal(split.col, c("KINASE", "KINASE", "PERT", "PATH"))
  expect_true(is.null(ht@left_annotation))
})

test_that("plot_ssgsea_heatmap: split_by_db_category=FALSE disables PTMSigDB-style splitting", {
  terms <- c("KINASE-PSP_CDK1", "PERT-PSP_EGF")
  scores <- matrix(rnorm(4), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(4), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms, split_by_db_category = FALSE)
  expect_null(ht@matrix_param$row_split)
})

test_that("plot_ssgsea_heatmap: a single shared prefix does not trigger a split (nothing to split by)", {
  terms <- c("KINASE-PSP_CDK1", "KINASE-PSP_CDK2", "KINASE-PSP_CDK3")
  scores <- matrix(rnorm(6), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(6), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms)
  expect_null(ht@matrix_param$row_split)
})

test_that("plot_ssgsea_heatmap: generic non-prefixed terms fall back to hierarchical clustering", {
  terms <- c("REACTOME_A", "REACTOME_B", "REACTOME_C")
  scores <- matrix(rnorm(6), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(6), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms)
  expect_null(ht@matrix_param$row_split)
  expect_true(ht@row_dend_param$cluster)
})

test_that("plot_ssgsea_heatmap: cluster_rows=FALSE is respected for the generic case", {
  terms <- c("REACTOME_A", "REACTOME_B")
  scores <- matrix(rnorm(4), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(4), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms, cluster_rows = FALSE)
  expect_false(ht@row_dend_param$cluster)
})

test_that("plot_ssgsea_heatmap: color scale is zero-centered but NOT forced symmetric, sized to the actual data range", {
  terms <- c("A", "B")
  scores <- matrix(c(2.3, -7.8, 1, -1), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(4), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms)
  ## breaks should span [floor(min), ceiling(max)] = [-8, 3] -- NOT
  ## [-8, 8] -- forcing symmetry around the larger-magnitude side wastes
  ## color resolution and misleadingly implies values that don't occur
  col_fun <- ht@matrix_color_mapping@col_fun
  breaks <- environment(col_fun)$breaks
  expect_equal(range(breaks), c(-8, 3))
})

test_that("plot_ssgsea_heatmap: legend name reflects the result's actual score type", {
  terms <- c("A", "B")
  scores <- matrix(c(1, -1, 2, -2), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(4), nrow = 2, dimnames = list(terms, c("S1", "S2")))

  result.nes <- make_result(terms, scores, fdr)
  ht.nes <- plot_ssgsea_heatmap(result.nes, terms = terms)
  expect_equal(ht.nes@matrix_color_mapping@name, "NES")

  result.es <- make_result(terms, scores, fdr)
  result.es$output.score.type <- "ES"
  ht.es <- plot_ssgsea_heatmap(result.es, terms = terms)
  expect_equal(ht.es@matrix_color_mapping@name, "ES")
})

test_that("plot_ssgsea_heatmap: cells are fixed-size (square), not stretched to fill the device", {
  terms <- c("A", "B", "C")
  scores <- matrix(rnorm(6), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(6), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms, cell.size = 12)
  expect_equal(ht@matrix_param$width, grid::unit(2 * 12, "points"))
  expect_equal(ht@matrix_param$height, grid::unit(3 * 12, "points"))
})

test_that("plot_ssgsea_heatmap: significance asterisks switch to white on dark cell fills", {
  terms <- c("A", "B", "C")
  scores <- matrix(rnorm(3), nrow = 3, dimnames = list(terms, "S1"))
  fdr <- matrix(0.001, nrow = 3, ncol = 1, dimnames = list(terms, "S1")) ## all significant
  result <- make_result(terms, scores, fdr, sample.names = "S1")
  ht <- plot_ssgsea_heatmap(result, terms = terms)
  cell_fun <- ht@matrix_param$cell_fun

  star.col <- character(0)
  fake.grid.text <- function(label, x, y, gp = grid::gpar()) {
    star.col[[length(star.col) + 1]] <<- gp$col
  }
  testthat::local_mocked_bindings(grid.text = fake.grid.text, .package = "grid")

  pdf(tempfile(fileext = ".pdf")) ## grid.rect() needs an open device
  on.exit(dev.off(), add = TRUE)
  grid::grid.newpage()
  cell_fun(1, 1, grid::unit(0.5, "npc"), grid::unit(0.5, "npc"),
           grid::unit(1, "npc"), grid::unit(1, "npc"), fill = "blue3") ## very dark
  cell_fun(1, 1, grid::unit(0.5, "npc"), grid::unit(0.5, "npc"),
           grid::unit(1, "npc"), grid::unit(1, "npc"), fill = "#FFEEEE") ## very light

  expect_equal(star.col, c("white", "black"))
})

test_that("plot_ssgsea_heatmap: default term selection defers to select_top_terms()", {
  terms <- c("A", "B", "C")
  scores <- matrix(c(5, -8, 1, 2, -3, 0.5), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(c(0.01, 0.01, 0.9, 0.2, 0.2, 0.9), nrow = 3, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, top_n = 10, fdr.max = 0.05)
  ## only A and B pass fdr.max=0.05 in >=1 sample; C should be excluded
  expect_setequal(rownames(ht@matrix), c("A", "B"))
})

test_that("plot_ssgsea_heatmap: errors clearly when nothing is selected or a named term is missing", {
  terms <- c("A", "B")
  scores <- matrix(c(1, 1, 1, 1), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(0.9, nrow = 2, ncol = 2, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  expect_error(plot_ssgsea_heatmap(result, fdr.max = 0.01), "No gene sets")
  expect_error(plot_ssgsea_heatmap(result, terms = "NOT_A_TERM"), "not found")
})

test_that("plot_ssgsea_heatmap: samples= restricts to a subset of columns", {
  terms <- c("A", "B")
  scores <- matrix(1:6, nrow = 2, dimnames = list(terms, c("S1", "S2", "S3")))
  fdr <- matrix(runif(6), nrow = 2, dimnames = list(terms, c("S1", "S2", "S3")))
  result <- make_result(terms, scores, fdr, sample.names = c("S1", "S2", "S3"))

  ht <- plot_ssgsea_heatmap(result, terms = terms, samples = c("S1", "S3"))
  expect_equal(colnames(ht@matrix), c("S1", "S3"))
})

test_that("plot_ssgsea_heatmap: show_row_names/show_category/show_heatmap_legend can be individually disabled", {
  terms <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA")
  scores <- matrix(c(5, -3, 2, 1), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(c(0.01, 0.02, 0.03, 0.04), nrow = 2, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht.default <- plot_ssgsea_heatmap(result, terms = terms)
  expect_true(ht.default@row_names_param$show)
  expect_true(ht.default@heatmap_param$show_heatmap_legend)

  ht.hidden <- plot_ssgsea_heatmap(result, terms = terms, show_row_names = FALSE,
                                    show_category = FALSE, show_heatmap_legend = FALSE)
  expect_false(ht.hidden@row_names_param$show)
  expect_false(ht.hidden@heatmap_param$show_heatmap_legend)
  expect_null(ht.hidden@left_annotation)
})

test_that("plot_ssgsea_heatmap: cluster_rows accepts a precomputed clustering function, giving identical row order across heatmaps built from different data (e.g. orig/fixed/v3 comparison panels)", {
  terms <- c("A", "B", "C", "D")
  set.seed(1)
  scores1 <- matrix(rnorm(8), nrow = 4, dimnames = list(terms, c("S1", "S2")))
  scores2 <- scores1 + matrix(rnorm(8, sd = 0.2), nrow = 4) ## different data, same terms
  fdr <- matrix(runif(8), nrow = 4, dimnames = list(terms, c("S1", "S2")))
  result1 <- make_result(terms, scores1, fdr)
  result2 <- make_result(terms, scores2, fdr)

  ## a clustering function that ignores whatever matrix it's asked to
  ## cluster and always ties back to one precomputed hclust -- this is the
  ## mechanism the manuscript's 3-panel (orig/fixed/v3) comparison figures
  ## rely on to keep rows aligned across panels built from different data
  precomputed <- hclust(dist(scores1))
  cluster_fn <- function(mat_slice) precomputed

  pdf(tempfile(fileext = ".pdf")) ## draw() needs an open device
  on.exit(dev.off(), add = TRUE)
  ht1 <- ComplexHeatmap::draw(plot_ssgsea_heatmap(result1, terms = terms, cluster_rows = cluster_fn))
  ht2 <- ComplexHeatmap::draw(plot_ssgsea_heatmap(result2, terms = terms, cluster_rows = cluster_fn))
  expect_equal(ComplexHeatmap::row_order(ht1), ComplexHeatmap::row_order(ht2))
})

test_that("plot_ssgsea_heatmap: show_row_dend=FALSE suppresses the dendrogram without changing row order", {
  terms <- c("A", "B", "C", "D")
  scores <- matrix(rnorm(8), nrow = 4, dimnames = list(terms, c("S1", "S2")))
  fdr <- matrix(runif(8), nrow = 4, dimnames = list(terms, c("S1", "S2")))
  result <- make_result(terms, scores, fdr)

  ht <- plot_ssgsea_heatmap(result, terms = terms, show_row_dend = FALSE)
  expect_false(ht@row_dend_param$show)
})
