#' ssGSEA/PTM-SEA summary heatmap
#'
#' Reworks the lab's existing `pw_hm()` heatmap (base `pheatmap`) onto
#' `ComplexHeatmap`. Selects gene sets via [select_top_terms()] by default
#' (FDR-significant, then top-N by direction), colors a symmetric blue-white-
#' red scale centered at zero, and marks FDR-significant cells with `*`.
#'
#' Row grouping is auto-detected from the selected terms' names:
#'
#' - If any term is a bundled [hallmark_categories] name, rows are split and
#'   colored by HALLMARK process category (matching `pw_hm()`'s behavior
#'   exactly), with the `HALLMARK_` prefix stripped from row labels.
#' - Otherwise, if `split_by_db_category` and the terms carry a shared
#'   `CATEGORY-name` prefix convention (as in PTMSigDB: `KINASE-`, `PERT-`,
#'   `PATH-`, `DISEASE-`) with more than one distinct prefix present, rows
#'   are split by that prefix.
#' - Otherwise, rows are hierarchically clustered (`cluster_rows`, default
#'   `TRUE`) using `ComplexHeatmap`'s own built-in clustering rather than
#'   `pw_hm()`'s `seriation::seriate()` -- avoids an extra dependency for
#'   what's typically a very similar row order.
#'
#' @param result An `ssgsea_result`, as returned by [ssGSEA3()].
#' @param terms Character vector of gene set names to plot. If `NULL`
#'   (default), gene sets are chosen via [select_top_terms()] using `top_n`/
#'   `direction`/`fdr.max`.
#' @param top_n,direction,fdr.max Passed to [select_top_terms()] when `terms`
#'   is `NULL`; `fdr.max` is also the threshold used for marking significant
#'   cells with `*`.
#' @param samples Character vector of sample names (columns) to include.
#'   Defaults to all samples in `result`.
#' @param split_by_db_category If `TRUE` (default), auto-detect and split by
#'   a shared `CATEGORY-name` prefix (e.g. PTMSigDB) when HALLMARK grouping
#'   doesn't apply.
#' @param cluster_rows Whether to hierarchically cluster rows when neither
#'   HALLMARK nor prefix-based grouping applies. Ignored (rows are left in
#'   `terms` order within each split) when either grouping does apply.
#' @param cell.size Fixed cell width/height, in points, matching `pw_hm()`'s
#'   own `cw`/`ch` default of 10 -- keeps cells square regardless of the
#'   number of rows/columns, rather than stretching to fill the plotting
#'   device (which is what a bare `ComplexHeatmap::Heatmap()` call does by
#'   default when `width`/`height` aren't set).
#' @param title Optional heatmap title.
#' @return A [ComplexHeatmap::Heatmap-class] object -- draw it with
#'   `ComplexHeatmap::draw()`, or save with e.g. `pdf()`/`png()` around that
#'   call.
#' @export
plot_ssgsea_heatmap <- function(result, terms = NULL, top_n = 10,
                                 direction = c("both", "pos", "neg"), fdr.max = 0.05,
                                 samples = NULL, split_by_db_category = TRUE,
                                 cluster_rows = TRUE, cell.size = 10, title = NULL) {
  if (is.null(terms)) {
    terms <- select_top_terms(result, top_n = top_n, direction = direction, fdr.max = fdr.max)
  }
  if (length(terms) == 0) {
    stop("No gene sets to plot (none selected/none survived filtering).")
  }
  missing.terms <- setdiff(terms, result$gs.names)
  if (length(missing.terms) > 0) {
    stop("Gene set(s) not found in `result`: ", paste(missing.terms, collapse = ", "))
  }

  samples <- if (is.null(samples)) result$sample.names else samples
  mat <- result$scores[terms, samples, drop = FALSE]
  fdr <- result$fdr[terms, samples, drop = FALSE]

  is.hallmark <- any(terms %in% names(hallmark_categories))
  prefixes <- sub("^(.*?)-.*", "\\1", terms)
  is.multi.category <- split_by_db_category && length(unique(prefixes)) > 1 && all(grepl("-", terms))

  row_split <- NULL
  row_labels <- terms
  left_annotation <- NULL
  do.cluster <- cluster_rows

  if (is.hallmark) {
    category <- unname(hallmark_categories[terms])
    category[is.na(category)] <- ""
    row_split <- factor(category, levels = unique(category))
    row_labels <- sub("^HALLMARK_", "", terms)
    ## the color mapping must cover every distinct value actually present,
    ## including "" for terms absent from the bundled lookup (mixed-database
    ## selections, or a HALLMARK set added upstream but not yet in our table)
    cat.colors <- hallmark_category_colors[intersect(names(hallmark_category_colors), unique(category))]
    if ("" %in% category) cat.colors[""] <- "grey90"
    left_annotation <- ComplexHeatmap::rowAnnotation(
      category = category,
      col = list(category = cat.colors),
      show_annotation_name = FALSE
    )
    do.cluster <- FALSE
  } else if (is.multi.category) {
    row_split <- factor(prefixes, levels = unique(prefixes))
    do.cluster <- FALSE
  }

  ## zero-centered but *not* forced symmetric -- an equal-magnitude scale
  ## wastes most of its resolution (and looks like a bug: a legend that goes
  ## to -40 when nothing in the data is more negative than -12) whenever the
  ## positive and negative extremes differ, which real scores usually do
  max.val <- ceiling(max(mat, na.rm = TRUE))
  min.val <- floor(min(mat, na.rm = TRUE))
  col_fun <- circlize::colorRamp2(c(min.val, 0, max.val), c("blue3", "white", "#DC3220"))

  ComplexHeatmap::Heatmap(
    mat,
    name = result$output.score.type,
    col = col_fun,
    na_col = "gray",
    width = grid::unit(ncol(mat) * cell.size, "points"),
    height = grid::unit(nrow(mat) * cell.size, "points"),
    row_split = row_split,
    cluster_row_slices = FALSE,
    cluster_rows = do.cluster,
    cluster_columns = FALSE,
    row_labels = row_labels,
    ## ComplexHeatmap's own default (6cm) silently truncates long gene set
    ## names (e.g. full HALLMARK_* names) regardless of device size -- this
    ## widens the cap rather than raising it per-plot from the caller side
    row_names_max_width = grid::unit(12, "cm"),
    left_annotation = left_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = "#555555", fill = fill))
      if (!is.na(fdr[i, j]) && fdr[i, j] < fdr.max) {
        ## switch to white on dark fills (e.g. deep blue) so the marker
        ## stays legible -- perceived-brightness (YIQ) formula, standard
        ## threshold of 0.5 on a 0-1 scale
        rgb <- grDevices::col2rgb(fill) / 255
        luminance <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
        star.col <- if (luminance < 0.5) "white" else "black"
        ## nudge down slightly -- "*" sits high within its own glyph box in
        ## most fonts, so an unadjusted y leaves it looking top-aligned
        grid::grid.text("*", x, y - grid::unit(0.15, "char"), gp = grid::gpar(col = star.col))
      }
    },
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    row_title_rot = 0,
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9)
  )
}
