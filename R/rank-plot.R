#' Select top-N FDR-significant gene sets by direction
#'
#' Picks gene sets that reach `fdr.max` in at least one sample, then keeps
#' the `top_n` largest by absolute score (its most extreme per-sample value),
#' separately for positive- and negative-scoring gene sets.
#'
#' @param result An `ssgsea_result`, as returned by [ssGSEA3()].
#' @param top_n Number of gene sets to keep per direction.
#' @param direction `"both"`, `"pos"`, or `"neg"`.
#' @param fdr.max Maximum (best-case, across samples) FDR-adjusted p-value
#'   for a gene set to be eligible.
#' @return Character vector of selected gene set names.
#' @export
select_top_terms <- function(result, top_n = 10, direction = c("both", "pos", "neg"), fdr.max = 0.05) {
  direction <- match.arg(direction)
  min.fdr <- apply(result$fdr, 1, min, na.rm = TRUE)
  best.score <- apply(result$scores, 1, function(x) x[which.max(abs(x))])

  sig <- names(min.fdr)[min.fdr <= fdr.max]
  if (length(sig) == 0) return(character(0))

  pick <- function(nms) {
    if (length(nms) == 0) return(character(0))
    nms[order(-abs(best.score[nms]))][seq_len(min(top_n, length(nms)))]
  }

  pos <- sig[best.score[sig] > 0]
  neg <- sig[best.score[sig] < 0]

  switch(direction,
    both = c(pick(pos), pick(neg)),
    pos = pick(pos),
    neg = pick(neg)
  )
}

#' @keywords internal
.alpha_color <- function(color, alpha = 80, maxColorValue = 255) {
  rgb <- grDevices::col2rgb(color)
  grDevices::rgb(rgb[1], rgb[2], rgb[3], alpha = alpha, maxColorValue = maxColorValue)
}

#' Rank ("mountain") plot diagnostic for ssGSEA/PTM-SEA results
#'
#' For each selected gene set, draws one page with one panel per sample:
#' the sample's original expression values, ranked high to low, with the
#' gene set's members marked along that ranking. This is the classic GSEA
#' enrichment-plot diagnostic, adapted to ssGSEA/PTM-SEA -- useful for
#' visually confirming that a significant score is actually driven by
#' members clustered toward one end of the ranking, not a fluke.
#'
#' Unlike the original implementation (which re-parsed the written score/
#' p-value/FDR GCT files back off disk for every plot), this draws directly
#' from the `result` object already in memory.
#'
#' @param result An `ssgsea_result`, as returned by [ssGSEA3()].
#' @param terms Character vector of gene set names to plot. If `NULL`
#'   (default), gene sets are chosen via [select_top_terms()] using `top_n`/
#'   `direction`/`fdr.max`.
#' @param top_n,direction,fdr.max Passed to [select_top_terms()] when `terms`
#'   is `NULL`.
#' @param samples Character vector of sample names to plot (as panels on each
#'   gene set's page). Defaults to all samples in `result`.
#' @param file Path to write a (possibly multi-page) PDF to. If `NULL`,
#'   plots to the current graphics device instead (one page per term).
#' @param mfrow Panel layout (rows, cols) per page; samples beyond
#'   `prod(mfrow)` for a single gene set spill onto additional pages.
#' @return (Invisibly) the character vector of gene sets plotted.
#' @export
plot_rank_diagnostics <- function(result, terms = NULL, top_n = 10,
                                   direction = c("both", "pos", "neg"), fdr.max = 0.05,
                                   samples = NULL, file = "rank-plots.pdf",
                                   mfrow = c(3, 3)) {
  if (is.null(result$data) || is.null(result$gene.names)) {
    stop("`result` doesn't include the original expression data -- ",
         "re-run ssGSEA3() to get an up-to-date result object.")
  }

  if (is.null(terms)) {
    terms <- select_top_terms(result, top_n = top_n, direction = direction, fdr.max = fdr.max)
  }
  if (length(terms) == 0) {
    warning("No gene sets to plot (none selected/none survived filtering).")
    return(invisible(character(0)))
  }
  missing.terms <- setdiff(terms, result$gs.names)
  if (length(missing.terms) > 0) {
    stop("Gene set(s) not found in `result`: ", paste(missing.terms, collapse = ", "))
  }

  samples <- if (is.null(samples)) result$sample.names else samples

  if (!is.null(file)) {
    grDevices::pdf(file, width = 9.5, height = 9.5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  for (gs.name in terms) {
    signature <- result$gs[[gs.name]]
    is.directional <- !is.null(result$gs.direction)
    if (is.directional) gs.direction <- result$gs.direction[[gs.name]]

    graphics::par(mfrow = mfrow)
    for (samp.name in samples) {
      score <- result$scores[gs.name, samp.name]
      fdr <- result$fdr[gs.name, samp.name]

      data.expr <- result$data[, samp.name]
      valid.idx <- which(!(is.na(data.expr) | is.infinite(data.expr)))
      data.expr <- data.expr[valid.idx]
      gn <- result$gene.names[valid.idx]

      ord.idx <- order(data.expr, decreasing = TRUE)
      gn <- gn[ord.idx]
      data.expr <- data.expr[ord.idx]

      graphics::plot(data.expr, pch = 20, col = "darkgrey", lwd = 4, type = "l",
                     xlab = "Rank", ylab = "Expression", main = paste(gs.name, samp.name, sep = "\n"),
                     ylim = range(data.expr), yaxs = "i")
      graphics::abline(h = 0, lty = "dashed", lwd = 2, col = "grey70")

      if (is.directional) {
        ## `signature`'s members have their ";u"/";d" suffix already
        ## stripped by match_gene_sets(); direction is tracked in the
        ## parallel `gs.direction` vector instead
        sig.u <- signature[gs.direction == "u"]
        sig.d <- signature[gs.direction == "d"]
        loc.u <- stats::na.omit(match(sig.u, gn))
        loc.d <- stats::na.omit(match(sig.d, gn))

        if (length(loc.u) > 0) {
          graphics::rug(loc.u, col = "darkred", side = 3, lwd = 3, ticksize = 0.02)
          graphics::points(loc.u, data.expr[loc.u], col = .alpha_color("darkred", 150), pch = 16, cex = 2)
        }
        if (length(loc.d) > 0) {
          graphics::rug(loc.d, col = "darkblue", side = 1, lwd = 3, ticksize = 0.02)
          graphics::points(loc.d, data.expr[loc.d], col = .alpha_color("darkblue", 150), pch = 16, cex = 2)
        }
        graphics::legend("bottom", legend = c(
          paste("No. down-regulated in signature:", length(sig.d)),
          paste("No. found in data set:", length(loc.d))
        ), inset = 0.05, bty = "n", text.col = "darkblue")
        graphics::legend("top", legend = c(
          paste("No. up-regulated in signature:", length(sig.u)),
          paste("No. found in data set:", length(loc.u))
        ), inset = 0.05, bty = "n", text.col = "darkred")
      } else {
        loc <- which(gn %in% signature)
        graphics::rug(loc, col = .alpha_color("darkred", 50), side = 3, lwd = 2, ticksize = 0.02)
        graphics::points(loc, data.expr[loc], col = .alpha_color("darkred", 150), pch = 16, cex = 2)

        loc.quart <- stats::quantile(loc)
        graphics::rug(loc.quart, col = "darkblue", side = 3, lwd = 2, ticksize = 0.03)
        y.span <- max(data.expr) - 0.04 * max(data.expr - min(data.expr))
        y.span2 <- max(data.expr) - 0.02 * max(data.expr - min(data.expr))
        graphics::rect(loc.quart[2], y.span, loc.quart[4], max(data.expr), border = "darkblue", lwd = 2, col = NA)
        graphics::rect(loc.quart[1], y.span2, loc.quart[2], y.span2, border = "darkblue", lwd = 2)
        graphics::rect(loc.quart[4], y.span2, loc.quart[5], y.span2, border = "darkblue", lwd = 2)

        graphics::legend("bottom", legend = c(
          paste("No. in signature:", length(signature)),
          paste("No. found in data set (non-redund.):", sum(signature %in% gn)),
          paste("No. found in data set (redundant):", length(loc))
        ), inset = 0.05, bty = "n", text.col = "darkred")
      }
      graphics::legend("right", legend = paste0("NES=", round(score, 3), " (p.adj=", round(fdr, 5), ")"),
                        bty = "n", inset = 0.2, cex = 1.5)
    }
    graphics::par(mfrow = c(1, 1))
  }

  invisible(terms)
}
