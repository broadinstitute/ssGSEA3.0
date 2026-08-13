#' Assemble a complete ssGSEA/PTM-SEA result object
#'
#' Combines the per-gene-set [project_geneset()] outputs into score/p-value/
#' FDR matrices and per-signature annotations, drops gene sets that scored in
#' zero samples, and bundles everything (including gene sets tested and
#' per-sample random-walk detail) into an `ssgsea_result` object -- enough to
#' both write the output GCT files ([write_ssgsea_gcts()]) and drive
#' downstream diagnostics (e.g. [plot_rank_diagnostics()]) without re-reading
#' anything from disk.
#'
#' @param tmp A named list (one entry per gene set, named by gene set) of
#'   [project_geneset()] results.
#' @param gs.names,gs.descs,gs.size,gs.direction Gene set names, descriptions,
#'   sizes, and (optionally) PTMsigDB directionality, as returned by
#'   [match_gene_sets()].
#' @param gs Matched gene sets (member lists), as returned by
#'   [match_gene_sets()] -- retained for downstream diagnostics.
#' @param sample.names Character vector of sample names.
#' @param sample.descs A data frame of per-sample annotations (may have zero
#'   rows/columns).
#' @param output.score.type `"ES"` or `"NES"` -- which per-gene-set vector to
#'   use as the reported score.
#' @param global.fdr If `TRUE`, FDR is computed once across the pooled
#'   p-values from all samples; otherwise FDR is computed per sample column.
#' @param extended.output If `TRUE`, per-signature overlap statistics are
#'   included as additional row annotations.
#' @param data,gene.names The original (untransformed) data matrix and its
#'   row names, retained so [plot_rank_diagnostics()] can redraw the ranked
#'   expression line for a sample without re-reading the input GCT. Ranks are
#'   always derived from *this* (not the internally rank/log-normalized
#'   scoring matrix), matching the original diagnostic plot's behavior and
#'   keeping the y-axis in interpretable expression units.
#' @return An object of class `ssgsea_result`.
#' @keywords internal
build_ssgsea_result <- function(tmp, gs.names, gs.descs, gs.size, gs.direction, gs,
                                 sample.names, sample.descs,
                                 output.score.type = c("NES", "ES"),
                                 global.fdr = FALSE, extended.output = TRUE,
                                 data = NULL, gene.names = NULL) {
  output.score.type <- match.arg(output.score.type)
  N.gs <- length(tmp)

  pval.matrix <- matrix(unlist(lapply(tmp, function(x) x$p.val.vector)), byrow = TRUE, nrow = N.gs)
  score.field <- if (output.score.type == "ES") "ES.vector" else "NES.vector"
  score.matrix <- matrix(unlist(lapply(tmp, function(x) x[[score.field]])), byrow = TRUE, nrow = N.gs)

  ol.matrix <- matrix(unlist(lapply(tmp, function(x) x$OL.vector)), byrow = TRUE, nrow = N.gs)
  ol.numb.matrix <- matrix(unlist(lapply(tmp, function(x) x$OL.numb.vector)), byrow = TRUE, nrow = N.gs)
  ol.perc.matrix <- matrix(unlist(lapply(tmp, function(x) x$OL.perc.vector)), byrow = TRUE, nrow = N.gs)

  random.walk <- lapply(tmp, function(x) x$random.walk)

  ## FDR p-values
  fdr.matrix <- pval.matrix
  if (global.fdr) {
    fdr.matrix <- matrix(stats::p.adjust(unlist(fdr.matrix), method = "fdr"), ncol = ncol(fdr.matrix))
  } else {
    for (i in seq_len(ncol(fdr.matrix))) fdr.matrix[, i] <- stats::p.adjust(fdr.matrix[, i], method = "fdr")
  }

  No.columns.scored <- apply(score.matrix, 1, function(x) sum(!is.na(x)))

  Signature.set.overlap <- data.frame(ol.matrix, stringsAsFactors = FALSE)
  colnames(Signature.set.overlap) <- paste("Signature.set.overlap", sample.names, sep = ".")
  Signature.set.overlap.percent <- data.frame(ol.perc.matrix, stringsAsFactors = FALSE)
  colnames(Signature.set.overlap.percent) <- paste("Signature.set.overlap.percent", sample.names, sep = ".")

  if (extended.output) {
    gs.descs.df <- data.frame(Signature.set.description = gs.descs, Signature.set.size = gs.size,
                               Signature.set.overlap.percent, Signature.set.overlap,
                               No.columns.scored, stringsAsFactors = FALSE)
  } else {
    gs.descs.df <- data.frame(Signature.set.description = gs.descs, Signature.set.size = gs.size,
                               stringsAsFactors = FALSE)
  }

  ## drop gene sets that scored in zero samples
  locs <- which(No.columns.scored > 0)
  if (length(locs) == 0) {
    stop("No gene set scored in any sample -- nothing to report.\n")
  }

  score.matrix <- matrix(score.matrix[locs, ], nrow = length(locs), dimnames = list(gs.names[locs], sample.names))
  pval.matrix <- matrix(pval.matrix[locs, ], nrow = length(locs), dimnames = list(gs.names[locs], sample.names))
  fdr.matrix <- matrix(fdr.matrix[locs, ], nrow = length(locs), dimnames = list(gs.names[locs], sample.names))
  gs.descs.df <- gs.descs.df[locs, , drop = FALSE]
  Signature.set.overlap <- Signature.set.overlap[locs, , drop = FALSE]

  gs.names <- gs.names[locs]
  rownames(gs.descs.df) <- rownames(Signature.set.overlap) <- gs.names
  gs <- gs[locs]
  if (!is.null(gs.direction)) gs.direction <- gs.direction[locs]
  random.walk <- random.walk[locs]

  structure(
    list(
      scores = score.matrix,
      pvalues = pval.matrix,
      fdr = fdr.matrix,
      gs.descs = gs.descs.df,
      overlap = Signature.set.overlap,
      gs.names = gs.names,
      gs = gs,
      gs.direction = gs.direction,
      sample.names = sample.names,
      sample.descs = sample.descs,
      random.walk = random.walk,
      output.score.type = output.score.type,
      data = data,
      gene.names = gene.names
    ),
    class = "ssgsea_result"
  )
}

#' @export
print.ssgsea_result <- function(x, ...) {
  cat("<ssgsea_result>", nrow(x$scores), "gene sets x", length(x$sample.names), "samples\n")
  cat("  score type:", x$output.score.type, "\n")
  invisible(x)
}

#' Write an `ssgsea_result` object out as GCT files
#'
#' Writes the scores/p-values/FDR matrices (each as its own GCT, plus one
#' combined GCT with p-values/FDR as additional row annotations), and,
#' optionally, one GCT per surviving signature containing the original
#' (untransformed) expression values of just that signature's members.
#'
#' @param result An `ssgsea_result`, as returned by [build_ssgsea_result()].
#' @param output.prefix File prefix for output GCT files.
#' @param m.org Original (untransformed) data matrix, genes/sites x samples.
#' @param gene.names Row names of `m.org`.
#' @param gene.descs Row annotation data frame for `m.org` (may have zero
#'   rows/columns).
#' @param gct.src,gct.version `src`/`version` slots to stamp on output GCTs.
#' @param export.signat.gct If `TRUE`, also write one GCT per signature to a
#'   `signature_gct/` subdirectory.
#' @param appenddim See [cmapR::write_gct()].
#' @return (Invisibly) a character vector of the file paths written.
#' @keywords internal
write_ssgsea_gcts <- function(result, output.prefix, m.org, gene.names, gene.descs,
                               gct.src = "", gct.version = "#1.3",
                               export.signat.gct = TRUE, appenddim = FALSE) {
  written <- character(0)
  sample.names <- result$sample.names
  sample.descs <- result$sample.descs
  cdesc.final <- if (nrow(sample.descs) > 0) {
    d <- data.frame(sample.descs, stringsAsFactors = FALSE)
    rownames(d) <- sample.names
    d
  } else {
    data.frame(row.names = sample.names)
  }

  make_gct <- function(mat, rdesc) {
    g <- cmapR::GCT()
    g@mat <- data.matrix(mat)
    g@rid <- rownames(mat)
    g@cid <- sample.names
    g@rdesc <- rdesc
    g@cdesc <- cdesc.final
    g@src <- gct.src
    g
  }

  V.GCT <- make_gct(result$scores, result$gs.descs)
  fn <- paste0(output.prefix, "-scores.gct")
  cmapR::write_gct(V.GCT, ofile = fn, appenddim = appenddim)
  written <- c(written, fn)

  P.GCT <- make_gct(result$pvalues, result$gs.descs)
  fn <- paste0(output.prefix, "-pvalues.gct")
  cmapR::write_gct(P.GCT, ofile = fn, appenddim = appenddim)
  written <- c(written, fn)

  F.GCT <- make_gct(result$fdr, result$gs.descs)
  fn <- paste0(output.prefix, "-fdr-pvalues.gct")
  cmapR::write_gct(F.GCT, ofile = fn, appenddim = appenddim)
  written <- c(written, fn)

  ## combined GCT: scores as the data matrix, pvalues/fdr as extra row annotations
  fdr.tmp <- result$fdr
  colnames(fdr.tmp) <- paste("fdr-pvalue", sample.names, sep = ".")
  pval.tmp <- result$pvalues
  colnames(pval.tmp) <- paste("pvalue", sample.names, sep = ".")
  all.descs <- data.frame(result$gs.descs, pval.tmp, fdr.tmp, stringsAsFactors = FALSE)
  ALL.GCT <- make_gct(result$scores, all.descs)
  ALL.GCT@version <- gct.version
  fn <- paste0(output.prefix, "-combined.gct")
  cmapR::write_gct(ALL.GCT, ofile = fn, appenddim = appenddim)
  written <- c(written, fn)

  if (export.signat.gct) {
    dir.create("signature_gct", showWarnings = FALSE)
    for (sig.name in result$gs.names) {
      signat <- result$overlap[sig.name, ]
      gene.names.tmp <- unique(unlist(strsplit(as.character(signat), "\\|")))
      gene.names.tmp <- gene.names.tmp[gene.names.tmp != ""]

      if (!is.null(result$gs.direction)) {
        gene.names.tmp.direction <- sub("^.*;(u|d)$", "\\1", gene.names.tmp)
        gene.names.tmp <- sub(";u$|;d$", "", gene.names.tmp)
        names(gene.names.tmp.direction) <- gene.names.tmp
      }

      gene.names.tmp.idx <- unique(unlist(lapply(gene.names.tmp, function(x) which(gene.names %in% gene.names.tmp))))
      if (length(gene.names.tmp.idx) == 0) next

      if (nrow(sample.descs) > 0) {
        sample.descs.tmp <- data.frame(
          sample.descs,
          signature.score = as.numeric(result$scores[sig.name, ]),
          signature.pvalue = as.numeric(result$pvalues[sig.name, ]),
          signature.fdr.pvalue = as.numeric(result$fdr[sig.name, ]),
          stringsAsFactors = FALSE
        )
      } else {
        sample.descs.tmp <- data.frame(
          signature.score = as.numeric(result$scores[sig.name, ]),
          signature.pvalue = as.numeric(result$pvalues[sig.name, ]),
          signature.fdr.pvalue = as.numeric(result$fdr[sig.name, ]),
          stringsAsFactors = FALSE
        )
      }
      rownames(sample.descs.tmp) <- sample.names

      if (nrow(gene.descs) > 0) {
        gene.descs.tmp <- data.frame(gene.descs[gene.names.tmp.idx, ])
        if (!is.null(result$gs.direction)) {
          signature.direction <- gene.names.tmp.direction[gene.names[gene.names.tmp.idx]]
          gene.descs.tmp <- data.frame(signature.direction, gene.descs.tmp, stringsAsFactors = FALSE)
        }
      } else if (!is.null(result$gs.direction)) {
        signature.direction <- gene.names.tmp.direction[gene.names[gene.names.tmp.idx]]
        gene.descs.tmp <- data.frame(signature.direction, stringsAsFactors = FALSE)
      } else {
        gene.descs.tmp <- NULL
      }

      gct.tmp <- cmapR::GCT()
      gct.tmp@mat <- data.matrix(m.org[gene.names.tmp.idx, , drop = FALSE])
      gct.tmp@rid <- make.unique(gene.names[gene.names.tmp.idx], sep = "_")
      gct.tmp@cid <- sample.names
      gct.tmp@cdesc <- sample.descs.tmp
      if (!is.null(gene.descs.tmp)) gct.tmp@rdesc <- gene.descs.tmp
      gct.tmp@src <- gct.src

      ## Windows OS supports a maximum of 259 characters in an absolute file
      ## path -- a long signature name (or a deeply nested working directory)
      ## can exceed this and fail to write. Truncate the filename if needed,
      ## leaving a margin of 15 characters for cmapR::write_gct()'s own
      ## appended dimension suffix (e.g. "_n4x50.gct" from appenddim=TRUE).
      max.nchar.file.path <- 259
      fn <- file.path("signature_gct", make.names(sig.name))
      abs.fn <- file.path(getwd(), fn)
      over <- nchar(abs.fn) + 15 - max.nchar.file.path
      if (over > 0) fn <- substr(fn, 1, nchar(fn) - over)
      cmapR::write_gct(gct.tmp, ofile = fn, appenddim = TRUE)
      written <- c(written, fn)
    }
  }

  invisible(written)
}
