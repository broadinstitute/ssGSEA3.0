#' Combine random-walk statistics into a final enrichment score
#'
#' @param max.ES Maximum running enrichment score.
#' @param min.ES Minimum running enrichment score.
#' @param RES Running enrichment score vector (values at gene-set hit positions).
#' @param gaps Gaps between consecutive gene-set hit positions.
#' @param valleys Running score just before each hit.
#' @param statistic `"Kolmogorov-Smirnov"` (peak deviation) or
#'   `"area.under.RES"` (signed area under the running score).
#' @return A list with `RES`, `ES`, and `arg.ES` (index of the extremum).
#' @keywords internal
.enrichment_statistic <- function(max.ES, min.ES, RES, gaps, valleys, statistic) {
  if (statistic == "Kolmogorov-Smirnov") {
    if (max.ES > -min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
  }
  if (statistic == "area.under.RES") {
    if (max.ES > -min.ES) {
      arg.ES <- which.max(RES)
    } else {
      arg.ES <- which.min(RES)
    }
    gaps <- gaps + 1
    RES <- c(valleys, 0) * gaps + 0.5 * (c(0, RES) - c(valleys, 0)) * gaps
    ES <- sum(RES)
  }
  list(RES = RES, ES = ES, arg.ES = arg.ES)
}

#' Calculate a (ss)GSEA enrichment score for one ranked gene list
#'
#' Computes the running/final enrichment score of `gene.set2` (locations of a
#' gene set's members within `ordered.gene.list`) against one sample's ranked
#' expression values, optionally splitting the calculation into up/down
#' components when `gene.set.direction` (PTMsigDB-style) is supplied.
#'
#' @param ordered.gene.list Row indices, ordered by decreasing expression.
#' @param gene.set2 Row indices (into the same space as `ordered.gene.list`)
#'   of the current gene set's members.
#' @param weight Exponent applied to (absolute) expression values when
#'   weighting steps of the running-sum statistic; `0` gives every hit equal
#'   weight.
#' @param n.rows Number of rows (genes/sites) in the current sample's data.
#' @param correl.type `"rank"`, `"z.score"`, or `"symm.rank"`.
#' @param gene.set.direction `NULL`, or a character vector of `"u"`/`"d"`
#'   aligned with the (pre-directionality-split) gene set.
#' @param data.expr Current sample's expression values, in original
#'   (unordered) row order.
#' @param statistic `"Kolmogorov-Smirnov"` or `"area.under.RES"`.
#' @param min.overlap Minimum number of up- or down-regulated members
#'   required before a lone (size-1) up/down component is allowed to
#'   contribute to a directional (PTMsigDB-style) score; see Details.
#' @return A list of enrichment-score results (`ES`, `RES`, `indicator`,
#'   `correl.vector`, `step.up`, `step.down`, and -- when directional --
#'   `ES.all`/`arg.ES`).
#' @keywords internal
gsea_score <- function(ordered.gene.list, gene.set2, weight, n.rows, correl.type,
                        gene.set.direction, data.expr, statistic = "area.under.RES",
                        min.overlap = 1) {
  ## #######################################
  ## weighting
  ## #######################################
  if (weight == 0) {
    correl.vector <- rep(1, n.rows)
  } else if (weight > 0) {
    ## NOTE: the original permutation-based code shuffled gene labels here,
    ## which breaks gene-gene correlation and understates the null's spread
    ## when weight != 0 -- see bootstrap_es_pvalue() for how this is corrected.
    if (correl.type == "rank") {
      correl.vector <- data.expr
    } else if (correl.type == "symm.rank") {
      correl.vector <- data.expr
      correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows / 2)],
                               correl.vector,
                               correl.vector + correl.vector - correl.vector[ceiling(n.rows / 2)])
    } else if (correl.type == "z.score") {
      x <- data.expr
      correl.vector <- (x - mean(x)) / stats::sd(x)
    }
  }

  N <- length(ordered.gene.list)

  if (!is.null(gene.set.direction)) {
    d.idx <- which(gene.set.direction == "d")
    Nh.d <- length(d.idx)
    Nm.d <- N - Nh.d

    u.idx <- which(gene.set.direction == "u")
    Nh.u <- length(u.idx)
    Nm.u <- N - Nh.u

    ## up-regulated part
    if (Nh.u > 1) {
      tag.u <- sign(match(ordered.gene.list, gene.set2[u.idx], nomatch = 0))
      ind.u <- which(tag.u == 1)

      correl.vector.u <- correl.vector[ind.u]
      correl.vector.u <- abs(correl.vector.u)^weight

      sum.correl.u <- sum(correl.vector.u)

      up.u <- correl.vector.u / sum.correl.u
      gaps.u <- (c(ind.u - 1, N) - c(0, ind.u))
      down.u <- gaps.u / Nm.u

      RES.u <- cumsum(c(up.u, up.u[Nh.u]) - down.u)
      valleys.u <- RES.u[1:Nh.u] - up.u

      max.ES.u <- max(RES.u)
      min.ES.u <- min(valleys.u)

      score.res <- .enrichment_statistic(max.ES.u, min.ES.u, RES.u[1:Nh.u], gaps.u, valleys.u, statistic)
      ES.u <- score.res$ES
      arg.ES.u <- score.res$arg.ES
      RES.u <- score.res$RES
    } else {
      correl.vector.u <- rep(0, N)
      ES.u <- 0
      RES.u <- 0
      ind.u <- NULL
      arg.ES.u <- NA
      up.u <- 0
      down.u <- 0
    }

    ## down-regulated part
    if (Nh.d > 1) {
      tag.d <- sign(match(ordered.gene.list, gene.set2[d.idx], nomatch = 0))
      ind.d <- which(tag.d == 1)

      correl.vector.d <- correl.vector[ind.d]
      correl.vector.d <- abs(correl.vector.d)^weight

      sum.correl.d <- sum(correl.vector.d)

      up.d <- correl.vector.d / sum.correl.d
      gaps.d <- (c(ind.d - 1, N) - c(0, ind.d))
      down.d <- gaps.d / Nm.d

      RES.d <- cumsum(c(up.d, up.d[Nh.d]) - down.d)
      valleys.d <- RES.d[1:Nh.d] - up.d

      max.ES.d <- max(RES.d)
      min.ES.d <- min(valleys.d)

      score.res <- .enrichment_statistic(max.ES.d, min.ES.d, RES.d[1:Nh.d], gaps.d, valleys.d, statistic)
      ES.d <- score.res$ES
      arg.ES.d <- score.res$arg.ES
      RES.d <- score.res$RES
    } else {
      correl.vector.d <- rep(0, N)
      ES.d <- 0
      RES.d <- 0
      ind.d <- NULL
      arg.ES.d <- NA
      up.d <- 0
      down.d <- 0
    }

    ## if only one of up/down has a single member, and the other side
    ## doesn't meet min.overlap on its own, drop both sides to zero
    ## (arg.ES.u/arg.ES.d are deliberately left untouched here -- they're
    ## fully overwritten by `arg.ES <- c(arg.ES.u, arg.ES.d)` below)
    if (Nh.d == 1 & Nh.u < min.overlap | Nh.u == 1 & Nh.d < min.overlap) {
      ES.u <- ES.d <- RES.u <- RES.d <- 0
      ind.u <- ind.d <- NULL
    }

    ES <- ES.u - ES.d
    RES <- list(u = RES.u, d = RES.d)
    arg.ES <- c(arg.ES.u, arg.ES.d)
    correl.vector <- list(u = correl.vector.u, d = correl.vector.d)
    ind <- list(u = ind.u, d = ind.d)
    step.up <- list(u = up.u, d = up.d)
    step.down <- list(u = 1 / Nm.u, d = 1 / Nm.d)

    return(list(ES = ES, ES.all = list(u = ES.u, d = ES.d), arg.ES = arg.ES, RES = RES,
                indicator = ind, correl.vector = correl.vector, step.up = step.up, step.down = step.down))
  }

  ## ##############################################################
  ## original ssGSEA code without directionality
  ## ##############################################################
  Nh <- length(gene.set2)
  Nm <- N - Nh

  tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch = 0))
  ind <- which(tag.indicator == 1)
  correl.vector <- abs(correl.vector[ind])^weight
  sum.correl <- sum(correl.vector)

  up <- correl.vector / sum.correl
  gaps <- (c(ind - 1, N) - c(0, ind))
  down <- gaps / Nm

  RES <- cumsum(c(up, up[Nh]) - down)
  valleys <- RES[1:Nh] - up

  max.ES <- max(RES)
  min.ES <- min(valleys)

  score.res <- .enrichment_statistic(max.ES, min.ES, RES[1:Nh], gaps, valleys, statistic)

  list(ES = score.res$ES, arg.ES = score.res$arg.ES, RES = score.res$RES,
       indicator = ind, correl.vector = correl.vector, step.up = up, step.down = 1 / Nm)
}

#' Map gene set members to row positions, keeping direction labels aligned
#'
#' `which(gene.names == x)` for each member of `gene.set` can return zero
#' matches (member absent from `gene.names`) or more than one (redundant/
#' multi-probe `gene.names`) -- naively `unlist()`-ing across members then
#' shrinks or grows the result out of step with `gene.set.direction` (one
#' label per *original* member). This keeps the two aligned by expanding
#' each member's direction label to match however many positions it actually
#' occupies (zero, one, or many).
#'
#' Without this, [gsea_score()]'s directional branch can index past the end
#' of a too-short position vector (returning `NA`s that later blow up a
#' `min()`/`max()` call) or silently pair positions with the wrong member's
#' direction. This is what made a resampled subset of `gene.names` (as in
#' [.resample_bootstrap_es()]) unsafe with the original, un-aligned
#' construction -- gene-level bootstrap resampling can drop individual gene
#' set members from a given replicate, which ordinary permutation (which
#' never drops a gene, only reshuffles all of them) cannot.
#'
#' @param gene.set Character vector of gene set member names.
#' @param gene.set.direction `NULL`, or a character vector of `"u"`/`"d"`
#'   labels parallel to `gene.set`.
#' @param gene.names Row (gene/site) names to match `gene.set` against.
#' @return A list with `gene.set2` (matched row positions) and
#'   `gene.set.direction` (`NULL`, or labels aligned with `gene.set2`).
#' @keywords internal
.align_geneset_to_direction <- function(gene.set, gene.set.direction, gene.names) {
  positions <- lapply(gene.set, function(x) which(gene.names == x))
  gene.set2 <- as.numeric(unlist(positions))
  if (is.null(gene.set.direction)) {
    return(list(gene.set2 = gene.set2, gene.set.direction = NULL))
  }
  list(gene.set2 = gene.set2, gene.set.direction = rep(gene.set.direction, lengths(positions)))
}

#' Gene-level bootstrap resampling for one sample/gene-set pair
#'
#' Resamples genes (rows) with replacement `nperm` times, recomputing the
#' enrichment score each time to build the null distribution used by
#' [bootstrap_es_pvalue()]. Each replicate retries up to 200 times if a given
#' resample doesn't reach `min.overlap`; a replicate that never does becomes
#' `NA` (see [bootstrap_es_pvalue()] for how that propagates).
#'
#' @inheritParams gsea_score
#' @param gene.names Row (gene/site) names for the current sample, in the
#'   same order as `data.expr`.
#' @param gene.set Character vector of gene set member names (not indices).
#' @param nperm Number of bootstrap replicates.
#' @param min.overlap Minimum required overlap for a resample to count.
#' @return Numeric vector of length `nperm` (bootstrap replicate ES values).
#' @keywords internal
.resample_bootstrap_es <- function(gene.names, gene.set, data.expr, weight, n.rows,
                                    correl.type, gene.set.direction, min.overlap,
                                    nperm, statistic) {
  vapply(seq_len(nperm), function(...) {
    for (attempt in seq_len(200)) {
      boot.ix <- sample.int(n.rows, n.rows, replace = TRUE)
      gene.names.b <- gene.names[boot.ix]
      if (sum(gene.names.b %in% gene.set) < min.overlap) next
      expr.boot <- data.expr[boot.ix]
      aligned <- .align_geneset_to_direction(gene.set, gene.set.direction, gene.names.b)
      gene.list.boot <- order(expr.boot, decreasing = TRUE)
      return(gsea_score(gene.list.boot, aligned$gene.set2, weight, n.rows, correl.type,
                         aligned$gene.set.direction, expr.boot, statistic, min.overlap)$ES)
    }
    NA_real_
  }, numeric(1))
}

#' Score one gene set against every sample in a dataset
#'
#' The per-gene-set unit of work: for each sample column, handles missing
#' values, checks minimum overlap, computes the observed enrichment score,
#' and (when `nperm > 0`) runs the gene-level bootstrap to obtain an NES and
#' p-value via [bootstrap_es_pvalue()].
#'
#' @inheritParams gsea_score
#' @param data.array Numeric matrix (genes/sites x samples).
#' @param gene.names Row (gene/site) names, aligned with `data.array`.
#' @param n.cols Number of samples (columns of `data.array`).
#' @param gene.set Character vector of gene set member names.
#' @param nperm Number of bootstrap replicates; `0` skips bootstrapping
#'   entirely (`ES` and `NES` are then identical, and `p.value` is `1`).
#' @param min.overlap Minimum overlap (in members) required to score a
#'   sample.
#' @param size.G.current Original (pre-overlap) size of this gene set, used
#'   only to report overlap percentage.
#' @return A list with `ES.vector`, `NES.vector`, `p.val.vector`,
#'   `OL.vector`/`OL.numb.vector`/`OL.perc.vector` (overlap diagnostics),
#'   `gene.set.size`, and `random.walk` (per-sample [gsea_score()] output, for
#'   the rank-plot diagnostic).
#' @keywords internal
project_geneset <- function(data.array, gene.names, n.cols, n.rows, weight = 0,
                             statistic = "Kolmogorov-Smirnov", gene.set, nperm = 200,
                             correl.type = "rank", gene.set.direction = NULL,
                             min.overlap, size.G.current) {
  ES.vector <- NES.vector <- p.val.vector <- rep(NA, n.cols)
  OL.vector <- OL.numb.vector <- OL.perc.vector <- rep(NA, n.cols)
  random.walk <- vector("list", n.cols)

  gene.names.all <- gene.names
  gene.set.all <- gene.set
  gene.set.size <- length(gene.set)
  gene.set.direction.all <- gene.set.direction
  n.rows.all <- n.rows

  for (sample.index in seq_len(n.cols)) {
    gene.names <- gene.names.all
    gene.set <- gene.set.all
    gene.set.direction <- gene.set.direction.all
    n.rows <- n.rows.all

    data.expr <- data.array[, sample.index]

    na.idx <- which(is.na(data.expr) | is.infinite(data.expr))
    if (length(na.idx) > 0) {
      data.expr <- data.expr[-na.idx]
      gene.names <- gene.names[-na.idx]
      n.rows <- length(data.expr)

      gene.set.no.na <- which(gene.set %in% gene.names)
      gene.set <- gene.set[gene.set.no.na]
      gene.set.direction <- gene.set.direction[gene.set.no.na]
    }

    if (sum(gene.names %in% gene.set) < min.overlap) {
      random.walk[[sample.index]] <- NA
      OL.numb.vector[sample.index] <- sum(gene.names %in% gene.set)

      if (OL.numb.vector[sample.index] > 0) {
        if (is.null(gene.set.direction)) {
          OL.vector[sample.index] <- paste(intersect(gene.names, gene.set), collapse = "|")
        } else {
          OL.vector[sample.index] <- paste(paste(gene.set, gene.set.direction, sep = ";")[
            stats::na.omit(match(gene.names, gene.set))
          ], collapse = "|")
        }
      }
      OL.perc.vector[sample.index] <- round(100 * (OL.numb.vector[sample.index] / size.G.current), 1)
    } else {
      aligned <- .align_geneset_to_direction(gene.set, gene.set.direction, gene.names)
      gene.set2 <- aligned$gene.set2
      gene.list <- order(data.expr, decreasing = TRUE)

      GSEA.results <- gsea_score(gene.list, gene.set2, weight, n.rows, correl.type,
                                  aligned$gene.set.direction, data.expr, statistic, min.overlap)
      ES.vector[sample.index] <- GSEA.results$ES

      if (is.null(gene.set.direction)) {
        OL.vector[sample.index] <- paste(unique(gene.names[gene.list][GSEA.results$indicator]), collapse = "|")
        OL.numb.vector[sample.index] <- length(unique(gene.names[gene.list][GSEA.results$indicator]))
        OL.perc.vector[sample.index] <- round(100 * (OL.numb.vector[sample.index] / size.G.current), 1)
      } else {
        ol.tmp <- c()
        if (length(GSEA.results$indicator$u) > 0)
          ol.tmp <- c(ol.tmp, paste(unique(gene.names[gene.list][GSEA.results$indicator$u]), "u", sep = ";"))
        if (length(GSEA.results$indicator$d) > 0)
          ol.tmp <- c(ol.tmp, paste(unique(gene.names[gene.list][GSEA.results$indicator$d]), "d", sep = ";"))

        OL.vector[sample.index] <- paste(ol.tmp, collapse = "|")
        OL.numb.vector[sample.index] <- length(ol.tmp)
        OL.perc.vector[sample.index] <- round(100 * (OL.numb.vector[sample.index] / size.G.current), 1)
      }

      random.walk[[sample.index]] <- GSEA.results

      if (nperm == 0) {
        NES.vector[sample.index] <- ES.vector[sample.index]
        p.val.vector[sample.index] <- 1
      } else {
        boot.es <- .resample_bootstrap_es(gene.names, gene.set, data.expr, weight, n.rows,
                                           correl.type, gene.set.direction, min.overlap,
                                           nperm, statistic)
        boot.res <- bootstrap_es_pvalue(ES.vector[sample.index], boot.es)
        NES.vector[sample.index] <- boot.res$NES
        p.val.vector[sample.index] <- boot.res$p.value
      }
    }
  }

  list(ES.vector = ES.vector, NES.vector = NES.vector, p.val.vector = p.val.vector,
       OL.vector = OL.vector, OL.numb.vector = OL.numb.vector, OL.perc.vector = OL.perc.vector,
       gene.set.size = gene.set.size, random.walk = random.walk)
}
