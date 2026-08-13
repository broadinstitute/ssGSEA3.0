#' Bootstrap-based p-value and normalized enrichment score
#'
#' Converts an observed enrichment score plus a vector of gene-level bootstrap
#' replicate enrichment scores into a normalized enrichment score (NES) and a
#' two-sided p-value, using a z-score approximation (see Details).
#'
#' The original ssGSEA/GSEA permutation approach shuffles gene labels to build
#' a null distribution, which breaks gene-gene correlation and produces an
#' artificially narrow null whenever `weight != 0` -- inflating significance
#' (even random gene sets that happen to include a few highly-expressed genes
#' get an inflated ES). This bootstrap approach instead resamples genes (rows)
#' with replacement, recenters the resulting bootstrap ES distribution around
#' zero (mimicking a permutation-style null), and computes a z-score of the
#' observed ES against the bootstrap standard deviation (Efron & Tibshirani
#' 1993, ch. 16). A full bootstrap-t correction is not tractable here since
#' estimating each bootstrap replicate's own standard error would require a
#' second level of resampling (O(n^2)); this z-score approximation is used
#' instead (see also MacKinnon 2002; Efron 1982).
#'
#' `boot.es` may contain `NA`s (a replicate that never reached `min.overlap`
#' after its resampling attempts, see [project_geneset()]) -- in keeping with
#' the original implementation, these are *not* removed before computing the
#' mean/SD, so a single failed replicate propagates to `NA` NES/p-value for
#' that sample. In practice this only happens for gene sets whose overlap is
#' very close to `min.overlap`.
#'
#' @param obs.es Observed enrichment score.
#' @param boot.es Numeric vector of bootstrap replicate enrichment scores.
#' @return A list with `NES` and `p.value`.
#' @keywords internal
bootstrap_es_pvalue <- function(obs.es, boot.es) {
  B <- length(boot.es)
  if (B == 0) return(list(NES = obs.es, p.value = 1))

  ## The naive two-sided bootstrap p-value formula is:
  ##   p = (1/B) * sum_i I( |ES*_i - ES_hat| >= |ES_hat - 0| )
  ## where ES*_i is the i-th bootstrap replicate's ES, ES_hat is the
  ## observed ES, and I(.) is the indicator function. This is biased
  ## (under-dispersed, too many p ~ 0) because it (i) assumes
  ## mean(ES*_i) == ES_hat and (ii) treats the null's SD as identical to
  ## the bootstrap distribution's SD. The gold-standard fix is a
  ## "studentized"/bootstrap-t correction, which isn't tractable here since
  ## it would need the standard error of each ES*_i, i.e. a second level of
  ## bootstrapping (O(n^2)). Instead we recenter the bootstrap ES
  ## distribution around zero (mimicking a permutation-style null) and
  ## compute a z-score of the observed ES against the bootstrap SD -- akin
  ## to a z-test with the SD estimated via bootstrap sampling.
  boot.es <- boot.es - mean(boot.es)
  se_boot <- stats::sd(boot.es)
  z_score <- obs.es / se_boot

  if (obs.es >= 0) {
    ## NES: same normalization convention as permutation-based ssGSEA (ES
    ## divided by the mean of same-sign null values), so NES stays
    ## comparable in magnitude/interpretation to the original method's.
    pos.phi <- boot.es[boot.es >= 0]
    pos.m <- if (length(pos.phi) == 0) 0.5 else mean(pos.phi)
    NES <- obs.es / pos.m
    p.value <- 2 * stats::pnorm(z_score, lower.tail = FALSE)
  } else {
    ## NES: same normalization convention as permutation-based ssGSEA (mean
    ## of same-sign bootstrap ES), mirrored for the negative-ES case.
    neg.phi <- boot.es[boot.es < 0]
    neg.m <- if (length(neg.phi) == 0) 0.5 else mean(neg.phi)
    NES <- obs.es / abs(neg.m)
    p.value <- 2 * stats::pnorm(z_score, lower.tail = TRUE)
  }

  list(NES = NES, p.value = p.value)
}
