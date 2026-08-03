#' Normalize a dataset's samples before ssGSEA scoring
#'
#' @param m Numeric matrix (genes/sites x samples).
#' @param sample.norm.type `"rank"` (rank-transform, scaled to 0-10000,
#'   each sample independently), `"log.rank"` (same, then log-transformed),
#'   `"log"` (values below 1 are clamped to 1, then log-transformed), or
#'   `"none"` (unchanged). Missing values are preserved (excluded from
#'   ranking, not imputed) in all cases.
#' @return `m`, normalized.
#' @keywords internal
normalize_samples <- function(m, sample.norm.type = c("rank", "log", "log.rank", "none")) {
  sample.norm.type <- match.arg(sample.norm.type)
  Ng <- nrow(m)

  rank_cols <- function(m) {
    apply(m, 2, function(x) {
      x.rank <- rep(NA, length(x))
      keep.idx <- which(!is.na(x))
      x.rank[keep.idx] <- rank(x[keep.idx], ties.method = "average")
      x.rank
    })
  }

  if (sample.norm.type == "rank") {
    m <- rank_cols(m)
    m <- 10000 * m / Ng
  } else if (sample.norm.type == "log.rank") {
    m <- rank_cols(m)
    m <- log(10000 * m / Ng + exp(1))
  } else if (sample.norm.type == "log") {
    m[m < 1] <- 1
    m <- log(m + exp(1))
  }
  m
}
