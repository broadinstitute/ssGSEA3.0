#' Parse a single-data-row GCT file by hand
#'
#' Works around a bug in `cmapR::parse_gctx()` (confirmed against cmapR
#' 1.18.0): when a `.gct` file has exactly one data row (`nrmat == 1`) and
#' more than one row-descriptor column (`nrhd > 1`), its internal `initialize`
#' method builds `m[, 2:(nrhd + 1)]` from a matrix that itself has only one
#' row -- since one of the two dimensions is 1, R's default `drop = TRUE`
#' collapses this to a plain vector, and `as.matrix()` on that vector then
#' builds an `nrhd x 1` column instead of the intended `1 x nrhd` row. The
#' following `dimnames(rdesc) <- list(rid, rhd)` then fails because `rid`
#' (length 1) doesn't match `nrow(rdesc)` (which is actually `nrhd`). Notably,
#' cmapR already has a mirror-image fix for the analogous column-descriptor
#' case (`if (nchd == 1) cdesc <- t(cdesc)`) -- this row case was just missed.
#'
#' This only handles the plain-text GCT v1.2/v1.3 format (not `.gctx`/HDF5,
#' which isn't affected) and is only ever invoked as a fallback once
#' `cmapR::parse_gctx()` has already failed with this specific error.
#'
#' @param gct.str Path to a GCT file with exactly one data row.
#' @return A [cmapR::GCT-class] object.
#' @keywords internal
.parse_single_row_gct <- function(gct.str) {
  lines <- readLines(gct.str)
  version <- lines[1]
  meta <- as.numeric(strsplit(lines[2], "\t")[[1]])
  nrhd <- meta[3]
  nchd <- meta[4]

  header <- strsplit(lines[3], "\t")[[1]]
  rhd <- if (nrhd > 0) header[2:(nrhd + 1)] else character(0)
  cid <- header[(nrhd + 2):length(header)]

  data.line <- strsplit(lines[3 + nchd + 1], "\t")[[1]]
  rid <- data.line[1]
  mat.vals <- as.numeric(data.line[(nrhd + 2):length(data.line)])
  mat <- matrix(mat.vals, nrow = 1, dimnames = list(rid, cid))

  if (nrhd > 0) {
    rdesc <- as.data.frame(as.list(data.line[2:(nrhd + 1)]), stringsAsFactors = FALSE)
    names(rdesc) <- rhd
    rdesc <- suppressWarnings(cmapR::fix.datatypes(rdesc))
  } else {
    rdesc <- data.frame(row.names = rid)
  }
  rownames(rdesc) <- rid
  rdesc$id <- rid

  g <- cmapR::GCT()
  g@mat <- mat
  g@rid <- rid
  g@cid <- cid
  g@rdesc <- rdesc
  g@cdesc <- data.frame(id = cid, row.names = cid)
  g@version <- version
  g@src <- gct.str
  g
}

#' Safely parse a GCT file, working around cmapR quirks
#'
#' `cmapR::parse_gctx()` refuses to parse a GCT file whose `rid` column is
#' not unique, and (separately, see [.parse_single_row_gct()]) mis-parses one
#' with exactly one data row and more than one row-descriptor column. This
#' wrapper falls back to reading the file as plain text and handling either
#' case by hand. The non-unique-`rid` fallback used to be duplicated between
#' the preprocessing and enrichment code paths; it now lives here once.
#'
#' @param gct.str Path to a GCT (v1.2/v1.3) file.
#' @return A [cmapR::GCT-class] object.
#' @keywords internal
read_gct_safe <- function(gct.str) {
  if (!file.exists(gct.str)) stop(glue::glue("File '{gct.str}' not found!"))

  gct <- try(cmapR::parse_gctx(gct.str), silent = TRUE)
  if (!methods::is(gct, "try-error")) return(gct)

  if (length(grep("dimnames.*array extent", gct)) > 0) {
    return(.parse_single_row_gct(gct.str))
  }

  if (length(grep("rid must be unique", gct)) == 0) {
    ## some other parse failure -- surface the original error rather than
    ## masking it behind a misleading "made ids unique" message
    stop(glue::glue(
      "\n\nError importing GCT file using cmapR::parse_gctx(). ",
      "Possible reasons: (1) an outdated cmapR version -- some releases ",
      "renamed 'parse.gctx()'; (2) the file isn't valid GCT format, see ",
      "https://clue.io/connectopedia/gct_format.\n\n",
      "Underlying error:\n{gct}"
    ))
  }

  warning("rids not unique! Making ids unique and re-exporting the GCT file...")
  gct.tmp <- readLines(gct.str)
  rid <- sub("\t.*", "", gct.tmp)
  meta <- as.numeric(unlist(strsplit(gct.tmp[2], "\t")))
  rid.idx <- (meta[4] + 3):length(rid)

  if (length(rid[rid.idx]) == length(unique(rid[rid.idx])))
    stop("parse.gctx() reported non-unique rids, but none were found on manual inspection -- please file a bug report.")

  rid[rid.idx] <- make.unique(rid[rid.idx], sep = "_")
  rest <- sub(".*?\t", "", gct.tmp)
  rest[1] <- ""
  gct.tmp2 <- paste(rid, rest, sep = "\t")
  gct.tmp2[1] <- sub("\t.*", "", gct.tmp2[1])

  gct.unique <- sub("\\.gct$", "_unique.gct", gct.str)
  writeLines(gct.tmp2, con = gct.unique)
  result <- cmapR::parse_gctx(fname = gct.unique)
  ## flag that `rid` was de-duplicated with a "_N" suffix, so callers that
  ## need the original (pre-dedup) identifiers -- e.g. for gene set matching
  ## -- know to strip it back off (see ssGSEA3())
  attr(result, "deduped") <- TRUE
  result
}
