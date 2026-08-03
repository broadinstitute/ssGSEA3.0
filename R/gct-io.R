#' Safely parse a GCT file, working around non-unique row IDs
#'
#' `cmapR::parse_gctx()` refuses to parse a GCT file whose `rid` column is
#' not unique. This wrapper falls back to reading the file as plain text,
#' de-duplicating `rid` with [make.unique()], writing a corrected copy next
#' to the original, and re-parsing that. This exact fallback used to be
#' duplicated between the preprocessing and enrichment code paths; it now
#' lives here once.
#'
#' @param gct.str Path to a GCT (v1.2/v1.3) file.
#' @return A [cmapR::GCT-class] object.
#' @keywords internal
read_gct_safe <- function(gct.str) {
  if (!file.exists(gct.str)) stop(glue::glue("File '{gct.str}' not found!"))

  gct <- try(cmapR::parse_gctx(gct.str), silent = TRUE)
  if (!methods::is(gct, "try-error")) return(gct)

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
  cmapR::parse_gctx(fname = gct.unique)
}
