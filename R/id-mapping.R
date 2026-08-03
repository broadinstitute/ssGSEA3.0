#' Map RefSeq accessions or gene symbols to UniProt accessions
#'
#' Uses the appropriate Bioconductor `org.*.eg.db` annotation package for the
#' given `organism` to map identifiers to UniProt accessions.
#'
#' This replaces the previous auto-detection behavior (which sampled a few
#' IDs and test-queried each of the human/mouse/rat/zebrafish annotation
#' packages in turn to guess the organism). Auto-detection meant every call
#' required all four `org.*.eg.db` packages to be installed, and could
#' silently guess wrong on a small or ambiguous ID sample. Callers must now
#' state the organism explicitly, and only the one relevant annotation
#' package needs to be installed.
#'
#' @param ids Character vector of accessions to map (RefSeq or gene symbols,
#'   per `keytype`). May contain version suffixes/site annotations after a
#'   `.` or `;`, which are stripped before querying.
#' @param keytype One of `"REFSEQ"` or `"SYMBOL"`, matching the type of `ids`.
#' @param organism One of `"human"`, `"mouse"`, `"rat"`, `"zebrafish"`.
#' @return A list with elements `keytype`, `organism`, and `id.map` (a
#'   data.frame with columns `id`, `id.query`, `id.mapped`, `id.concat`).
#'   Entries with no match have `id.mapped = NA`.
#' @export
RefSeq2UniProt <- function(ids,
                            keytype = c("REFSEQ", "SYMBOL"),
                            organism = c("human", "mouse", "rat", "zebrafish")) {
  keytype <- match.arg(keytype)
  organism <- match.arg(organism)

  orgdb_pkg <- switch(organism,
    human = "org.Hs.eg.db", mouse = "org.Mm.eg.db",
    rat = "org.Rn.eg.db", zebrafish = "org.Dr.eg.db"
  )

  if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
    stop(glue::glue(
      "Package '{orgdb_pkg}' is required to map organism = '{organism}' IDs ",
      "to UniProt, but is not installed. Install it with:\n",
      "  BiocManager::install('{orgdb_pkg}')"
    ))
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Package 'AnnotationDbi' is required for ID mapping but is not installed.")

  orgdb <- getExportedValue(orgdb_pkg, orgdb_pkg)

  ## strip version suffixes / site annotations, e.g. "NP_001.2;S100" -> "NP_001"
  id.query <- sub("(\\.|;).*", "", ids)
  names(id.query) <- ids

  id.mapped <- AnnotationDbi::mapIds(orgdb,
    keys = unname(id.query), column = "UNIPROT",
    keytype = keytype, multiVals = "first"
  )

  id.map <- data.frame(
    id = names(id.query), id.query = id.query,
    id.mapped = unname(id.mapped), id.concat = paste(ids, id.mapped, sep = "_"),
    stringsAsFactors = FALSE
  )

  list(keytype = keytype, organism = organism, id.map = id.map)
}
