#' Import one gene set database (GMT file)
#'
#' Reads a single GMT file, filters gene sets by size, and de-duplicates
#' members within each gene set.
#'
#' @param gs.db Path to a GMT file.
#' @param thres.min Minimum gene set size to keep.
#' @param thres.max Maximum gene set size to keep.
#' @return A list with `gs` (named list of member vectors), `gs.names`,
#'   `gs.desc`, and `size.G` (post-dedup set sizes).
#' @keywords internal
.read_gmt <- function(gs.db, thres.min = 2, thres.max = 2000) {
  temp <- readLines(gs.db)
  temp <- strsplit(temp, "\t")
  temp.size.G <- sapply(temp, function(x) length(x) - 2)

  rm.idx <- which(temp.size.G < thres.min | temp.size.G > thres.max)
  if (length(rm.idx) > 0) {
    temp <- temp[-rm.idx]
    temp.size.G <- temp.size.G[-rm.idx]
  }

  gs <- lapply(temp, function(x) x[3:length(x)])
  gs.names <- sapply(temp, function(x) x[1])
  gs.desc <- sapply(temp, function(x) x[2])

  ## de-duplicate members within each gene set
  gs.unique <- lapply(gs, unique)
  gs.unique.size.G <- sapply(gs.unique, length)
  gs.not.unique.idx <- which(gs.unique.size.G < temp.size.G)
  if (length(gs.not.unique.idx) > 0) {
    warning(
      "Duplicated gene set members detected. Removing redundant members from:\n",
      paste(gs.names[gs.not.unique.idx], collapse = "\n")
    )
    gs <- gs.unique
    temp.size.G <- gs.unique.size.G
  }
  size.G <- temp.size.G
  names(gs) <- names(gs.names) <- names(gs.desc) <- names(size.G) <- gs.names

  list(gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G)
}

#' Import and merge one or more gene set databases
#'
#' Reads one or more GMT files, filters each by set size, and merges them into
#' a single collection. If any gene set member carries a PTMsigDB `;u`/`;d`
#' directionality suffix, the whole collection is treated as PTMsigDB-style:
#' directionality is split out into `gs.direction` and stripped from members.
#'
#' @param gene.set.databases Character vector of paths to GMT files.
#' @param thres.min Minimum gene set size to keep.
#' @param thres.max Maximum gene set size to keep.
#' @return A list with `gs`, `gs.names`, `gs.descs`, `size.G`, and
#'   `gs.direction` (`NULL` for ordinary, non-directional gene sets).
#' @export
read_gene_sets <- function(gene.set.databases, thres.min = 2, thres.max = 2000) {
  GSDB <- lapply(gene.set.databases, .read_gmt, thres.min = thres.min, thres.max = thres.max)

  gs <- GSDB[[1]]$gs
  gs.names <- GSDB[[1]]$gs.names
  gs.descs <- GSDB[[1]]$gs.desc
  size.G <- GSDB[[1]]$size.G
  if (length(GSDB) > 1) {
    for (i in 2:length(GSDB)) {
      gs <- append(gs, GSDB[[i]]$gs)
      gs.names <- append(gs.names, GSDB[[i]]$gs.names)
      gs.descs <- append(gs.descs, GSDB[[i]]$gs.desc)
      size.G <- append(size.G, GSDB[[i]]$size.G)
    }
  }

  ## PTMsigDB-style directionality (';u'/';d' suffix on members)
  if (length(grep(";u$|;d$", gs[[1]])) > 0) {
    gs.direction <- lapply(gs, function(x) sub("^.*;(d|u)$", "\\1", x))
    gs <- lapply(gs, function(x) sub(";u$|;d$", "", x))
  } else {
    gs.direction <- NULL
  }

  list(gs = gs, gs.names = gs.names, gs.descs = gs.descs, size.G = size.G, gs.direction = gs.direction)
}

#' Match imported gene sets against a dataset's row (gene/site) names
#'
#' Intersects each gene set's members with `gene.names`, drops gene sets that
#' fall below `min.overlap`, and collapses duplicate gene set names (keeping
#' the first occurrence).
#'
#' @param gsdb A gene set collection as returned by [read_gene_sets()].
#' @param gene.names Row (gene/site) names of the dataset being tested.
#' @param min.overlap Minimum number of overlapping members required to keep a
#'   gene set.
#' @return `gsdb`, subset/filtered in place, plus `size.ol.G` (post-overlap
#'   set sizes).
#' @keywords internal
match_gene_sets <- function(gsdb, gene.names, min.overlap) {
  gs <- gsdb$gs
  gs.direction <- gsdb$gs.direction

  if (!is.null(gs.direction)) {
    ## filter members AND their parallel direction labels together, so the
    ## two stay index-aligned (a mismatch here previously caused gs.direction
    ## to retain stale, pre-filter length/content)
    for (nm in names(gs)) {
      keep <- gs[[nm]] %in% gene.names
      gs[[nm]] <- gs[[nm]][keep]
      gs.direction[[nm]] <- gs.direction[[nm]][keep]
    }
  } else {
    gs <- lapply(gs, function(x) intersect(x, gene.names))
  }
  size.ol.G <- sapply(gs, length)

  if (length(which(size.ol.G > 0)) == 0) {
    stop("No overlap to any gene sets found in your data set!\nPossible reasons: 1) organism other than human; 2) wrong format of gene names/site ids!\n")
  }

  keep.idx <- which(size.ol.G >= min.overlap)
  if (length(keep.idx) == 0) {
    ## NOTE: this message's wording is checked elsewhere (e.g. calling CLI code)
    ## to identify/tolerate minimum-overlap errors -- keep it stable if edited.
    stop("Dataset does not meet minimum-overlap (", min.overlap, ") to any gene sets in the provided database!\n")
  }

  gs.names <- gsdb$gs.names[keep.idx]
  gs <- gs[keep.idx]
  gs.descs <- gsdb$gs.descs[keep.idx]
  size.G <- gsdb$size.G[keep.idx]
  size.ol.G <- size.ol.G[keep.idx]
  if (!is.null(gs.direction)) gs.direction <- gs.direction[keep.idx]

  ## drop gene sets with a duplicated name, keeping the first occurrence
  locs <- match(unique(gs.names), gs.names)
  if (length(locs) < length(gs.names)) {
    gs <- gs[locs]
    gs.names <- gs.names[locs]
    gs.descs <- gs.descs[locs]
    size.G <- size.G[locs]
    size.ol.G <- size.ol.G[locs]
    if (!is.null(gs.direction)) gs.direction <- gs.direction[locs]
  }

  list(gs = gs, gs.names = gs.names, gs.descs = gs.descs, size.G = size.G,
       size.ol.G = size.ol.G, gs.direction = gs.direction)
}
