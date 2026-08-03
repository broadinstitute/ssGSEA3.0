#' Preprocess a GCT file for ssGSEA or PTM-SEA
#'
#' Converts a raw GCT file into one of three report levels used downstream by
#' `ssGSEA3()` (not yet ported):
#' \itemize{
#'   \item `"gcr"` (gene-centric-redundant): relabels rows with gene symbols,
#'     keeping every row (deduplicating names only).
#'   \item `"gc"` (gene-centric): collapses multiple rows per gene into one,
#'     per `mode`.
#'   \item `"ssc"` (single-site-centric): the PTM-SEA path. Explodes
#'     multiply-modified peptides into one row per site and builds a
#'     PTM-SEA-compatible site ID.
#' }
#'
#' @param gct.str Path to the input GCT (v1.2/v1.3) file.
#' @param level One of `"ssc"`, `"gc"`, `"gcr"` (see above).
#' @param id.type Input site-ID notation: `"sm"` (Spectrum Mill) or `"ph"`
#'   (Philosopher). Only used when `level = "ssc"`.
#' @param id.type.out Output site-ID type: `"uniprot"`, `"refseq"`, or
#'   `"seqwin"`. Only used when `level = "ssc"`.
#' @param acc.type Accession type of the `rid` in the input GCT: `"uniprot"`,
#'   `"refseq"`, or `"symbol"`. Only used when `level = "ssc"`.
#' @param organism One of `"human"`, `"mouse"`, `"rat"`, `"zebrafish"`.
#'   Required only when mapping to UniProt (`level = "ssc"`,
#'   `id.type.out = "uniprot"`, `acc.type` in `c("refseq", "symbol")`).
#' @param seqwin.col Column containing `|`-separated flanking sequences; only
#'   used when `id.type.out = "seqwin"`.
#' @param gene.col Column listing gene symbols; used for gene-centric levels
#'   and by `fix.gene.names`.
#' @param fix.gene.names If `TRUE`, replace `gene.col` with gene names
#'   extracted from `fix.gene.names.column` via `fix.gene.names.regexpr` --
#'   a workaround for gene symbols corrupted by spreadsheet
#'   date-autoformatting (e.g. "SEPT10" becoming a date).
#' @param fix.gene.names.column Column to extract corrected gene names from.
#' @param fix.gene.names.regexpr Regular expression identifying the gene name
#'   within `fix.gene.names.column`.
#' @param humanize.gene.names If `TRUE`, upper-case gene symbols (for mouse or
#'   rat data). Only used when `level = "gc"`.
#' @param loc If `TRUE` (recommended), only fully localized PTM sites are
#'   kept. Only used when `level = "ssc"`.
#' @param mode How multiple rows per gene are combined when `level = "gc"`,
#'   and how multiply-detected duplicate sites are combined when
#'   `level = "ssc"` (only `"mean"`, `"median"`, `"abs.max"` apply to
#'   `"ssc"`): `"mean"`, `"median"`, `"sd"` (keep the single row with the
#'   highest cross-sample standard deviation -- i.e. the most
#'   variable/informative probe for that gene; `"gc"` only), `"SGT"`
#'   (subgroup-top: first Spectrum Mill subgroup per protein group; `"gc"`
#'   only), or `"abs.max"` (per sample, the most extreme value -- for
#'   signed, log-transformed p-values).
#' @param SGT.col Column used to collapse to subgroup-top (`mode = "SGT"`).
#' @param mod.res Modified residue(s), e.g. `"S|T|Y"` or `"K"`.
#' @param mod.type Modification type, e.g. `"p"` or `"ac"`.
#' @param appenddim See [cmapR::write_gct()].
#'
#' @return (Invisibly) the path to the written output GCT file.
#' @export
preprocessGCT <- function(gct.str,
                           level = c("ssc", "gc", "gcr"),
                           id.type = c("sm", "ph"),
                           id.type.out = c("uniprot", "refseq", "seqwin"),
                           acc.type = c("uniprot", "refseq", "symbol"),
                           organism = NULL,
                           seqwin.col = "VMsiteFlanks",
                           gene.col = "geneSymbol",
                           fix.gene.names = FALSE,
                           fix.gene.names.column = "id.description",
                           fix.gene.names.regexpr = ".* GN\\=(.*)$",
                           humanize.gene.names = FALSE,
                           loc = TRUE,
                           mode = c("mean", "median", "sd", "SGT", "abs.max"),
                           SGT.col = "subgroupNum",
                           mod.res = c("S|T|Y", "K"),
                           mod.type = c("p", "ac", "ub"),
                           appenddim = TRUE) {
  level <- match.arg(level)
  mode <- match.arg(mode)
  mod.res <- match.arg(mod.res)
  mod.type <- match.arg(mod.type)
  acc.type <- match.arg(acc.type)
  id.type <- match.arg(id.type)
  id.type.out <- match.arg(id.type.out)

  cat("importing gct file: ", gct.str, " ...\n")
  gct <- read_gct_safe(gct.str)
  if (is.null(gct@rid)) gct@rid <- rownames(gct@mat)
  if (is.null(gct@cid)) gct@cid <- colnames(gct@mat)
  gct@rdesc$original.id <- gct@rid

  if (fix.gene.names) {
    if (!gene.col %in% colnames(gct@rdesc)) stop(glue::glue("Column {gene.col} not found!"))
    if (!fix.gene.names.column %in% colnames(gct@rdesc)) stop(glue::glue("Column {fix.gene.names.column} not found!"))
    gct@rdesc[[gene.col]] <- sub(fix.gene.names.regexpr, "\\1", gct@rdesc[[fix.gene.names.column]])
  }

  gct <- switch(level,
    gcr = .preprocess_gcr(gct, gene.col = gene.col),
    gc  = .preprocess_gc(gct, gene.col = gene.col, mode = mode, SGT.col = SGT.col, humanize.gene.names = humanize.gene.names),
    ssc = .preprocess_ssc(gct, id.type = id.type, id.type.out = id.type.out, acc.type = acc.type,
                           organism = organism, seqwin.col = seqwin.col, loc = loc,
                           mod.res = mod.res, mod.type = mod.type, mode = mode)
  )

  ## ###################################
  ## export GCT
  fn <- sub(".*/", "", gct.str)
  sr <- level %in% c("ssc", "gcr")
  glue.str <- glue::glue("_{level}-{if (sr) id.type.out else ''}{if (sr && loc) '_loc' else ''}")
  if (appenddim) {
    fn <- if (length(grep("_n[0-9]*x[0-9]*\\.gct$", fn)) > 0) {
      sub("_n[0-9]*x[0-9]*\\.gct$", glue.str, fn)
    } else {
      sub("\\.gct$", glue.str, fn)
    }
  } else {
    fn <- sub("\\.gct$", paste0(glue.str, ".gct"), fn)
  }
  fn <- gsub("-\\.", ".", fn)
  fn <- gsub("-$", "", fn)

  cmapR::write_gct(gct, ofile = fn, appenddim = appenddim)
  invisible(fn)
}

## #####################################################################
##                     gene-centric-redundant (gcr)
## #####################################################################
#' @keywords internal
.preprocess_gcr <- function(gct, gene.col) {
  mat <- gct@mat; rid <- gct@rid; cid <- gct@cid; rdesc <- gct@rdesc

  if (!gene.col %in% colnames(rdesc)) stop(glue::glue("Column {gene.col} not found!"))
  genes <- rdesc[[gene.col]]

  rm.idx <- which(is.na(genes) | nchar(genes) == 0)
  if (length(rm.idx) > 0) {
    genes <- genes[-rm.idx]
    rid <- rid[-rm.idx]
    mat <- mat[-rm.idx, , drop = FALSE]
    rdesc <- rdesc[-rm.idx, , drop = FALSE]
    warning(glue::glue("Removed {length(rm.idx)} rows due to missing gene symbol in column '{gene.col}'"))
  }

  gct@mat <- data.matrix(mat)
  gct@rid <- make.unique(genes, sep = "_")
  gct@cid <- cid
  gct@rdesc <- rdesc
  gct
}

## #####################################################################
##                          gene-centric (gc)
## #####################################################################
#' @keywords internal
.preprocess_gc <- function(gct, gene.col, mode, SGT.col, humanize.gene.names) {
  mat <- gct@mat; rid <- gct@rid; cid <- gct@cid; rdesc <- gct@rdesc

  if (!gene.col %in% colnames(rdesc)) stop(glue::glue("Column {gene.col} not found!"))
  genes <- rdesc[[gene.col]]
  if (humanize.gene.names) genes <- toupper(genes)

  rm.idx <- which(is.na(genes) | nchar(genes) == 0)
  if (length(rm.idx) > 0) {
    genes <- genes[-rm.idx]
    mat <- mat[-rm.idx, , drop = FALSE]
    rdesc <- rdesc[-rm.idx, , drop = FALSE]
    warning(glue::glue("Removed {length(rm.idx)} rows due to missing gene symbol in column '{gene.col}'"))
  }

  if (mode == "mean") {
    mat.gc <- stats::aggregate(mat, FUN = function(x) mean(x, na.rm = TRUE), by = list(genes))
  } else if (mode == "median") {
    mat.gc <- stats::aggregate(mat, FUN = function(x) stats::median(x, na.rm = TRUE), by = list(genes))
  } else if (mode == "abs.max") {
    mat.gc <- stats::aggregate(mat, FUN = function(x) x[which.max(abs(x))], by = list(genes))
  } else if (mode == "sd") {
    ## keep, per gene, the single row with the highest cross-sample standard
    ## deviation -- i.e. the most variable/informative probe for that gene --
    ## rather than aggregating across rows. Mirrors the SGT approach below.
    row.sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
    keep.idx <- unlist(tapply(seq_along(genes), genes, function(idx) idx[which.max(row.sd[idx])]))
    genes <- genes[keep.idx]
    mat.gc <- data.frame(id = genes, mat[keep.idx, , drop = FALSE], check.names = FALSE)
    rdesc <- rdesc[keep.idx, , drop = FALSE]
  } else if (mode == "SGT") {
    if (!SGT.col %in% colnames(rdesc)) stop(glue::glue("Column {SGT.col} not found!"))
    sgt <- rdesc[[SGT.col]]
    ## keep lowest subgroup per gene
    sgt.idx <- tapply(sgt, genes, function(x) x[which.min(as.numeric(sub(".*\\.", "", x)))])
    ## drop duplicated subgroup numbers introduced when there are >9
    ## subgroups (a PANOPLY parsing-module artifact), e.g. 234.10 -> 234.1
    sgt.idx <- unique(sgt.idx)
    keep.idx <- match(sgt.idx, sgt)
    genes <- names(sgt.idx)
    mat.gc <- data.frame(id = genes, mat[keep.idx, , drop = FALSE], check.names = FALSE)
    rdesc <- rdesc[keep.idx, , drop = FALSE]
  }

  rid <- as.character(mat.gc[, 1])
  mat.gc <- mat.gc[, -1, drop = FALSE]

  rdesc.gc <- stats::aggregate(rdesc, FUN = function(x) paste(unique(x), collapse = "|"), by = list(genes))
  rdesc.gc <- rdesc.gc[, -1, drop = FALSE]

  rownames(mat.gc) <- rownames(rdesc.gc) <- rid
  gct@mat <- data.matrix(mat.gc)
  gct@rid <- rid
  gct@cid <- cid
  gct@rdesc <- rdesc.gc
  gct
}

## #####################################################################
##                  single-site-centric (ssc, PTM-SEA)
## #####################################################################
#' @keywords internal
.preprocess_ssc <- function(gct, id.type, id.type.out, acc.type, organism,
                             seqwin.col, loc, mod.res, mod.type, mode) {
  mat <- gct@mat; rid <- gct@rid; cid <- gct@cid; rdesc <- gct@rdesc
  n <- ncol(mat)

  site.ids <- gsub(" ", "", rid)
  names(site.ids) <- rid

  ## #####################################
  ##   make ids UniProt-centric
  ## #####################################
  if (id.type.out == "uniprot") {
    if (acc.type == "refseq") {
      ## RefSeq only (drop e.g. UniProt accessions of lab contaminants), then map to UniProt
      idx <- grep("^(NP_|XP_|YP_|ZP_)", site.ids)
      site.ids <- site.ids[idx]
      rid <- rid[idx]
      mat <- mat[idx, , drop = FALSE]
      rdesc <- rdesc[idx, , drop = FALSE]

      if (is.null(organism))
        stop("`organism` must be specified (one of 'human', 'mouse', 'rat', 'zebrafish') when mapping to UniProt (id.type.out = 'uniprot', acc.type = 'refseq').")
      up <- RefSeq2UniProt(sub("^(.*?_.*?)_.*", "\\1", site.ids), keytype = "REFSEQ", organism = organism)$id.map$id.mapped
      mapped.idx <- which(!is.na(up))

      up <- up[mapped.idx]
      site.ids <- site.ids[mapped.idx]
      rid <- rid[mapped.idx]
      mat <- mat[mapped.idx, , drop = FALSE]
      rdesc <- rdesc[mapped.idx, , drop = FALSE]

      site.ids <- paste0(up, sub("^(.*?_.*?)(_.*)", "\\2", site.ids))
      names(site.ids) <- rid
    } else if (acc.type == "symbol") {
      ## gene symbol + site, e.g. EIF2S2.S2
      if (is.null(organism))
        stop("`organism` must be specified (one of 'human', 'mouse', 'rat', 'zebrafish') when mapping to UniProt (id.type.out = 'uniprot', acc.type = 'symbol').")
      up <- RefSeq2UniProt(sub("^(.*?)-.*", "\\1", site.ids), keytype = "SYMBOL", organism = organism)$id.map$id.mapped
      mapped.idx <- which(!is.na(up))

      up <- up[mapped.idx]
      rid <- rid[mapped.idx]
      mat <- mat[mapped.idx, , drop = FALSE]
      rdesc <- rdesc[mapped.idx, , drop = FALSE]

      site.ids <- glue::glue("{up};{sub('^(.*?)-(.*)', '\\\\2', site.ids[mapped.idx])}-{mod.type}")
      names(site.ids) <- rid
    }
    ## acc.type == "uniprot": already UniProt-centric, nothing to do
  }

  ## ####################################
  ## extract all sites per entry
  ## ####################################
  if (id.type == "sm") {
    if (loc) {
      loc.idx <- which(sapply(strsplit(sub("^.*_([1-9]_[0-9])_.*", "\\1", site.ids), "_"), function(x) length(unique(x))) == 1)
      rid <- rid[loc.idx]
      site.ids <- site.ids[loc.idx]
      mat <- mat[loc.idx, , drop = FALSE]
      rdesc <- rdesc[loc.idx, , drop = FALSE]
    }

    if (acc.type == "refseq" && id.type.out == "refseq") {
      idx <- grep("^(NP_|XP_|YP_|ZP_)", site.ids)
      site.ids <- site.ids[idx]
      rdesc <- rdesc[idx, , drop = FALSE]

      site.var <- sub("^(.{2}_.*?_.*?)_.*", "\\1", site.ids)
      names(site.var) <- names(site.ids)

      all.sites <- lapply(strsplit(site.var, tolower(mod.res)), function(x) {
        prot.acc <- sub("^(.{2}_.*?)_.*", "\\1", x[1])
        x <- sub(paste0(prot.acc, "_"), "", x)
        paste(prot.acc, grep(toupper(mod.res), x, value = TRUE), sep = ";")
      })
    } else if (id.type.out == "seqwin") {
      if (!(seqwin.col %in% colnames(rdesc))) stop(glue::glue("column '{seqwin.col}' not found!"))
      site.var <- gsub("-", "_", toupper(rdesc[[seqwin.col]]))
      names(site.var) <- names(site.ids)
      all.sites <- strsplit(site.var, "\\|")
    } else {
      site.var <- sub("^(.*?_.*?)_.*", "\\1", site.ids)
      names(site.var) <- names(site.ids)
      all.sites <- lapply(strsplit(site.var, tolower(mod.res)), function(x) {
        prot.acc <- sub("_.*", "", x[1])
        x <- sub(".*?_", "", x)
        paste(prot.acc, grep(toupper(mod.res), x, value = TRUE), sep = ";")
      })
    }
  } else if (id.type == "ph") {
    loc.idx <- which(sapply(strsplit(sub("^.*_([1-9]_[0-9])_{0,1}.*", "\\1", site.ids), "_"), function(x) x[2] > 0))
    rid <- rid[loc.idx]
    site.ids <- site.ids[loc.idx]
    mat <- mat[loc.idx, , drop = FALSE]
    rdesc <- rdesc[loc.idx, , drop = FALSE]

    if (acc.type == "refseq" && id.type.out == "refseq") {
      idx <- grep("^(NP_|XP_|YP_|ZP_)", site.ids)
      site.ids <- site.ids[idx]
      rdesc <- rdesc[idx, , drop = FALSE]

      site.var <- gsub("([0-9])([A-Z])", "\\1 \\2", sub(".*_", "", site.ids))
      prot.acc <- sub("^(.{2}_.*?\\.[0-9]{1,2})_.*", "\\1", site.ids)
      site.var <- paste(prot.acc, site.var, collapse = "_")
      names(site.var) <- names(site.ids)

      all.sites <- lapply(strsplit(site.var, " "), function(x) {
        prot.acc <- sub("^(.{2}_.*?)_.*", "\\1", x[1])
        x <- sub(paste0(prot.acc, "_"), "", x)
        paste(prot.acc, grep(toupper(mod.res), x, value = TRUE), sep = ";")
      })
    } else if (id.type.out == "seqwin") {
      if (!seqwin.col %in% colnames(rdesc)) stop(glue::glue("column '{seqwin.col}' not found!"))
      site.var <- gsub("-", "_", toupper(rdesc[[seqwin.col]]))
      names(site.var) <- names(site.ids)
      all.sites <- strsplit(site.var, "\\|")
    } else {
      site.var <- sub("^(.*?_.*?)_.*", "\\1", site.ids)
      names(site.var) <- names(site.ids)
      all.sites <- lapply(strsplit(site.var, tolower(mod.res)), function(x) {
        prot.acc <- sub("_.*", "", x[1])
        x <- sub(".*?_", "", x)
        paste(prot.acc, grep(toupper(mod.res), x, value = TRUE), sep = ";")
      })
    }
  }

  ## ###################################################################
  ##   create single-site-centric report
  ##
  ##   redundant sites on multiply modified peptides: exclude sites on
  ##   multiply-phosphorylated peptides for which we've also detected the
  ##   singly phosphorylated version.
  ## ###################################################################
  n.sites <- sapply(all.sites, length)
  all.sites.mult <- all.sites[which(n.sites > 1)]
  all.sites.sing <- all.sites[which(n.sites == 1)]
  all.sites.mult <- lapply(all.sites.mult, function(x) {
    rm.idx <- which(x %in% unlist(all.sites.sing))
    if (length(rm.idx) > 0) x <- x[-rm.idx]
    x
  })
  all.sites.mult <- all.sites.mult[sapply(all.sites.mult, length) > 0]
  all.sites <- append(all.sites.mult, all.sites.sing)
  all.sites.unique <- unique(unlist(all.sites))

  sites.dup <- unlist(all.sites)[duplicated(unlist(all.sites))]
  sites.nondup <- setdiff(all.sites.unique, sites.dup)

  mat.ss <- matrix(NA, nrow = length(all.sites.unique), ncol = ncol(mat), dimnames = list(all.sites.unique, cid))

  map.idx <- lapply(rownames(mat.ss), function(x) names(all.sites)[grep(paste0(x, "($|\\\")"), all.sites)])
  names(map.idx) <- rownames(mat.ss)

  ## combine expression of duplicated sites, e.g. sites detected on multiple
  ## versions of multiply-phosphorylated peptides. Uses the same `mode` as
  ## the gene-centric path, but only mean/median/abs.max make sense here --
  ## 'sd'/'SGT' are gene-centric-only concepts (the original code let these
  ## reach an opaque eval(parse(text=mode)) call that failed with a cryptic
  ## "could not find function" error; this now fails with a clear message).
  if (!mode %in% c("mean", "median", "abs.max"))
    stop(glue::glue("mode = '{mode}' is not supported for level = 'ssc' (only 'mean', 'median', 'abs.max' apply here; 'sd' and 'SGT' are gene-centric-only concepts)."))

  for (s in sites.dup) {
    mat.tmp <- mat[map.idx[[s]], , drop = FALSE]
    mat.tmp2 <- if (mode == "abs.max") {
      apply(mat.tmp, 2, function(x) x[which.max(abs(x))])
    } else {
      apply(mat.tmp, 2, get(mode), na.rm = TRUE)
    }
    mat.tmp2[lengths(mat.tmp2) == 0] <- NA
    mat.ss[s, ] <- unlist(mat.tmp2)
  }
  mat.ss[names(map.idx[sites.nondup]), ] <- data.matrix(mat[unlist(map.idx[sites.nondup]), , drop = FALSE])

  ## multiple sites (rids) have been combined -- pick a single rid to represent them
  idx <- sapply(map.idx, function(x) x[1])
  rid <- rid[idx]
  names(rid) <- names(map.idx)
  rdesc <- rdesc[idx, , drop = FALSE]
  ## `%in%` (not `==`) since `x` may list >1 peptide name for a duplicated
  ## site -- `==` against a >1-length vector triggers an (unintended, harmless
  ## in practice) recycling warning without changing the result for length-1
  ## `x`, but `%in%` is the actually-correct membership check either way.
  VMsitesAll <- sapply(map.idx, function(x) paste(names(all.sites)[names(all.sites) %in% x], collapse = "|"))
  rdesc <- data.frame(rdesc, VMsitesAll)

  ptm.site.ids <- paste(names(map.idx), mod.type, sep = ifelse(grepl("^-", mod.type), "", "-"))

  rownames(mat.ss) <- rownames(rdesc) <- rid <- ptm.site.ids

  gct@mat <- mat.ss
  gct@rid <- rid
  gct@cid <- cid
  gct@rdesc <- rdesc
  gct
}
