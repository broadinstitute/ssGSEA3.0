## 20161013 modified by Karsten Krug
## adapt the ssGSEA code to:
##   1) work with site specific signature sets
##   2) take directionality of regulation into account
##   3) multi-threaded using 'doParallel'
##   4) handle missing values
##   5) cmapR functions for data import/export
##   6) produce more extensive output, e.g. dataset/database overlaps, etc.
## 20260415 edited by D R Mani
##   change permutation testing to bootstrap p-value estimation
##   see bootstrap_es_pvalue() (R/ssgsea-bootstrap.R) and the "NOTE"s in
##   gsea_score() (R/ssgsea-core.R) for details

#' Run one gene set through every sample, dispatching parallel or serial
#'
#' @importFrom foreach %dopar%
#' @keywords internal
.score_all_genesets <- function(m, gene.names, matched, weight, statistic, nperm,
                                 correl.type, min.overlap, par, spare.cores, seed,
                                 verbose, log.file) {
  Ns <- ncol(m)
  Ng <- nrow(m)
  N.gs <- length(matched$gs)

  run_one <- function(gs.i) {
    if (!is.null(log.file)) cat(matched$gs.names[gs.i], "\n", file = log.file, append = TRUE)
    gene.set.direction <- if (!is.null(matched$gs.direction)) matched$gs.direction[[gs.i]] else NULL
    project_geneset(
      data.array = m, gene.names = gene.names, n.cols = Ns, n.rows = Ng,
      weight = weight, statistic = statistic, gene.set = matched$gs[[gs.i]],
      nperm = nperm, correl.type = correl.type, gene.set.direction = gene.set.direction,
      min.overlap = min.overlap, size.G.current = matched$size.G[gs.i]
    )
  }

  if (par) {
    cl <- parallel::makeCluster(max(1, parallel::detectCores() - spare.cores))
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)

    run_task <- run_one
    if (!is.null(seed)) {
      ## foreach's default scheduling assigns gene sets to workers
      ## dynamically (whichever worker is free next runs the next gene
      ## set), so clusterSetRNGStream()'s per-*worker* substreams would
      ## make results depend on scheduling timing, not just the seed.
      ## Instead, generate one independent L'Ecuyer-CMRG substream per
      ## gene set up front, keyed by gs.i -- so which worker happens to
      ## run a given gene set no longer affects its result.
      old.kind <- RNGkind()[1]
      old.seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        RNGkind(old.kind)
        if (!is.null(old.seed)) assign(".Random.seed", old.seed, envir = .GlobalEnv)
      }, add = TRUE)
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
      gs.stream <- vector("list", N.gs)
      state <- .Random.seed
      for (i in seq_len(N.gs)) {
        gs.stream[[i]] <- state
        state <- parallel::nextRNGStream(state)
      }
      run_task <- function(gs.i) {
        assign(".Random.seed", gs.stream[[gs.i]], envir = .GlobalEnv)
        run_one(gs.i)
      }
    }

    gs.i <- NULL ## appease R CMD check re: foreach's NSE binding
    tmp <- foreach::foreach(gs.i = seq_len(N.gs), .packages = "doParallel") %dopar% run_task(gs.i)
  } else {
    if (!is.null(seed)) set.seed(seed)
    if (verbose) message("Scoring ", N.gs, " gene sets...")
    tmp <- lapply(seq_len(N.gs), run_one)
  }

  names(tmp) <- matched$gs.names
  tmp
}

#' Single-sample GSEA / PTM-SEA enrichment, with bootstrap significance
#'
#' Scores every sample in `input.ds` against every gene set in
#' `gene.set.databases`, using a per-gene-set running-sum enrichment
#' statistic (Subramanian et al. 2005; Abazeed et al. 2013) and, for
#' `output.score.type = "NES"`, a gene-level bootstrap to obtain a
#' normalized enrichment score and p-value (see [bootstrap_es_pvalue()] for
#' why this replaced the original permutation-based approach).
#'
#' @param input.ds Path to a GCT file (first column: gene symbols or site
#'   ids), or an already-parsed [cmapR::GCT-class] object -- e.g. the result
#'   of [preprocessGCT()], to avoid re-reading it from disk.
#' @param output.prefix File prefix for output GCT files.
#' @param gene.set.databases Character vector of paths to GMT gene set
#'   database files.
#' @param sample.norm.type See [normalize_samples()].
#' @param weight When `0`, all genes are weighted equally; when `> 0`, actual
#'   (normalized) expression values are used and can shift the score --
#'   see `correl.type`.
#' @param statistic `"area.under.RES"` or `"Kolmogorov-Smirnov"`, see
#'   [gsea_score()].
#' @param output.score.type `"NES"` (bootstrap normalized enrichment score
#'   and p-value) or `"ES"` (raw enrichment score only, no bootstrapping,
#'   `p.value` is always `1`).
#' @param nperm Number of gene-level bootstrap replicates (ignored when
#'   `output.score.type = "ES"`).
#' @param min.overlap Minimum number of a gene set's members that must be
#'   present in the data for it to be scored.
#' @param correl.type `"rank"`, `"z.score"`, or `"symm.rank"` -- only
#'   relevant when `weight > 0`, see [gsea_score()].
#' @param global.fdr If `TRUE`, FDR is computed once across all samples'
#'   pooled p-values; otherwise per sample.
#' @param extended.output If `TRUE`, include per-signature overlap
#'   statistics as additional row annotations in the output GCTs.
#' @param par If `TRUE`, score gene sets in parallel across
#'   `parallel::detectCores() - spare.cores` workers.
#' @param spare.cores Number of CPU cores to leave idle when `par = TRUE`.
#' @param export.signat.gct If `TRUE`, also write one GCT per signature (to a
#'   `signature_gct/` subdirectory) containing that signature's members'
#'   original expression values.
#' @param param.file If `TRUE`, write a `<output.prefix>_parameters.txt`
#'   record of the parameters used.
#' @param seed Optional random seed, for reproducible bootstrap resampling.
#' @param verbose If `FALSE`, suppress progress messages (the log file, if
#'   `log.file` is not `NULL`, is still written).
#' @param log.file Path to a log file recording basic progress/timing;
#'   `NULL` to skip logging entirely.
#' @return An object of class `ssgsea_result` (see [build_ssgsea_result()]);
#'   the same information is also written to `<output.prefix>-scores.gct`,
#'   `-pvalues.gct`, `-fdr-pvalues.gct`, and `-combined.gct` (plus, if
#'   `export.signat.gct`, one GCT per signature).
#'
#' @section Fixed limitation (directional/PTMsigDB gene sets): earlier
#'   versions of this port could desync a directional (PTM-SEA/PTMsigDB
#'   `;u`/`;d`) gene set's position/direction bookkeeping inside
#'   [gsea_score()] whenever a bootstrap replicate happened to resample only
#'   *some* of a gene set's members -- sometimes erroring out, sometimes
#'   silently contributing a wrong enrichment score to that replicate's
#'   bootstrap null. This is now fixed via [.align_geneset_to_direction()],
#'   which keeps the two in lockstep regardless of which members a given
#'   replicate happens to include. The vulnerable construction traces back to
#'   code essentially unchanged since directionality support was added to the
#'   original ssGSEA2.0 in 2016, but it could never actually manifest there:
#'   ssGSEA2.0's gene-label *permutation* reorders every gene without ever
#'   dropping one, so its (fixed, computed-once) position/direction vectors
#'   never desync. It's specifically ssGSEA3.0's gene-level bootstrap
#'   resampling *with replacement* -- which can and does drop individual
#'   genes from a given replicate -- that exposes it; this is a real
#'   consequence of the bootstrap rework, not a latent defect ssGSEA2.0 use
#'   was ever exposed to.
#' @export
ssGSEA3 <- function(input.ds,
                     output.prefix,
                     gene.set.databases,
                     sample.norm.type = c("rank", "log", "log.rank", "none"),
                     weight = 0,
                     statistic = c("area.under.RES", "Kolmogorov-Smirnov"),
                     output.score.type = c("NES", "ES"),
                     nperm = 1000,
                     min.overlap = 10,
                     correl.type = c("rank", "z.score", "symm.rank"),
                     global.fdr = FALSE,
                     extended.output = TRUE,
                     par = FALSE,
                     spare.cores = 1,
                     export.signat.gct = TRUE,
                     param.file = TRUE,
                     seed = NULL,
                     verbose = TRUE,
                     log.file = "run.log") {
  ## For results similar to the Java version, use weight=0; when weight==0,
  ## sample.norm.type and correl.type do not matter. When weight > 0, the
  ## combination of sample.norm.type and correl.type dictates how the
  ## expression values in input.ds are transformed to obtain the score --
  ## use this setting with care, since the transformations can skew scores
  ## towards +ve or -ve values:
  ##   - sample.norm.type=="none" uses actual expression values; combined
  ##     with correl.type=="rank", genes are weighted by actual values
  ##   - sample.norm.type=="rank" weights genes proportional to rank
  ##   - sample.norm.type=="log" can be used for log-transforming input data
  ##   - correl.type=="z.score" standardizes the (normalized) input values
  ##     before using them to calculate scores
  sample.norm.type <- match.arg(sample.norm.type)
  statistic <- match.arg(statistic)
  output.score.type <- match.arg(output.score.type)
  correl.type <- match.arg(correl.type)

  if (!is.null(log.file)) cat("##", format(Sys.time()), "\n", file = log.file)

  if (param.file) {
    input.ds.str <- if (methods::is(input.ds, "GCT")) {
      path <- attr(input.ds, "path")
      if (!is.null(path)) path else if (nzchar(input.ds@src)) input.ds@src else "<in-memory GCT object>"
    } else {
      input.ds
    }
    param.str <- c(
      paste("##", Sys.time()),
      paste("input gct:", input.ds.str, sep = "\t"),
      paste("gene.set.database:", paste(gene.set.databases, collapse = ", "), sep = "\t"),
      paste("sample.norm.type:", sample.norm.type, sep = "\t"),
      paste("weight:", weight, sep = "\t"),
      paste("statistic:", statistic, sep = "\t"),
      paste("output.score.type", output.score.type, sep = "\t"),
      paste("nperm (gene-level bootstrap replicates):", nperm, sep = "\t"),
      paste("global.fdr:", global.fdr, sep = "\t"),
      paste("min.overlap:", min.overlap, sep = "\t"),
      paste("correl.type:", correl.type, sep = "\t"),
      paste("export.signat.gct:", export.signat.gct, sep = "\t"),
      paste("run.parallel:", par, sep = "\t"),
      paste("seed:", if (is.null(seed)) "NA" else seed, sep = "\t")
    )
    writeLines(param.str, con = paste(output.prefix, "parameters.txt", sep = "_"))
  }

  ## ##################################################
  ## import dataset
  dataset <- read_gct_safe(input.ds)
  m <- dataset@mat
  gene.names <- dataset@rid
  if (isTRUE(attr(dataset, "deduped"))) gene.names <- sub("_[0-9]{1,5}$", "", gene.names)
  gene.descs <- dataset@rdesc
  sample.names <- dataset@cid
  sample.descs <- dataset@cdesc
  m.org <- m
  gct.version <- dataset@version
  gct.src <- dataset@src

  ## ##################################################
  ## import + match gene set databases
  gsdb <- read_gene_sets(gene.set.databases, thres.min = min.overlap, thres.max = 2000)
  matched <- match_gene_sets(gsdb, gene.names, min.overlap)

  ## ##################################################
  ## sample normalization
  m <- normalize_samples(m, sample.norm.type)

  ## ES-only mode never bootstraps -- see bootstrap_es_pvalue()'s docs for
  ## why nperm=0 gives ES==NES and p.value==1 without any resampling
  nperm.use <- if (output.score.type == "ES") 0 else nperm

  tt <- Sys.time()
  tmp <- .score_all_genesets(m, gene.names, matched, weight, statistic, nperm.use,
                              correl.type, min.overlap, par, spare.cores, seed,
                              verbose, log.file)
  if (verbose) message("main loop: ", format(Sys.time() - tt))

  result <- build_ssgsea_result(tmp, matched$gs.names, matched$gs.descs, matched$size.G,
                                 matched$gs.direction, matched$gs, sample.names, sample.descs,
                                 output.score.type = output.score.type, global.fdr = global.fdr,
                                 extended.output = extended.output,
                                 data = m.org, gene.names = gene.names)

  write_ssgsea_gcts(result, output.prefix, m.org, gene.names, gene.descs,
                     gct.src = gct.src, gct.version = gct.version,
                     export.signat.gct = export.signat.gct, appenddim = FALSE)

  result
}
