#' HALLMARK gene set process categories
#'
#' Named character vector mapping each of the 50 MSigDB HALLMARK gene set
#' names to a broad process category (signaling, immune, metabolic, etc.),
#' used to color-code and group rows in summary heatmaps.
#'
#' @format A named character vector of length 50; names are HALLMARK gene set
#'   IDs, values are one of "signaling", "immune", "development",
#'   "proliferation", "cellular component", "metabolic", "DNA damage", or
#'   "pathway".
#' @source Broad Institute MSigDB HALLMARK collection; category assignments
#'   originate from the lab's ssGSEA reporting scripts (`pw_hm()`).
"hallmark_categories"

#' Color palette for HALLMARK process categories
#'
#' @format A named character vector of hex color codes, one per category
#'   value in [hallmark_categories].
"hallmark_category_colors"

utils::globalVariables(c("hallmark_categories", "hallmark_category_colors"))
