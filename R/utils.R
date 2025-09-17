#' Pull mouse equivalents of human gene sets
#'
#' Orthology map mouse to human gene sets using biomart
#' @param geneset A list of genes to orthology map
#' @param mart a biomart object
#' @param keytype keytypes to search mart
#' @param value value to return
#'
#' @importFrom biomaRt select
#' @importFrom dplyr select filter
#' @export
ortho_map <- function(geneset, mart, keytype='ensembl_gene_id',
                      value='mmusculus_homolog_ensembl_gene'){
  return_cols <- c(value,
                   keytype,
                   "mmusculus_homolog_orthology_confidence",
                   "mmusculus_homolog_perc_id")
  homologs <- biomaRt::select(mart, keys=geneset,
                              keytype=keytype, columns=return_cols)
  mouse_genes <- homologs %>%
    dplyr::filter(value!="" & mmusculus_homolog_orthology_confidence == 1) %>%
    dplyr::select(all_of(value))
  return(mouse_genes)
}

#' Strip version from ensenmbl geneID
#'
#' @param genes Vector of ensembl geneid versions
#'
#' @export
strip_version <- function(genes){
  return(gsub("\\..*", "", genes))
}

#' Returns whether a string is a color
#' Credit from https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
#' @noRd

#'
#' @param x string to test for color or not

isColor <- function(x)
{
  res <- try(col2rgb(x), silent = TRUE)
  return(!"try-error" %in% class(res))
}
