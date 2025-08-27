#' Plot GSEA results
#'
#' @param results dataframe of gene expression
#' @param p column with adjusted p values, default PValue
#' @param logFC column with logFC values
#' @param geneset set of genes to test, in the format GSNAME,GENES
#' @param name column with gene name
#' @param ranking ranking metric for genes. "signp" for signlogFC*Pval, "FC" for logFC, logpfc for -logp\*FC
#' @param pcutoff P value cutoff for returned genes
#' @param min.size min gs size for GSEA
#' @param max.size max gs size for gsea
#' @param eps boundary for calculating p value
#'
#' @import ClusterProfiler
#' @export
run_gsea <- function(results, geneset, p="PValue",
                     logFC="logFC", name="gene_id",
                     min.size=10, max.size=500,
                     pcutoff=0.05, eps=1e-10){
  results<-results %>% dplyr::distinct({{ name }}) %>%
    dplyr::filter({{ name }}!="" & !is.na({{ name }}))
  if(ranking=="signp"){
    genes=sign(results[[logFC]])*results[[PValue]]
  } else if(ranking=="FC"){
    genes=results[[logFC]]
  }else if(ranking=="logpfc"){
    genes=-log10(results[[p]])*results[[logFC]]
  }else{
    stop("Unsupported ranking metric!")
  }
  names(genes) <- results[[name]]
  genes <- genes[is.finite(genes)]
  genes <- sort(genes, decreasing=T)
  res <- clusterProfiler::GSEA(genes, TERM2GENE = geneset,
                               pvalueCutoff = pcutoff, minGSSize = min.size,
                               maxGSSize = max.size, eps=eps)
  return(res)
}

