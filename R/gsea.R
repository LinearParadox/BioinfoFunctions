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
#' @importFrom clusterProfiler GSEA
#' @importFrom dplyr filter distinct
#' @importFrom rlang .data
#' @export
run_gsea <- function(results, geneset, p="PValue",
                     logFC="logFC", name="gene_id",
                     min.size=10, max.size=500,
                     ranking=c("logpfc", "FC", "signp"),
                     pcutoff=0.05, eps=1e-10){
  ranking <- match.arg(ranking)
  results<-results %>% dplyr::distinct(.data[[name]], .keep_all = T) %>%
    dplyr::filter(.data[[name]] !="" & !is.na(.data[[name]]))
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

#' Wrapper function around various barplots for gsea
#'
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 labs
#' @param gseaRes GSEA results from run_gsea
#' @param n number of pathways to plot
#' @param direction which direction to plot. all plots all,
#'  both plots up and down in one graph
#'  up and down plot only up and down respectively
#' @param fill what color to use for fill. Can be a vector of length 2,
#'  a color string which will be reused or a column in GSEA res
#' @param return.type whether to return a patchwork object or a list of plots
#' @param sort.by  what to sort by in the case of the both option
#'
gsea_barplot <- function(gseaRes, n=15, direction=c("all", "up", "down"),
                         fill=NA, return.type=c("patchwork", "list")){
  direction=match.arg(direction)
  return.type = match.arg(return.type)
  if(direction=="up"){
    return(.gsea_barplot_up(gseaRes, n=n, fill=fill))
  } else if(direction=="down"){
    return(.gsea_barplot_down(gseaRes, n=n, fill=fill))
  }
  else if(direction == "all"){
    if(!is.list(fill)){
      fill <- rep(fill, 2)
    }
    up <- .gsea_barplot_up(gseaRes, n=n, fill=fill[[1]])+ggplot2::labs(title=paste0("Top ", n, " Significantly up Pathways"))
    down <- .gsea_barplot_down(gseaRes, n=n, fill=fill[[2]]) + ggplot2::labs(title=paste0("Top ", n, " Significantly down Pathways"))
    if(return.type=="patchwork"){
      return(patchwork::wrap_plots(up, down, ncol = 2))
    }
    else{
      return(list(up=up, down=down))
    }
  }
  #TODO: IMPLEMENT SAME GRAPH GSEA PLOT
}

#' Plot barplot to plot GSEA barplot up
#' @importFrom dplyr filter slice_max
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme theme_minimal element_blank theme_void
#' @importFrom forcats fct_reorder
#' @importFrom rlang .data
#' @param gseaRes GSEA results from run_gsea
#' @param n number of pathways to plot
#' @param fill fill color of bars
#' @noRd
#'
#'

.gsea_barplot_up <- function(gseaRes, n = 15, fill = NA) {
  if (is.na(fill)) fill <- "#FA8072"

  res <- gseaRes %>%
    dplyr::filter(NES > 0) %>%
    dplyr::slice_max(NES, n = n, with_ties = FALSE)

  if (nrow(res) < 1) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  plt <- ggplot2::ggplot(
    res,
    ggplot2::aes(
      y = -log10(.data[["p.adjust"]]),
      x = forcats::fct_reorder(.data[["Description"]], -log10(.data[["p.adjust"]]))
    )
  ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

  if (isColor(fill)) {
    plt <- plt + ggplot2::geom_col(fill = fill) +
      ggplot2::theme(legend.position = "none")
  } else {
    plt <- plt + ggplot2::geom_col(ggplot2::aes(fill = .data[[fill]]))
  }

  plt
}



#' Plot barplot to plot GSEA barplot down
#' @importFrom dplyr filter slice_max
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme theme_minimal element_blank theme_void
#' @importFrom forcats fct_reorder
#' @importFrom rlang .data
#' @param gseaRes GSEA results from run_gsea
#' @param n number of pathways to plot
#' @param fill fill color of bars
#' @noRd
#'
#'
#'
.gsea_barplot_down <- function(gseaRes, n = 15, fill = NA) {
  if (is.na(fill)) fill <- "#87CEFA"

  res <- gseaRes %>%
    dplyr::filter(NES < 0) %>%
    dplyr::slice_min(NES, n = n, with_ties = FALSE)

  if (nrow(res) < 1) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  plt <- ggplot2::ggplot(
    res,
    ggplot2::aes(
      y = -log10(.data[["p.adjust"]]),
      x = forcats::fct_reorder(.data[["Description"]], -log10(.data[["p.adjust"]]))
    )
  ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

  if (isColor(fill)) {
    plt <- plt + ggplot2::geom_col(fill = fill) +
      ggplot2::theme(legend.position = "none")
  } else {
    plt <- plt + ggplot2::geom_col(ggplot2::aes(fill = .data[[fill]]))
  }

  plt
}




