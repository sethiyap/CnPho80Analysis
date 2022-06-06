#' plot_qPCR_foldchange
#'
#' @param fc_tibble tibble of fold-change values determined as 2^-∆∆Ct, can be
#'   output of compute_qPCR_foldchange function. With rows as sample name ie.,
#'   WT_R1, WT_R2, WT_R3, treatment_R1,... and columns are genes. Replicates
#'   should always be in the format "_R1" or "_1" or "_set1"
#' @param display_pval Logical, to display pvalue on the barplot. default: NULL
#' @param nrow numeric, number of rows in which the genes should be displayed.
#'   default: 2
#' @param control_group character vector, sample or group which should be
#'   considered as reference to compute pvalue from t.test.
#'
#' @return barplot of the fold-change
#' @export
#' @import tidyr
#' @import dplyr
#' @import ggpubr
#' @import rstatix
#' @import forcats
#' @import stats
#' @import utils
#' @import tibble
#'
#' @details pvalue calculation
#'   source:\url{https://stats.oarc.ucla.edu/other/mult-pkg/faq/pvalue-htm/}
#' @examples
#' library(magrittr)
#' gene_1 <- c(1.08, 1.01 , 0.917, 0.401, 0.356, 0.450)
#' gene_2 <- c(1.02, 1.03, 1.002,0.8, 0.76, 0.79)
#' rname <- c(paste("WT", seq(1:3), sep = "_"), paste("WT+drug", seq(1:3), sep = "_"))
#' qPCR_fc <- tidyr::tibble(gene_1, gene_2)
#' rownames(qPCR_fc) <- rname
#'
#' qPCR_fc <- qPCR_fc %>% tibble::rownames_to_column(var = "Sample")
#' plot_qPCR_foldchange(fc_tibble = qPCR_fc, display_pval = FALSE, nrow=1)
#'
plot_qPCR_foldchange <- function(fc_tibble, display_pval=NULL, control_group="WT",nrow=2){

  library(ggpubr)
  # First column name, always sample
  names(fc_tibble)[1] = "sample"

  fc_tibble_group <- fc_tibble %>%
                      tidyr::separate(col = sample, into = c("sample", "replicate"), sep = "_") %>%
                      dplyr::select(-c(replicate)) %>%
                      tidyr::gather(gene, value, -sample)

  pval_tibble <- fc_tibble_group %>%
                      dplyr::mutate(gene=forcats::as_factor(gene), sample=forcats::as_factor(sample)) %>%
                      dplyr::group_by(gene) %>%
                      rstatix::t_test(value ~ sample, paired = FALSE, ref.group = control_group) %>%
                      dplyr::mutate(p=p/2,
                                    # left tailed: treatment < control, right tailed: treatment > control
                                    # However, in population of genes if some are control > treatment or control < treatment.
                                    # then converting the two-tailed pvalue into one-tailed useful (by dividing the two-tailed pvalue into half)
                                    # as when we don't know the direction of population

                                    pval_significance= dplyr::if_else(p < 0.05 & p >= 0.01,"*",
                                                                      dplyr::if_else(p <= 0.01 & p >= 0.001, "**", # adding pvalue significance
                                                                                     dplyr::if_else(p <= 0.001 & p >= 0.0001, "***", "ns")))) %>%

                      rstatix::adjust_pvalue(method = "fdr") %>%
                      rstatix::add_significance("p.adj")



  # Add p-values onto the bar plots
  pval_tibble <- pval_tibble %>%
                    rstatix::add_xy_position(fun = "mean_sd", x = "sample", scales="free") %>%
                    dplyr::mutate(y.position=y.position-0.06)


  # color palette
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # plot barplot
  if(is.null(display_pval)==TRUE){
    bp <- ggpubr::ggbarplot(fc_tibble_group, x = "sample", y = "value", add = "mean_sd",
                            fill= "sample", palette = cbp1,facet.by = "gene",
                            position = position_dodge(0.8), scales="free")


  }else if(display_pval == TRUE ){
    bp <- ggpubr::ggbarplot(fc_tibble_group, x = "sample", y = "value", add = "mean_sd",
                            fill= "sample", palette = cbp1,facet.by = "gene",
                            position = position_dodge(0.8), scales="free") +
      ggpubr::stat_pvalue_manual(pval_tibble,  label = "p", tip.length = 0.01)


  }else {
    bp <- ggpubr::ggbarplot(fc_tibble_group, x = "sample", y = "value", add = "mean_sd",
                            fill= "sample", palette = cbp1,facet.by = "gene",
                            position = position_dodge(0.8), scales="free") +
      ggpubr::stat_pvalue_manual(pval_tibble,  label = "pval_significance", tip.length = 0.01)
  }

  bp <- bp +
    ggplot2::ylab(expression(bold(2^-("\u0394\u0394C"[T]))))+
    ggplot2::facet_wrap(~gene, nrow=nrow, scales="free")+
    ggplot2::theme(axis.text.y = ggplot2::element_text(color="black", size=12, face="bold"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(color="black", size=12, face="bold"),
                   strip.text = ggplot2::element_text(color="black", size=12, face="bold"),
                   strip.background = ggplot2::element_rect(fill = "white"))

  return(bp)
}

