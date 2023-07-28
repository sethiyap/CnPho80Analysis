#' plot_pvales
#'
#' @param dat_tble tibble or data frame with rows=sample and column=conditions (+Pi or -Pi)
#' @param display_pval logical, to display pvalue or not. TRUE: Pvalue determined by t-test, FALSE: asterisk, NULL: no pvalue displayed
#' @param control_group character vector, which sample to be considered as control to determine pvalue. eg. "WT" *special characters will through error
#' @param scales character vector, to keep free scales for each group on y-axis or not. eg. "free" or NULL
#' @param y_label character vector, label displayed on y-axis.
#' @param nrow numeric, number of rows in which groups should be displayed
#'
#' @return barplot with pvalues displayed
#' @export
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import rlang
#' @import ggrepel
#' @import forcats
#' @import rstatix
#' @examples
#'
#' library(magrittr)
#' gene_1 <- c(1.08, 1.01 , 0.917, 0.401, 0.356, 0.450)
#' gene_2 <- c(1.02, 1.03, 1.002,0.8, 0.76, 0.79)
#' rname <- c(paste("WT", seq(1:3), sep = "_"), paste("WT+drug", seq(1:3), sep = "_"))
#' tibble_fc <- tidyr::tibble(gene_1, gene_2)
#' rownames(tibble_fc) <- rname
#'
#' tibble_fc <- tibble_fc %>% tibble::rownames_to_column(var = "sample")
#'
#' plot_pvalues(dat_tble = tibble_fc, display_pval = FALSE, control_group = "WT", scales = NULL, y_label = "acid phosphatase activity")
#'
plot_pvalues <- function(dat_tble, display_pval=NULL, control_group="WT", scales="free", y_label="fold_change", nrow=1){

  library(ggpubr)
  library(magrittr)
  # First column name, always sample
  names(dat_tble)[1] = "sample"

  dat_tble_group <- dat_tble %>%
                      tidyr::separate(col = sample, into = c("sample", "replicate"), sep = "_") %>%
                      dplyr::select(-c(replicate)) %>%
                      tidyr::gather(condition, value, -sample) %>%
                      dplyr::mutate(condition=forcats::as_factor(condition))

  pval_tibble <- dat_tble_group %>%
                    dplyr::mutate(sample=forcats::as_factor(sample)) %>%
                    dplyr::group_by(condition) %>%
                    rstatix::t_test(value ~ sample, paired = FALSE, ref.group = control_group) %>%
                    dplyr::mutate(p=p/2,
                                  # left tailed: treatment < control, right tailed: treatment > control
                                  # However, in population of genes if some are control > treatment or control < treatment.
                                  # then converting the two-tailed pvalue into one-tailed useful (by dividing the two-tailed pvalue into half)
                                  # as when we don't know the direction of population

                                  pval_significance= dplyr::if_else(p < 0.05 & p > 0.01,"*",
                                                                    dplyr::if_else(p < 0.01 & p > 0.001, "**", # adding pvalue significance
                                                                                   dplyr::if_else(p < 0.001 & p > 0.0001, "***",
                                                                                                  dplyr::if_else(p< 0.0001, "****", "ns"))))) %>%

                    rstatix::adjust_pvalue(method = "fdr") %>%
                    rstatix::add_significance("p.adj")

  cbp1 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

  # Add p-values onto the bar plots
  pval_tibble <- pval_tibble %>%
                      rstatix::add_xy_position(fun = "mean_sd", x = "sample", scales=scales) %>%
                      dplyr::mutate(y.position=y.position-0.06)

  if(is.null(display_pval)==TRUE){
    bp <- ggpubr::ggbarplot(dat_tble_group, x = "sample", y = "value", add = "mean_sd",
                            fill= "sample", palette = cbp1,facet.by = "condition",
                            position = position_dodge(0.8), scales=scales)


  }else if(display_pval == TRUE ){
    bp <- ggpubr::ggbarplot(dat_tble_group, x = "sample", y = "value", add = "mean_sd",
                            fill= "sample", palette = cbp1,facet.by = "condition",
                            position = position_dodge(0.8), scales=scales) +
      ggpubr::stat_pvalue_manual(pval_tibble,  label = "p", tip.length = 0.01)


  }else {
    bp <- ggpubr::ggbarplot(dat_tble_group, x = "sample", y = "value", add = "mean_sd",
                            fill= "sample", palette = cbp1,facet.by = "condition",
                            position = position_dodge(0.8), scales=scales) +
      ggpubr::stat_pvalue_manual(pval_tibble,  label = "pval_significance", tip.length = 0.01)
  }

  bp <- bp +
    ggplot2::ylab(label = y_label)+
    ggplot2::facet_wrap(~condition, nrow=nrow, scales=scales)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(color="black", size=12, face="bold"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(color="black", size=12, face="bold"),
                   strip.text = ggplot2::element_text(color="black", size=12, face="bold"),
                   strip.background = ggplot2::element_rect(fill = "white"))

  return(bp)
}
