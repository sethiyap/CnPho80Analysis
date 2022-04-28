
#' qPCR_foldchange_heatmap
#'
#' @param fc_tibble tibble of fold-change values determined as 2^-∆∆Ct, can be
#'   output of compute_qPCR_foldchange function. With rows as sample name ie.,
#'   WT_R1, WT_R2, WT_R3, treatment_R1,... and columns are genes. Replicates
#'   should always be in the format "_R1" or "_1" or "_set1"
#' @param control character vector, determining name of the control samples eg.
#'   WT
#' @param treatment character vector, determining name of the control samples
#' @param threshold logical, to limit the color scale of heatmap. default: TRUE
#'   eg. treatment
#' @import tidyr
#' @import dplyr
#' @import ggpubr
#' @import rstatix
#' @import forcats
#' @import stats
#' @import utils
#' @import tibble
#'
#' @return heatmap of control and treatment pvalues
#' @export
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
#' qPCR_foldchange_heatmap(fc_tibble = qPCR_fc)
#'
qPCR_foldchange_heatmap <- function(fc_tibble, control="WT", treatment="WT+drug", threshold=TRUE){

  # First column name, always sample
  names(fc_tibble)[1] = "sample"

  fc_tibble_group <- fc_tibble %>%
                        tidyr::separate(col = sample, into = c("sample", "replicate"), sep = "_") %>%
                        dplyr::select(-c(replicate)) %>%
                        tidyr::gather(gene, value, -sample)

  fc_tibble_mean <- fc_tibble_group %>%
                        dplyr::group_by(sample, gene) %>%
                        dplyr::summarise(mean=mean(value))

  pval_tibble <- fc_tibble_group %>%
    dplyr::mutate(gene=forcats::as_factor(gene)) %>%
    dplyr::group_by(gene) %>%
    rstatix::t_test(value ~ sample, paired = FALSE) %>%
    dplyr::mutate(p=p/2,
                  # left tailed: treatment < control, right tailed: treatment > control
                  # However, in population of genes if some are control > treatment or control < treatment.
                  # then converting the two-tailed pvalue into one-tailed useful (by dividing the two-tailed pvalue into half)
                  # as when we don't know the direction of population

                  p_significance= dplyr::if_else(p < 0.05 & p > 0.01,"*",
                                                    dplyr::if_else(p < 0.01 & p > 0.001, "**", # adding pvalue significance
                                                                   dplyr::if_else(p < 0.001 & p > 0.0001, "***", "ns")))) %>%

    rstatix::adjust_pvalue(method = "fdr") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(fun = "mean_sd", x = "sample", scales="free") %>%
    dplyr::mutate(y.position=y.position-0.06)

  fc_tibble_mean <- fc_tibble_mean %>%
                      dplyr::inner_join(pval_tibble, by="gene") %>%
                      dplyr::select(sample, gene, mean, p, p_significance) %>%
                      dplyr::mutate(p=replace(p, sample==control, ""),
                                    p_significance=replace(p_significance, sample==control, ""))

  if(threshold==TRUE){

          readnumber <- function()
          {
            n <- readline(prompt="Enter threshold to limit the color scale: ")
            n <- as.double(n)
            if (is.na(n)){
              n <- readnumber()
            }
            return(n)
          }

    read_num <- readnumber()
    fc_tibble_mean <- fc_tibble_mean %>%
      dplyr::mutate(value=dplyr::if_else(mean >= read_num,read_num,mean))
  }
  else{
    fc_tibble_mean <- fc_tibble_mean
  }

    gg_tile <- ggplot2::ggplot(fc_tibble_mean, ggplot2::aes(sample, gene, fill=value, label=p))+
                        ggplot2::geom_tile(color="white", size=1)+
                        ggplot2::scale_fill_gradient(low = "#00AFBB", high ="#E7B800", trans="log2")+
                        ggplot2::geom_text(color="black", size=5)+
                        ggplot2::theme_bw()

    gg_tile <- gg_tile + ggplot2::theme(axis.text = ggplot2::element_text(color="black", size=12, face="bold"))+
      ggplot2::guides(fill = ggplot2::guide_colourbar(title.position="top", title.hjust = 0.5,title =expression(bold(2^-("\u0394\u0394C"[T])))))

    return(gg_tile)
}
