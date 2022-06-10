#' compare_two_rnaseq
#'
#' @param rnaseq_dat tibble, table of multiple rnaseq data log2 fold-change.
#'   first column should be gene-d, then data and then label (optional)
#' @param data_1 string, first data to be plotted in comparison
#' @param data_2 string, second data to be plotted in comparison
#'
#' @return a scatterplot highlighting commonly up- and down- regulated genes
#'
#' @export
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import rlang
#' @import ggrepel
#'
#'
#' @examples
#'
#'
compare_two_rnaseq <- function(rnaseq_dat, data_1, data_2){

  names(rnaseq_dat)[1] = "gene"

  if(isTRUE(tidyr::contains(colnames(rnaseq_dat),vars = "label")==1)){
    dat <- rnaseq_dat %>%
                        dplyr::select(gene,!!rlang::sym(data_1),!!rlang::sym(data_2), label) %>%
                        tidyr::gather(sample, value, -gene, -label)


  }else{
    dat <- rnaseq_dat %>%
                        dplyr::select(gene,!!rlang::sym(data_1),!!rlang::sym(data_2)) %>%
                        tidyr::gather(sample, value, -gene)

  }
  scatterplot_dat <- dat %>%
                        dplyr::group_by(gene) %>%
                        dplyr::filter(sum(abs(value)) > 0) %>%
                        tidyr::spread(sample, value) %>%
                        dplyr::mutate(
                          DEG=dplyr::if_else(!!rlang::sym(data_1) >= 1 & !!rlang::sym(data_2) >= 1, "up in both",
                                             dplyr::if_else(!!rlang::sym(data_1) <= -1 & !!rlang::sym(data_2) <= -1,"down in both", "no-change")))

  if(isTRUE(tidyr::contains(colnames(rnaseq_dat),vars = "label")==1)){
  gg= ggplot2::ggplot(scatterplot_dat, ggplot2::aes(!!rlang::sym(data_1), !!rlang::sym(data_2),  color=DEG, label=label))  +
    ggplot2::geom_point()+
    ggrepel::geom_text_repel(color="black", box.padding = 0.5, max.overlaps = 12)
  }else{
    gg= ggplot2::ggplot(scatterplot_dat, ggplot2::aes(!!rlang::sym(data_1), !!rlang::sym(data_2),  color=DEG))  +
      ggplot2::geom_point()
  }

gg = gg+
    ggplot2::scale_color_manual(values = c("blue", "grey", "red"))+
    ggplot2::theme_bw()+
    ggplot2::lims(x=c(-9, 9), y=c(-9,9))+
    ggplot2::geom_hline(yintercept = c(-1,1), linetype=2, color="gray45")+
    ggplot2::geom_vline(xintercept = c(-1,1), linetype=2, color="gray45")+
    ggplot2::theme(axis.text = ggplot2::element_text(color="black", size = 12),
                   legend.text = ggplot2::element_text(color="black", size = 12))

return(gg)

}
