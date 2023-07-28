
#' compute_qPCR_foldchange
#'
#' @param dat_qPCR tibble, a tibble of qPCR output. With rows as sample name
#'   ie., WT_R1, WT_R2, WT_R3, treatment_R1,... and columns are genes.
#'   Replicates should always be in the format "_R1" or "_1" or "_set1"
#' @param house_keeping character vector, containing name of the house-keeping
#'   control eg. "TUB2" for tubulin gene.
#' @param control character vector, determining name of the control samples eg.
#'   WT
#' @param treatment character vector, determining name of the treatment samples
#'   eg. treatment
#' @param return_pval logical, to return pvalue or tibble of fold-change
#' @return ggplot of the fold change and tibble of foldchange, pvalue
#'   (determined by t-test)
#' @export
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @import purrr
#' @import stats
#' @import utils
#' @import tibble
#' @import forcats
#'
#' @details pvalue calculation
#'   source:\url{https://stats.oarc.ucla.edu/other/mult-pkg/faq/pvalue-htm/}
#' @examples
#' library(magrittr)
#' gene_of_interest <- c(22.1, 22.3, 22.4, 8, 8.9, 8.3)
#' house_keeping_gene <- c(20, 20.1, 20.3,20, 20.4, 19.99)
#' rname <- c(paste("WT", seq(1:3), sep = "_"), paste("WT+drug", seq(1:3), sep = "_"))
#'
#' dt_qPCR <- tidyr::tibble(gene_of_interest, house_keeping_gene)
#' rownames(dt_qPCR) <- rname
#'
#' dt_qPCR <- dt_qPCR %>% tibble::rownames_to_column(var = "Sample")
#'
#' compute_qPCR_foldchange(dat_qPCR = dt_qPCR, house_keeping = "house_keeping_gene",
#' control = "WT", treatment = "WT+drug", return_pval=TRUE)
#'
compute_qPCR_foldchange <- function(dat_qPCR=dat_qPCR, house_keeping="TUB2", control="WT", treatment="WT+drug", return_pval=TRUE){

  # press cmd+shift+alt+R to insert roxygen
  # devtools::document() to roxygenize into html

  # First column name, always sample
  names(dat_qPCR)[1] = "sample"

  target <- c(control, treatment)
  dCt_df <- dat_qPCR %>%
    tidyr::gather(gene, value,-sample, -house_keeping) %>% # add house-keeping gene Ct value into new column
    dplyr::mutate(sample=forcats::as_factor(sample)) %>%
    dplyr::mutate(dCt = value - !!rlang::sym(house_keeping)) # calculate ∆Ct=Ct(gene of interest) - Ct(house-keeping)

  ddCt_df <- dCt_df %>%
    dplyr::select(c(sample, gene, dCt)) %>%
    tidyr::separate(sample, into = c("sample", "replicate"), sep = "_") %>%
    tidyr::spread(sample, dCt) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(avg_control = mean(!!rlang::sym(control))) %>% # average control ∆Ct value
    tidyr::gather(sample, dCt, -avg_control, -gene, -replicate) %>%
    dplyr::mutate(ddCt = dCt-avg_control, fold_change = 2^-ddCt) %>% # calculate ∆∆Ct = ∆Ct (gene of interest) - Control ∆Ct value average
    dplyr::arrange(factor(sample, levels = target)) %>%
    dplyr::mutate(sample=forcats::as_factor(sample))

  # return the foldchange value
  fc_val <- ddCt_df %>%
    dplyr::select(c(gene, sample, replicate, fold_change)) %>%
    tidyr::spread(gene, fold_change) %>%
    tidyr::unite(col = "sample",sample:replicate, sep = "_")

  if(return_pval==TRUE){

    foldchange_pvalue <- ddCt_df %>%
      dplyr::select(c(gene, sample, replicate, fold_change)) %>%
      tidyr::spread(sample, fold_change) %>%
      dplyr::group_by(gene) %>%
      tidyr::nest() %>%
      dplyr::mutate(pvalue= purrr::map(data, function(x){ # compute pvalue
            cntrl <- x %>% dplyr::pull(!!rlang::sym(control))
            treat <- x %>% dplyr::pull(!!rlang::sym(treatment))

        # paired t-test is used to determine significance
        # as we are testing the expression of same population across different conditions
        # also divide by 2 to get one-tail significance

        tt <- t.test(cntrl,treat,exact=FALSE)
        tt_p <- tt$p.value

        # left tailed: treatment < control, right tailed: treatment > control
        # However, in population of genes if some are control > treatment or control < treatment.
        # then converting the two-tailed pvalue into one-tailed useful (by dividing the two-tailed pvalue into half)
        # as when we don't know the direction of population

        return(tt_p/2)

      })) %>%
      tidyr::unnest_wider(pvalue, simplify = TRUE, names_sep = "_") %>%
      dplyr::mutate(pvalue=pvalue_1, significance=dplyr::if_else(pvalue < 0.05 & pvalue > 0.01,"*",
                                                dplyr::if_else(pvalue < 0.01 & pvalue > 0.001, "**", # adding pvalue significance
                                                               dplyr::if_else(pvalue < 0.001 & pvalue > 0.0001, "***",
                                                                              dplyr::if_else(pvalue< 0.0001, "****", "ns"))))) %>%
      dplyr::select(c(gene, pvalue, significance))

    return(foldchange_pvalue)
  }else{
    return(fc_val)
  }
}
