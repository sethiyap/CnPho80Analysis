#' plot_sruvival_curve
#'
#' @param survival_data status should be 1 or 0, 1=dead, 0=alive
#' @param palette character vector to determine a specific color palette or NULL for default palette
#' @param display_pval logical, to display pvalue or not
#'
#' @return survival plot
#' @export
#'
#' @examples
#' @import dplyr
#' @import survival
#' @import ggplot2
#' @import rlang
#' @import survminer
#'
#'
#'
plot_survival_curve <- function(survival_data, palette=NULL, display_pval=TRUE){

  library(magrittr)
  # column names
  names(survival_data) = c("sample", "day", "status")

  # factor sample
  survival_data <- survival_data %>% dplyr::mutate(sample=forcats::as_factor(sample))

  if(is.null(palette)==TRUE){
    palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  }else{
    palette = palette
  }

  # calculate survival probability

  fit <- survival::survfit(survival::Surv(day, status) ~ sample, data = survival_data)

  # plot survival data
  if(display_pval==TRUE){
    survival_plot <- survminer::ggsurvplot(fit, data = survival_data, pval = TRUE, ggtheme = ggplot2::theme_bw(),
                                           palette = palette,
                                           font.main = c(16, "bold", "darkblue"),
                                           font.x = c(14, "bold", "black"),
                                           font.y = c(14, "bold", "black"),
                                           font.tickslab = c(12, "plain", "black"),
                                           font.legend = c(12, "bold", "black"),
                                           legend="right",
                                           legend.labs=gsub(pattern = "sample=", replacement = "", x = names(fit$strata)))
  }else{
    survival_plot <- survminer::ggsurvplot(fit, data = survival_data, pval = FALSE, ggtheme = ggplot2::theme_bw(),
                                           palette = palette,
                                           font.main = c(16, "bold", "darkblue"),
                                           font.x = c(14, "bold", "black"),
                                           font.y = c(14, "bold", "black"),
                                           font.tickslab = c(12, "plain", "black"),
                                           font.legend = c(12, "bold", "black"),
                                           legend="right",
                                           legend.labs=gsub(pattern = "sample=", replacement = "", x = names(fit$strata)))
  }



  return(survival_plot)


}

