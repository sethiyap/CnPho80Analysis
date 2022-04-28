
#' plot_growth_rate_bars
#'
#' @param dat_gr tibble of growth data in the format such that first column is
#' Time and then all the samples containing OD. If replicates add after '_'. example: WT_Set1 or WT-drug_Set1
#' @param is_replicates logical to indicate if the data contains replicates
#' @param is_facet logical to facet the data according to samples
#'
#' @return ggplot of OD600 and percent growth. Tibble of percent growth rates
#' @export
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import forcats
#' @import stats
#' @import utils
#' @import tibble
#' @examples
#' library(magrittr)
#' WT_R1 <- c(0.14, 59.2, 41.3)
#' WT_R2 <- c(0.14, 54.8, 37.7)
#' time <- c("0h", "24h", "48h")
#'
#' OD_dat <- tidyr::tibble(WT_R1, WT_R2)
#' rownames(OD_dat) <- time
#'
#' OD_dat <- OD_dat %>% tibble::rownames_to_column(var = "time")
#'
#'
#' plot_growth_rate_bars(dat_gr = OD_dat, is_replicates = TRUE, is_facet = FALSE)
#'
plot_growth_rate_bars <- function(dat_gr=dat_gr, is_replicates=FALSE, is_facet=FALSE){

  # First column name, always Time
  names(dat_gr)[1] = "Time"

  # color palette
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # number of colors
  legend_order <-  length(unique(dat_gr$Time))+1

  od_aes = function(gg = plot, color_pal = cbp1, is_facet=is_facet){

    if(isFALSE(is_facet)){
      gg_od <- gg + ggplot2::theme_bw() +
        ggplot2::ylab(expression(OD[600]))+
        ggplot2::theme(axis.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(color="black", size=12, face="bold"))+
        ggplot2::scale_fill_manual(values = color_pal)
    }else{
      gg_od <- gg + ggplot2::theme_bw()+
        ggplot2::facet_wrap(~condition, nrow = 1, scales="free_x")+
        ggplot2::ylab(expression(OD[600]))+
        ggplot2::theme(axis.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       strip.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       strip.background = ggplot2::element_rect(fill = "white"))+
        ggplot2::scale_fill_manual(values = color_pal)

    }

    return(gg_od)

  }

  percent_aes <- function(gg=plot, color_pal = cbp1, is_facet=is_facet, legend_order= legened_order){

    if(isFALSE(is_facet)){
      gg_od <- gg + ggplot2::theme_bw() +
        ggplot2::ylab("percent_growth_rate wrt 0h")+
        ggplot2::theme(axis.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(color="black", size=12, face="bold"))+
        ggplot2::scale_fill_manual(values=color_pal[2:legend_order])
    }else{
      gg_od <- gg + ggplot2::theme_bw()+
        ggplot2::facet_wrap(~condition, nrow = 1, scales="free_x")+
        ggplot2::ylab("percent_growth_rate wrt 0h")+
        ggplot2::theme(axis.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       strip.text = ggplot2::element_text(color="black", size=12, face="bold"),
                       strip.background = ggplot2::element_rect(fill = "white"))+
        ggplot2::scale_fill_manual(values=color_pal[2:legend_order])

    }

  }

  if(is_replicates==FALSE){


    # data grouping
    gr_group <- dat_gr %>%
      tidyr::gather(condition, values, -Time) %>%
      dplyr::mutate(condition=forcats::as_factor(condition),
                    Time = forcats::as_factor(Time)) %>%
      dplyr::group_by(condition)


    # plot OD600
    gg_od <- gr_group %>%
      ggplot2::ggplot(ggplot2::aes(condition, values, fill=Time, label=round(values,2)))+
      ggplot2::geom_col(position=ggplot2::position_dodge(preserve = "total"), color="black")+
      ggplot2::geom_text(position = ggplot2::position_dodge(width = .9), vjust=-0.3)

    gg_od <- od_aes(gg = gg_od, color_pal = cbp1, is_facet = is_facet)

    # percent growth rate calculation
    gr_calculator <-    gr_group %>%
      dplyr::mutate(time = stringr::str_replace_all(string = Time, pattern = "[^0-9.-]", replacement = ""),
                    time=as.numeric(time)) %>%
      dplyr::mutate(gr = (values - dplyr::first(values))/dplyr::first(values)) %>% # subtracts from 0hr time-point
      dplyr::mutate(gr_time = gr/time, percent_growth = gr_time*100) %>%
      tidyr::drop_na() %>%
      dplyr::ungroup() %>%
      dplyr::select(c(Time, condition, percent_growth))

    gg_percent <- gr_calculator %>%
      ggplot2::ggplot(ggplot2::aes(condition, percent_growth, fill=Time, label=round(percent_growth,2)))+
      ggplot2::geom_col(position=ggplot2::position_dodge(preserve = "total"), color="black")

    percent_plot <- percent_aes(gg = gg_percent, color_pal = cbp1, is_facet = is_facet,legend_order = legend_order)

  } else{

    # data grouping
    dat_gr_rep <- dat_gr %>%
      tidyr::gather(replicate, values, -Time) %>%
      tidyr::separate(col = replicate, into = c("condition", "replicate"), sep = "_") %>%
      dplyr::mutate(condition=forcats::as_factor(condition), Time=forcats::as_factor(Time))

    # calculate SD and mean for OD
    dat_od_sd <- dat_gr_rep %>%
      dplyr::group_by(condition, Time) %>%
      dplyr::mutate(mean_od= mean(values), sd_od = sd(values)) %>%
      dplyr::filter(replicate==max(replicate)) %>%
      dplyr::select(-c(replicate, values))

    # plot OD
    gg_od <- dat_od_sd %>%
      ggplot2::ggplot(ggplot2::aes(condition, mean_od, fill=Time))+
      ggplot2::geom_col(position=ggplot2::position_dodge(preserve = "total"))+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_od-sd_od, ymax=mean_od+sd_od),
                             width=.2,size=0.6, position=ggplot2::position_dodge(0.9))

    gg_od <- od_aes(gg = gg_od, color_pal = cbp1, is_facet = is_facet)


    # calculate SD and mean for percent
    gr_calculator <- dat_gr_rep %>%
      dplyr::group_by(condition, replicate) %>%
      dplyr::mutate(time = stringr::str_replace_all(string = Time, pattern = "[^0-9.-]", replacement = ""),
                    time=as.numeric(time)) %>%
      dplyr::mutate(gr = (values - dplyr::first(values))/dplyr::first(values)) %>% # subtracts from 0hr time-point
      dplyr::mutate(gr_time = gr/time, percent_growth = gr_time*100) %>%
      tidyr::drop_na() %>%
      dplyr::ungroup() %>%
      dplyr::group_by(condition, Time) %>%
      dplyr::mutate(mean_percent=mean(percent_growth), sd_percent = sd(percent_growth)) %>%
      dplyr::filter(replicate==max(replicate)) %>%
      dplyr::select(-c(replicate, values, gr, gr_time))

    # plot growth calculator
    gg_pr <- gr_calculator %>%
      ggplot2::ggplot(ggplot2::aes(condition,mean_percent , fill=Time))+
      ggplot2::geom_col(position=ggplot2::position_dodge(preserve = "total"))+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_percent-sd_percent, ymax=mean_percent+sd_percent),
                             width=.2,size=0.6, position=ggplot2::position_dodge(0.9))

    percent_plot <- percent_aes(gg = gg_pr, color_pal = cbp1, is_facet = is_facet,legend_order = legend_order)
  }

  gr_calculator <- gr_calculator %>% dplyr::select(condition, percent_growth, Time) %>% tidyr::spread(condition, percent_growth)
  message ("plotting OD...")
  print(gg_od)

  message ("plotting percent growth rate...")
  print(percent_plot)

  return(gr_calculator)
}

