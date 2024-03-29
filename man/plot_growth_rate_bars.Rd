% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_growth_rate_bars.R
\name{plot_growth_rate_bars}
\alias{plot_growth_rate_bars}
\title{plot_growth_rate_bars}
\usage{
plot_growth_rate_bars(dat_gr = dat_gr, is_replicates = FALSE, is_facet = FALSE)
}
\arguments{
\item{dat_gr}{tibble of growth data in the format such that first column is
Time and then all the samples containing OD. If replicates add after '_'. example: WT_Set1 or WT-drug_Set1}

\item{is_replicates}{logical to indicate if the data contains replicates}

\item{is_facet}{logical to facet the data according to samples}
}
\value{
ggplot of OD600 and percent growth. Tibble of percent growth rates
}
\description{
plot_growth_rate_bars
}
\examples{
library(magrittr)
WT_R1 <- c(0.14, 59.2, 41.3)
WT_R2 <- c(0.14, 54.8, 37.7)
time <- c("0h", "24h", "48h")

OD_dat <- tidyr::tibble(WT_R1, WT_R2)
rownames(OD_dat) <- time

OD_dat <- OD_dat \%>\% tibble::rownames_to_column(var = "time")


plot_growth_rate_bars(dat_gr = OD_dat, is_replicates = TRUE, is_facet = FALSE)

}
