#' Sample concentration-response data set from invitrodb
#'
#' A data set containing 100 chemicals worth of data for the Tox21 assay
#' TOX21_ERa_BLA_Agonist_ratio, which measures response to estrogen receptor agonists.
#' The data can be accessed further through the Comptox Chemicals Dashboard
#' (\url{https://comptox.epa.gov/dashboard}).
#'
#' This data is extracted from the released version of the ToxCast database,
#' invitrodb, at level 3 (mc3) and contains the concentration-response information.
#'
#' A data frame with 32175 rows and 7 variables:
#'   \itemize{
#'   \item dtxsid - DSSTox generic substance ID
#'   \item casrn - Chemical Abstracts Registry Number (CASRN)
#'   \item name - chemical name
#'   \item spid - sample ID - there can be multiple samples per chemical
#'   \item logc - log10(concentration), micromolar (uM)
#'   \item resp - response in \%
#'   \item assay - name of the assay / assay component endpoint name
#' }
#'
#' @source \url{https://doi.org/10.23645/epacomptox.6062623.v5}
"mc3"

#' Sample concentration-response data set from HTTR
#'
#' A data set containing 6 of the active transcriptional signatures after
#' perturbation of MCF7 cells with Clomiphene citrate (1:1).
#'
#' A data frame with 6 rows and 8 variables:
#' \itemize{
#'   \item sample_id - experimental sample ID
#'   \item dtxsid - DSSTox generic substance ID
#'   \item name - chemical name
#'   \item signature - transcriptional signature name
#'   \item cutoff - the 95\% confidence interval from the baseline response (2 lowest concentrations)
#'   \item onesd - one standard deviation of the baseline response
#'   \item conc - experimental concentrations, micromolar (uM)
#'   \item resp - transcriptional signature response for each experimental concentrations, ssGSEA score
#' }
#'
#' @references Joshua A. Harrill, Logan J. Everett, Derik E. Haggard,
#'   Thomas Sheffield, Joseph L. Bundy, Clinton M. Willis, Russell S. Thomas,
#'   Imran Shah, Richard S. Judson, High-Throughput Transcriptomics Platform for
#'   Screening Environmental Chemicals, Toxicological Sciences, Volume 181,
#'   Issue 1, May 2021, Pages 68 - 89, \url{https://doi.org/10.1093/toxsci/kfab009}.
"signatures"
