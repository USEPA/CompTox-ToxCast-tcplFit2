#' Sample concentration-response data set from invitrodb
#'
#' A data set containing 100 chemicals worth of data for the Tox21 assay
#' TOX21_ERa_BLA_Agonist_ratio, whcih measures response to  estrogen receptor agonists.
#' The data can be accessed further through the Comptox Chemicals Dashboard:
#' https://comptox.epa.gov/dashboard
#'
#' This data is extracted from the database invitrodb, at level 3 (conc-response data)
#'
#' A data frame with 32175 rows and 6 variables:
#'   \itemize{
#'   \item dtxsid - DSSTox generic substance ID
#'   \item casrn - Chemical Abstracts Registry Number (CASRN)
#'   \item name- chemical name
#'   \item spid - sample ID - ther can be multiple samples per chemical
#'   \item logc - log10(concentraiotn uM)
#'   \item resp - response in %
#'   \item assay - name of the assay / assay component endpoint
#'   ...
#' }
"mc3"

