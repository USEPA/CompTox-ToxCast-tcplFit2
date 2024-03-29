% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mc3}
\alias{mc3}
\title{Sample concentration-response data set from invitrodb}
\format{
An object of class \code{data.frame} with 32175 rows and 7 columns.
}
\source{
\doi{https://doi.org/10.23645/epacomptox.6062623.v5}
}
\usage{
mc3
}
\description{
A data set containing 100 chemicals worth of data for the Tox21 assay
TOX21_ERa_BLA_Agonist_ratio, which measures response to estrogen receptor agonists.
The data can be accessed further through the Comptox Chemicals Dashboard
(\url{https://comptox.epa.gov/dashboard}).
}
\details{
This data is extracted from the released version of the ToxCast database,
invitrodb, at level 3 (mc3) and contains the concentration-response information.

A data frame with 32175 rows and 7 variables:
  \itemize{
  \item dtxsid - DSSTox generic substance ID
  \item casrn - Chemical Abstracts Registry Number (CASRN)
  \item name - chemical name
  \item spid - sample ID - there can be multiple samples per chemical
  \item logc - log10(concentration), micromolar (uM)
  \item resp - response in \%
  \item assay - name of the assay / assay component endpoint name
}
}
\keyword{datasets}
