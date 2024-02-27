#' Dichotomous Model Functions
#'
#' Implementations of dichotomous model. Models include: Multistage, Weibull,
#' Gamma, Logistic, Log-Logistic, Probit, Log-Probit, and Dichotomous Hill.
#'
#' To stay inline with the current assumption that baseline response is 0,
#' g (background) is assumed to be 0 and is dropped from all implementations.
#' However it is to be noted that in BMDS User Guide pg. 91 it warns the users
#' to not assume 0 baseline if they have more than zero responses in a dose
#' group that has a dose value of zero. We will likely to revisit this part
#' in the future. But as of now, g = 0.
#'
#' parameter vector ps <- c(a, b, v)
#'
#'

#' Multistage Model
#'
#' @param beta Vector of dose coefficients
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#'
multistage <- function(beta, x){
  # beta is a vector of dose coefficients

  k <- seq(1, length(beta))

  r <- lapply(x, function(i) sum(beta*(i^k)))
  return ( 1 - exp(-unlist(r))  )

}

#' Weibull Model
#'
#' @param ps Vector of parameters: a,b,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#'
weibull <- function(ps, x){
  # a = ps[1], b = ps[2]
  return( 1 - exp(-ps[2]*(x^ps[1]))  )
}


#' Integral Function side Gamma model
#'
#' @param ps Vector of parameters: a,b
#' @param x Vector of concentrations (regular units)
#'
#' @return Function value
#' @export
#'
gammainner <- function(a, t){
  # a = ps[1], b = ps[2]
  t^(a-1)*exp(-t)
  }

#' Gamma Model
#'
#' @param ps Vector of parameters: a, b
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#'
gamma_d <- function(ps, x){
  # a = ps[1], b = ps[2]
  #r <- lapply(x, function(i) integrate(gammainner,
  #                                     lower=0,
  #                                     upper=ps[2]*i,
  #                                     a = ps[1])$value
  #            )
  r <- pgamma(x, ps[1], ps[2])
  return( r )
}


#' Logistic Model
#'
#' @param ps Vector of parameters: a, b
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#'
logistic <- function(ps, x){
  # a = ps[1], b = ps[2]
  return( 1 / (1 + exp(-ps[1]-ps[2]*x) ) )
}


#' Log-Logistic Model
#'
#' @param ps Vector of parameters: a, b
#' @param x Vector of concentrations (regular units, convert to natrual-log scale in function)
#'
#' @return Vector of model responses
#' @export
#'
llogistic <- function(ps, x){
  # a = ps[1], b = ps[2]
  return( 1 / (1+exp(-ps[1]-ps[2]*log(x)))  )

}

#' probitinner
#'
#' @param t
#'
#' @return function value
#' @export
#'
probitinner <- function(t) {(1/sqrt(2*pi))*exp(-(t^2)/2)}


#' Probit Model
#'
#' @param ps Vector of parameters: a, b
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#'
probit <- function(ps, x) {
  # a = ps[1], b = ps[2]
  x <- ps[1] + ps[2]*x
  #r <- lapply(x, function(i) integrate(probitinner,
  #                                     lower=-Inf,
  #                                     upper=i)$value
  #)
  r <- pnorm(x)
  return(r)
}


#' Log-Probit Model
#'
#' @param ps Vector of parameters: a, b
#' @param x Vector of concentrations (regular units, convert to natural-log scale in function)
#'
#' @return Vector of model responses
#' @export
#'

lprobit <- function(ps, x) {
  # a = ps[1], b = ps[2]
  x <- ps[1] + ps[2]*log(x)
  #r <- lapply(x, function(i) integrate(probitinner,
  #                                     lower=-Inf,
  #                                     upper=i)$value
  #)
  r <- pnorm(x)
  return(r)
}


#' Dichotomous Hill Model
#'
#' @param ps Vector of parameters: a, b, v
#' @param x Vector of concentrations (regular units, convert to natural-log scale in function)
#'
#' @return Vector of model responses
#' @export
#'
d_hill <- function(ps, x) {
  # a = ps[1], b = ps[2], v = ps[3]

  ps[3] / (1 + exp(-ps[1]-ps[2]*log(x)))

}


