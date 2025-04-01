#' Top Likelihood
#'
#' Probability of top being above cutoff.
#'
#' Should only be called by hitcontinner. Uses profile likelihood, similar
#' to bmdbounds. Here, the y-scale type parameter is substituted in such a
#' way that the top equals the cutoff. Then the log-likelihood is compared to
#' the maximum log-likelihood using chisq function to retrieve probability.
#'
#' @param fname Model function name (equal to model name except hill which
#'   uses "hillfn")
#' @param cutoff Desired cutoff.
#' @param conc Vector of concentrations.
#' @param resp Vector of responses.
#' @param ps Vector of parameters, must be in order: a, tp, b, ga, p, la, q, er
#' @param top Model top.
#' @param mll Winning model maximum log-likelihood.
#' @param errfun Which error distribution to assume for each point, defaults to
#'   "dt4". "dt4" is the original 4 degrees of freedom t-distribution. Another
#'   supported distribution is "dnorm", the normal distribution.
#' @param poly2.biphasic Which fitting method to use for poly2. If poly2.biphasic = TRUE, allows for biphasic polynomial 2
#'   model fits (i.e. both monotonic and non-monotonic). (Defaults to TRUE.)
#'
#' @importFrom stats pchisq
#'
#' @return Probability of top being above cutoff.
#' @export
#'
#' @examples
#' fname = "hillfn"
#' conc = c(.03,.1,.3,1,3,10,30,100)
#' resp = c(0,.1,0,.2,.6,.9,1.1,1)
#' ps = c(1.033239, 2.453014, 1.592714, er = -3.295307)
#' top = 1.023239
#' mll = 12.71495
#' toplikelihood(fname, cutoff = .8, conc, resp, ps, top, mll)
#' toplikelihood(fname, cutoff = 1, conc, resp, ps, top, mll)
#' toplikelihood(fname, cutoff = 1.2, conc, resp, ps, top, mll)
toplikelihood = function(fname, cutoff, conc, resp, ps, top, mll, errfun = "dt4", poly2.biphasic = TRUE){
  #cutoff needs to account for sign otherwise reparameterization will flip the model
  cutoff = cutoff*sign(top)

  #reparameterize so that top is exactly at cutoff
  if(fname == "exp2"){
    x_top = acy(y = top, modpars = list(a=ps[1],b=ps[2],er=ps[3]),type=fname)
    ps[1] = cutoff/( exp(x_top/ps[2]) - 1 )
  } else if(fname == "exp3"){
    x_top = acy(y = top, modpars = list(a=ps[1],b=ps[2],p=ps[3],er=ps[4]),type=fname)
    ps[1] = cutoff/( exp((x_top/ps[2])^ps[3]) - 1 )
  } else if(fname == "exp4"){
    if (top == ps[1]) {
      ps[1] = cutoff
    } else {
      x_top = acy(y = top, modpars = list(tp=ps[1],ga=ps[2],er=ps[3]),type=fname)
      ps[1] = cutoff/( 1 - 2^(-x_top/ps[2]))
    }
  } else if(fname == "exp5"){
    if (top == ps[1]) {
      ps[1] = cutoff
    } else{
      x_top = acy(y = top, modpars = list(tp=ps[1], ga=ps[2], p=ps[3], er=ps[4]), type = fname)
      ps[1]  = cutoff/ (1-2^(-(x_top/ps[2])^ps[3]))
    }
  } else if(fname == "hillfn"){
    if (top == ps[1]){
      ps[1] = cutoff
    } else {
      x_top = acy(y = top, modpars = list(tp=ps[1], ga=ps[2], p=ps[3], er=ps[4]), type = "hill")
      ps[1] = cutoff * ( 1 + (ps[2]/x_top)^ps[3])
    }
  } else if(fname == "gnls"){
    if (top == ps[1]) {
      ps[1] = cutoff
    } else {
      x_top = acy(y = top, modpars = list(tp=ps[1],ga=ps[2],p=ps[3],la=ps[4],q=ps[5],er=ps[6]), type = fname)
      ps[1] = cutoff * (( 1 + (ps[2]/x_top)^ps[3])*( 1 + (x_top/ps[4])^ps[5]))
    }
  } else if(fname == "poly1"){
    x_top = acy(y = top, modpars = list(a=ps[1],er=ps[2]),type = fname)
    ps[1] = cutoff/x_top
  } else if(fname == "poly2"){
    x_top = acy(y = top, modpars = list(a=ps[1], b=ps[2], er=ps[3]), type = fname, poly2.biphasic = poly2.biphasic)
    ps[1] = cutoff/(x_top/ps[2] + (x_top/ps[2])^2 )
  } else if(fname == "pow"){
    x_top = acy(y = top, modpars = list(a=ps[1], p=ps[2], er=ps[3]), type = fname)
    ps[1] = cutoff/(x_top^ps[2])
  }
  #get loglikelihood of top exactly at cutoff, use likelihood profile test
  # to calculate probability of being above cutoff
  loglik = tcplObj(p = ps, conc = conc, resp = resp, fname = fname, errfun = errfun)
  if(abs(top) >= abs(cutoff)) out = (1 + pchisq(2*(mll - loglik), 1))/2
  if(abs(top) < abs(cutoff)) out = (1 - pchisq(2*(mll - loglik), 1))/2

  return(out)

}
