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
toplikelihood = function(fname, cutoff, conc, resp, ps, top, mll){

  #reparameterize so that top is exactly at cutoff
  if(fname == "exp2"){
    ps[1] = cutoff/( exp(max(conc)/ps[2]) - 1 )
  } else if(fname == "exp3"){
    ps[1] = cutoff/( exp((max(conc)/ps[2])^ps[3]) - 1 )
  } else if(fname == "exp4"){
    ps[1] = cutoff
  } else if(fname == "exp5"){
    ps[1] = cutoff
  } else if(fname == "hillfn"){
    ps[1] = cutoff
  } else if(fname == "gnls"){
    #approximating actual top with theoretical top for convenience.
    ps[1] = cutoff
  } else if(fname == "poly1"){
    ps[1] = cutoff/max(conc)
  } else if(fname == "poly2"){
    ps[1] = cutoff/(max(conc)/ps[2] + (max(conc)/ps[2])^2 )
  } else if(fname == "pow"){
    ps[1] = cutoff/(max(conc)^ps[2])
  }
  #get loglikelihood of top exactly at cutoff, use likelihood profile test
  # to calculate probability of being above cutoff
  loglik = tcplObj(p = ps, conc = conc, resp = resp, fname = fname)
  if(abs(top) >= cutoff) out = (1 + pchisq(2*(mll - loglik), 1))/2
  if(abs(top) < cutoff) out = (1 - pchisq(2*(mll - loglik), 1))/2

  return(out)

}
