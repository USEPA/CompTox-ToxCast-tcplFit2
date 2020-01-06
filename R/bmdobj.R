#' BMD Objective Function
#'
#' Utility function for bmdbounds
#'
#' @param bmd Benchmark dose.
#' @param fname Function name: "exp2", "exp3", "exp4", "exp5", "hillfn", "gnls",
#'   "poly1", "poly2", or "pow".
#' @param bmr Benchmark response.
#' @param conc Vector of concentrations NOT in log units.
#' @param resp Vector of corresponding responses.
#' @param ps Named list of paramters.
#' @param mll Maximum log-likelihood of winning model.
#' @param onesp One-sided p-value.
#' @param partype Number for parameter type. Type 1 is y-scaling: a or tp.
#'   Type 2 is x-scaling: b or ga, when available, a otherwise. Type 3 is
#'   power scaling: p when available, then b or ga, then a if no others.
#'   Since bmd is linked to the x-scale, type 2 should always be used. Other
#'   types can also be vulnerable to underflow/overflow.
#'
#' @importFrom stats qchisq
#'
#' @return Objective function value to find the zero of.
#' @export
bmdobj= function(bmd, fname, bmr, conc, resp, ps, mll, onesp, partype = 2){

  #implements the BMD substitutions in Appendix A of the Technical Report.
  #Changes one of the existing parameters to an explicit bmd parameter through
  # the magic of algebra.
  if(fname == "exp2"){
    if(partype == 1) ps["a"] = bmr/( exp(bmd/ps["b"]) - 1 )
    if(partype == 2) ps["b"] = bmd/( log(bmr/ps["a"] + 1) )
    if(partype == 3) ps["b"] = bmd/( log(bmr/ps["a"] + 1) )
  } else if(fname == "exp3"){
    if(partype == 1) ps["a"] = bmr/( exp((bmd/ps["b"])^ps["p"]) - 1 )
    if(partype == 2) ps["b"] = bmd/( log(bmr/ps["a"] + 1) )^(1/ps["p"])
    if(partype == 3) ps["p"] = log(log(bmr/ps["a"] + 1))/log(bmd/ps["b"])
  } else if(fname == "exp4"){
    if(partype == 1) ps["tp"] = bmr/( 1-2^(-bmd/ps["ga"]))
    if(partype == 2) ps["ga"] = bmd/( -log2(1-bmr/ps["tp"]) )
    if(partype == 3) ps["ga"] = bmd/( -log2(1-bmr/ps["tp"]) )
  } else if(fname == "exp5"){
    if(partype == 1) ps["tp"] = bmr/( 1-2^(-(bmd/ps["ga"])^ps["p"] ))
    if(partype == 2) ps["ga"] = bmd/(( -log2(1-bmr/ps["tp"]) )^(1/ps["p"]))
    if(partype == 3) ps["p"] = log( -log2( 1 - bmr/ps["tp"]) )/log(bmd/ps["ga"])
  } else if(fname == "hillfn"){
    if(partype == 1) ps["tp"] = bmr*( 1 + (ps["ga"]/bmd)^ps["p"])
    if(partype == 2) ps["ga"] = bmd* (ps["tp"]/bmr - 1)^(1/ps["p"])
    if(partype == 3) ps["p"] = log(ps["tp"]/bmr - 1)/log(ps["ga"]/bmd)
  } else if(fname == "gnls"){
    #gnls bounds don't include loss part parameterization
    if(partype == 1) ps["tp"] = bmr * ( (1 + (ps["ga"]/bmd)^ps["p"]) * ( 1 + (bmd/ps["la"])^ps["q"]) )
    if(partype == 2) ps["ga"] = bmd * ( ps["tp"]/(bmr*(1 + (bmd/ps["la"])^ps["q"])) - 1)^(1/ps["p"])
    if(partype == 3) ps["p"] = log(ps["tp"]/(bmr*(1 + (bmd/ps["la"])^ps["q"])) - 1) / log(ps["ga"]/bmd)
  } else if(fname == "poly1"){
    if(partype == 1) ps["a"] = bmr/bmd
    if(partype == 2) ps["a"] = bmr/bmd
    if(partype == 3) ps["a"] = bmr/bmd
  } else if(fname == "poly2"){
    if(partype == 1) ps["a"] = bmr/(bmd/ps["b"] + (bmd/ps["b"])^2 )
    if(partype == 2) ps["b"] = 2*bmd/(sqrt(1 + 4*bmr/ps["a"]) - 1)
    if(partype == 3) ps["b"] = 2*bmd/(sqrt(1 + 4*bmr/ps["a"]) - 1)
  } else if(fname == "pow"){
    if(partype == 1) ps["a"] = bmr/(bmd^ps["p"])
    if(partype == 2) ps["a"] = bmr/(bmd^ps["p"])
    if(partype == 3) ps["p"] = log(bmr/ps["a"])/log(bmd)
  }

  loglik = tcplObj(p = ps, conc = conc, resp = resp, fname = fname)
  #for bmd bounds, we want the difference between the max log-likelihood and the
  #bounds log-likelihood to be equal to chi-squared at 1-2*onesp (typically .9)
  #with one degree of freedom divided by two.
  return(mll - loglik - qchisq(1-2*onesp,1)/2)

}

