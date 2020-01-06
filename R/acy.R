#' Activity Concentration y
#'
#' Returns concentration at which model equals y.
#'
#' Mathematically inverts model functions of the given type, except for gnls,
#' which is numerically inverted. gnls returns NA when y > tp. Other options
#' return the actual top (as opposed to theoretical tp) and top location for
#' gnls model. gnls model defaults to giving concentration on gain side. Only
#' one of getloss, returntop, and returntoploc should be TRUE at a time.  If
#' top location solution fails for gnls, top is set to tp. Returns NA if gnls
#' numerical solver fails.
#'
#' @param y Activity value at which the concentration is desired. y
#'   should be less than the model's top, if there is one, and greater
#'   than zero.
#' @param modpars List of named model parameters. Model parameters can
#'   include: "a", "b", "ga", "la", "p", "q", "tp". ga and la should NOT
#'   be in log units.
#' @param type Model type; must be one of: "exp1", "exp2", "exp3", "exp4",
#'   "gnls", "hill",  "poly1", "poly2", "pow".
#' @param returntop When TRUE, returns actual top value for gnls. Has no
#'   effect for other models.
#' @param returntoploc When TRUE, returns concentration of top for gnls.
#'   Has no effect for other models. If top location can't be found,
#'   NA is returned.
#' @param getloss When TRUE, returns value on loss side of curve for gnls.
#'   Has no effect for other models.
#' @param verbose When TRUE, shows warnings.
#'
#' @return Ouputs concentration at activity y, or gnls top or top concentration,
#'   when applicable.
#'
#' @importFrom stats uniroot
#'
#'
#' @export
#'
#' @examples
#' acy(1, list(ga = 10, tp = 2, p = 3), type = "hill")
#' acy(1, list(ga = .1, tp = 2, p = 3, q = 3,la = 10), type = "gnls")
#' acy(1, list(ga = .1, tp = 2, p = 3, q = 3,la = 10), type = "gnls", getloss = TRUE)
#' acy(1, list(ga = .1, tp = 2, p = 3, q = 3,la = 10), type = "gnls", returntop = TRUE)
#' acy(1, list(ga = .1, tp = 2, p = 3, q = 3,la = 10), type = "gnls", returntoploc = TRUE)
#'
acy <- function(y, modpars, type = "hill", returntop = F, returntoploc = F, getloss =F, verbose = F) {
  #variable binding to pass cmd checks
  a <- b <- tp <- ga <- p <- q <- la <- NULL
  #Put model parameters in environment: a,b,tp,ga,p,q,la,er
  list2env(modpars, envir = environment())

  #warnings
  if(!returntop){
    if(!is.null(modpars$tp) && abs(y) >= abs(tp)) {
      if(verbose) warning("y is greater than top in function acy, returning NA")
      return(NA)
    }
    if(!is.null(modpars$tp) && sign(y) != sign(tp)) {
      if(verbose) warning("y is wrong sign in function acy, returning NA")
      return(NA)
    }
  }

  #Invert most models analytically, gnls numerically.
  if(type == "poly1"){
    return(y/a)
  } else if(type == "poly2"){
    return( b*(-1 + sqrt(1 + 4*y/a))/2)
  } else if(type == "pow"){
    return((y/a)^(1/p))
  } else if(type == "exp2"){
    return(b*log(y/a + 1))
  } else if(type == "exp3"){
    return(b*(log(y/a + 1))^(1/p) )
  } else if(type == "exp4"){
    return(-ga*log2(1-y/tp))
  } else if(type == "exp5"){
    return(ga*(-log2(1-y/tp))^(1/p))
  } else if(type =="hill"){
    return(ga/(tp/y-1)^(1/p))
  } else if(type == "gnls"){
    #gnls top can be much lower than tp, so first find top location by setting derivative to zero
    #gnls outputs fraction of actual top, not theoretical top tp
    toploc = try(uniroot(gnlsderivobj, c(ga,la), tp = tp, ga = ga, p = p, la = la, q = q, tol = 1e-8)$root)

    #If toploc fails, set topval to tp, set toploc to NA
    if(class(toploc) == "try-error"){
      if(verbose) warning("toploc could not be found numerically")
      topval = tp
      toploc = NA_real_
    } else {
      #get actual top, as opposed to theoretical tp.
      topval = gnls(c(tp, ga, p, la, q), toploc)
    }
    if(returntoploc) return(toploc)
    if(returntop) return(topval)

    #If y >= top, don't try to solve.
    if(abs(y) > abs(topval)) {
      if(verbose) warning("y is greater than gnls top in function acy, returning NA")
      return(NA)
    }
    if(y == topval) return(toploc)

    #solve for acy
    if(getloss) {
      output = try(uniroot(acgnlsobj, c(toploc, 1e5), y = y, tp = tp, ga = ga, p = p, la = la, q = q, tol = 1e-8)$root)
    } else {
      output = try(uniroot(acgnlsobj, c(1e-8, toploc), y = y, tp = tp, ga = ga, p = p, la = la, q = q, tol = 1e-8)$root)
    }
    if(class(output) == "try-error") return(NA_real_) else return(output)
  }

  return(NA)
}

#' GNLS Derivative Objective Function
#'
#' Derivative of the gnls function set to zero for top location solver.
#'
#' @param x Concentration.
#' @param tp Top.
#' @param ga Gain AC50.
#' @param p Gain power.
#' @param la Loss AC50.
#' @param q Loss power.
#'
#' @return Value of gnls derivative at x.
#' @export
gnlsderivobj = function(x,tp,ga,p,la,q){
  a = ga^p
  b = la^(-q)
  return(b*q*x^(q+p) + a*b*(q-p)*x^q -a*p)

}


#' AC GNLS Objective Function
#'
#' GNLS objective function set to y for gnls solver.
#'
#' @param x Concentration.
#' @param y Desired activity level.
#' @param tp Top.
#' @param ga Gain AC50.
#' @param p Gain power.
#' @param la Loss AC50.
#' @param q Loss power.
#'
#' @return Difference between GNLS model repsone at x and y.
#' @export
acgnlsobj = function(x,y,tp,ga,p,la,q){
  #y is desired y value
  return(gnls(c(tp, ga, p, la, q), x)-y)
}
