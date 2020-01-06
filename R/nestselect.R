#' Nest Select
#'
#' Chooses between nested models.
#'
#' @param aics Named vector of model aics (can include extra models).
#' @param mod1 Name of model 1, the model with fewer degrees of freedom.
#' @param mod2 Name of model 2, the model with more degrees of freedom.
#' @param dfdiff Absolute difference in number of degrees of freedom
#'   (i.e. the difference in parameters).
#' @param pval P-value for nested model test.
#' @importFrom stats pchisq
#'
#' @return Named aic vector with losing model removed.
#' @export
#'
#' @examples
#' aics = c(-5,-6,-3)
#' names(aics) = c("poly1", "poly2", "hill")
#' nestselect(aics, "poly1", "poly2", 1)
#'
#' aics = c(-5,-7,-3)
#' names(aics) = c("poly1", "poly2", "hill")
#' nestselect(aics, "poly1", "poly2", 1)
nestselect = function(aics, mod1, mod2, dfdiff, pval = .05){
  if(is.na(aics[mod1])){
    loser = mod1 #if model 1 AIC is NA it is the loser
  } else if(isTRUE(aics[mod2] <= aics[mod1])){
    #if both aics exist and model 2 aic is lower and model 2 passes the ratio test,
    #model 2 wins and model 1 loses; otherwise, model 2 loses.
    chisq = aics[mod1] - aics[mod2] + 2*dfdiff #2*loglikelihood(poly2) - 2*loglikelihood(poly1)
    ptest =  1-pchisq(chisq, dfdiff)
    if(ptest < pval) loser = mod1 else loser = mod2
  } else loser = mod2

  return(aics[names(aics) != loser])
}
