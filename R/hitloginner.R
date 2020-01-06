#' Hit Logic Inner (Discrete)
#'
#' Contains hit logic, called directly during CR fitting or later through "hitlogic".
#'
#' The purpose of this function is to keep the actual hit rules in one
#' location so it can be called during CR fitting, and then again after the fact
#' for a variety of cutoffs. Curves fit with constant winning should have
#' top = NA, generating a miss.
#'
#' @param conc Vector of concentrations (No longer necessary).
#' @param resp Vector of responses.
#' @param top Model top.
#' @param cutoff Desired cutoff.
#' @param ac50 Model AC50 (No longer necessary).
#'
#' @return Outputs 1 for hit, 0 for miss.
#' @export
#'
#' @examples
#' hitloginner(resp = 1:8, top = 7, cutoff = 5) #hit
#' hitloginner(resp = 1:8, top = 7, cutoff = 7.5) #miss: top too low
#' hitloginner(resp = 1:8, top = 9, cutoff = 8.5) #miss: no response> cutoff
#' hitloginner(resp = 1:8, top = NA, cutoff = 5) #miss: no top (constant)
hitloginner = function(conc = NULL, resp, top, cutoff, ac50 = NULL){


  n_gt_cutoff = sum(abs(resp)>cutoff)

  #hitlogic - hit must have: at least one point above abs cutoff,
  # a defined top (implying there is a winning non-constant model),
  #and an abs. top greater than the cutoff
  hitcall = 0
  if(n_gt_cutoff>0 && !is.na(top) && abs(top)>cutoff) hitcall <- 1

  return(hitcall)
}
