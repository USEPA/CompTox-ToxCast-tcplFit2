#' Hit Logic (Discrete)
#'
#' Wrapper that computes discrete hitcalls for a provided concRespCore dataframe.
#'
#' @param indf Dataframe similar to concRespCore inpu Must contain "conc" and "resp"
#'   columns if xs and ys are not provided. Must contain "cutoff" and "bmad_factor"
#'   columns if newbmad is not NULL. Must contain "top" and "ac50" columns. "conc"
#'   and "resp" entries should be a single string with values separated by |.
#' @param newbmad (Deprecated) New number of bmads to use for the cutoff.
#' @param xs List of concentration vectors that can be provided for speed.
#' @param ys List of response vectors that can be provided for speed.
#' @param newcutoff Vector of new cutoff values to use. Length should be equal
#'   to rows in indf.
#'
#' @return Vector of hitcalls with length equal to number of rows in indf.
#' @export
#'
#' @examples
#' conc = rep(".03|.1|.3|1|3|10|30|100",2)
#' resp = rep("0|0|.1|.1|.5|.5|1|1",2)
#' indf = data.frame(top = c(1,1), ac50 = c(3,4), conc = conc, resp = resp,
#'   stringsAsFactors = FALSE)
#' hitlogic(indf, newcutoff = c(.8,1.2))
hitlogic = function(indf, newbmad = NULL, xs = NULL, ys = NULL, newcutoff = NULL){

  #extract cutoff from newbmad or newcutoff
  if(!is.null(newbmad)) cutoff = indf$cutoff/indf$bmad_factor*newbmad
  if(!is.null(newcutoff)) cutoff = newcutoff

  #reformat concs and resps, if necessary
  if(is.null(xs)){
    xs = strsplit(indf$conc, "\\|")
    xs = lapply(xs, as.numeric)
  }
  if(is.null(ys)){
    ys = strsplit(indf$resp, "\\|")
    ys = lapply(ys, as.numeric)
  }

  #run hitoginner for each row of indf
  hitcall = mapply(hitloginner, conc= xs, resp = ys, top = indf$top, cutoff = cutoff, ac50 = indf$ac50)

  return(hitcall)
}





