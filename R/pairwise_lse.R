#' Pairwise logSumExp
#'
#' Implements logSumExp trick on pairs.
#'
#' @param la log(a).
#' @param lb log(b).
#'
#' @return log(a + b)
#' @export
#'
#' @examples
#' pairwise_lse(la = logq + epi_prior, lb = log(1-q) + epi_data)
pairwise_lse <- function(la, lb) {
  #log(a+b)=[la=log(a),lb=log(b)]=log(exp(la)+exp(lb))=la+log(1+exp(lb−la))
  lc <- pmax(la, lb)
  ld <- pmin(la, lb)

  return(lc + log(1 + exp(ld - lc)))
}
