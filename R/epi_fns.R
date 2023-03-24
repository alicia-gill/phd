#' Epidemic functions
#'
#' A collection of functions for the birth rate in used to generate epidemics
#'
#' @param t time
#' @param stop_time length of the epidemic
#' @param a start rate
#' @param b end rate
#'
#' @return function
#' @export

constant <- function(t, stop_time, a=0.3) {
  return(a)
}

change_pt <- function(t, stop_time, cp=0.5, a=0.5, b=0.1) {
  if (t < stop_time*cp) {
    return(a)
  } else {
    return(b)
  }
}

linear_dec <- function(t, stop_time, a=0.5, b=0.1) {
  return(a + t*(b-a)/stop_time)
}

linear_inc <- function(t, stop_time, a=0.1, b=0.5) {
  return(a + t*(b-a)/stop_time)
}
