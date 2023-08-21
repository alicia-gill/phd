#' Sample from Infected
#'
#' Samples and keeps leaves with probability pi.
#'
#' @param ptree object of class phylo.
#' @param day how many days in the past the most recent leaf is.
#' @param pi0 probability of sampling a leaf at the present day.
#' @param pi1 probability of sampling a leaf before the present day.
#'
#' @return object of class phylo; the tree keeping only the sampled leaves.
#' @export
#'
#' @examples
#' sample_het(ptree = full_tree, pi0 = 0, pi1 = 0.01)
sample_het <- function(ptree, day=0, pi0, pi1) {
  if (pi0 < 0 | pi0 > 1 | pi1 < 0 | pi1 > 1) {
    warning("pi must be between 0 and 1")
    break
  }

  #number of leaves
  n_leaves <- length(ptree$tip.label)

  distance <- distToRoot(ptree)
  max <- max(distance)

  #calculate how much to add to edge lengths to make leaves end on particular days
  change <- (max - distance) %% 1
  vector <- (change > 0.99999999 | change < 0.00000001)
  change[vector] <- 0

  #add on those differences
  edges <- (ptree$edge[,2] <= n_leaves)
  ptree$edge.length[edges] <- ptree$edge.length[edges] + change

  new_distance <- distToRoot(ptree)
  new_max <- max(new_distance)

  if (day > 0) {
    day1 <- 1:n_leaves
    day0 <- NULL
  } else {
    day1 <- (1:n_leaves)[round(new_max - new_distance, 0) > 0]
    day0 <- (1:n_leaves)[-day1]
  }
  #sample leaves
  u0 <- runif(length(day0))
  keep0 <- day0[u0 <= pi0]
  u1 <- runif(length(day1))
  keep1 <- day1[u1 <= pi1]
  keep <- c(keep0, keep1)

  if (!length(keep)) {
    new_tree <- NULL
    day <- NULL
  } else {
    if (!length(keep0)) {
      day <- min((new_max - new_distance)[keep1]) + day
    } else {
      day <- day
    }
    new_tree <- ape::keep.tip(ptree, keep)
    class(new_tree) <- "phylo"
  }

  output <- list("ptree"=new_tree, "day"=day)
  return(output)
}
