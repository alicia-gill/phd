#' Sample from Infected
#'
#' Samples and keeps leaves with probability pi.
#'
#' @param ptree object of class phylo.
#' @param pi0 probability of sampling a leaf at the present day.
#' @param pi1 probability of sampling a leaf before the present day.
#'
#' @return object of class phylo; the tree keeping only the sampled leaves.
#' @export
#'
#' @examples
#' sample_het(ptree = full_tree, pi0 = 0, pi1 = 0.01)
sample_het <- function(ptree, pi0, pi1) {
  if (pi0 < 0 | pi0 > 1 | pi1 < 0 | pi1 > 1) {
    warning("pi must be between 0 and 1")
    break
  }

  #number of leaves
  n_leaves <- length(ptree$tip.label)

  #distance from root node (always coded as n_leaves+1) and leaves (always codes from 1 to n_leaves)
  distance <- ape::dist.nodes(ptree)[1:n_leaves, n_leaves+1]

  #calculate how much to add to edge lengths to make leaves end on particular days
  max <- max(distance)
  change <- (max - distance) %% 1

  #add on those differences
  edges <- which(ptree$edge[,2] <= n_leaves)
  ptree$edge.length[edges] <- ptree$edge.length[edges] + change

  distance <- ape::dist.nodes(ptree)[1:n_leaves, n_leaves+1]
  max <- max(distance)

  #option 1 - split into present/not present and sample separately
  #label leaves
  day1 <- (1:n_leaves)[which(distance < max)]
  day0 <- (1:n_leaves)[-day1]
  #sample leaves
  u0 <- runif(length(day0))
  keep0 <- day0[u0 <= pi0]
  u1 <- runif(length(day1))
  keep1 <- day1[u1 <= pi1]
  keep <- c(keep0, keep1)

  new_tree <- ape::keep.tip(ptree, keep)
  return(new_tree)
}