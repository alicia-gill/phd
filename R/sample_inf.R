#' Sample from Infected
#'
#' Samples and keeps those infected at the present day (i.e. the leaves) with probability pi.
#'
#' @param ptree object of class phylo.
#' @param pi probability of sampling a leaf.
#'
#' @return object of class phylo; the tree keeping only the sampled leaves.
#' @export
#'
#' @examples
#' sample_inf(ptree = full_tree, pi = 0.1)
sample_inf <- function(ptree, pi) {
  if (pi < 0 | pi > 1) {
    warning("pi must be between 0 and 1")
    break
  }

  #number of leaves
  n_leaves <- length(ptree$tip.label)

  #distance from root node (always coded as n_leaves+1) and leaves (always codes from 1 to n_leaves)
  distance <- ape::dist.nodes(ptree)[1:n_leaves, n_leaves+1]

  #only keep those with maximum distance
  max <- round(max(distance), digits=10)
  keep <- which(distance == max)

  n <- length(keep)
  u <- runif(n)
  keep <- keep[u <= pi]

  new_tree <- ape::keep.tip(ptree, keep)
  return(new_tree)
}
