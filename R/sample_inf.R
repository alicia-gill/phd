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
