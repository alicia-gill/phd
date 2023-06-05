distToRoot <- function (phy) {
  e <- rep(0,nrow(phy$edge))
  e[phy$edge[,2]] <- 1:nrow(phy$edge)
  ntip <- length(phy$tip.label)
  d <- rep(0,ntip)
  for (i in 1:ntip) {
    w <- i
    while (1) {
      r <- e[w]
      if (r==0) {
        break
      }
      d[i] <- d[i] + phy$edge.length[r]
      w <- phy$edge[r,1]
    }
  }
  return(d)
}
