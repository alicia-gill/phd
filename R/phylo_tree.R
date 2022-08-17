#' Phylogenetic Tree
#'
#' Represents the epidemic as a phylogenetic tree.
#'
#' @param epidemic output from epidemic().
#' @param stop_time number of days the epidemic simulation was run for.
#'
#' @return object of class phylo.
#' @export
#'
#' @examples
#' phylo_tree(epidemic = epi, stop_time = 50)
phylo_tree <- function(epidemic, stop_time) {
  n_total <- nrow(epidemic)

  #sets recovery time as stop time (t_stop) for those still infectious at the end
  n_inf <- epidemic$index[is.na(epidemic$rem_time) == T]
  epidemic$rem_time[n_inf] <- stop_time

  #initial empty dataset
  tree_data <- list("start_node","end_node","start_time","end_time","edge_length")

  #edge matrix
  for (i in 1:n_total) {
    #c is a set of internal nodes which map to each other depending on who their parent
    c <- c(n_total + i, n_total + epidemic$index[epidemic$parent == i])
    l <- length(c)
    if (l > 1) {
      for (j in 1:(l-1)) {
        tree_data$start_node <- c(tree_data$start_node, c[j])
        tree_data$end_node <- c(tree_data$end_node, c[j+1])
      }
    }
    tree_data$start_node <- c(tree_data$start_node, c[l])
    tree_data$end_node <- c(tree_data$end_node, i)
  }

  E <- length(tree_data$start_node)
  if (E != 2 * n_total - 1) {
    print("Wrong number of edges!")
    break
  }

  for (i in 1:E) {
    #start time is infection time of whoever generated that node
    #eg person 1 generates node n_total+1, person 2 generates node n_total+2, etc
    tree_data$start_time[i] <- epidemic$inf_time[tree_data$start_node[i] - n_total]
    #if the end node is a leaf, then end time is recovery time of that leaf
    if (tree_data$end_node[i] <= n_total) {
      tree_data$end_time[i] <- epidemic$rem_time[tree_data$end_node[i]]
    } else { #otherwise, it's an infection time
      tree_data$end_time[i] <- epidemic$inf_time[tree_data$end_node[i] - n_total]
    }
    #length is end_time-start_time
    tree_data$edge_length[i] <- tree_data$end_time[i] - tree_data$start_time[i]
  }

  #turn into a tree
  tree <- list()
  class(tree) <- "phylo"

  tree$edge <- matrix(c(tree_data$start_node,tree_data$end_node),ncol=2,byrow=F)
  tree$Nnode <- n_total
  tree$tip.label <- 1:n_total
  tree$edge.length <- tree_data$edge_length
  tree$node.label <- (n_total+1):(E+1)

  return(tree)
}
