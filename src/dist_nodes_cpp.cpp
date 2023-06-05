#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]
//' Distances from Leaves to Root of Phylogenetic Tree
//'
//' Computes the distance from the leaves to the root of a phylogenetic tree.
//'
//' @param ptree_nodes matrix of nodes and edge lengths.
//' @param root_node label of root node.
//'
//' @return distances from leaves to root
//' @export
//'
//' @examples
//' dist_nodes_cpp <- function(ptree_nodes, root_node)
//'
//' @useDynLib phd, .registration = TRUE
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericVector dist_nodes_cpp(IntegerMatrix ptree_nodes, NumericVector ptree_lengths, int root_node) {
  int n_leaves = root_node - 1;
  NumericVector dist (n_leaves);

  IntegerVector v = ptree_nodes( _, 1);

  IntegerVector j(1);

  for (int i = 0; i < n_leaves; i++) {
    //int j = i;
    j[0] = i;
    while (j[0] != root_node) {
      IntegerVector k = match(j, v);
      int l = k[0] - 1;
      dist[i] += ptree_lengths[l];
      j[0] = ptree_nodes(l, 0);
    }
  }

  return dist;
}
