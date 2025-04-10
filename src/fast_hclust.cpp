#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List fast_hclust_cpp(NumericMatrix dist_matrix, String method) {
  int n = dist_matrix.nrow();
  if (dist_matrix.ncol() != n)
    stop("La matrice de distance doit être carrée.");
  
  // Initialisation
  std::vector<int> cluster_ids(n);
  std::vector<bool> active(n, true);
  for (int i = 0; i < n; ++i)
    cluster_ids[i] = -(i + 1);
  
  IntegerMatrix merge(n - 1, 2);
  NumericVector height(n - 1);
  int next_cluster_id = 1;
  
  // Liste des indices actifs (évite boucle complète)
  std::vector<int> active_indices(n);
  for (int i = 0; i < n; ++i) active_indices[i] = i;
  
  for (int step = 0; step < n - 1; ++step) {
    double min_dist = R_PosInf;
    int i_min = -1, j_min = -1;
    
    // Rechercher la paire la plus proche parmi les actifs
    for (size_t a = 0; a < active_indices.size(); ++a) {
      int i = active_indices[a];
      for (size_t b = a + 1; b < active_indices.size(); ++b) {
        int j = active_indices[b];
        double d = dist_matrix(i, j);
        if (d < min_dist) {
          min_dist = d;
          i_min = i;
          j_min = j;
        }
      }
    }
    
    // Fusion des clusters
    merge(step, 0) = cluster_ids[i_min];
    merge(step, 1) = cluster_ids[j_min];
    height[step] = min_dist;
    cluster_ids[i_min] = next_cluster_id++;
    active[j_min] = false;
    
    // Mise à jour des distances
    for (int k : active_indices) {
      if (k == i_min || k == j_min || !active[k]) continue;
      
      double dist_i = dist_matrix(i_min, k);
      double dist_j = dist_matrix(j_min, k);
      double new_dist;
      
      if (method == "single") {
        new_dist = std::min(dist_i, dist_j);
      } else if (method == "complete") {
        new_dist = std::max(dist_i, dist_j);
      } else if (method == "average") {
        new_dist = (dist_i + dist_j) / 2.0;
      } else {
        stop("Méthode non supportée : single, complete, average");
      }
      
      dist_matrix(i_min, k) = new_dist;
      dist_matrix(k, i_min) = new_dist;
    }
    
    // Mise à jour de la liste des indices actifs
    std::vector<int> new_active_indices;
    for (int idx : active_indices) {
      if (active[idx]) new_active_indices.push_back(idx);
    }
    active_indices = new_active_indices;
  }
  
  return List::create(
    Named("merge") = merge,
    Named("height") = height
  );
}
