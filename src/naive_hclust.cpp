#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List naive_hclust_cpp(NumericMatrix dist_matrix, String method) {
  int n = dist_matrix.nrow();
  
  // Vérifications
  if (dist_matrix.ncol() != n) {
    stop("dist_matrix must be a square matrix.");
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (dist_matrix(i, j) != dist_matrix(j, i)) {
        stop("dist_matrix must be symmetric.");
      }
    }
  }
  
  // Initialisation des clusters
  std::vector<std::vector<int>> clusters(n);
  for (int i = 0; i < n; i++) {
    clusters[i].push_back(i + 1); // 1-based indexing for R
  }
  
  // Matrices pour stocker les résultats
  IntegerMatrix merge(n - 1, 2);
  NumericVector height(n - 1);
  std::vector<int> cluster_ids(n);
  for (int i = 0; i < n; i++) {
    cluster_ids[i] = -(i + 1); // Étiquettes initiales négatives
  }
  int next_cluster_id = 1;
  
  // Boucle principale
  for (int step = 0; step < n - 1; step++) {
    double min_dist = R_PosInf;
    int to_merge_i = -1, to_merge_j = -1;
    int idx_i = -1, idx_j = -1;
    
    // Recherche de la paire de clusters la plus proche
    for (int i = 0; i < clusters.size() - 1; i++) {
      for (int j = i + 1; j < clusters.size(); j++) {
        double d = 0.0;
        if (method == "single") {
          d = R_PosInf;
          for (int ei : clusters[i]) {
            for (int ej : clusters[j]) {
              double dist = dist_matrix(ei - 1, ej - 1);
              if (dist < d) d = dist;
            }
          }
        } else if (method == "complete") {
          d = -R_PosInf;
          for (int ei : clusters[i]) {
            for (int ej : clusters[j]) {
              double dist = dist_matrix(ei - 1, ej - 1);
              if (dist > d) d = dist;
            }
          }
        } else if (method == "average") {
          double sum = 0.0;
          int count = 0;
          for (int ei : clusters[i]) {
            for (int ej : clusters[j]) {
              sum += dist_matrix(ei - 1, ej - 1);
              count++;
            }
          }
          d = sum / count;
        } else {
          stop("Unsupported method: choose 'single', 'complete', or 'average'.");
        }
        
        if (d < min_dist) {
          min_dist = d;
          to_merge_i = cluster_ids[i];
          to_merge_j = cluster_ids[j];
          idx_i = i;
          idx_j = j;
        }
      }
    }
    
    // Fusion des clusters
    merge(step, 0) = to_merge_i;
    merge(step, 1) = to_merge_j;
    height[step] = min_dist;
    
    // Mettre à jour les clusters
    clusters[idx_i].insert(clusters[idx_i].end(), clusters[idx_j].begin(), clusters[idx_j].end());
    clusters.erase(clusters.begin() + idx_j);
    cluster_ids[idx_i] = next_cluster_id;
    cluster_ids.erase(cluster_ids.begin() + idx_j);
    next_cluster_id++;
  }
  
  // Retourner les résultats
  return List::create(
    Named("merge") = merge,
    Named("height") = height
  );
}




