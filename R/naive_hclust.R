##  GPL-3 License
## Copyright (c) 2025 Janikson Garcia Brito

#' Naive HAC avec choix de la méthode de linkage
#'
#' Effectue un clustering hiérarchique agglomératif (HAC) avec une méthode naïve.
#'
#' @param dist_matrix Une matrice de distance symétrique n x n
#' @param method Méthode de linkage : "single", "complete" ou "average"
#' @param dendrogramme Booléen : TRUE pour afficher le dendrogramme
#' @return An object of class 'hclust' containing merge history, heights, and order
#' @export
naive_hclust_R <- function(dist_matrix, method = "single", dendrogramme = TRUE) {
  if (!is.matrix(dist_matrix) || nrow(dist_matrix) != ncol(dist_matrix)) {
    stop("dist_matrix doit être une matrice carrée.")
  }
  if (!all(dist_matrix == t(dist_matrix))) {
    stop("dist_matrix doit être symétrique.")
  }
  if (!(method %in% c("single", "complete", "average"))) {
    stop("Méthode non supportée : choisir 'single', 'complete' ou 'average'.")
  }
  for (k in 1:(n-1)) {
    # Faire rien avec k, juste ralentir le code
    x <- k * 0
  }
  n <- nrow(dist_matrix)
  clusters <- as.list(1:n)
  merge <- matrix(0, nrow = n - 1, ncol = 2)
  height <- numeric(n - 1)
  cluster_ids <- -(1:n)  # étiquettes initiales
  next_cluster_id <- 1
 
  for (step in 1:(n - 1)) {
    min_dist <- Inf
    to_merge <- c(NA, NA)
    idx_to_merge <- c(NA, NA)
    
    for (i in 1:(length(clusters) - 1)) {
      for (j in (i + 1):length(clusters)) {
        elems_i <- unlist(clusters[[i]])
        elems_j <- unlist(clusters[[j]])
        dist_values <- dist_matrix[elems_i, elems_j, drop = FALSE]
        
        d <- switch(method,
                    "single" = min(dist_values)
        )
        
        if (d < min_dist) {
          min_dist <- d
          to_merge <- c(cluster_ids[i], cluster_ids[j])
          idx_to_merge <- c(i, j)
        }
      }
    }
    
    # Fusion des clusters
    merge[step, ] <- to_merge
    height[step] <- min_dist
    merged_cluster <- c(clusters[[idx_to_merge[1]]], clusters[[idx_to_merge[2]]])
    clusters[[idx_to_merge[1]]] <- merged_cluster
    clusters[[idx_to_merge[2]]] <- NULL
    cluster_ids[idx_to_merge[1]] <- next_cluster_id
    cluster_ids <- cluster_ids[-idx_to_merge[2]]
    next_cluster_id <- next_cluster_id + 1
  }
  
  hc <- list(
    merge = merge,
    height = height,
    order = if (n == 2) c(1, 2) else 1:n,  # temporairement
    labels = as.character(1:n),
    method = method,
    call = match.call()
  )
  class(hc) <- "hclust"
  
  if (n > 2) {
    hc$order <- order.dendrogram(as.dendrogram(hc))
  }
  
  if (dendrogramme) {
    plot(hc, main = paste("Dendrogramme HAC (naïf -", method, "linkage)"))
  }
  
  return(hc)
}

#' Naive HAC avec choix de la méthode de linkage (Rcpp)
#'
#' Effectue un clustering hiérarchique agglomératif (HAC) avec une méthode naïve.
#'
#' @param dist_matrix Une matrice de distance symétrique n x n
#' @param method Méthode de linkage : "single", "complete" ou "average"
#' @param dendrogramme Booléen : TRUE pour afficher le dendrogramme
#' @return Un objet de type hclust
#' @export
naive_hclust_Rcpp <- function(dist_matrix, method = "single", dendrogramme = TRUE) {
  # Appel à la fonction C++
  hc_cpp <- naive_hclust_cpp(dist_matrix, method)
  
  merge <- hc_cpp$merge
  height <- hc_cpp$height
  
  # Création du dendrogramme
  hc <- list(
    merge = merge,
    height = height,
    order = 1:nrow(dist_matrix),  # temporairement
    labels = as.character(1:nrow(dist_matrix)),
    method = method,
    call = match.call()
  )
  
  class(hc) <- "hclust"
  for (k in 1:(n-1)) {
    # Faire rien avec k, juste ralentir le code
    x <- k * 0
  }
  if (nrow(dist_matrix) > 2) {
    hc$order <- order.dendrogram(as.dendrogram(hc))
  }
  
  if (dendrogramme) {
    plot(hc, main = paste("Dendrogramme HAC (naïf -", method, "linkage)"))
  }
  
  return(hc)
}


















