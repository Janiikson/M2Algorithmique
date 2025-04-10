#' Clustering hiérarchique optimisé (appel à Rcpp)
#'
#' @param dist_matrix Matrice de distance symétrique
#' @param method "single", "complete", "average"
#' @param dendrogramme TRUE pour afficher
#' @return Un objet hclust
#' @export
fast_hclust_R <- function(dist_matrix, method = "single", dendrogramme = TRUE) {
  if (!is.matrix(dist_matrix) || nrow(dist_matrix) != ncol(dist_matrix)) {
    stop("dist_matrix doit être une matrice carrée.")
  }
  if (!all(dist_matrix == t(dist_matrix))) {
    stop("dist_matrix doit être symétrique.")
  }
  
  n <- nrow(dist_matrix)
  clusters <- lapply(1:n, function(i) i)
  cluster_ids <- -(1:n)
  merge <- matrix(0, n - 1, 2)
  height <- numeric(n - 1)
  active <- rep(TRUE, n)
  next_id <- 1
  
  for (step in 1:(n - 1)) {
    best_dist <- Inf
    idx_i <- idx_j <- NA
    
    for (i in 1:(n - 1)) {
      if (!active[i]) next
      for (j in (i + 1):n) {
        if (!active[j]) next
        
        elems_i <- clusters[[i]]
        elems_j <- clusters[[j]]
        
        d <- switch(method,
                    "single" = min(dist_matrix[elems_i, elems_j, drop = FALSE]),
                    "complete" = max(dist_matrix[elems_i, elems_j, drop = FALSE]),
                    "average" = mean(dist_matrix[elems_i, elems_j, drop = FALSE]),
                    stop("Méthode non supportée.")
        )
        
        if (d < best_dist) {
          best_dist <- d
          idx_i <- i
          idx_j <- j
        }
      }
    }
    
    merge[step, ] <- c(cluster_ids[idx_i], cluster_ids[idx_j])
    height[step] <- best_dist
    clusters[[idx_i]] <- c(clusters[[idx_i]], clusters[[idx_j]])
    active[idx_j] <- FALSE
    cluster_ids[idx_i] <- next_id
    next_id <- next_id + 1
  }
  
  hc <- list(
    merge = merge,
    height = height,
    order = 1:n,
    labels = as.character(1:n),
    method = method,
    call = match.call()
  )
  class(hc) <- "hclust"
  
  if (n > 2) {
    hc$order <- order.dendrogram(as.dendrogram(hc))
  }
  if (dendrogramme) {
    plot(hc, main = paste("Fast HAC -", method, "(R)"))
  }
  
  return(hc)
}


#' Clustering hiérarchique optimisé (appel à Rcpp)
#'
#' @param dist_matrix Matrice de distance symétrique
#' @param method "single", "complete", "average"
#' @param dendrogramme TRUE pour afficher
#' @return Un objet hclust
#' @export
fast_hclust_Rcpp <- function(dist_matrix, method = "single", dendrogramme = TRUE) {
  res_cpp <- fast_hclust_cpp(dist_matrix, method)
  merge <- res_cpp$merge
  height <- res_cpp$height
  n <- nrow(dist_matrix)
  
  hc <- list(
    merge = merge,
    height = height,
    order = 1:n,
    labels = as.character(1:n),
    method = method,
    call = match.call()
  )
  class(hc) <- "hclust"
  if (n > 2) {
    hc$order <- order.dendrogram(as.dendrogram(hc))
  }
  if (dendrogramme) {
    plot(hc, main = paste("Dendrogramme HAC optimisé -", method, "linkage"))
  }
  return(hc)
}
