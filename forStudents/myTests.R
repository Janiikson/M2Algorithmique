library(M2AlgoHAC)
library(dplyr)
library(microbenchmark)
library(ggplot2)
library(usethis)
library(devtools)
devtools::document()
devtools::load_all()


# Générer des données de test
set.seed(123)
n <- 10
data <- matrix(runif(n * 2), nrow = n, ncol = 2)
dist_matrix <- as.matrix(dist(data, method = "euclidean"))

library(FNN)

# Exécuter les quatre méthodes
hc_naive_R <- naive_hclust_R(dist_matrix, method = "single", dendrogramme = FALSE)
hc_naive_Rcpp <- naive_hclust_Rcpp(dist_matrix, method = "single", dendrogramme = FALSE)
# Fast HAC (R)
hc_fast_R <- fast_hclust_R(dist_matrix, method = "single")
# Fast HAC (Rcpp)
hc_fast_Rcpp <- fast_hclust_Rcpp(dist_matrix, method = "single")

# Comparer les matrices merge
cat("Comparaison des matrices merge :\n")
print(all.equal(hc_naive_R$merge, hc_naive_Rcpp$merge))
print(all.equal(hc_naive_R$merge, hc_fast_R$merge))
print(all.equal(hc_naive_R$merge, hc_fast_Rcpp$merge))

# Comparer les hauteurs
cat("\nComparaison des hauteurs :\n")
print(all.equal(hc_naive_R$height, hc_naive_Rcpp$height))
print(all.equal(hc_naive_R$height, hc_fast_R$height))
print(all.equal(hc_naive_R$height, hc_fast_Rcpp$height))

# Visualiser les dendrogrammes pour une vérification visuelle
par(mfrow = c(2, 2))
plot(hc_naive_R, main = "Naive R")
plot(hc_naive_Rcpp, main = "Naive Rcpp")
plot(hc_fast_R, main = "Fast R")
plot(hc_fast_Rcpp, main = "Fast Rcpp")



# Générer des données de test
set.seed(123)
n <- 10
data <- matrix(runif(n * 2), nrow = n, ncol = 2)

# Tester fast_hclust_Rcpp
hc_fast_Rcpp <- fast_hclust_Rcpp(dist_matrix, method = "complete", dendrogram = TRUE)
print(hc_fast_Rcpp)









