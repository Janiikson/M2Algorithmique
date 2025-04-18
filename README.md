[![Build Status](https://travis-ci.com/Janiikson/M2AlgoHAC.svg?branch=main)](https://travis-ci.com/Janiikson/M2AlgoHAC)

# M2AlgoHAC Vignette

### Janikson Garcia Brito  
#### Université d'Évry - Paris-Saclay  
📅 Avril 3, 2025

---

## Sommaire

- [Quick Start](#quick-start)
- [The 4 Algorithms at Fixed Data Length](#the-4-algorithms-at-fixed-data-length)
- [Microbenchmark](#microbenchmark)
- [Time Complexity](#time-complexity)

---

## Quick Start

Le package `M2algorithmique` est un **exemple de package R** développé dans le cadre du cours d’algorithmique du Master 2 Data Science à l’Université d'Évry - Paris-Saclay. Il inclut plusieurs stratégies algorithmiques (notamment `HAC_naive` et `HAC_improved`) implémentées en R et en Rcpp.

---

# Hierarchical Agglomerative Clustering (HAC)

## Introduction

**Hierarchical Agglomerative Clustering (HAC)** est une technique de regroupement qui construit une hiérarchie de clusters en fusionnant itérativement les paires les plus proches selon une métrique de similarité (par exemple : la distance euclidienne). Le résultat est un **dendrogramme** illustrant cette hiérarchie.

### Algorithmic Challenge

Le principal défi algorithmique réside dans la complexité computationnelle. Une implémentation naïve construit une matrice de distance \( n \times n \) (**O(n²)**), mise à jour à chaque itération (**n - 1** fusions), ce qui donne une complexité totale de **O(n³)**.

---

## Problem Statement

Ce coût algorithmique devient vite prohibitif avec des jeux de données volumineux. Il est donc essentiel d’étudier des implémentations optimisées pour rendre cette méthode scalable.

---

## Objectives

1. Implémenter l’algorithme **HAC naïf** en **R** et **Rcpp**.
2. Comparer les performances entre les versions naïve et améliorée.
3. Étudier l’impact de la mise à jour des matrices de distance sur la complexité.
4. Explorer l’utilisation de structures de données efficaces (par exemple : **k-d trees**).

---

## Next Steps

- Tester et profiler les performances des différentes versions de l’algorithme.
- Identifier les goulots d’étranglement.
- Optimiser la complexité via des structures plus efficaces et une meilleure stratégie de fusion.

---

## Installation du package

Installe d'abord `devtools`, puis le package via GitHub :

```r
# install.packages("devtools")
devtools::install_github("Janiikson/M2Algorithmique")
library(FastHierarchicalClust)
```
 
 
