
################################################################
################################################################
#################### Data and Parameters #######################
################################################################
################################################################

rm(list = ls())
## Preliminaries
library(ivpack)
library(igraph)
library(RColorBrewer)
set.seed(1234) # set seed

#####################################
########### Loads Dataset ###########
#####################################

setwd("~/Dropbox/local_summer_camp/structural_group_conformism_clean") # set working directory
load("estimWave2withcontrols.RData")
source("fcts_gmm.R") # load functions


#####################################
### Make network graphs examples ####
#####################################
mancolors <- brewer.pal(n=3, name="Blues")
# classroom 33 (size=18, q=0.33)
# classroom 29 (size=29, q=0.43)
i <- 33
g <- graph_from_adjacency_matrix(G[[i]])
g$palette <- mancolors[c(1,3)]
V(g)$color <- (1-(X[[i]][,3]-1))+1
V(g)$label <- NA
plot(g)
print(mean(X[[i]][,3]==1))
print(nrow(X[[i]]))

i <- 29
g <- graph_from_adjacency_matrix(G[[i]])
g$palette <- mancolors[c(1,3)]
V(g)$color <- (1-(X[[i]][,3]-1)) + 1
V(g)$label <- NA
plot(g)
print(mean(X[[i]][,3]==1))
print(nrow(X[[i]]))


for (i in 1:36){
  print("**")
  print(i)
  g <- graph_from_adjacency_matrix(G[[i]])
  V(g)$color <- (1-(X[[i]][,3]-1))*2
  V(g)$label <- NA
  plot(g)
  print(mean(X[[i]][,3]==1))
  print(nrow(X[[i]]))
}

