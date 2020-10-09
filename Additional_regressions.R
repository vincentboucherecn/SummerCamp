
################################################################
################################################################
#################### Data and Parameters #######################
################################################################
################################################################

rm(list = ls())
## Preliminaries
library(ivpack)
library(igraph)
set.seed(1234) # set seed

#####################################
########### Loads Dataset ###########
#####################################

setwd("~/Dropbox/local_summer_camp/structural_group_conformism") # set working directory
load("estimWave2withcontrols.RData")
source("fcts_gmm.R") # load functions

#####################################
########### Build dataset ###########
#####################################

dtapairs <- buildpairs() # build dataset
# initialize vars
dtapairs$SS <- 0
dtapairs$ST <- 0
dtapairs$TS <- 0
dtapairs$TT <- 0
## pairs caracts
dtapairs[(dtapairs$nati==1 & dtapairs$natj==1),"SS"] <- 1 # Syrian-Syrian pair
dtapairs[(dtapairs$nati==1 & dtapairs$natj==2),"ST"] <- 1 # Syrian named Turkish pair
dtapairs[(dtapairs$nati==2 & dtapairs$natj==1),"TS"] <- 1 # Turkish named Syrian pair
dtapairs[(dtapairs$nati==2 & dtapairs$natj==2),"TT"] <- 1 # Turkish-Turkish pair
dtapairs$samegender <- 0
dtapairs[(dtapairs$genderi==dtapairs$genderj),"samegender"] <- 1
dtapairs$school <- factor(dtapairs$school)

## Syrian-Syrian
regSS <- lm(gij2 ~ qSyrians + gij1 + samegender + school, data=dtapairs[dtapairs$SS==1,])
cluster.robust.se(regSS, dtapairs[dtapairs$SS==1,"class"])

## Syrian-Turkish
regST <- lm(gij2 ~ qSyrians + gij1 + samegender + school, data=dtapairs[dtapairs$ST==1,])
cluster.robust.se(regST, dtapairs[dtapairs$ST==1,"class"])

## Turkish-Turkish
regTT <- lm(gij2 ~ qSyrians + gij1 + samegender + school, data=dtapairs[dtapairs$TT==1,])
cluster.robust.se(regTT, dtapairs[dtapairs$TT==1,"class"])

## Turkish-Syrian
regTS <- lm(gij2 ~ qSyrians + gij1 + samegender + school, data=dtapairs[dtapairs$TS==1,])
cluster.robust.se(regTS, dtapairs[dtapairs$TS==1,"class"])

## fraction of "success"
mean(dtapairs[dtapairs$SS==1,"gij2"])
mean(dtapairs[dtapairs$ST==1,"gij2"])
mean(dtapairs[dtapairs$TT==1,"gij2"])
mean(dtapairs[dtapairs$TS==1,"gij2"])


#####################################
### Make network graphs examples ####
#####################################

# classroom 33 (size=18, q=)
# classroom 29 (size=29, q=0.43)
i <- 33
g <- graph_from_adjacency_matrix(G[[i]])
V(g)$color <- (1-(X[[i]][,3]-1))*2
V(g)$label <- NA
plot(g)
print(mean(X[[i]][,3]==1))
print(nrow(X[[i]]))

i <- 29
g <- graph_from_adjacency_matrix(G[[i]])
V(g)$color <- (1-(X[[i]][,3]-1))*2
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

