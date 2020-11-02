
################################################################
################################################################
#################### Data and Parameters #######################
################################################################
################################################################


rm(list = ls())
## Preliminaries
library(ggplot2)
library(mvtnorm)
library(RColorBrewer)
library(ivpack)
set.seed(1234) # set seed

#####################################
########### Loads Dataset ###########
#####################################

setwd("~/Dropbox/local_summer_camp/structural_group_conformism_clean") # set working directory
load("estimWave2withcontrols.RData")
source("fcts_gmm.R") # load functions

nsim <- 500 # overwrite the number of simulations
# generates list of errors (epsilon using the paper's notation)
E <- vector("list",length(X))
for (c in 1:length(X)){
  nt <- nrow(X[[c]])
  E[[c]] <- matrix(rnorm(nt*nsim),nt,nsim)
}

################################################################
################################################################
########################## Model Fit ###########################
################################################################
################################################################

nbsim <- 500
resall <- as.data.frame(matrix(0,nbsim,6))
colnames(resall) <- c("pss","pst","ptt","pts","ss","st")
dtalevels <- resall

for (bsim in 1:nbsim){
  thetat <- rmvnorm(n=1,theta,VC) # random draw of theta
  allout <- gendata(thetat,1)
  d1 <- allout[[1]]
  rd1 <- lm(s ~ female + factor(school) + factor(ayear) + factor(region) + lang + friends + neighbours + q, data=d1[d1$syrian==1,])
  resall[bsim,"ss"] <- rd1$coefficients["q"]
  rd1 <- lm(s ~ female + factor(school) + friends + neighbours + q, data=d1[d1$syrian==0,])
  resall[bsim,"st"] <- rd1$coefficients["q"]
  
  d2 <- allout[[2]]
  rd2 <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==1,])
  resall[bsim,"pss"] <- rd2$coefficients["q"]
  
  rd2 <- lm(g ~ sgender + factor(school) + parents + skill + wave1 + q, data=d2[d2$type==2,])
  resall[bsim,"pst"] <- rd2$coefficients["q"]
  
  rd2 <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==3,])
  resall[bsim,"ptt"] <- rd2$coefficients["q"]
  
  rd2 <- lm(g ~ sgender + factor(school) + parents + skill + wave1 + q, data=d2[d2$type==4,])
  resall[bsim,"pts"] <- rd2$coefficients["q"]
  
  dtalevels[bsim,"ss"] <- mean(d1[d1$syrian==1,"s"])
  dtalevels[bsim,"st"] <- mean(d1[d1$syrian==0,"s"])
  dtalevels[bsim,"pss"] <- mean(d2[d2$type==1,"g"])
  dtalevels[bsim,"pst"] <- mean(d2[d2$type==2,"g"])
  dtalevels[bsim,"ptt"] <- mean(d2[d2$type==3,"g"])
  dtalevels[bsim,"pts"] <- mean(d2[d2$type==4,"g"])
}

allout <- gendata(thetat,0)
d1 <- allout[[1]]
rdss <- lm(s ~ female + factor(school) + factor(ayear) + factor(region) + lang + friends + neighbours + q, data=d1[d1$syrian==1,])
sdss <- cluster.robust.se(rdss, d1[d1$syrian==1,"class"])

rdst <- lm(s ~ female + factor(school) + friends + neighbours + q, data=d1[d1$syrian==0,])
sdst <- cluster.robust.se(rdst, d1[d1$syrian==0,"class"])

d2 <- allout[[2]]
rd2ss <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==1,])
sd2ss <- cluster.robust.se(rd2ss, d2[d2$type==1,"class"])

rd2st <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==2,])
sd2st <- cluster.robust.se(rd2st, d2[d2$type==2,"class"])

rd2tt <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==3,])
sd2tt <-cluster.robust.se(rd2tt, d2[d2$type==3,"class"])

rd2ts <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==4,])
sd2ts <- cluster.robust.se(rd2ts, d2[d2$type==4,"class"])

reallevels_mean <- (c(mean(d2[d2$type==1,"g"]),mean(d2[d2$type==2,"g"]),mean(d2[d2$type==3,"g"]),mean(d2[d2$type==4,"g"]),mean(d1[d1$syrian==1,"s"]),mean(d1[d1$syrian==0,"s"]) ))



## list of variables names
varnames <- c("$p^{SS}$","$p^{ST}$","$p^{TT}$","$p^{TS}$","$s^{S}$","$s^{T}$")
simav <- as.numeric(colMeans(resall))*100
simse <- as.numeric(apply(resall,2,sd))*100
dtaav <- c(rd2ss$coefficients["q"],rd2st$coefficients["q"],rd2tt$coefficients["q"],rd2ts$coefficients["q"],rdss$coefficients["q"],rdst$coefficients["q"])*100
dtase <- c(sd2ss[14,2],sd2st[14,2],sd2tt[14,2],sd2ts[14,2],sdss[29,2],sdst[15,2])*100
vnobs <- c(nrow(d2[d2$type==1,]),nrow(d2[d2$type==2,]),nrow(d2[d2$type==3,]),nrow(d2[d2$type==4,]),nrow(d1[d1$syrian==1,]),nrow(d1[d1$syrian==0,]) )
### print as a LaTeX table

print("**** baseline, refugee shares ****")

for (i in 1:length(varnames)){
  cat(paste(varnames[i], " & ", vnobs[i], " & ", format(round(simav[i], 3), nsmall = 3), " & (", format(round(simse[i], 3), nsmall = 3),") & ",
            format(round(dtaav[i], 3), nsmall = 3), " & (", format(round(dtase[i], 3), nsmall = 3),") \\\\ \n",sep=""))
}

## list of variables names
varnames <- c("$g^{SS}$","$g^{ST}$","$g^{TT}$","$g^{TS}$","$s^{S}$","$s^{T}$")
lsimav <- as.numeric(colMeans(dtalevels))
lsimse <- as.numeric(apply(dtalevels,2,sd))

print("**** baseline, levels ****")
### print as a LaTeX table
for (i in 1:length(varnames)){
  cat(paste(varnames[i], " & ", format(round(lsimav[i], 3), nsmall = 3), " & (", format(round(lsimse[i], 3), nsmall = 3),") & ",
            format(round(reallevels_mean[i], 3), nsmall = 3), " \\\\ \n",sep=""))
}


#########################################################################
#########################################################################
################### counterfactuals turkish skill #######################
#########################################################################
#########################################################################
# 
# resall_sk <- as.data.frame(matrix(0,nbsim,6))
# colnames(resall_sk) <- c("pss","pst","ptt","pts","ss","st")
# dtalevels_sk <- resall_sk
# 
# for (bsim in 1:nbsim){
#   thetat <- rmvnorm(n=1,theta,VC) # random draw of theta
#   allout <- gendata_highskill(thetat,1)
#   d1 <- allout[[1]]
#   rd1 <- lm(s ~ female + factor(school) + factor(ayear) + factor(region) + lang + friends + neighbours + q, data=d1[d1$syrian==1,])
#   resall_sk[bsim,"ss"] <- rd1$coefficients["q"]
#   rd1 <- lm(s ~ female + factor(school) + friends + neighbours + q, data=d1[d1$syrian==0,])
#   resall_sk[bsim,"st"] <- rd1$coefficients["q"]
#   
#   d2 <- allout[[2]]
#   rd2 <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==1,])
#   resall_sk[bsim,"pss"] <- rd2$coefficients["q"]
#   
#   rd2 <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==2,])
#   resall_sk[bsim,"pst"] <- rd2$coefficients["q"]
#   
#   rd2 <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==3,])
#   resall_sk[bsim,"ptt"] <- rd2$coefficients["q"]
#   
#   rd2 <- lm(g ~ sgender + factor(school) + wave1 + q, data=d2[d2$type==4,])
#   resall_sk[bsim,"pts"] <- rd2$coefficients["q"]
#   
#   dtalevels_sk[bsim,"ss"] <- mean(d1[d1$syrian==1,"s"])
#   dtalevels_sk[bsim,"st"] <- mean(d1[d1$syrian==0,"s"])
#   dtalevels_sk[bsim,"pss"] <- mean(d2[d2$type==1,"g"])
#   dtalevels_sk[bsim,"pst"] <- mean(d2[d2$type==2,"g"])
#   dtalevels_sk[bsim,"ptt"] <- mean(d2[d2$type==3,"g"])
#   dtalevels_sk[bsim,"pts"] <- mean(d2[d2$type==4,"g"])
# }
# 
# 
# 
# simav1 <- as.numeric(colMeans(resall_sk))*100
# simse1 <- as.numeric(apply(resall_sk,2,sd))*100
# levav1 <- as.numeric(colMeans(dtalevels_sk))
# levsd1 <- as.numeric(apply(dtalevels_sk,2,sd))
# levav0 <- as.numeric(colMeans(dtalevels))
# levsd0 <- as.numeric(apply(dtalevels,2,sd))
# 
# print("**** High skill, all ****")
# ### print as a LaTeX table
# for (i in 1:length(varnames)){
#   cat(paste(varnames[i], " & ", format(round(simav[i], 3), nsmall = 3), " & (", format(round(simse[i], 3), nsmall = 3),") & ",
#             format(round(simav1[i], 3), nsmall = 3), " & (", format(round(simse1[i], 3), nsmall = 3),") & ",
#             format(round(levav0[i], 3), nsmall = 3), " & (", format(round(levsd0[i], 3), nsmall = 3),") & ",
#             format(round(levav1[i], 3), nsmall = 3), " & (", format(round(levsd1[i], 3), nsmall = 3), ") \\\\ \n",sep=""))
# }


#########################################################################
#########################################################################
########################## matched moments ##############################
#########################################################################
#########################################################################

# 
# mmatch <- as.data.frame(matrix(0,nbsim,(29+ 2*18+length(gamhat))))
# colnames(mmatch) <- c("Syrian-Turkish pair",
#                          "Turkish-Syrian pair", "Same Gender", "Parents live in diverse neighborhood$^\\dagger$","Fluent in Turkish (Syrians)", "Linked in wave 1", "School 1", "School 2", "School 3", "School 4", "School 5",
#                          "School 7", "School 8", "School 9", "School 10", "School 11",
#                          
#                          "Syrian", "Male", "Parents have diverse friends$^\\dagger$", "Parents live in diverse Neighborhood$^\\dagger$",
#                          "Arrived in 2012 (Syrians)", "Arrived in 2013 (Syrians)", "Arrived in 2014 (Syrians)", "Arrived in 2015 (Syrians)",
#                          "Arrived in 2016 (Syrians)", "Arrived in 2017 (Syrians)", "Arrived in 2018 (Syrians)", "Syrian region 3 (Syrians)",
#                          "Syrian region 4 (Syrians)", "Syrian region 5 (Syrians)", "Syrian region 6 (Syrians)", "Syrian region 7 (Syrians)",
#                          "Syrian region 8 (Syrians)", "Fluent in Turkish (Syrians)", "School 1", "School 2", "School 3",
#                          "School 4", "School 5", "School 6", "School 7", "School 8", "School 9", "School 10", "School 11",
#                          
#                       "($\\phi^S$) Syrian", "($\\phi^S$) Male", "($\\phi^S$) Parents have diverse friends$^\\dagger$", "($\\phi^S$) Parents live in diverse Neighborhood$^\\dagger$",
#                       "($\\phi^S$) Arrived in 2012 (Syrians)", "($\\phi^S$) Arrived in 2013 (Syrians)", "($\\phi^S$) Arrived in 2014 (Syrians)", "($\\phi^S$) Arrived in 2015 (Syrians)",
#                       "($\\phi^S$) Arrived in 2016 (Syrians)", "($\\phi^S$) Arrived in 2017 (Syrians)", "($\\phi^S$) Arrived in 2018 (Syrians)", "($\\phi^S$) Syrian region 3 (Syrians)",
#                       "($\\phi^S$) Syrian region 4 (Syrians)", "($\\phi^S$) Syrian region 5 (Syrians)", "($\\phi^S$) Syrian region 6 (Syrians)", "($\\phi^S$) Syrian region 7 (Syrians)",
#                       "($\\phi^S$) Syrian region 8 (Syrians)", "($\\phi^S$) Fluent in Turkish (Syrians)",
#                       
#                       "($\\phi^T$) Syrian", "($\\phi^T$) Male", "($\\phi^T$) Parents have diverse friends$^\\dagger$", "($\\phi^T$) Parents live in diverse Neighborhood$^\\dagger$",
#                       "($\\phi^T$) Arrived in 2012 (Syrians)", "($\\phi^T$) Arrived in 2013 (Syrians)", "($\\phi^T$) Arrived in 2014 (Syrians)", "($\\phi^T$) Arrived in 2015 (Syrians)",
#                       "($\\phi^T$) Arrived in 2016 (Syrians)", "($\\phi^T$) Arrived in 2017 (Syrians)", "($\\phi^T$) Arrived in 2018 (Syrians)", "($\\phi^T$) Syrian region 3 (Syrians)",
#                       "($\\phi^T$) Syrian region 4 (Syrians)", "($\\phi^T$) Syrian region 5 (Syrians)", "($\\phi^T$) Syrian region 6 (Syrians)", "($\\phi^T$) Syrian region 7 (Syrians)",
#                       "($\\phi^T$) Syrian region 8 (Syrians)", "($\\phi^T$) Fluent in Turkish (Syrians)")
# 
# for (bsim in 1:nbsim){
#   thetat <- rmvnorm(n=1,theta,VC) # random draw of theta
#   mmatch[bsim,] <- c(GMMjacob1(thetat),GMMjacob2(thetat))
# }
# 
# mommean <- colMeans(mmatch)
# momstdv <- apply(mmatch,2,sd)
# 
# for (i in 1:ncol(mmatch)){
#   cat(paste(colnames(mmatch)[i], " & ", format(round(mommean[i], 3), nsmall = 3), " & (", format(round(momstdv[i], 3), nsmall = 3), ") \\\\ \n",sep=""))
# }
# 
# 
# simdata <- gendata(theta,1)
# rdata <- gendata(theta,0)
# gsimdata <- simdata[[2]]
# grdata <- rdata[[2]]
# mytable <- table(grdata$g,gsimdata$g)
# library(gmodels)
# CrossTable(grdata$g,gsimdata$g)
# prop.table(mytable) # row percentages
# 

save.image("estimWave2withcontrols_modelfit.RData")

