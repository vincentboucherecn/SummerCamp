
################################################################
################################################################
#################### Data and Parameters #######################
################################################################
################################################################


rm(list = ls())
## Preliminaries
library(numDeriv)
library(Matrix)
library(maxLik)
library(ggplot2)
library(tidyverse)
library(ggridges)
library(mvtnorm)
library(wesanderson)
library(RColorBrewer)
set.seed(1234) # set seed

#####################################
########### Loads Dataset ###########
#####################################

setwd("~/Dropbox/local_summer_camp/structural_group_conformism") # set working directory
load("estimWave2withcontrols.RData")
source("fcts_gmm.R") # load functions

gdata <- graphdata()
gdata <- as.data.frame(gdata)
colnames(gdata) <- c("q","hs","ht","ss","st","dss","dtt","dst","dts","pss","ptt","pst","pts","ihs","iht","sfrac","tfrac","bs","bt") # name columns
gdata$fq <- with(gdata, cut(q, breaks=quantile(q, probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Syrians (deciles)
levels(gdata$fq) <- as.character(1:10)

################################################################
################################################################
##################### Post-estimation ## #######################
################################################################
################################################################


###
### Marginal bias:
###
mbias <-  marginal_bias(theta)[[1]]
mgamma <- marginal_bias(theta)[[2]]

## simulate SE
 bmgamma <- matrix(NA,length(mgamma),500)
 for (boot in 1:500){
   btheta <- rmvnorm(1,theta,VC)
   bmgamma[,boot] <- marginal_bias(btheta)[[2]]
 }
 bsegamma <- apply(bmgamma,1,sd)
 gammaname <- c("Syrian-Turkish pair", "Turkish-Syrian pair", "Same Gender", "Linked in wave 1", "School 1", "School 2", "School 3", "School 4", "School 5",
                "School 7", "School 8", "School 9", "School 10", "School 11")
 for (i in 1:length(mgamma)){
   cat(paste(gammaname[i], " & ", format(round(mgamma[i], 3), nsmall = 3), " & (", format(round(bsegamma[i], 3), nsmall = 3),") \\\\ \n",sep=""))
 }

######################################################
######################################################
############## policy simulations
######################################################
######################################################

killhet <- 1 # keeps only gender, nationality and wave1 outcomes

if (killhet == 1){
  #theta[(length(bethat)-10):length(bethat)] <- mean(theta[(length(bethat)-10):length(bethat)])
  #theta[(length(bethat)+length(gamhat)-9):(length(bethat)+length(gamhat))] <- mean(theta[(length(bethat)+length(gamhat)-9):(length(bethat)+length(gamhat))])
  theta[3:18] <- 0  
}
 
ncounter <- 100 # number of simulations
thetat <- theta # takes the estimated value of theta
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),22)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","ss","st","dss","dtt","dst","dts","pss","ptt","pst","pts","ihs","iht","sfrac","tfrac","bs","bt","mix","simu","class") # name columns
pos <- 1 # position index
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:19){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),20] <- mx
    results[pos:(pos+length(X)-1),21] <- sim
    results[pos:(pos+length(X)-1),22] <- 1:length(X)
    pos <- pos + length(X)
  }
}

results1 <- results # baseline results

thetat <- theta # takes the estimated value of theta
thetat[(length(bethat)+1):(length(bethat)+2)] <- 0 # remove ethnic bias
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),22)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","ss","st","dss","dtt","dst","dts","pss","ptt","pst","pts","ihs","iht","sfrac","tfrac","bs","bt","mix","simu","class") # name columns
pos <- 1 # position index
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:19){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),20] <- mx
    results[pos:(pos+length(X)-1),21] <- sim
    results[pos:(pos+length(X)-1),22] <- 1:length(X)
    pos <- pos + length(X)
  }
}

results0 <- results # no bias results

thetat <- theta # takes the estimated value of theta
thetat[(length(bethat)+1):(length(bethat)+2)] <- thetat[(length(bethat)+1):(length(bethat)+2)]*2  # double ethnic bias
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),22)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","ss","st","dss","dtt","dst","dts","pss","ptt","pst","pts","ihs","iht","sfrac","tfrac","bs","bt","mix","simu","class") # name columns
pos <- 1 # position index
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:19){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),20] <- mx
    results[pos:(pos+length(X)-1),21] <- sim
    results[pos:(pos+length(X)-1),22] <- 1:length(X)
    pos <- pos + length(X)
  }
}

results2 <- results # high biases results

results0$delta <- 0 # code no bias
results1$delta <- 1 # code baseline
results2$delta <- 2 # code high bias

results <- rbind(results0,results1,results2) # combine datasets
results$fq <- with(results, cut(q, breaks=quantile(q, probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Syrians (deciles)
results$f1mq <- with(results, cut((100-q), breaks=quantile((100-q), probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Turkish (deciles)
levels(results$fq) <- levels(results$f1mq) <- as.character(1:10)

theme_set(theme_minimal() + theme( axis.title = element_text(size=15), legend.position ="top",legend.direction = "horizontal")) # set theme for ggplot2

######################################################
######################################################
############## Print and save many graphs ############
######################################################
######################################################

p <- ggplot(results, aes(x=factor(fq), y=ss, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Syrians (deciles)") + ylab("Socialization of Syrians") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(0.1,0.75) + ggsave("social_syrian.pdf")

p <- ggplot(results, aes(x=factor(f1mq), y=st, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Turkish (deciles)") + ylab("Socialization of Turkish") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(0.1,0.75) + ggsave("social_turkish.pdf")

p <- ggplot(results, aes(x=factor(fq), y=pss, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Syrians (deciles)") + ylab("Probability of a Syrian-Syrian link") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(0,0.4) + ggsave("pss.pdf")

p <- ggplot(results, aes(x=factor(fq), y=pst, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Syrians (deciles)") + ylab("Probability of a Syrian-Turkish link") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(0,0.4) + ggsave("pst.pdf")

p <- ggplot(results, aes(x=factor(f1mq), y=ptt, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Turkish (deciles)") + ylab("Probability of a Turkish-Turkish link") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(0,0.4) + ggsave("ptt.pdf")

p <- ggplot(results, aes(x=factor(f1mq), y=pts, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Turkish (deciles)") + ylab("Probability of a Turkish-Syrian link") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(0,0.4) + ggsave("pts.pdf")


p <- ggplot(results, aes(x=factor(fq), y=ihs, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Syrians (deciles)") + ylab("Inbreeding Homophily (Syrian kids)") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(-1,1) + ggsave("ihs.pdf")

p <- ggplot(results, aes(x=factor(f1mq), y=iht, fill=factor(delta,labels=c('No ethnic bias','Baseline ethnic bias','High ethnic bias')))) + geom_boxplot() +
  xlab("Share of Turkish (deciles)") + ylab("Inbreeding Homophily (Turkish kids)") + labs(fill="Ethnic biases: ") + scale_fill_manual(values=brewer.pal(n=3, name="Paired")) + ylim(-1,1) + ggsave("iht.pdf")

reghs <- lm(ihs ~ q + I(q^2), data=results1)
summary(reghs)

reght <- lm(iht ~ I(100-q) + I((100-q)^2), data=results1)
summary(reght)

results1$fq <- with(results1, cut(q, breaks=quantile(q, probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Syrians (deciles)
results1$f1mq <- with(results1, cut((100-q), breaks=quantile((100-q), probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Turkish (deciles)
levels(results1$fq) <- levels(results1$f1mq) <- as.character(1:10)
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=pss)) + xlab("Refugee share (deciles)") + ylab("Probability of a Syrian-Syrian link") +
  geom_smooth(aes(x=q/9, y=pss), method=lm, formula = y ~ poly(x,1)) + ylim(0,0.4) + ggsave("pss-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=pst)) + xlab("Refugee share (deciles)") + ylab("Probability of a Syrian-Turkish link") + 
  geom_smooth(aes(x=q/9, y=pst), method=lm, formula = y ~ poly(x,1)) + ylim(0,0.4) + ggsave("pst-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=pts)) + xlab("Refugee share (deciles)") + ylab("Probability of a Turkish-Syrian link") +
  geom_smooth(aes(x=q/9, y=pts), method=lm, formula = y ~ poly(x,1)) + ylim(0,0.4) + ggsave("pts-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ptt)) + xlab("Refugee share (deciles)") + ylab("Probability of a Turkish-Turkish link") +
  geom_smooth(aes(x=q/9, y=ptt), method=lm, formula = y ~ poly(x,1)) + ylim(0,0.4) + ggsave("ptt-baseline.pdf")

p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ihs)) + xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=ihs), method=lm, formula = y ~ poly(x,2)) + ylim(-1,1) + ggsave("ihs-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=iht)) + xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=iht), method=lm, formula = y ~ poly(x,2)) + ylim(-1,1) + ggsave("iht-baseline.pdf")

p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ss)) + xlab("Refugee share (deciles)") + ylab("Socialization of Syrians") +
  geom_smooth(aes(x=q/9, y=ss), method=lm, formula = y ~ poly(x,1)) + ylim(0.1,0.75) + ggsave("social_syrian-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=st)) + xlab("Refugee share (deciles)") + ylab("Socialization of Turkish") +
  geom_smooth(aes(x=q/9, y=st), method=lm, formula = y ~ poly(x,1), color = "indianred3") + ylim(0.1,0.75) + ggsave("social_turkish-baseline.pdf")

ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=st), width=0.5) + xlab("Refugee share (deciles)") + ylab("Socialization of Turkish") +
  geom_smooth(aes(x=q/9, y=st), method=lm, formula = y ~ poly(x,1)) + ylim(0.1,0.75)

ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ss)) + xlab("Refugee share (deciles)") + ylab("Socialization of Syrians") +
  geom_smooth(aes(x=q/9, y=ss), method=lm, formula = y ~ poly(x,1)) + ylim(0.1,0.75)

p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=hs)) + xlab("Refugee share (deciles)") + ylab("Homophily index  (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=hs), method=lm, formula = y ~ poly(x,2)) + ylim(0,1) + ggsave("shom-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ht)) + xlab("Refugee share (deciles)") + ylab("Homophily index (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=ht), method=lm, formula = y ~ poly(x,2)) + ylim(0,1) + ggsave("thom-baseline.pdf")

p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=sfrac)) + xlab("Refugee share (deciles)") + ylab("Number of Turkish friends  (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=sfrac), method=lm, formula = y ~ poly(x,1)) + ylim(0,1) + ggsave("sntf-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=tfrac)) + xlab("Refugee share (deciles)") + ylab("Number of Syrian friends (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=tfrac), method=lm, formula = y ~ poly(x,1)) + ylim(0,1) + ggsave("tnsf-baseline.pdf")

nlist <- rep(0,length(X))
for (i in 1:length(X)){
  nlist[i] <- nrow(X[[i]])
}

p <- ggplot() + xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=ihs), method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(-1,1) + scale_x_discrete(limits = as.character(1:10)) + geom_point(aes(x=q/9, y=ihs), data=gdata) + ggsave("ihs-data.pdf")
p <- ggplot() + xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=iht), method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(-1,1) + scale_x_discrete(limits = as.character(1:10)) + geom_point(aes(x=q/9, y=iht), data=gdata) + ggsave("iht-data.pdf")

p <- ggplot() + xlab("Refugee share (deciles)") + ylab("Homophily index (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=hs), method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(0,1) + scale_x_discrete(limits = as.character(1:10)) + geom_point(aes(x=q/9, y=hs), data=gdata) + ggsave("hs-data.pdf")
p <- ggplot() + xlab("Refugee share (deciles)") + ylab("Homophily index (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=ht), method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(0,1) + scale_x_discrete(limits = as.character(1:10)) + geom_point(aes(x=q/9, y=ht), data=gdata) + ggsave("ht-data.pdf")


save.image("estimWave2withcontrolsandgraphs.RData")
