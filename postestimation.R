
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
set.seed(1234) # set seed

#####################################
########### Loads Dataset ###########
#####################################

setwd("~/Dropbox/local_summer_camp/structural_group_conformism_clean") # set working directory
load("estimWave2withcontrols.RData")
source("fcts_gmm.R") # load functions

gdata <- graphdata()
gdata <- as.data.frame(gdata)
colnames(gdata) <- c("q","hs","ht","ihs","iht") # name columns
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
 gammaname <- c("Syrian-Turkish pair", "Turkish-Syrian pair", "Same Gender","Parents live in diverse neighborhood$^\\dagger$","Syrian-Syrian pair","Fluent in Turkish (Syrians)", "Linked in wave 1", "School 1", "School 2", "School 3", "School 4", "School 5",
                "School 7", "School 8", "School 9", "School 10", "School 11")
 for (i in 1:length(mgamma)){
   cat(paste(gammaname[i], " & ", format(round(mgamma[i], 3), nsmall = 3), " & (", format(round(bsegamma[i], 3), nsmall = 3),") \\\\ \n",sep=""))
 }

######################################################
######################################################
############## policy simulations
######################################################
######################################################

ncounter <- 500 # number of simulations
thetat <- theta # takes the estimated value of theta
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),14)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","pss","ptt","pst","pts","ihs","iht","ss","st","mix","simu","class") # name columns
pos <- 1 # position index
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:11){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),12] <- mx
    results[pos:(pos+length(X)-1),13] <- sim
    results[pos:(pos+length(X)-1),14] <- 1:length(X)
    pos <- pos + length(X)
  }
}

results1 <- results # baseline results

thetat <- theta # takes the estimated value of theta
thetat[(length(bethat)+1):(length(bethat)+2)] <- 0 # remove ethnic bias
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),14)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","pss","ptt","pst","pts","ihs","iht","ss","st","mix","simu","class") # name columns
pos <- 1 # position index
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:11){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),12] <- mx
    results[pos:(pos+length(X)-1),13] <- sim
    results[pos:(pos+length(X)-1),14] <- 1:length(X)
    pos <- pos + length(X)
  }
}

results0 <- results # no bias results

thetat <- theta # takes the estimated value of theta
thetat[(length(bethat)+1):(length(bethat)+2)] <- thetat[(length(bethat)+1):(length(bethat)+2)]*2  # double ethnic bias
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),14)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","pss","ptt","pst","pts","ihs","iht","ss","st","mix","simu","class") # name columns
pos <- 1 # position index
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:11){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),12] <- mx
    results[pos:(pos+length(X)-1),13] <- sim
    results[pos:(pos+length(X)-1),14] <- 1:length(X)
    pos <- pos + length(X)
  }
}

results2 <- results # high biases results

thetat <- theta # takes the estimated value of theta
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),14)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","pss","ptt","pst","pts","ihs","iht","ss","st","mix","simu","class") # name columns
pos <- 1 # position index
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype_no_congestion(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:11){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),12] <- mx
    results[pos:(pos+length(X)-1),13] <- sim
    results[pos:(pos+length(X)-1),14] <- 1:length(X)
    pos <- pos + length(X)
  }
}

resultsnc <- results # high biases results

thetat <- theta # takes the estimated value of theta
mixes <- seq(from=0.2,to=0.8,by=0.2) # looks at weights between 0.2 and 0.8 (we do not learn anything from perfectly segregated classrooms)
results <- as.data.frame(matrix(NA,(length(X)*length(mixes)*ncounter),14)) # initialize results dataframe
colnames(results) <- c("q","hs","ht","pss","ptt","pst","pts","ihs","iht","ss","st","mix","simu","class") # name columns
pos <- 1 # position index
for (c in 1:length(X)){
  X[[c]][(X[[c]][,14]==0),14] <- 1 #Every Syrian is fluent
}
for (mx in mixes){ # for each mixing value
  for (sim in 1:ncounter){ # for ncounter simulations
    lstswap <- shufflepeople(mx) # shuffle students
    X <- swappeople(lstswap) # overwrite X
    E <- swaperrors(lstswap) # overwrite E
    Z <- buildZ() #  overwrite Z
    fracdata <- sametype(thetat) # simulate data
    
    ## store results in dataframe
    for (kk in 1:11){
      results[pos:(pos+length(X)-1),kk] <- fracdata[[kk]]
    }
    results[pos:(pos+length(X)-1),12] <- mx
    results[pos:(pos+length(X)-1),13] <- sim
    results[pos:(pos+length(X)-1),14] <- 1:length(X)
    pos <- pos + length(X)
  }
}

resultshs <- results # highskill results



results0$delta <- 0 # code no bias
results1$delta <- 1 # code baseline
results2$delta <- 2 # code high bias
resultsnc$delta <- 4 # code no congestion
resultshs$delta <- 5 # code high skill

### for preference biases comparisons
results <- rbind(results0,results1,results2) # combine datasets
results$fq <- with(results, cut(q, breaks=quantile(q, probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Syrians (deciles)
results$f1mq <- with(results, cut((100-q), breaks=quantile((100-q), probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Turkish (deciles)
results$f1mq <- results$fq # overwrite
levels(results$fq) <- levels(results$f1mq) <- as.character(1:10)

### for congestion comparisons
results_nc <- rbind(results1,resultsnc) # combine datasets
results_nc$fq <- with(results_nc, cut(q, breaks=quantile(q, probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Syrians (deciles)
results_nc$f1mq <- with(results_nc, cut((100-q), breaks=quantile((100-q), probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Turkish (deciles)
results_nc$f1mq <- results_nc$fq # overwrite
levels(results_nc$fq) <- levels(results_nc$f1mq) <- as.character(1:10)

### for highskill comparisons
results_hs <- rbind(results0,results1,resultshs) # combine datasets
results_hs$fq <- with(results_hs, cut(q, breaks=quantile(q, probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Syrians (deciles)
results_hs$f1mq <- with(results_hs, cut((100-q), breaks=quantile((100-q), probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Turkish (deciles)
results_hs$f1mq <- results_hs$fq # overwrite
levels(results_hs$fq) <- levels(results_hs$f1mq) <- as.character(1:10)

### for baseline graphs
results1$fq <- with(results1, cut(q, breaks=quantile(q, probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Syrians (deciles)
results1$f1mq <- with(results1, cut((100-q), breaks=quantile((100-q), probs=seq(0,1, by=0.1), na.rm=TRUE), include.lowest=TRUE)) # fraction of Turkish (deciles)
levels(results1$fq) <- levels(results1$f1mq) <- as.character(1:10)

theme_set(theme_minimal() + theme( axis.title = element_text(size=15), legend.position ="top",legend.direction = "horizontal" ,legend.text=element_text(size=12) ) ) # set theme for ggplot2

######################################################
######################################################
############## Print and save many graphs ############
######################################################
######################################################

mancolors <- brewer.pal(n=3, name="Blues")
preflevels <- c('No ethnic bias','Baseline ethnic bias','High ethnic bias')
conglevels <- c('Baseline','No congestion')
fluencylevels <- c('No ethnic bias', 'Baseline','No language barrier')

######### Baseline graphs ##############


# homophily indices
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=hs)) + xlab("Refugee share (deciles)") + ylab("Homophily index  (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=hs), colour=mancolors[3], method=lm, formula = y ~ poly(x,2)) + ylim(0,1) + ggsave("shom-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ht)) + xlab("Refugee share (deciles)") + ylab("Homophily index (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=ht), colour=mancolors[3], method=lm, formula = y ~ poly(x,2)) + ylim(0,1) + ggsave("thom-baseline.pdf")

# inbreeding homophily indices
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ihs)) + xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=ihs), colour=mancolors[3], method=lm, formula = y ~ poly(x,2)) + ylim(-1,1) + ggsave("ihs-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=iht)) + xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=iht), colour=mancolors[3], method=lm, formula = y ~ poly(x,2)) + ylim(-1,1) + ggsave("iht-baseline.pdf")


ffit <-lm(ihs ~ I(q/9) + I((q/9)^2),data=results1)$coefficients

# Socialization
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=ss)) + xlab("Refugee share (deciles)") + ylab("Socialization (Syrian kids)") +
  geom_smooth(aes(x=q/9, y=ss), colour=mancolors[3], method=lm, formula = y ~ poly(x,2)) + ylim(0,1) + ggsave("ss-baseline.pdf")
p <- ggplot(results1) + geom_boxplot(aes(x=factor(fq), y=st)) + xlab("Refugee share (deciles)") + ylab("Socialization (Turkish kids)") +
  geom_smooth(aes(x=q/9, y=st), colour=mancolors[3], method=lm, formula = y ~ poly(x,2)) + ylim(0,1) + ggsave("st-baseline.pdf")


########## Changes in preference biases #############

p <- ggplot(results, aes(x=factor(fq), y=ihs, fill=factor(delta,labels=preflevels)),color=mancolors) + geom_boxplot() +
  xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Syrian kids)") + labs(fill="Ethnic biases: ") + 
  geom_smooth(aes(x=q/9, y=ihs, fill=factor(delta,labels=preflevels)), colour="black", method=lm, formula = y ~ poly(x,2)) + scale_fill_manual(values=mancolors) + ylim(-1,1) +
  ggsave("ihs.pdf")

p <- ggplot(results, aes(x=factor(f1mq), y=iht, fill=factor(delta,labels=preflevels))) + geom_boxplot() +
  xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Turkish kids)") + labs(fill="Ethnic biases: ") + 
  geom_smooth(aes(x=q/9, y=iht, fill=factor(delta,labels=preflevels)), colour="black", method=lm, formula = y ~ poly(x,2)) + scale_fill_manual(values=mancolors) + ylim(-1,1) +
  ggsave("iht.pdf")

########## Changes in congestion #############

p <- ggplot(results_nc, aes(x=factor(fq), y=ihs, fill=factor(delta,labels=conglevels))) + geom_boxplot() +
  xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Syrian kids)") + labs(fill="Congestion: ") +
  geom_smooth(aes(x=q/9, y=ihs, fill=factor(delta,labels=conglevels)), colour="black", method=lm, formula = y ~ poly(x,2)) + scale_fill_manual(values=mancolors[2:3]) + ylim(-1,1) +
  ggsave("ihs-congestion.pdf")

p <- ggplot(results_nc, aes(x=factor(f1mq), y=iht, fill=factor(delta,labels=conglevels))) + geom_boxplot() +
  xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Turkish kids)") + labs(fill="Congestion: ") +
  geom_smooth(aes(x=q/9, y=iht, fill=factor(delta,labels=conglevels)), colour="black", method=lm, formula = y ~ poly(x,2)) + scale_fill_manual(values=mancolors[2:3]) + ylim(-1,1) +
  ggsave("iht-congestion.pdf")

########## Changes in fluency #############
flcolors <- c(mancolors[1],mancolors[3],mancolors[2])

p <- ggplot(results_hs, aes(x=factor(fq), y=ihs, fill=factor(delta,labels=fluencylevels))) + geom_boxplot() +
  xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Syrian kids)") + labs(fill="Fluency: ") +
  geom_smooth(aes(x=q/9, y=ihs, fill=factor(delta,labels=fluencylevels)), colour="black", method=lm, formula = y ~ poly(x,2)) + scale_fill_manual(values=flcolors) + ylim(-1,1) +
  ggsave("ihs-fluent.pdf")

p <- ggplot(results_hs, aes(x=factor(f1mq), y=iht, fill=factor(delta,labels=fluencylevels))) + geom_boxplot() +
  xlab("Refugee share (deciles)") + ylab("Inbreeding Homophily (Turkish kids)") + labs(fill="Fluency: ") +
  geom_smooth(aes(x=q/9, y=iht, fill=factor(delta,labels=fluencylevels)), colour="black", method=lm, formula = y ~ poly(x,2)) + scale_fill_manual(values=flcolors) + ylim(-1,1) +
  ggsave("iht-fluent.pdf")


########## data graphs ############

nlist <- rep(0,length(X))
for (i in 1:length(X)){
  nlist[i] <- nrow(X[[i]])
}

p <- ggplot() + xlab("Refugee share (%)") + ylab("Inbreeding Homophily (Syrian kids)") +
  geom_smooth(aes(x=q, y=ihs), colour=mancolors[3], method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(-1,1) + geom_point(aes(x=q, y=ihs), data=gdata) + ggsave("ihs-data.pdf")
p <- ggplot() + xlab("Refugee share (%)") + ylab("Inbreeding Homophily (Turkish kids)") +
  geom_smooth(aes(x=q, y=iht), colour=mancolors[3], method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(-1,1) + geom_point(aes(x=q, y=iht), data=gdata) + ggsave("iht-data.pdf")

p <- ggplot() + xlab("Refugee share (%)") + ylab("Homophily index (Syrian kids)") +
  geom_smooth(aes(x=q, y=hs), colour=mancolors[3], method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(0,1) + geom_point(aes(x=q, y=hs), data=gdata) + ggsave("hs-data.pdf")
p <- ggplot() + xlab("Refugee share (%)") + ylab("Homophily index (Turkish kids)") +
  geom_smooth(aes(x=q, y=ht), colour=mancolors[3], method=lm, formula = y ~ poly(x,2),data=gdata) + ylim(0,1) + geom_point(aes(x=q, y=ht), data=gdata) + ggsave("ht-data.pdf")


save.image("estimWave2withcontrolsandgraphs.RData")
