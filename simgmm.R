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
set.seed(1234) # set seed

#####################################
############ Parameters #############
#####################################

nsim <- 100 # number of simulations for numerical integration

#####################################
########### Loads Dataset ###########
#####################################

setwd("~/Dropbox/local_summer_camp/structural_group_conformism_clean") # set working directory
dta <- read.csv("./Syria_data.csv") # load data (csv file directly exported from the excel file)
dta0 <- dta # backup dataset
dta$st_id <- as.numeric(dta$Name) # numerical id for students
dta$sch_id <- as.numeric(dta$School) # numerical id for schools
dta$numnat <- as.numeric(dta$Nationality) # numerical id for nationality
dta <- dta[dta$sch_id != 9,] # remove school Mevlana Sosyal Tesis
dta$memp <- 1 # mother employment variable
dta[(as.numeric(dta$Mother.current.occ)==3 |as.numeric(dta$Mother.current.occ)==5),"memp"] <- 0 # mother employment variable
dta$femp <- 1 # father employment variable
dta[(as.numeric(dta$Father.current.occ)==3 |as.numeric(dta$Father.current.occ)==5),"femp"] <- 0  # father employment variable
dta$ayear <- as.numeric(substr(dta$Arrival.date,4,8)) # year of arrival
dta[which(dta$ayear==2104),"ayear"] <- 2014 # recode typo
dta$Gender <- as.numeric(dta$Gender) # gender to numeric: FEMALE=1, MALE=2
dta[dta$Region.in.Syria==".","Region.in.Syria"] <- NA # recode region in Syria
dta$Region.in.Syria <- as.numeric(dta$Region.in.Syria) # recode region in Syria
dta[dta$Language==".","Language"] <- NA  # recode language
dta$Language <- as.numeric(dta$Language) # language as numeric
dta[dta$Language %in% c(6,7),"Language"] <- 1 # recode language 1=fluent
dta[dta$Language %in% c(2:5),"Language"] <- 0 # recode language 0= not fluent
dta$Sticker.left <- as.numeric(dta$Sticker.left)/10 # fraction of stickers given
dta$sch_id_old <- dta$sch_id # backup school id
dta[dta$sch_id>9,"sch_id"] <- dta[dta$sch_id>9,"sch_id"]-1 # recode school ids
for (var in c("Name","F1","F2","F3","F4","F5")){
  dta[,var] <- as.character(dta[,var])
}
source("fcts_gmm.R") # load functions

#####################################
########## Builds Dataset ###########
#####################################

biglist <- builddata() # format data
X1 <- biglist[[1]] # Wave 1 individual variables
X2 <- biglist[[2]] # Wave 2 individual variables
G1 <- biglist[[3]] # Wave 1 network
G2 <- biglist[[4]] # Wave 2 network
rm(biglist)

## estimates based on second wave
G <- G2
X <- X2 #(stickers and language are wave 1)

# generates list of errors (epsilon using the paper's notation)
E <- vector("list",length(X))
for (c in 1:length(X)){
  nt <- nrow(X[[c]])
  E[[c]] <- matrix(rnorm(nt*nsim),nt,nsim)
}

## computes the number of pairs
n1 <- 0 # number of pairs
n2 <- 0 # number of individuals
for (c in 1:length(X)){
  nt <- nrow(X[[c]])
  n1 <- n1 + nt*(nt-1)
  n2 <- n2 + nt
}

Z <- buildZ() # builds list of pair characteristics


################################################################
################################################################
################# Calibrate starting values ####################
################################################################
################################################################
dta2SLS <- TSLS_comp()

###############
######## Syrian 2SLS
###############

expls <- c("Male", "ParentsF", "ParentsN", "A2012", "A2013", "A2014", "A2015", "A2016", "A2017", "A2018",
          "Region3", "Region4", "Region5", "Region6", "Region7", "Region8", "Fluent", "S1", "S2", "S3", "S4", "S5",
          "S6", "S7", "S8", "S9", "S10", "S11")
instr <- c("GSyrian", "GMale", "GParentsF", "GParentsN", "GA2012", "GA2013", "GA2014", "GA2015", "GA2016", "GA2017", "GA2018",
           "GRegion3", "GRegion4", "GRegion5", "GRegion6", "GRegion7", "GRegion8", "GFluent")

formus <- paste(c("Y ~ 0 +", paste(expls,collapse=" + ")," + GY | . -GY + "),collapse = "")
formus <- paste(formus, paste(instr,collapse=" + ") )
outs <- ivreg(as.formula(formus),data=dta2SLS[dta2SLS$Syrian==1,] )

###############
######## Turkish 2SLS
###############

explt <- c("Male", "ParentsF", "ParentsN", "S1", "S2", "S3", "S4", "S5",
           "S6", "S7", "S8", "S9", "S10", "S11")
instr <- c("GSyrian", "GMale", "GParentsF", "GParentsN", "GA2012", "GA2013", "GA2014", "GA2015", "GA2016", "GA2017", "GA2018",
           "GRegion3", "GRegion4", "GRegion5", "GRegion6", "GRegion7", "GRegion8", "GFluent")

formut <- paste(c("Y ~ 0 +", paste(explt,collapse=" + ")," + GY | . -GY + "),collapse = "")
formut <- paste(formut, paste(instr,collapse=" + ") )
outt <- ivreg(as.formula(formut),data=dta2SLS[dta2SLS$Syrian==0,] )

###############
######## Total 2SLS
###############

expl <- c("Syrian", "Male", "ParentsF", "ParentsN", "A2012", "A2013", "A2014", "A2015", "A2016", "A2017", "A2018",
          "Region3", "Region4", "Region5", "Region6", "Region7", "Region8", "Fluent", "S1", "S2", "S3", "S4", "S5",
          "S6", "S7", "S8", "S9", "S10", "S11")
instr <- c("GSyrian", "GMale", "GParentsF", "GParentsN", "GA2012", "GA2013", "GA2014", "GA2015", "GA2016", "GA2017", "GA2018",
           "GRegion3", "GRegion4", "GRegion5", "GRegion6", "GRegion7", "GRegion8", "GFluent")

formu <- paste(c("Y ~ 0 +", paste(expl,collapse=" + ")," + GY | . -GY + "),collapse = "")
formu <- paste(formu, paste(instr,collapse=" + ") )
out_tsls <- ivreg(as.formula(formu),data=dta2SLS )

#########
#### Starting values
#########

bethat <- out_tsls$coefficients[1:(length(out_tsls$coefficients)-1)]/(1-out_tsls$coefficients[length(out_tsls$coefficients)])
phihat <- c(outs$coefficients[length(outs$coefficients)],outt$coefficients[length(outt$coefficients)])

lstvars <- c(1:5,7,8) # variables to use in Z

## Initial values for gamma are computed using the conditional likelihood
pr <- maxLik(likprob,start=runif((length(lstvars)+11)),method="BHHH")
print(summary(pr))
gamhat <- pr$estimate
sighat <- gamhat[length(gamhat)] # point estimate for sigma
gamhat <- gamhat[1:(length(gamhat)-1)]  # keep point estimate for gamma

################################################################
################################################################
################### Launch the estimation ######################
################################################################
################################################################

theta0 <- c(bethat,gamhat,log(phihat/(1-phihat)),log(sighat)) # initial values (including variable change for positive values)
out <- optim(theta0,GMM) # minimize GMM objective (unweighted)
theta <- out$par # theta is estimator
theta[(length(theta)-2):(length(theta)-1)] <- exp(theta[(length(theta)-2):(length(theta)-1)])/(1+exp(theta[(length(theta)-2):(length(theta)-1)])) # redo the variable change for phi
theta[length(theta)] <- exp(theta[length(theta)]) # redo variable change for sigma

allW <- compW(theta) # compute optimal weighted matrices
W1 <- allW[[1]] # first set of moments
W2 <- allW[[2]] # second set of moments

theta0 <- c(bethat,gamhat,log(phihat/(1-phihat)),log(sighat)) # initial values (including variable change for positive values)

out <- optim(theta0,WGMM) # minimize GMM objective (optimal)
theta <- out$par # theta is estimator
theta[(length(theta)-2):(length(theta)-1)] <- exp(theta[(length(theta)-2):(length(theta)-1)])/(1+exp(theta[(length(theta)-2):(length(theta)-1)])) # redo the variable change for phi
theta[length(theta)] <- exp(theta[length(theta)]) # redo variable change for sigma

VC <- GMMvarcov(theta) # variance-covariance matrix (optimal)
stheta <- sqrt(diag(VC)) # standard errors

## list of variables names
varnames <- c("Syrian", "Male", "Parents have diverse friends$^\\dagger$", "Parents live in diverse Neighborhood$^\\dagger$",
              "Arrived in 2012 (Syrians)", "Arrived in 2013 (Syrians)", "Arrived in 2014 (Syrians)", "Arrived in 2015 (Syrians)",
              "Arrived in 2016 (Syrians)", "Arrived in 2017 (Syrians)", "Arrived in 2018 (Syrians)", "Syrian region 3 (Syrians)",
              "Syrian region 4 (Syrians)", "Syrian region 5 (Syrians)", "Syrian region 6 (Syrians)", "Syrian region 7 (Syrians)",
              "Syrian region 8 (Syrians)", "Fluent in Turkish (Syrians)", "School 1", "School 2", "School 3",
              "School 4", "School 5", "School 6", "School 7", "School 8", "School 9", "School 10", "School 11", "Syrian-Turkish pair",
              "Turkish-Syrian pair", "Same Gender", "Parents live in diverse neighborhood$^\\dagger$","Syrian-Syrian pair","Fluent in Turkish (Syrians)", "Linked in wave 1", "School 1", "School 2", "School 3", "School 4", "School 5",
              "School 7", "School 8", "School 9", "School 10", "School 11", "$\\phi^S$", "$\\phi^T$", "$\\sigma$")

### print as a LaTeX table
for (i in 1:length(theta)){
  cat(paste(varnames[i], " & ", format(round(theta[i], 3), nsmall = 3), " & (", format(round(stheta[i], 3), nsmall = 3),") \\\\ \n",sep=""))
}


save.image("estimWave2withcontrols.RData")

