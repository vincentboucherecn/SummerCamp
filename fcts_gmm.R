
#############################################################################
#############################################################################
#############################################################################
#######################     Format data     #################################
#############################################################################
#############################################################################
#############################################################################

builddata <- function(){
  
  ###
  ### This function formats the original data into lists to be used by the estimators
  ###
  
  classes <- unique(dta$Class) # list of classes
  ## computes number of pairs
  X1 <- vector("list",length(classes))
  X2 <- vector("list",length(classes))
  G1 <- vector("list",length(classes))
  G2 <- vector("list",length(classes))
  
  for (c in 1:length(classes)){
    print(c)
    tdta1 <- dta[(dta$Class==classes[c] & dta$Wave==1),] # wave 1 school c
    tdta2 <- dta[(dta$Class==classes[c] & dta$Wave==2),] # wave 2 school c
    nt <- nrow(tdta1) # size of the class
    Gt2 <- Gt1 <- matrix(0,nt,nt)
    Xt1 <- Xt2 <- matrix(NA,nt,15)
    
    ### Build matrices of individual variables
    # wave 1
    Xt1[,1] <- tdta1$st_id # student id
    Xt1[,2] <- tdta1$sch_id # school id
    Xt1[,3] <- tdta1$numnat # nationality
    Xt1[,4] <- tdta1$Gender # gender
    Xt1[,5] <- tdta1$Siblings # number of siblings
    Xt1[,6] <- tdta1$Mother.ed # mother education
    Xt1[,7] <- tdta1$Father.ed # father education
    Xt1[,8] <- tdta1$memp # mother employment status
    Xt1[,9] <- tdta1$femp # father employment status
    Xt1[,10] <- tdta1$Parents..friends # parents' friends composition
    Xt1[,11] <- tdta1$Parents..neighbors # parents' neighborhood composition
    Xt1[,12] <- tdta1$ayear # arrival year
    Xt1[,13] <- tdta1$Region.in.Syria # region in Syria
    Xt1[,14] <- tdta1$Language # turkish skill
    Xt1[,15] <- tdta1$Sticker.left # number of stickers left
    
    ## idem for wave2
    Xt2[,1] <- tdta2$st_id
    Xt2[,2] <- tdta2$sch_id
    Xt2[,3] <- tdta2$numnat
    Xt2[,4] <- tdta2$Gender
    Xt2[,5] <- tdta2$Siblings
    Xt2[,6] <- tdta2$Mother.ed
    Xt2[,7] <- tdta2$Father.ed
    Xt2[,8] <- tdta2$memp
    Xt2[,9] <- tdta2$femp
    Xt2[,10] <- tdta2$Parents..friends
    Xt2[,11] <- tdta2$Parents..neighbors
    Xt2[,12] <- tdta2$ayear
    Xt2[,13] <- tdta2$Region.in.Syria
    Xt2[,14] <- tdta1$Language #### overwrite we do not use wave2 language skill
    Xt2[,15] <- tdta1$Sticker.left #### overwrite we do not use wave2 stickers
    
    #for all pairs, construct network structures
    for (i in 1:nt){
      for (j in 1:nt){
        if (j != i){
          gt1 <- gt2 <- 0
          for (k in c("F1","F2","F3","F4","F5")){
            if (tdta1[j,"Name"]==tdta1[i,k]){
              gt1 <- gt1 + 1 # i named j at friendship position k in wave 1
            }
            if (tdta2[j,"Name"]==tdta2[i,k]){
              gt2 <- gt2 + 1 # i named j at friendship position k in wave 2
            }
          }
          Gt1[i,j] <- gt1
          Gt2[i,j] <- gt2
        }
      }
    }
    G1[[c]] <- Gt1
    G2[[c]] <- Gt2
    X1[[c]] <- Xt1
    X2[[c]] <- Xt2
  }
  return(list(X1,X2,G1,G2))
}

buildZ <- function(){

  ###
  ### This function generates the list of pairwise characteristics
  ###

  Z <- vector("list",length(X))
  for (c in 1:length(X)){
    Zt <- vector("list",8)
    Xt <- X[[c]]
    nt <- nrow(Xt)
    Zt[[1]] <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # syrian-turkish
    Zt[[2]] <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # turkish-syrian
    Zt[[3]] <- matrix(as.numeric(matrix(Xt[,4],nt,nt)==t(matrix(Xt[,4],nt,nt))),nt,nt) # same gender
    divn <- as.numeric(Xt[,11] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of neighbours are other than ones type
    Zt[[4]] <- matrix(divn,nt,nt)*t(matrix(divn,nt,nt))*(Zt[[1]]+Zt[[2]]) # only TS and ST links
    Zt[[5]] <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # syrian-syrian
    #Zt[[6]] <- position 6 is not used
    Xt14 <- Xt[,14]
    Xt14[is.na(Xt14)] <- 0
    Zt[[7]] <- (Zt[[1]] + Zt[[2]])*(matrix(Xt14,nt,nt) + t(matrix(Xt14,nt,nt))) # language for syrian-turkish links
    Zt[[8]] <- G1[[c]] # linked in wave 1
    
    Z[[c]] <- Zt
  }
  return(Z)
}


buildpairs <- function(){
  ## this function builds a large dataset for all DIRECTED pairs in the sample for preliminary analysis
  
  classes <- unique(dta$Class) # list of classes
  ## computes number of pairs
  Npairs <- 0
  for (i in classes){
    Npairs <- Npairs + sum(as.numeric((dta$Class==i) & dta$Wave==1))*(sum(as.numeric((dta$Class==i) & dta$Wave==1))-1)
  }
  
  ## create empty dataset
  cnames <- c("idi","idj","gij1","gij2","school","class","nati","natj","genderi","genderj","langi1","langj1","langi2","langj2",
              "sibi","sibj","medi","medj","fedi","fedj","pfri","pfrj","pnei","pnej","plani","planj","sticki1","stickj1","sticki2","stickj2","ayeari","ayearj","qSyrians")
  dtapairs <- as.data.frame(matrix(NA,Npairs,length(cnames)))
  colnames(dtapairs) <- cnames
  
  
  pos <- 1 # current pair position
  for (c in 1:length(classes)){
    print(c)
    tdta1 <- dta[(dta$Class==classes[c] & dta$Wave==1),] # wave 1 school c
    tdta2 <- dta[(dta$Class==classes[c] & dta$Wave==2),] # wave 2 school c
    tdta1$qSyrians <- mean(as.numeric(tdta1$numnat==1))*100 # % of Syrians
    nt <- nrow(tdta1) # size of the class
    
    #for all pairs, fillin variables
    for (i in 1:nt){
      for (j in 1:nt){
        if (j != i){
          dtapairs[pos,"qSyrians"] <- tdta1[i,"qSyrians"]
          dtapairs[pos,"idi"] <- tdta1[i,"st_id"]
          dtapairs[pos,"idj"] <- tdta1[j,"st_id"]
          dtapairs[pos,"class"] <- tdta1[i,"Class"]
          dtapairs[pos,"school"] <- tdta1[i,"sch_id"]
          dtapairs[pos,"nati"] <- tdta1[i,"numnat"]
          dtapairs[pos,"natj"] <- tdta1[j,"numnat"]
          dtapairs[pos,"genderi"] <- tdta1[i,"Gender"]
          dtapairs[pos,"genderj"] <- tdta1[j,"Gender"]
          dtapairs[pos,"langi1"] <- tdta1[i,"Language"]
          dtapairs[pos,"langj1"] <- tdta1[j,"Language"]
          dtapairs[pos,"langi2"] <- tdta2[i,"Language"]
          dtapairs[pos,"langj2"] <- tdta2[j,"Language"]
          dtapairs[pos,"sibi"] <- tdta1[i,"Siblings"]
          dtapairs[pos,"sibj"] <- tdta1[j,"Siblings"]
          dtapairs[pos,"medi"] <- tdta1[i,"Mother.ed"]
          dtapairs[pos,"medj"] <- tdta1[j,"Mother.ed"]
          dtapairs[pos,"fedi"] <- tdta1[i,"Father.ed"]
          dtapairs[pos,"fedj"] <- tdta1[j,"Father.ed"]
          dtapairs[pos,"pfri"] <- tdta1[i,"Parents..friends"]
          dtapairs[pos,"pfrj"] <- tdta1[j,"Parents..friends"]
          dtapairs[pos,"pnei"] <- tdta1[i,"Parents..neighbors"]
          dtapairs[pos,"pnej"] <- tdta2[j,"Parents..neighbors"]
          dtapairs[pos,"plani"] <- tdta1[i,"Parental.language"]
          dtapairs[pos,"planj"] <- tdta1[j,"Parental.language"]
          dtapairs[pos,"ayeari"] <- tdta1[i,"ayear"]
          dtapairs[pos,"ayearj"] <- tdta1[j,"ayear"]
          dtapairs[pos,"regioni"] <- tdta1[i,"Region.in.Syria"]
          dtapairs[pos,"regionj"] <- tdta1[j,"Region.in.Syria"]
          dtapairs[pos,"sticki1"] <- tdta1[i,"Sticker.left"]
          dtapairs[pos,"stickj1"] <- tdta1[j,"Sticker.left"]
          dtapairs[pos,"sticki2"] <- tdta2[i,"Sticker.left"]
          dtapairs[pos,"stickj2"] <- tdta2[j,"Sticker.left"]
          
          ## compute networks structures
          gt1 <- gt2 <- 0
          for (k in c("F1","F2","F3","F4","F5")){
            if (tdta1[j,"Name"]==tdta1[i,k]){
              gt1 <- gt1 + 1
            }
            if (tdta2[j,"Name"]==tdta2[i,k]){
              gt2 <- gt2 + 1
            }
          }
          dtapairs[pos,"gij1"] <- gt1
          dtapairs[pos,"gij2"] <- gt2
          pos <- pos + 1
        }    
      }
    }
  }
  return(dtapairs)
}


#############################################################################
#############################################################################
#############################################################################
###########################     Model   #####################################
#############################################################################
#############################################################################
#############################################################################

typef <- function(Xt){

  ###
  ### This function gives the type of each student (ethnicity x gender)
  ###
  
  typeout <- rep(0,nrow(Xt))
  type1 <- Xt[,3]==1 & Xt[,4]==2
  type2 <- Xt[,3]==1 & Xt[,4]==1
  type3 <- Xt[,3]==2 & Xt[,4]==2
  type4 <- Xt[,3]==2 & Xt[,4]==1
  typeout[type1] <- 1
  typeout[type2] <- 2
  typeout[type3] <- 3
  typeout[type4] <- 4
  return(typeout)
}

Dmat <- function(gamma,c){

  ###
  ### This function computes the matrix of preference biases for class c, given gamma
  ###
  
  Zt <- Z[[c]]
  Xt <- X[[c]]
  nt <- nrow(Xt)
  D <- matrix(0,nt,nt)
  pos <- 1
  for (i in lstvars){
    D <- D + Zt[[i]]*gamma[pos] # each observed pair characteristic
    pos <- pos + 1
  }
  if (Xt[1,2]<6){ # school 6 dropped due to colinearity
    D <- D + gamma[(pos+Xt[1,2]-1)] # school dummy
  }
  if (Xt[1,2]>6){
    D <- D + gamma[(pos+Xt[1,2]-2)] # school dummy
  }
  D <- pnorm(D)
  diag(D) <- 0
  return(D)
}

bvec <- function(bet,c){
  
  ###
  ### This function computes the vector (b - epsilon) for class c, given a value for bet(a)
  ###
  
  Xt <- X[[c]]
  nt <- nrow(Xt)
  b <- matrix(as.numeric(Xt[,3]==1),nt,1)*bet[1] # syrian
  b <- b + matrix(as.numeric(Xt[,4]==2),nt,1)*bet[2] # male
  divf <- as.numeric(Xt[,10] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,10] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of friends are of the other nationality
  divn <- as.numeric(Xt[,11] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of neighbours are of the other nationality
  b <- b + matrix(divf,nt,1)*bet[3]
  b <- b + matrix(divn,nt,1)*bet[4]
  ay <- Xt[,12]
  ay[is.na(ay)] <- 0
  b <- b + as.numeric(ay==2012)*bet[5] # arrival year
  b <- b + as.numeric(ay==2013)*bet[6] # arrival year
  b <- b + as.numeric(ay==2014)*bet[7] # arrival year
  b <- b + as.numeric(ay==2015)*bet[8] # arrival year
  b <- b + as.numeric(ay==2016)*bet[9] # arrival year
  b <- b + as.numeric(ay==2017)*bet[10] # arrival year
  b <- b + as.numeric(ay==2018)*bet[11] # arrival year
  rs <- Xt[,13]
  rs[is.na(rs)] <- 0
  b <- b + as.numeric(rs==3)*bet[12] # region syria
  b <- b + as.numeric(rs==4)*bet[13] # region syria
  b <- b + as.numeric(rs==5)*bet[14] # region syria
  b <- b + as.numeric(rs==6)*bet[15] # region syria
  b <- b + as.numeric(rs==7)*bet[16] # region syria
  b <- b + as.numeric(rs==8)*bet[17] # region syria
  
  l <- Xt[,14]
  l[is.na(l)] <- 0
  b <- b + l*bet[18] # language (if Syrian)
  
  b <- b + bet[(18+Xt[1,2])]
  return(b)
}



Seq <- function(bet,gamma,phi,sig,c,sim){
  
  ###
  ### This function computes the equilibrium value for s in class c, given all of the parameters and the error for simulation sim
  ###
  
  Xt <- X[[c]]
  Et <- E[[c]][,sim]
  nt <- nrow(Xt)
  b <- bvec(bet,c)
  syr <- matrix(as.numeric(Xt[,3]==1),nt,1)
  b <- b + Et*sig # random shock
  b[c(syr)==1] <- b[c(syr)==1]*(1-phi[1])
  b[c(syr)==0] <- b[c(syr)==0]*(1-phi[2])
  D <- matrix(1/(nt-1),nt,nt)
  diag(D) <- 0
  S0 <- matrix(0,nt,1) # initial socialization guess
  sc <- 0 # stoping flag
  while (sc==0){
    S1 <- b + (phi[1]*syr + phi[2]*(1-syr))*(D%*%S0) # best response (unconstrained)
    S1 <- pmax(0,pmin(1,S1)) # best response (constrained)
    if (sum(abs(S1-S0))<1e-6){ # convergence criterion
      sc <- 1
    }
    S0 <- S1
  }
  return(S1)
}

Pmat <- function(bet,gamma,phi,sig,c,sim){
  
  ###
  ### This function computes the probability matrix P
  ###
  
  D <- Dmat(gamma,c)
  st <- Seq(bet,gamma,phi,sig,c,sim) # compute equilibrium S without assuming interior solution
  D <- D*(st%*%t(st))
  Xt <- X[[c]]
  nt <- nrow(Xt)
  typel <- typef(Xt)
  for (j in 1:nt){
    rowi <- as.numeric(typel==typel[j])
    Gtype <- t(matrix(rep(rowi,nt),nt,nt))
    Gtype2 <- matrix(as.numeric(G1[[c]]==G1[[c]][,j]),nt,nt)
    Gtype <- Gtype*Gtype2
    tvar <- D[,j]/(Gtype%*%st)
    tvar[is.nan(tvar)] <- 0
    D[,j] <- tvar
  }
  
  return(D)
}


#############################################################################
#############################################################################
#############################################################################
#########################     Estimator   ###################################
#############################################################################
#############################################################################
#############################################################################


M1i <- function(bet,gamma,phi,sig){

  ###
  ### This function computes the matrix of simulated unconditional moments for the network
  ###
  
  moments <- matrix(0,n1,length(gamma))
  for (sim in 1:nsim){ # monte carlo integration
    pos <- 1
    for (c in 1:length(X)){ # for all classes
      nt <- nrow(X[[c]])
      P <- Pmat(bet,gamma,phi,sig,c,sim) # linking probabilities
      Gt <- G[[c]] # network
      for (i in 1:nt){ #for all pairs
        for (j in 1:nt){
          if (i != j){
            d <- (Gt[i,j]-P[i,j]) # errors
            z <- rep(0,length(gamma)) # initialize vector of zij
            pos2 <- 1
            for (k in lstvars){ # explanatory variables
              z[pos2] <- Z[[c]][[k]][i,j]
              pos2 <- pos2 + 1
            }
            if (X[[c]][1,2]<6){ # school 6 dropped due to colinearity
              z[(pos2-1+X[[c]][1,2])] <- 1 # school dummy
            }
            if (X[[c]][1,2]>6){
              z[(pos2-2+X[[c]][1,2])] <- 1 # school dummy
            }
            moments[pos,] <- moments[pos,] + d*z/nsim # moment condition
            pos <- pos + 1
          }
        }
        
      }
    }
  }
  return(moments)
}


M2i <- function(bet,gamma,phi){
  
  ###
  ### This function computes the matrix of moments for the outcome
  ###
  
  moments <- matrix(0,n2,(29+ 2*18))
  for (c in 1:length(X)){ # for all classes
    Xt <- X[[c]]
    nt <- nrow(Xt)
    # position in the moment matrix
    if (c == 1){
      p1 <- 1
      p2 <- nt
    } else {
      p1 <- p2 + 1
      p2 <- p1 + nt -1
    }
    D <- matrix(1/(nt-1),nt,nt)
    diag(D) <- 0
    syr <- matrix(as.numeric(Xt[,3]==1),nt,1)
    b <- bvec(bet,c) # b vector
    b[c(syr)==1] <- b[c(syr)==1]*(1-phi[1])
    b[c(syr)==0] <- b[c(syr)==0]*(1-phi[2])
    e <- Xt[,15] - b - phi[1]*(D%*%Xt[,15])*matrix(as.numeric(Xt[,3]==1),nt,1) - phi[2]*(D%*%Xt[,15])*matrix(as.numeric(Xt[,3]==2),nt,1) # errors
    bz <- matrix(0,nt,(29+ 2*18)) # initialize matrix of moments
    bz[,1] <- as.numeric(Xt[,3]==1) # syrian
    bz[,2] <- as.numeric(Xt[,4]==2) # male
    bz[,3] <- as.numeric(Xt[,10] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,10] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of friends are other than ones type
    bz[,4] <- as.numeric(Xt[,11] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of neighbours are other than ones type
    
    ay <- Xt[,12]
    ay[is.na(ay)] <- 0
    bz[,5] <- as.numeric(ay==2012) # arrival year
    bz[,6] <- as.numeric(ay==2013) # arrival year
    bz[,7] <- as.numeric(ay==2014) # arrival year
    bz[,8] <- as.numeric(ay==2015) # arrival year
    bz[,9] <- as.numeric(ay==2016) # arrival year
    bz[,10] <- as.numeric(ay==2017) # arrival year
    bz[,11] <- as.numeric(ay==2018) # arrival year
    rs <- Xt[,13]
    rs[is.na(rs)] <- 0
    bz[,12] <- as.numeric(rs==3) # region syria
    bz[,13] <- as.numeric(rs==4) # region syria
    bz[,14] <- as.numeric(rs==5) # region syria
    bz[,15] <- as.numeric(rs==6) # region syria
    bz[,16] <- as.numeric(rs==7) # region syria
    bz[,17] <- as.numeric(rs==8) # region syria
    
    l <- Xt[,14]
    l[is.na(l)] <- 0
    bz[,18] <- l # language (if Syrian)
    
    bz[,(18+Xt[1,2])] <- 1 # school dummy
    Smask <- matrix(rep(as.numeric(Xt[,3]==1),18),nt,18)
    bz[,30:(30+18-1)] <- (D%*%bz[,1:18])*Smask
    bz[,(30+18):ncol(bz)] <- (D%*%bz[,1:18])*(1-Smask)
    mmat <- matrix(rep(c(e),ncol(bz)),nt,ncol(bz))
    mmat <- mmat*bz
    moments[p1:p2,] <- mmat
  }
  return(moments)
}

GMM <- function(theta){
  
  ###
  ### This function computes the estimator's objective function
  ###
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  t3 <- exp(t3)/(1+exp(t3))
  t4 <- exp(t4)
  
  MM1 <- M1i(t1,t2,t3,t4) # network moments
  MM2 <- M2i(t1,t2,t3) # socialization moments
  Q1 <- sum(colMeans(MM1)^2) # unweighted GMM
  Q2 <- sum(colMeans(MM2)^2) # unewighted GMM
  obj <- Q1 + Q2
  print(obj)
  return(obj)
}

WGMM <- function(theta){
  
  ###
  ### This function computes the estimator's objective function
  ###
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  t3 <- exp(t3)/(1+exp(t3))
  t4 <- exp(t4)
  
  MM1 <- GMMjacob1(theta) # network moments
  MM2 <- GMMjacob2(theta) # socialization moments
  Q1 <- c(t(MM1)%*%W1%*%MM1)
  Q2 <- c(t(MM2)%*%W2%*%MM2)
  obj <- Q1 + Q2
  print(obj)
  return(obj)
}


GMMjacob1 <- function(theta){

  ###
  ### Function wrapper for the expected moments for the network
  ###
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  MM <- M1i(t1,t2,t3,t4)
  return(c(colMeans(MM)))
}

GMMjacob2 <- function(theta){
  
  ###
  ### Function wrapper for the expected moments for the outcome
  ###
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  MM <- M2i(t1,t2,t3)
  return(c(colMeans(MM)))
}

compW <- function(theta){
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  t3 <- exp(t3)/(1+exp(t3))
  t4 <- exp(t4)
  
  MM1 <- M1i(t1,t2,t3,t4)
  MM2 <- M2i(t1,t2,t3)
  H1 <- matrix(0,ncol(MM1),ncol(MM1))
  H2 <- matrix(0,ncol(MM2),ncol(MM2))
  
  for (i in 1:nrow(MM1)){
    H1 <- H1 + matrix(MM1[i,],ncol(MM1),1)%*%matrix(MM1[i,],1,ncol(MM1))/n1 # outer-product Q1
  }
  print("H1")
  for (i in 1:nrow(MM2)){
    H2 <- H2 + matrix(MM2[i,],ncol(MM2),1)%*%matrix(MM2[i,],1,ncol(MM2))/n2 # outer-product Q2
  }
  return(list(solve(H1),solve(H2)))
}


GMMvarcov <- function(theta){
  
  ###
  ### This function computes the variance-covariance matrix of the estimator
  ###
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  
  GG1 <- jacobian(GMMjacob1,theta)
  print("Jacobian 1")
  GG2 <- jacobian(GMMjacob2,theta)
  print("Jacobian 2")
  A1 <- t(GG1)%*%W1%*%GG1 # "hessian" Q1
  A2 <- t(GG2)%*%W2%*%GG2 # "hessian" Q2
  
  MM1 <- M1i(t1,t2,t3,t4)
  MM2 <- M2i(t1,t2,t3)
  H1 <- matrix(0,ncol(MM1),ncol(MM1))
  H2 <- matrix(0,ncol(MM2),ncol(MM2))
  
  for (i in 1:nrow(MM1)){
    H1 <- H1 + matrix(MM1[i,],ncol(MM1),1)%*%matrix(MM1[i,],1,ncol(MM1))/n1 # outer-product Q1
  }
  print("H1")
  for (i in 1:nrow(MM2)){
    H2 <- H2 + matrix(MM2[i,],ncol(MM2),1)%*%matrix(MM2[i,],1,ncol(MM2))/n2 # outer-product Q2
  }
  B1 <- t(GG1)%*%W1%*%H1%*%W1%*%GG1
  B2 <- t(GG2)%*%W2%*%H2%*%W2%*%GG2
  print("H2")
  V <- chol2inv(n1*A1 + n2*A2) # the solve function produces produces numerical errors
  VC <- V%*%( n1*B1 + n2*B2 )%*%V # variance-covariance matrix for theta
  return(VC)
}

#############################################################################
#############################################################################
#############################################################################
###########################   Post estimation   #############################
#############################################################################
#############################################################################
#############################################################################

marginal_bias <- function(theta){
  
  ###
  ### This function computes the marginal effects on the preference bias delta
  ###
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  
  mbias <- 0 # initialize mean bias
  fmgamma <- rep(0,length(t2)) # for all regressors
  
  for (c in 1:length(X)){ # for all classrooms
    D <- Dmat(t2,c) # get D
    D2 <- pnorm(qnorm(D)) # first derivative of the normal cdf=normal pdf
    diag(D) <- NA # ignores diagnonal elements
    diag(D2) <- NA # ignores diagnonal elements
    mbias <- mbias + sum(D,na.rm=T)/n1 # average bias
    
    for (par in 1:length(t2)){ # forall binary variables
      if (par<=length(lstvars)){ # binary variables that are not school dummies
        
        ## =1
        Zt <- Z[[c]] # get pairwise vars
        Zt[[lstvars[par]]] <- 1 # set current variable to 1
        Xt <- X[[c]] # get indiv variables
        nt <- nrow(Xt) # classroom size
        
        ## compute preference biasses under current variable =1
        D <- matrix(0,nt,nt)
        pos <- 1
        for (i in lstvars){
          D <- D + Zt[[i]]*t2[pos] # each observed pair characteristic
          pos <- pos + 1
        }
        if (Xt[1,2]<6){ # school 6 dropped due to colinearity
          D <- D + t2[(pos-1+Xt[1,2])] # school dummy
        }
        if (Xt[1,2]>6){
          D <- D + t2[(pos-2+Xt[1,2])] # school dummy
        }
        D <- pnorm(D)
        diag(D) <- 0
        
        D1  <- D # keep matrix
        
        ## =0
        Zt <- Z[[c]] # get pairwise vars
        Zt[[lstvars[par]]] <- 0 # set current var to 0
        Xt <- X[[c]] # get indiv. vars
        nt <- nrow(Xt) # classroom size
        
        ## compute preference biasses under current variable =0
        D <- matrix(0,nt,nt)
        pos <- 1
        for (i in lstvars){
          D <- D + Zt[[i]]*t2[pos] # each observed pair characteristic
          pos <- pos + 1
        }
        if (Xt[1,2]<6){ # school 6 dropped due to colinearity
          D <- D + t2[(pos-1+Xt[1,2])] # school dummy
        }
        if (Xt[1,2]>6){
          D <- D + t2[(pos-2+Xt[1,2])] # school dummy
        }
        D <- pnorm(D)
        diag(D) <- 0
        
        D0  <- D
      }
      
      if (par>length(lstvars)){ # for school dummies
        ## =1
        Zt <- Z[[c]] # get pairwise vars
        Xt <- X[[c]] # get indiv vars
        nt <- nrow(Xt) # clasroom size
        
        ## compute bias matrix
        D <- matrix(0,nt,nt)
        pos <- 1
        for (i in lstvars){
          D <- D + Zt[[i]]*t2[pos] # each observed pair characteristic
          pos <- pos + 1
        }
        ## if school is not current school, add current school dummy
        if (Xt[1,2]<6){ # school 6 dropped due to colinearity
          if (par != (pos-1+Xt[1,2])){
            D <- D + t2[par]
          }
          D <- D + t2[(pos-1+Xt[1,2])] # school dummy
        }
        if (Xt[1,2]>6){
          if (par != (pos-2+Xt[1,2])){
            D <- D + t2[par]
          }
          D <- D + t2[(pos-2+Xt[1,2])] # school dummy
        }
        D <- pnorm(D)
        diag(D) <- 0
        
        D1  <- D # store matrix
        
        ## =0
        Zt <- Z[[c]] # get pairwise vars
        Xt <- X[[c]] # get indiv vars
        nt <- nrow(Xt) # classroom size
        
        ## compute bias matrix
        D <- matrix(0,nt,nt)
        pos <- 1
        for (i in lstvars){
          D <- D + Zt[[i]]*t2[pos] # each observed pair characteristic
          pos <- pos + 1
        }
        ## if school is not current school, add current school dummy
        if (Xt[1,2]<6){
          if (par != (pos-1+Xt[1,2])){
            D <- D + t2[(pos-1+Xt[1,2])]
          }
        }
        if (Xt[1,2]>6){
          if (par != (pos-2+Xt[1,2])){
            D <- D + t2[(pos-2+Xt[1,2])]
          }
        }
        D <- pnorm(D)
        diag(D) <- 0
        
        D0  <- D # keep matrix
      }
      deltaD <- D1-D0 # marginal variation in preference bias
      diag(deltaD) <- NA # ignore diagonal elements
      fmgamma[par] <- fmgamma[par] + sum(deltaD,na.rm=T)/n1 # average marginal effect
    }
  }

  return(list(mbias,fmgamma)) # return average bias and average marginal effects
}

shufflepeople <- function(mix){
  
  ###
  ### shuffles students: create list of students (type, class0, position0, class1, position1)
  ### 0 is original position and 1 is shuffled position. mix is level of mixing. 0.5 is random, 1 is (almost) perfect segregation
  ###
  
  ## 
  liststd <- matrix(NA,n2,5)
  csize <- rep(NA,length(X))
  pos <- 1
  for (c in 1:length(X)){
    Xt <- X[[c]]
    nt <- nrow(Xt)
    csize[c] <- nt
    for (i in 1:nt){
      liststd[pos,1] <- Xt[i,3] # type
      liststd[pos,2] <- c # class
      liststd[pos,3] <- i # position in class
      pos <- pos + 1
    }
  }
  #class types
  nsyrians <- sum(as.numeric(liststd[,1]==1))
  fillsyrians <- 0
  ctype <- rep(1,length(X)) # all classes are "turkish"
  nclist <- sample(1:length(X),length(X),replace=F)
  for (c in nclist){
    if (fillsyrians<nsyrians){
      ctype[c] <- 0
      fillsyrians <- fillsyrians + csize[c]
    }
  }
  cleft <- csize
  nilist <- sample(1:n2,n2,replace=F)
  for (i in nilist){ # shuffle order of individuals
    if (liststd[i,1]==1){ # if syrian
      proba <- ((1-ctype)*mix + (1-mix)*ctype)*as.numeric(cleft>0)
    }
    if (liststd[i,1]==2){ # if turkish
      proba <- (ctype*mix + (1-mix)*(1-ctype))*as.numeric(cleft>0)
    }
    proba <- proba/sum(proba)
    wclass <- sample(1:length(X),1,replace=F,proba)
    liststd[i,4] <- wclass
    liststd[i,5] <- cleft[wclass]
    cleft[wclass] <- cleft[wclass] - 1
  }
  return(liststd)
}

swappeople <- function(listin){
  Xnew <- X
  for (i in 1:n2){
    Xnew[[listin[i,4]]][listin[i,5],] <- X[[listin[i,2]]][listin[i,3],]
  }
  return(Xnew)
}

swaperrors <- function(listin){
  Enew <- E
  for (i in 1:n2){
    Enew[[listin[i,4]]][listin[i,5],] <- E[[listin[i,2]]][listin[i,3],]
  }
  return(Enew)
}

sametype <- function(theta){
  
  ###
  ### This function computes equilibrium quantities (for graphs in counterfactual analysis)
  ###
  
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  
  
  q <- rep(NA,length(X)) # share of syrians
  shom <- rep(NA,length(X)) # homophily: syrians
  sihom <- rep(NA,length(X)) # inbreeding homophily: syrians
  thom <- rep(NA,length(X)) # homophily: turkish
  tihom <- rep(NA,length(X)) # inbreeding homophily: turkish
  pss <- rep(NA,length(X)) # average pSS
  ptt <- rep(NA,length(X)) # average pTT
  pst <- rep(NA,length(X)) # average pST
  pts <- rep(NA,length(X)) # average pTS
  st <- rep(NA,length(X)) # average st
  ss <- rep(NA,length(X)) # average ss
  
  for (c in 1:length(X)){ # for all classrooms
    D <- Dmat(t2,c) # bias matrix
    Xt <- X[[c]] # indiv. vars
    nt <- nrow(Xt)
    q[c] <- mean(as.numeric(Xt[,3]==1))*100 # share of syrians
    St <- Seq(t1,t2,t3,t4,c,sample(1:nsim,1)) # simulate socialization levels
    P <- D*(St%*%t(St)) # probability of linking
    ss[c] <- mean(St[Xt[,3]==1])
    st[c] <- mean(St[Xt[,3]==2])
    typel <- typef(Xt)
    for (j in 1:nt){
      rowi <- as.numeric(typel==typel[j])
      Gtype <- t(matrix(rep(rowi,nt),nt,nt))
      oldG <- matrix(rbinom(nt*nt,1,0.15),nt,nt)
      Gtype2 <- matrix(as.numeric(oldG==oldG[,j]),nt,nt)
      Gtype <- Gtype*Gtype2
      P[,j] <- P[,j]/(Gtype%*%St)
    }
    
    nt <- nrow(Xt)
    iST <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # syrian-turkish indicator
    diag(iST) <- 0
    iTS <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # turkish-syrian indicator
    diag(iTS) <- 0
    iSS <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # syrian-syrian indicator
    diag(iSS) <- 0
    iTT <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # turkish-turkish indicator
    diag(iTT) <- 0
    ## draw network using P
    Gt <- matrix(runif((nt*nt)),nt,nt)
    Gt <- matrix(as.numeric(P>=Gt),nt,nt)
    diag(Gt) <- 0
    if (sum(as.numeric(Xt[,3]==1))>1 & sum(as.numeric(Xt[,3]==1))<(nt-1)){
      shom[c] <- sum(Gt*iSS)/max(sum(Gt*iSS) + sum(Gt*iST) ,1) # homophily index Syrians
      sihom[c] <- (shom[c]-q[c]/100)/(1-q[c]/100)
      thom[c] <- sum(Gt*iTT)/max(sum(Gt*iTT) + sum(Gt*iTS) ,1) # homophily index Syrians
      tihom[c] <- (thom[c]-(1-q[c]/100))/(q[c]/100)
      pss[c] <- mean(P[iSS==1]) # average pSS
      ptt[c] <- mean(P[iTT==1]) # average pTT
      pst[c] <- mean(P[iST==1]) # average pST
      pts[c] <- mean(P[iTS==1]) # average pTS
    }
  }
  return(list(q,shom,thom,pss,ptt,pst,pts,sihom,tihom,ss,st))
}

sametype_no_congestion <- function(theta){
  
  ###
  ### This function computes equilibrium quantities without congestion (for graphs in counterfactual analysis)
  ###
  
  
  t1 <- theta[1:length(bethat)] # beta
  t2 <- theta[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- theta[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- theta[(length(bethat)+length(gamhat)+3)] # sigma
  
  
  q <- rep(NA,length(X)) # share of syrians
  shom <- rep(NA,length(X)) # homophily: syrians
  sihom <- rep(NA,length(X)) # inbreeding homophily: syrians
  thom <- rep(NA,length(X)) # homophily: turkish
  tihom <- rep(NA,length(X)) # inbreeding homophily: turkish
  pss <- rep(NA,length(X)) # average pSS
  ptt <- rep(NA,length(X)) # average pTT
  pst <- rep(NA,length(X)) # average pST
  pts <- rep(NA,length(X)) # average pTS
  st <- rep(NA,length(X)) # average st
  ss <- rep(NA,length(X)) # average ss
  for (c in 1:length(X)){ # for all classrooms
    D <- Dmat(t2,c) # bias matrix
    Xt <- X[[c]] # indiv. vars
    nt <- nrow(Xt)
    q[c] <- mean(as.numeric(Xt[,3]==1))*100 # share of syrians
    St <- Seq(t1,t2,t3,t4,c,sample(1:nsim,1)) # simulate socialization levels
    P <- D*(St%*%matrix(1,1,nt)) # probability of linking
    ss[c] <- mean(St[Xt[,3]==1])
    st[c] <- mean(St[Xt[,3]==2])
    nt <- nrow(Xt)
    iST <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # syrian-turkish indicator
    diag(iST) <- 0
    iTS <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # turkish-syrian indicator
    diag(iTS) <- 0
    iSS <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # syrian-syrian indicator
    diag(iSS) <- 0
    iTT <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # turkish-turkish indicator
    diag(iTT) <- 0
    ## draw network using P
    Gt <- matrix(runif((nt*nt)),nt,nt)
    Gt <- matrix(as.numeric(P>=Gt),nt,nt)
    diag(Gt) <- 0
    btemp <- bvec(t1,c)
    if (sum(as.numeric(Xt[,3]==1))>1 & sum(as.numeric(Xt[,3]==1))<(nt-1)){
      shom[c] <- sum(Gt*iSS)/max(sum(Gt*iSS) + sum(Gt*iST) ,1) # homophily index Syrians
      sihom[c] <- (shom[c]-q[c]/100)/(1-q[c]/100)
      thom[c] <- sum(Gt*iTT)/max(sum(Gt*iTT) + sum(Gt*iTS) ,1) # homophily index Syrians
      tihom[c] <- (thom[c]-(1-q[c]/100))/(q[c]/100)
      pss[c] <- mean(P[iSS==1]) # average pSS
      ptt[c] <- mean(P[iTT==1]) # average pTT
      pst[c] <- mean(P[iST==1]) # average pST
      pts[c] <- mean(P[iTS==1]) # average pTS
    }
  }
  return(list(q,shom,thom,pss,ptt,pst,pts,sihom,tihom,ss,st))
}

graphdata <- function(){
  
  ###
  ### This function computes equilibrium quantities on DATA (not for counterfactual analysis)
  ###
  
  
  q <- rep(NA,length(X)) # share of syrians
  shom <- rep(NA,length(X)) # homophily: syrians
  sihom <- rep(NA,length(X)) # inbreeding homophily: syrians
  thom <- rep(NA,length(X)) # homophily: turkish
  tihom <- rep(NA,length(X)) # inbreeding homophily: turkish

  for (c in 1:length(X)){ # for all classrooms
    Xt <- X[[c]] # indiv. vars
    nt <- nrow(Xt)
    q[c] <- mean(as.numeric(Xt[,3]==1))*100 # share of syrians
    St <- Xt[,15] # simulate socialization levels

    nt <- nrow(Xt)
    iST <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # syrian-turkish indicator
    diag(iST) <- 0
    iTS <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # turkish-syrian indicator
    diag(iTS) <- 0
    iSS <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # syrian-syrian indicator
    diag(iSS) <- 0
    iTT <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # turkish-turkish indicator
    diag(iTT) <- 0
    ## draw network using P
    Gt <- G[[c]]
    if (sum(as.numeric(Xt[,3]==1))>1 & sum(as.numeric(Xt[,3]==1))<(nt-1)){
      shom[c] <- sum(Gt*iSS)/max(sum(Gt*iSS) + sum(Gt*iST) ,1) # homophily index Syrians
      sihom[c] <- (shom[c]-q[c]/100)/(1-q[c]/100)
      thom[c] <- sum(Gt*iTT)/max(sum(Gt*iTT) + sum(Gt*iTS) ,1) # homophily index Syrians
      tihom[c] <- (thom[c]-(1-q[c]/100))/(q[c]/100)
    }
  }
  return(list(q,shom,thom,sihom,tihom))
}

gendata <- function(thetat,sim){

  ###
  ### This function creates datesets of simulated data (for model fit)
  ###
    
  t1 <- thetat[1:length(bethat)] # beta
  t2 <- thetat[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- thetat[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- thetat[(length(bethat)+length(gamhat)+3)] # sigma
  
  indreg <- as.data.frame(matrix(NA,n2,11)) # initialize results dataframe
  colnames(indreg) <- c("s","syrian","female","school","q","friends","neighbours","ayear","region","lang","class") # name columns
  pairreg <- as.data.frame(matrix(NA,n1,9)) # initialize results dataframe
  colnames(pairreg) <- c("g","type","sgender","parents","skill","school","wave1","q","class") # name columns
  pos1 <- 1
  pos2 <- 1
  for (c in 1:length(X)){
    Xt <- X[[c]]
    nt <- nrow(Xt)
    Zt <-Z[[c]]
    if (sim==1){
      St <- Seq(t1,t2,t3,t4,c,sample(1:nsim,1)) # simulate socialization levels
      St <- pmax(St,1e-8)
    } else {
      St <- Xt[,15]
    }
    qt <- mean(as.numeric(Xt[,3]==1))*100
    if (sim==1){
      D <- Dmat(t2,c) # bias matrix
      P <- D*(St%*%t(St)) # probability of linking
      typel <- typef(Xt)
      for (j in 1:nt){
        rowi <- as.numeric(typel==typel[j])
        Gtype <- t(matrix(rep(rowi,nt),nt,nt))
        oldG <- G1[[c]]
        Gtype2 <- matrix(as.numeric(oldG==oldG[,j]),nt,nt)
        Gtype <- Gtype*Gtype2
        P[,j] <- P[,j]/(Gtype%*%St)
      }
      Gt <- matrix(as.numeric(P>=matrix(runif(nt*nt),nt,nt)),nt,nt)
    } else {
      Gt <- G[[c]]
    }
    
    diag(Gt) <- NA
    
    indreg[pos1:(pos1+nt-1),"s"] <- St
    indreg[pos1:(pos1+nt-1),"syrian"] <- as.numeric(Xt[,3]==1)
    indreg[pos1:(pos1+nt-1),"female"] <- as.numeric(Xt[,4]==1)
    indreg[pos1:(pos1+nt-1),"q"] <- qt
    indreg[pos1:(pos1+nt-1),"school"] <- rep(Xt[1,2],nt)
    indreg[pos1:(pos1+nt-1),"friends"] <- as.numeric(Xt[,10] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,10] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of friends are other than ones type
    indreg[pos1:(pos1+nt-1),"neighbours"] <- as.numeric(Xt[,11] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of neighbours are other than ones type
    indreg[pos1:(pos1+nt-1),"ayear"] <- Xt[,12]
    indreg[pos1:(pos1+nt-1),"region"] <- Xt[,13]
    indreg[pos1:(pos1+nt-1),"lang"] <- Xt[,14]
    indreg[pos1:(pos1+nt-1),"class"] <- c
    pos1 <- pos1 + nt
    
    tvec <- c(Gt)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"g"] <- tvec
    
    tmat1 <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # SS
    tmat2 <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # ST
    tmat3 <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # TT
    tmat4 <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # TS
    tmat <- tmat1 + 2*tmat2 + 3*tmat3 + 4*tmat4
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"type"] <- tvec
    
    tmat <- matrix(as.numeric(matrix(Xt[,4],nt,nt)==t(matrix(Xt[,4],nt,nt))),nt,nt) #same gender
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"sgender"] <- tvec
    
    tmat <- G1[[c]]
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"wave1"] <- tvec
    
    
    divn <- as.numeric(Xt[,11] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of neighbours are other than ones type
    tmat <- matrix(divn,nt,nt)*t(matrix(divn,nt,nt))*(tmat2+tmat4) # only TS and ST links
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"parents"] <- tvec
    
    Xt14 <- Xt[,14]
    Xt14[is.na(Xt14)] <- 0
    tmat <- (tmat2 + tmat4)*(matrix(Xt14,nt,nt) + t(matrix(Xt14,nt,nt))) # language for syrian-turkish links
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"skill"] <- tvec
    
    
    pairreg[pos2:(pos2+nt*(nt-1)-1),"q"] <- qt
    pairreg[pos2:(pos2+nt*(nt-1)-1),"school"] <- rep(Xt[1,2],(nt*(nt-1)))
    pairreg[pos2:(pos2+nt*(nt-1)-1),"class"] <- c
    
    pos2 <- pos2 + nt*(nt-1)
  }
  return(list(indreg,pairreg))
}

gendata_highskill <- function(thetat,sim){
  
  ###
  ### This function creates datesets of simulated data (for model fit)
  ###
  
  t1 <- thetat[1:length(bethat)] # beta
  t2 <- thetat[(length(bethat)+1):(length(bethat)+length(gamhat))] # gamma
  t3 <- thetat[(length(bethat)+length(gamhat)+1):(length(bethat)+length(gamhat)+2)] # phi
  t4 <- thetat[(length(bethat)+length(gamhat)+3)] # sigma
  
  indreg <- as.data.frame(matrix(NA,n2,12)) # initialize results dataframe
  colnames(indreg) <- c("s","syrian","female","school","q","friends","neighbours","ayear","region","lang","class") # name columns
  pairreg <- as.data.frame(matrix(NA,n1,7)) # initialize results dataframe
  colnames(pairreg) <- c("g","type","sgender","school","wave1","q","class") # name columns
  pos1 <- 1
  pos2 <- 1
  for (c in 1:length(X)){
    Xt <- X[[c]]
    Xt[Xt[,14]==0,14] <- 1 # all Syrian kids are fluent.
    nt <- nrow(Xt)
    Zt <-Z[[c]]
    if (sim==1){
      St <- SeqHS(t1,t2,t3,t4,c,sample(1:nsim,1),Xt) # simulate socialization levels
      St <- pmax(St,1e-8)
    } else {
      St <- Xt[,15]
    }
    qt <- mean(as.numeric(Xt[,3]==1))*100
    if (sim==1){
      D <- Dmat(t2,c) # bias matrix
      P <- D*(St%*%t(St)) # probability of linking
      typel <- typef(Xt)
      for (j in 1:nt){
        rowi <- as.numeric(typel==typel[j])
        Gtype <- t(matrix(rep(rowi,nt),nt,nt))
        oldG <- G1[[c]]
        Gtype2 <- matrix(as.numeric(oldG==oldG[,j]),nt,nt)
        Gtype <- Gtype*Gtype2
        P[,j] <- P[,j]/(Gtype%*%St)
      }
      Gt <- matrix(as.numeric(P>=matrix(runif(nt*nt),nt,nt)),nt,nt)
    } else {
      Gt <- G[[c]]
    }
    
    diag(Gt) <- NA
    
    indreg[pos1:(pos1+nt-1),"s"] <- St
    indreg[pos1:(pos1+nt-1),"syrian"] <- as.numeric(Xt[,3]==1)
    indreg[pos1:(pos1+nt-1),"female"] <- as.numeric(Xt[,4]==1)
    indreg[pos1:(pos1+nt-1),"q"] <- qt
    indreg[pos1:(pos1+nt-1),"school"] <- rep(Xt[1,2],nt)
    indreg[pos1:(pos1+nt-1),"friends"] <- as.numeric(Xt[,10] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,10] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of friends are other than ones type
    indreg[pos1:(pos1+nt-1),"neighbours"] <- as.numeric(Xt[,11] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of neighbours are other than ones type
    indreg[pos1:(pos1+nt-1),"ayear"] <- Xt[,12]
    indreg[pos1:(pos1+nt-1),"region"] <- Xt[,13]
    indreg[pos1:(pos1+nt-1),"lang"] <- Xt[,14]
    indreg[pos1:(pos1+nt-1),"class"] <- c
    pos1 <- pos1 + nt
    
    tvec <- c(Gt)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"g"] <- tvec
    
    tmat1 <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # SS
    tmat2 <- matrix(as.numeric(Xt[,3]==1),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # ST
    tmat3 <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==2),nt,nt)) # TT
    tmat4 <- matrix(as.numeric(Xt[,3]==2),nt,nt)*t(matrix(as.numeric(Xt[,3]==1),nt,nt)) # TS
    tmat <- tmat1 + 2*tmat2 + 3*tmat3 + 4*tmat4
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"type"] <- tvec
    
    tmat <- matrix(as.numeric(matrix(Xt[,4],nt,nt)==t(matrix(Xt[,4],nt,nt))),nt,nt) #same gender
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"sgender"] <- tvec
    
    tmat <- G1[[c]]
    diag(tmat) <- NA
    tvec <- c(tmat)
    tvec <- tvec[is.na(tvec)==F]
    pairreg[pos2:(pos2+nt*(nt-1)-1),"wave1"] <- tvec
    
    pairreg[pos2:(pos2+nt*(nt-1)-1),"q"] <- qt
    pairreg[pos2:(pos2+nt*(nt-1)-1),"school"] <- rep(Xt[1,2],(nt*(nt-1)))
    pairreg[pos2:(pos2+nt*(nt-1)-1),"class"] <- c
    
    pos2 <- pos2 + nt*(nt-1)
  }
  return(list(indreg,pairreg))
}


#############################################################################
#############################################################################
#############################################################################
#####################     Auxiliary function   ##############################
#############################################################################
#############################################################################
#############################################################################

likprob <- function(gammasig){
  gamma <- gammasig[1:(length(gammasig)-1)]
  sig <- gammasig[length(gammasig)]
  ###
  ### This function computes the conditional likelihood of the network
  ###
  obj <- rep(NA,n1)
  pos <- 1
  for (c in 1:length(X)){
    Xt <- X[[c]]
    nt <- nrow(Xt)
    #st <- pmax(X[[c]][,15],1e-6)
    st <- matrix(0,nt,1)
    for (sim in 1:nsim){
      st <- st + Seq(bethat,gamma,phihat,sig,c,sim)/nsim
    }
    D <- Dmat(gamma,c)*(st%*%t(st))
    typel <- typef(Xt)
    for (j in 1:nt){
     rowi <- as.numeric(typel==typel[j])
     Gtype <- t(matrix(rep(rowi,nt),nt,nt))
     Gtype2 <- matrix(as.numeric(G1[[c]]==G1[[c]][,j]),nt,nt)
     Gtype <- Gtype*Gtype2
     tvar <- D[,j]/(Gtype%*%st)
     tvar[is.nan(tvar)] <- 0
     D[,j] <- tvar
    }
    Gt <- G[[c]]
    diag(D) <- 0.5 # temporary value
    P <- Gt*log(D) + (1-Gt)*log(1-D)
    diag(P) <- NA
    tobj <- c(P)
    tobj <- tobj[is.na(tobj)==F]
    obj[pos:(pos-1+length(tobj))] <- tobj
    pos <- pos + length(tobj)
  }
  return(obj)
}


TSLS <- function(){
  
  ###
  ### This function computes the 2SLS estimator and variance-covariance matrix for the outcome (it is consistent)
  ###
  
  
  BX <- NULL
  BZ <- NULL
  BY <- NULL
  for (c in 1:length(X)){
    Xt <- X[[c]]
    nt <- nrow(Xt)
    D <- matrix(1/(nt-1),nt,nt)
    diag(D) <- 0
    bz <- matrix(0,nt,(29+ 18))
    bz[,1] <- as.numeric(Xt[,3]==1) # syrian
    bz[,2] <- as.numeric(Xt[,4]==2) # male
    bz[,3] <- as.numeric(Xt[,10] %in% c(2,4,5))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,10] %in% c(1,3,5))*as.numeric(Xt[,3]==2) # majority of friends are other than ones type
    bz[,4] <- as.numeric(Xt[,11] %in% c(2,4,5))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3,5))*as.numeric(Xt[,3]==2) # majority of neighbours are other than ones type
    
    ay <- Xt[,12]
    ay[is.na(ay)] <- 0
    bz[,5] <- as.numeric(ay==2012) # arrival year
    bz[,6] <- as.numeric(ay==2013) # arrival year
    bz[,7] <- as.numeric(ay==2014) # arrival year
    bz[,8] <- as.numeric(ay==2015) # arrival year
    bz[,9] <- as.numeric(ay==2016) # arrival year
    bz[,10] <- as.numeric(ay==2017) # arrival year
    bz[,11] <- as.numeric(ay==2018) # arrival year
    rs <- Xt[,13]
    rs[is.na(rs)] <- 0
    bz[,12] <- as.numeric(rs==3) # region syria
    bz[,13] <- as.numeric(rs==4) # region syria
    bz[,14] <- as.numeric(rs==5) # region syria
    bz[,15] <- as.numeric(rs==6) # region syria
    bz[,16] <- as.numeric(rs==7) # region syria
    bz[,17] <- as.numeric(rs==8) # region syria
    
    l <- Xt[,14]
    l[is.na(l)] <- 0
    bz[,18] <- l # language (if Syrian)
    bz[,(18+Xt[1,2])] <- 1
    bz[,30:(30+18-1)] <- (D%*%bz[,1:18])
    BX <- rbind(BX, cbind(bz[,1:29],(D%*%Xt[,15])))
    BY <- rbind(BY,matrix(Xt[,15],nt,1))
    BZ <- rbind(BZ, bz )
  }
  ZZ <- solve(t(BZ)%*%BZ)
  XX <- (t(BX)%*%BX)
  P1 <- solve(t(BX)%*%BZ%*%ZZ%*%t(BZ)%*%BX)
  betahat <- P1%*%t(BX)%*%BZ%*%ZZ%*%t(BZ)%*%BY
  err <- c(BY-BX%*%betahat)
  S <- matrix(0,ncol(BZ),ncol(BZ))
  for (i in 1:nrow(BZ)){
    S <- S + (err[i]^2)*t(matrix(BZ[i,],1,ncol(BZ)))%*%matrix(BZ[i,],1,ncol(BZ))
  }
  S <- S/nrow(BZ)
  varhat <- P1%*%t(BX)%*%BZ%*%ZZ%*%S%*%ZZ%*%t(BZ)%*%BX%*%P1
  olshat <- solve(t(BX[,1:(ncol(BX)-2)])%*%BX[,1:(ncol(BX)-2)])%*%t(BX[,1:(ncol(BX)-2)])%*%BY
  return(list(betahat,varhat,olshat))
}

bvecHS <- function(bet,c,Xt){
  
  ###
  ### This function computes the vector (b - epsilon) for class c, given a value for bet(a)
  ### ALTERNATIVE VERSION FOR FLUENCY COUNTERFACTUALS ------ NOT USED
  ###
  
  nt <- nrow(Xt)
  b <- matrix(as.numeric(Xt[,3]==1),nt,1)*bet[1] # syrian
  b <- b + matrix(as.numeric(Xt[,4]==2),nt,1)*bet[2] # male
  divf <- as.numeric(Xt[,10] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,10] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of friends are of the other nationality
  divn <- as.numeric(Xt[,11] %in% c(2,4))*as.numeric(Xt[,3]==1) + as.numeric(Xt[,11] %in% c(1,3))*as.numeric(Xt[,3]==2) # majority of neighbours are of the other nationality
  b <- b + matrix(divf,nt,1)*bet[3]
  b <- b + matrix(divn,nt,1)*bet[4]
  ay <- Xt[,12]
  ay[is.na(ay)] <- 0
  b <- b + as.numeric(ay==2012)*bet[5] # arrival year
  b <- b + as.numeric(ay==2013)*bet[6] # arrival year
  b <- b + as.numeric(ay==2014)*bet[7] # arrival year
  b <- b + as.numeric(ay==2015)*bet[8] # arrival year
  b <- b + as.numeric(ay==2016)*bet[9] # arrival year
  b <- b + as.numeric(ay==2017)*bet[10] # arrival year
  b <- b + as.numeric(ay==2018)*bet[11] # arrival year
  rs <- Xt[,13]
  rs[is.na(rs)] <- 0
  b <- b + as.numeric(rs==3)*bet[12] # region syria
  b <- b + as.numeric(rs==4)*bet[13] # region syria
  b <- b + as.numeric(rs==5)*bet[14] # region syria
  b <- b + as.numeric(rs==6)*bet[15] # region syria
  b <- b + as.numeric(rs==7)*bet[16] # region syria
  b <- b + as.numeric(rs==8)*bet[17] # region syria
  
  l <- Xt[,14]
  l[is.na(l)] <- 0
  b <- b + l*bet[18] # language (if Syrian)
  
  b <- b + bet[(18+Xt[1,2])]
  return(b)
}

SeqHS <- function(bet,gamma,phi,sig,c,sim,Xt){
  
  ###
  ### This function computes the equilibrium value for s in class c, given all of the parameters and the error for simulation sim
  ### ALTERNATIVE VERSION FOR FLUENCY COUNTERFACTUALS ------ NOT USED
  ###
  
  
  Et <- E[[c]][,sim]
  nt <- nrow(Xt)
  b <- bvecHS(bet,c,Xt)
  b <- b + Et*sig # random shock
  syr <- matrix(as.numeric(Xt[,3]==1),nt,1)
  b[c(syr)==1] <- b[c(syr)==1]*(1-phi[1])
  b[c(syr)==0] <- b[c(syr)==0]*(1-phi[2])
  
  D <- matrix(1/(nt-1),nt,nt)
  diag(D) <- 0
  S0 <- matrix(0,nt,1) # initial socialization guess
  sc <- 0 # stoping flag
  while (sc==0){
    S1 <- b + (phi[1]*syr + phi[2]*(1-syr))*(D%*%S0) # best response (unconstrained)
    S1 <- pmax(0,pmin(1,S1)) # best response (constrained)
    if (sum(abs(S1-S0))<1e-6){ # convergence criterion
      sc <- 1
    }
    S0 <- S1
  }
  return(S1)
}

