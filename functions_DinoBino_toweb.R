## Functions:

# Drawing a number/index from an empirical probability distribution.
Emprand <- function(x) {
  j <- runif(1) ## drawing random uniform number
  which(cumsum(x)/sum(x)>j)[1] #returning first entry in the cumsum/sum [empirical density function]
}

## Likelihood function of the poisson intensity conditioned on more than one occurrence.
likepoissint <- function(lambda,dt,nobs) {
  (((lambda*dt)^(nobs))/(factorial(nobs))*(exp(-lambda*dt)))/(1- exp(-lambda*dt))
}

## Negative Log likelihood for many observations
# Make more general, i.e. date multiple nobs, but one dt value (then repeat it)
# 
loglikepoissint <- function(lambda,dtin,nobs) {
  if (length(nobs)!=length(dtin)){
    dt = rep(dt,length(nobs))
  }else{
    dt = dtin
  } # if not equal lengths, assume dt has 1 value and repeat
  
  ll <- 0
  for (ii in (1:length(dt))){
    ll <- ll - log(likepoissint(lambda,dt[ii],nobs[ii]))
  }
  return(ll)
}

# Function to return phat and pci
estimatePoiss <- function(dTs,Occs){
  # Negative log likelihood of the observations.
  nllnow <- function(lambda){ loglikepoissint(lambda,dTs,Occs)}
  # Getting maximum likelihood estimate of the poisson rate.
  fit1 <- mle(nllnow,start = list(lambda=1),nobs=length(Occs),method="Brent",lower = 1e-8,upper=10)
  # Defining the likelihood ratio function to estimate confidence intervals
  mycis2 <- function(l_x) ((2*(-nllnow(coef(fit1))+nllnow(l_x)))-qchisq(0.95,1))^2
  
  # lower CI
  ci1<-optim(0.9*coef(fit1),mycis2,method="Brent",lower=1e-6,upper=coef(fit1))
  # upper CI
  ci2<-optim(1.1*coef(fit1),mycis2,method="Brent",lower=coef(fit1),upper=10)
  return(c(coef(fit1),ci1$par,ci2$par))
}

## Function to collect occurrances and durations from PBDB output
# This is custom made for dinosaur data, it will only extract occurrences from
# pbdb intervals 112:138
createDataArrs <- function(dinos){
  
  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Data has [species by interval] with number of occurrences per species in each interval.
  Data = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Times are the durations in a matrix of same size for ease of computation.
  Times = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Times[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        # OR after Ladinian (WE include Ladinian here)
        
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]			
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Data[tix,binow-(Bins[1,4]-1)] <- Data[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  results <- list(Data = Data, Times=Times)
  return(results)
}


createDataArrs_v2 <- function(dinos,removedoubles=FALSE){
  
  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  # This function adds a third list [Isbird] which is true/false if the taxon is a bird.
  # This is to 'unify' the datamatrices with dinos and theropods with and without birds.
  # if removedoubles is set to TRUE, each occurrence that spans more than one interval is ignored.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Data has [species by interval] with number of occurrences per species in each interval.
  Data = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  Isbird = array(NA,c(length(uniqspec),1))
  # Times are the durations in a matrix of same size for ease of computation.
  Times = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Times[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter this is for each unique species
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # TESTING TO INCLUDE LADINIAN_ 040615
        if (j1[jj]>j2[jj]) {
          if (removedoubles==TRUE){
            binow <- -100 #make it not count.
          }else{
            countdoubles=countdoubles+1
            bix = seq(j1[jj],j2[jj])
            x = Bins[bix-(Bins[1,4]-1),3]  		
            binow <- bix[Emprand(x)]
          }
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Data[tix,binow-(Bins[1,4]-1)] <- Data[tix,binow-(Bins[1,4]-1)]+1
          
          
        }
      }
    }    
    Isbird[tix] = any(dinos[dinos$mid==ii,]$cln==36616)
    Isbird[tix]
    tix <-tix+1
  }
  
  results <- list(Data = Data, Times=Times,Isbird=Isbird)
  return(results)
}


createDataArrs_v3 <- function(dinos){
  
  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Dats has [species by interval] with number of occurrences per species in each interval.
  Dats = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  uniqnam<-unique(dinos$mna[dinos$mra==3])
  rownames(Dats)<-uniqnam;
  colnames(Dats)<-interval.names;
  # Trying to not defined Data, but do it iteratively to accurately get the ones included. earlier (bf 070815)
  # this code included rows with no occurrences used.
  # The prealloc is needed for computation, how about removing the zeros
  # Times are the durations in a matrix of same size for ease of computation.
  Tims = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Tims[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]			
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Dats[tix,binow-(Bins[1,4]-1)] <- Dats[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  tmp<-rowSums(Dats)>0;
  Data = Dats[tmp,]
  Times = Tims[tmp,]
  results <- list(Data = Data, Times=Times)
  return(results)
}


createDataArrs_v4 <- function(dinos){
  
  ## Counting occurrences inside each bin for each uniqe GENUS
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$gnn)
  uniqspec = uniqspec[!is.na(uniqspec)] # if there's a NA
  
  # Dats has [species by interval] with number of occurrences per species in each interval.
  Dats = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  dinos[which(dinos$gnn==uniqspec[2]),]$gnl
  uniqnam = array(NA,dim=c(length(uniqspec),1))
  for (ii in 1:length(uniqspec)){
    uniqnam[ii] = unique(dinos[which(dinos$gnn==uniqspec[ii]),]$gnl)
  }
  # uniqnam<-unique(dinos$gnl)
  rownames(Dats)<-uniqnam;
  colnames(Dats)<-interval.names;
  # Trying to not defined Data, but do it iteratively to accurately get the ones included. earlier (bf 070815)
  # this code included rows with no occurrences used.
  # The prealloc is needed for computation, how about removing the zeros
  # Times are the durations in a matrix of same size for ease of computation.
  Tims = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Tims[,ii] <- Bins[ii,3]
    
  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[which(dinos$gnn==ii)]#dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[which(dinos$gnn==ii)]#dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]			
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Dats[tix,binow-(Bins[1,4]-1)] <- Dats[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  tmp<-rowSums(Dats)>0;
  Data = Dats[tmp,]
  Times = Tims[tmp,]
  results <- list(Data = Data, Times=Times)
  return(results)
}

# To estimate maximum likelihood of the true number of species given observed number and 
# binomial sampling probability.
# pdf = 
estimatetrue <- function(nobs,binomprob) {
  if (!is.na(binomprob)){
    n <- seq(0,nobs/binomprob * 4+10)
    liks <- log(dbinom(nobs,size=n,prob=binomprob))
    tmp <- n[which((2*(max(liks)-liks))<qchisq(0.95,1))]
    
    return(c(n[which.max(liks)],min(tmp),max(tmp)))
  } else {
    return(c(NA,NA,NA))
  }
}

