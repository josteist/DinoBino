## loading palobiodb

# Changing this line for github testing.

rm(list=ls())

library(paleobioDB)
library(stats4)

# source("C:/Users/josteist/Documents/DinoBino/R_code/dinbinfun.R")

# source("http://www.mn.uio.no/cees/english/people/researcher-postdoc/josteist/r_code/dinbinfun.r")
source("functions_DinoBino.R")
# Hmmm, this didn't seem to load as a set of functions...

## loading occurrences of theropods
## We also want to have early_int_no and late_int_no [ein and lin] to place occs inside intervals
## achieved by adding "time" to the show entry.


dinos  <- pbdb_occurrences(limit="all", base_name="Dinosauria", show=c("phylo", "ident", "time"))
ornits <- pbdb_occurrences(limit="all", base_name="Ornithischia", show=c("phylo", "ident", "time"))
sauros <- pbdb_occurrences(limit="all", base_name="Sauropoda", show=c("phylo", "ident", "time"))
theros <- pbdb_occurrences(limit="all", base_name="Theropoda", show=c("phylo", "ident", "time"))
theroswobird <- pbdb_occurrences(limit="all", base_name="Theropoda",exclude_id=36616, show=c("phylo", "ident", "time"))
dinoswobird  <- pbdb_occurrences(limit="all", base_name="Dinosauria",exclude_id=36616, show=c("phylo", "ident", "time"))

# We could group them nestedly, as they are in a phylogeny;

# Now made new crateDataArrs_v2 which can remove doubles, and returns an array of whether or no
# the entries are in class Aves.
M <- createDataArrs_v2(dinos)

# Dinosauria [All]
# Split into Ornithischia vs Sauropodomorpha&Theropoda
# Split into Ornithischia vs Sauropodomorpha vs Theropoda
# dinos <- pbdb_occurrences(limit="all", base_name="Ornithischia", show=c("phylo", "ident", "time"))

# NOTE: Should use the confidence intervals from the estimated sampling rate to estimate
# true diversities. If these are used then [under a variety of assumptions with varying
# speciation, extinction and different sampling rates] the true value (really true) is 
# inside this span in approx 95% of the simulated cases. Using only the mle for the sampling 
# rate makes the interval span the true value in approx 50% of the cases. The most likely
# true value varies quite some under the simulated model, but the ci's almost always (again)
# about 5%) spans the true simulated value.
# These simulations were done in Matlab (see sep.code), and simulated a BDS process (birth, death,
# fossilize) with rates spanning; spec_rate = 0...0.1, ext_rate =0...0.1, samplingPI = 0.01...0.35
# simulations were initiated with species numbers from 20 to 120 and duration was fixed at
# 5 (which will only scale the other rates). It seems like the MLE for true richness is usually
# biased downward, i.e. the estimated true richness is lower than the actual true richness)

## THis code works (for now) only for dinosaur data with specific bins and time intervals. It read
## in data from PBDB and generates data matrices necessary for the computation of richness curves
## both using the raw species counts and the estimated true diversity using the sampling rate
## estimations. 

## There is a lot of noise and some typos in the age limits in the data from PBDB. Inputting the
## consensus numbers from my own xlsheet
Bins = matrix(data=NA,nrow=27,ncol=6)

Bins[,1] = c(72.1 , 83.6, 86.3, 89.8, 93.9,
             100.5, 113, 125, 129.4, 132.9, 139.8, 145, 
             152.1, 157.3, 163.5, 166.1, 168.3, 
             170.3, 174.1, 182.7, 190.8, 199.3, 
             201.3, 208.5, 227, 237, 242)
Bins[1,2] = 66
Bins[2:27,2] = Bins[1:26,1]
Bins[,3] = Bins[,1]-Bins[,2]
Bins[,4] = seq(112,138)
## Bins are [Start End Duration intervalinx]
## Names for dinosaur bins. PBDB intervals 112:138
interval.names =c("Maastrichtian","Campanian","Santonian","Coniacian","Turonian","Cenomanian","Albian","Aptian","Barremian","Hauterivian","Valanginian","Berriasian","Tithonian","Kimmeridgian","Oxfordian","Callovian","Bathonian","Bajacian","Aalenian","Toarcian","Pliensbachian","Sinemurian","Hettangian","Rhaetian","Norian","Carnian","Ladinian")
interval2.names = c("Triassic","Jurassic","Cretaceaous")
# Periods Triassic, Jurassic, Cretaceous
Bins[,5] = c(3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1)
# Epoch, Triassic, Early Jurassic, Mid Jurassic, Late Jurassic, Early Cretaceaous, Late Cretaceous
Bins[,6] = c(6,6,6,6,6,6,5,5,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,1,1,1,1)
midpoints=(Bins[,1]+Bins[,2])/2


J1 = createDataArrs(dinos)
J2 = createDataArrs(dinoswobird)
J3 = createDataArrs(ornits)
J4 = createDataArrs(sauros)
J5 = createDataArrs(theros)
J6 = createDataArrs(theroswobird)
# Now Data and Times are made with the function above.

# Extracting number of species.
Nospecies = array(NA,c(27,6)); #interval by dinosaurgroup
J = J1;
jj =  1
Nospecies[,jj] <- getNospec(J$Data)

J = J2;
jj =  2
Nospecies[,jj] <- getNospec(J$Data)

J = J3;
jj =  3
Nospecies[,jj] <- getNospec(J$Data)

J = J4;
jj =  4
Nospecies[,jj] <- getNospec(J$Data)

J = J5;
jj =  5
Nospecies[,jj] <- getNospec(J$Data)

J = J6;
jj =  6
Nospecies[,jj] <- getNospec(J$Data)

## Testing, estimating the sampling probability using the functions defined below for each 
## interval

# Estimating global sampling rate, i.e. assuming a fixed fossil/sampling intensity for the 
# entire span of time.
# Comparing the different groups. THe are numbered as the arrays defined in lines 77-82

poisrates = matrix(NA,6,3)
# All Dinos
J = J1;
Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates[1,] <- estimatePoiss(dTs,Occs)

# Dinos wo birds
J = J2;
Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates[2,] <- estimatePoiss(dTs,Occs)

# Ornits
J = J3;
Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates[3,] <- estimatePoiss(dTs,Occs)

#Sauropods
J = J4;
Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates[4,] <- estimatePoiss(dTs,Occs)

#Theropods
J = J5;
Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates[5,] <- estimatePoiss(dTs,Occs)


#Theropods wo birds
J = J6;
Occs <- J$Data[J$Data>0] # Occurrences 
dTs  <- J$Times[J$Data>0] # List of durations for each of these occurrences
poisrates[6,] <- estimatePoiss(dTs,Occs)


# Making them be binomial probs.
p_glob = array(NA,c(27,3,6))
for (jj in 1:6){
for (ii in 1:27){
  p_glob[ii,1,jj] <- 1-exp(-poisrates[jj,1]*Bins[ii,3])
  p_glob[ii,2,jj] <- 1-exp(-poisrates[jj,2]*Bins[ii,3])
  p_glob[ii,3,jj] <- 1-exp(-poisrates[jj,3]*Bins[ii,3])
#     1-exp(-coef(fit1)*Bins[ii,3])
  
}
}


plot(midpoints,p_glob[,1,1],type = "o",lty=1,xlim = rev(range(midpoints)))
for (jj in 2:6){
  lines(midpoints,p_glob[,1,jj],type = "o",lty=jj)
}


# Estimating period specific sampling rates for each of the datasets
poisrates_period = array(NA,c(3,3,6)) # storing them here, period by [mle,ci1,ci2] by dinogroup
J = J1; jj=1;  
for (ii in 1:3){
  tmp = J$Data[,Bins[,5]==ii]
  tmp2 = J$Time[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_period[ii,,jj] <- estimatePoiss(dTs,Occs)
}

J = J2; jj=2;  
for (ii in 1:3){
  tmp = J$Data[,Bins[,5]==ii]
  tmp2 = J$Time[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_period[ii,,jj] <- estimatePoiss(dTs,Occs)
}

J = J3; jj=3;  
for (ii in 1:3){
  tmp = J$Data[,Bins[,5]==ii]
  tmp2 = J$Time[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  if (sum(Occs>1)>1){ # if enough data to estimate.
  poisrates_period[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

J = J4; jj=4;  
for (ii in 1:3){
  tmp = J$Data[,Bins[,5]==ii]
  tmp2 = J$Time[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_period[ii,,jj] <- estimatePoiss(dTs,Occs)
}

J = J5; jj=5;  
for (ii in 1:3){
  tmp = J$Data[,Bins[,5]==ii]
  tmp2 = J$Time[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_period[ii,,jj] <- estimatePoiss(dTs,Occs)
}

J = J6; jj=6;  
for (ii in 1:3){
  tmp = J$Data[,Bins[,5]==ii]
  tmp2 = J$Time[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_period[ii,,jj] <- estimatePoiss(dTs,Occs)
}

# Period specific Binomial rates
p_period = array(NA,c(27,3,6))
# Interval by [mle,ci1,ci2] by group
for (jj in 1:6){
  for (ii in 1:27){
    p_period[ii,1,jj] <- 1-exp(-poisrates_period[Bins[ii,5],1,jj]*Bins[ii,3])
    p_period[ii,2,jj] <- 1-exp(-poisrates_period[Bins[ii,5],2,jj]*Bins[ii,3])
    p_period[ii,3,jj] <- 1-exp(-poisrates_period[Bins[ii,5],3,jj]*Bins[ii,3])
    #     1-exp(-coef(fit1)*Bins[ii,3])
    
  }
}



## Estimating epoch specific sampling rates.
poisrates_epoch = array(NA,c(6,3,6));
# epoch by [mle,ci1,ci2] by group
J = J1
jj = 1;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_epoch[ii,,jj] <- estimatePoiss(dTs,Occs)
}

# Period by [mle,ci1,ci2] by group
J = J2
jj = 2;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_epoch[ii,,jj] <- estimatePoiss(dTs,Occs)
}

# Period by [mle,ci1,ci2] by group
J = J3
jj = 3;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  if (sum(Occs>1)>1){
  poisrates_epoch[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

# Period by [mle,ci1,ci2] by group
J = J4
jj = 4;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_epoch[ii,,jj] <- estimatePoiss(dTs,Occs)
}

# Period by [mle,ci1,ci2] by group
J = J5
jj = 5;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_epoch[ii,,jj] <- estimatePoiss(dTs,Occs)
}

# Period by [mle,ci1,ci2] by group
J = J6
jj = 6;
for (ii in 1:6){
  tmp = J$Data[,Bins[,6]==ii]
  tmp2 = J$Time[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poisrates_epoch[ii,,jj] <- estimatePoiss(dTs,Occs)
}



# Changing into binomial rates for comparison across time.
p_epoch = array(NA,c(27,3,6))
# binom prob interval by [mle,ci1,ci2] by dinosaurgroup with sampling rate 
# estimated per epoch.
for (jj in 1:6){
  for (ii in 1:27){
    p_epoch[ii,1,jj] <- 1-exp(-poisrates_epoch[Bins[ii,6],1,jj]*Bins[ii,3])
    p_epoch[ii,2,jj] <- 1-exp(-poisrates_epoch[Bins[ii,6],2,jj]*Bins[ii,3])
    p_epoch[ii,3,jj] <- 1-exp(-poisrates_epoch[Bins[ii,6],3,jj]*Bins[ii,3])
    #     1-exp(-coef(fit1)*Bins[ii,3])
    
  }
}



## Estimating interval specific sampling rates ----
poisrates_interval = array(NA,c(27,3,6)); #interval by [mle,ci1,ci2] by dinosaurgroup
# For each group
J = J1;
jj= 1;
for (ii in 1:27) {
  tmp = J$Data[,ii];
  tmp2 = J$Times[,ii];
  Occs = tmp[tmp>0]
  dTs  = tmp2[tmp>0];
  if (sum(Occs>1)>0) {
  poisrates_interval[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

# For each group
J = J2;
jj= 2;
for (ii in 1:27) {
  tmp = J$Data[,ii];
  tmp2 = J$Times[,ii];
  Occs = tmp[tmp>0]
  dTs  = tmp2[tmp>0];
  if (sum(Occs>1)>0) {
    poisrates_interval[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

# For each group
J = J3;
jj= 3;
for (ii in 1:27) {
  tmp = J$Data[,ii];
  tmp2 = J$Times[,ii];
  Occs = tmp[tmp>0]
  dTs  = tmp2[tmp>0];
  if (sum(Occs>1)>0) {
    poisrates_interval[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

# For each group
J = J4;
jj= 4;
for (ii in 1:27) {
  tmp = J$Data[,ii];
  tmp2 = J$Times[,ii];
  Occs = tmp[tmp>0]
  dTs  = tmp2[tmp>0];
  if (sum(Occs>1)>0) {
    poisrates_interval[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

# For each group
J = J5;
jj= 5;
for (ii in 1:27) {
  tmp = J$Data[,ii];
  tmp2 = J$Times[,ii];
  Occs = tmp[tmp>0]
  dTs  = tmp2[tmp>0];
  if (sum(Occs>1)>0) {
    poisrates_interval[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}

# For each group
J = J6;
jj= 6;
for (ii in 1:27) {
  tmp = J$Data[,ii];
  tmp2 = J$Times[,ii];
  Occs = tmp[tmp>0]
  dTs  = tmp2[tmp>0];
  if (sum(Occs>1)>0) {
    poisrates_interval[ii,,jj] <- estimatePoiss(dTs,Occs)
  }
}


# Converting the interval specific rates to binomial probabilities.

p_interval = array(NA,c(27,3,6)); #interval by [mle,ci1,ci2] by group
for (jj in 1:6){
  for (ii in 1:27){
    p_interval[ii,1,jj] <- 1-exp(-poisrates_interval[ii,1,jj]*Bins[ii,3])
    p_interval[ii,2,jj] <- 1-exp(-poisrates_interval[ii,2,jj]*Bins[ii,3])
    p_interval[ii,3,jj] <- 1-exp(-poisrates_interval[ii,3,jj]*Bins[ii,3])
    #     1-exp(-coef(fit1)*Bins[ii,3])
    
  }
}



# Temp code to export the numbers
# Number of species in the different groups
writeClipboard(as.character(Nospecies))
# Sampling probabilities as binomials. MLEs
writeClipboard(as.character(p_glob[,1,]))
writeClipboard(as.character(p_period[,1,]))
writeClipboard(as.character(p_epoch[,1,]))
writeClipboard(as.character(p_interval[,1,]))

# Late Cretaceous with cis
writeClipboard(as.character(p_interval[1:6,1,]))
writeClipboard(as.character(p_interval[1:6,2,]))
writeClipboard(as.character(p_interval[1:6,3,]))
               ## Converting the Poisson rates to binomial probabilities. 
p_binomial_all = array(NA,c(27,4,3))



# 
# Would want to get confidence bounds on the estimated true richnesses in the late cretaceous, which
# are the last 6 intervals.
colnames(Nospecies) <- c("Alldinos","Dinos_wo_birds","Ornits","Sauros","Theros","Theros_wo_birds")

# 

# Dirty scatterplots

p_matrix = p_interval[,1,]
colnames(p_matrix) <- 
  Dinogroups <-c("Alldinos","Dinos_wo_birds","Ornits","Sauros","Theros","Theros_wo_birds")


# Using interval estimated rates for all dinosaur w/o birds ----
plot(midpoints[1:6],p_interval[1:6,1,2],type="o",ylim=c(0,1),xlim=rev(range(midpoints[1:6])))
lines(midpoints[1:6],p_interval[1:6,2,2],type="l")
lines(midpoints[1:6],p_interval[1:6,3,2],type="l")

CretaRichness1 <- array(NA,c(6,3,6)) #using the upper and lower ci's from binom p
CretaRichness2 <- array(NA,c(6,3,6)) #using the upper and lower ci's from min(binomp) and max(N_true)
use_sampling<-c(2,2,2,2,2,2);
for (ii in 1:3) # for mle,ci1,ci2
{
  for (jj in 1:6) # for the cretaceous intervals
  {
    for (oo in 1:6) # for each species
    {
     if (!is.na(p_interval[jj,ii,use_sampling[oo]])){
      CretaRichness1[jj,ii,oo] = floor(Nospecies[jj,oo]/p_interval[jj,ii,use_sampling[oo]])
      tmp = estimatetrue(Nospecies[jj,oo],p_interval[jj,ii,use_sampling[oo]])
      if (ii==1){ CretaRichness2[jj,ii,oo] = tmp[1]}
      
      if (ii==2){ CretaRichness2[jj,ii,oo] = max(tmp)}
      if (ii==3){ CretaRichness2[jj,ii,oo] = min(tmp)}
    }
  }
  }
}

par(mfrow=c(2,3))
for (ii in 1:6){
plot(midpoints[1:6],CretaRichness2[1:6,1,ii],type="o",ylim=c(5,max(CretaRichness2[1:6,2,ii],na.rm=TRUE)),xlim=rev(range(midpoints[1:6])),xlab="Mya",ylab="Species richness",log="y")
lines(midpoints[1:6],CretaRichness2[1:6,2,ii],type="l")
lines(midpoints[1:6],CretaRichness2[1:6,3,ii],type="l")
title(Dinogroups[ii])
}


cor(p_matrix,use="pairwise.complete.obs")
## Below are old cuts -----
# Interval by sampling period estimation by [mle, lower ci, upper ci]
# Global, Period, Interval
for (ii in 1:27){
  # Global rates
  dt = Bins[ii,3] # Duration of this bin
  for (jj in 1:3){ # for mle, ci1, ci2
  p_binomial_all[ii,1,jj] <- 1-exp(-dt*poisrates[jj])
  p_binomial_all[ii,2,jj] <- 1-exp(-dt*poisrates_period[Bins[ii,5],jj])
  p_binomial_all[ii,3,jj] <- 1- exp(-dt*poisrates_epoch[Bins[ii,6],jj])
  p_binomial_all[ii,4,jj] <- 1-exp(-dt*poisrates_interval[ii,jj])
  }
  
  
}

## Plotting the binomial sampling rates over time for different lumped data.
ymin = 0
ytx  = ymin - 0.125*(ymax-ymin)
plot(midpoints,p_binomial_all[,1,1],type="o",xlab="M",ylab="Binomial sampling probability",xlim=rev(range(midpoints)),ylim=c(0,1))
lines(midpoints,p_binomial_all[,2,1],type="o",col = 2)
lines(midpoints,p_binomial_all[,3,1],type="o",col = 3)
lines(midpoints,p_binomial_all[,4,1],type="o",col = 4)
legend("bottomright",c("Global sampling rate","Period sampling rate","Epoch sampling rate","Interval sampling rate"),col=c(1:4),pch=c(1,1,1))
# abline(h=0)
# text(midpoints,ytx,interval.names,srt=90)
for (ii in 1:27){
  abline(v=Bins[ii,1],col="lightgray")
  #ep(Bins[ii,1],2),c(yplmin,0))
}
abline(v=Bins[1,2],col="lightgray")


N_true = matrix(NA,27,3)
N_true_all = array(NA,c(27,3,5))
# As sampling probs; interval by sampling method by [mle, ci1,ci2, ci3 ,ci4]
# Where ci1,ci2 uses the CIS from the binomial distribution using the MLE p_binom
# and ci3,ci4 is the ci's using the min and max sampling prob from the ci's of the binomial probability.

for (ii in 1:27){ #for each interval
  for (jj in 1:3){ # for each sampling rate - type [global, period, interval]
    # N_true[ii,] <- estimatetrue(Nospec[ii],p_glob[ii])
    for (oo in 1:3){ # for each estimate [MLE, p-binoci1 p-binoci2]
      # But for each sampling rate we have N[mle], N[lowci],N[highci] from the binomial
      # for ii==1 we store all of these, since we want these ci's.
      # for ii=2,3 we take the
      if (is.na(p_binomial_all[ii,jj,oo])){
        # don't compute since rate is not defined
      } else {
        n_tmp <- estimatetrue(Nospec[ii],p_binomial_all[ii,jj,oo])
        
        if (oo==1){
          N_true_all[ii,jj,1:3]<-n_tmp
        }
        else if (oo==2){
          N_true_all[ii,jj,5] <- max(n_tmp)
        }
        else if (oo==3){
          N_true_all[ii,jj,4] <- min(n_tmp)
        }
        
        
      }
    }
  }
}
# par(mfrow=c(2,2))



##
# THinking about what to plot;
N_true_all[1,,]
# SO for any interval we have 15 richnesses estimated;
# 3 by 5
# Method [global, period, interval]
# Estimate [N_MLE & P_MLE, N_ci1 & P_mle, N_ci2 & P_mle, N_ci3 & max(P_ci], N_ci & min(P_ci))]

# We want to compare MLE's
plot(midpoints,N_true_all[,1,1],type="o",lty=2,col=2,xlim = rev(range(midpoints)))
lines(midpoints,N_true_all[,2,1],type="o",lty=3,col=3)
lines(midpoints,N_true_all[,3,1],type="o",lty=4,col=4)
lines(midpoints,Nospec,type="o",lty=4,col=1)
legend("topleft",c("Observed richness","Global sampling rate","Period sampling rate","Interval sampling rate"),col=c(1:4),pch=c(1,1,1,1))

# Observed richness for sure
# MLE richness N_true_all[,1,1] - global sampling rate
# MLE richness N_true_all[,
# plot(midpoints,N_true_all[,1,1],type="o",xlim= rev(range(midpoints)),ylab='Sauropodomorphs richness',xlab="Myr")
#   lines(midpoints,Nospec)

# Global sampling with max/min (i.e. ci's using min/max Prob)
par(mfrow=c(2,2))
plot(midpoints,N_true_all[,1,1],type="o",xlim=rev(range(midpoints)))
lines(midpoints,N_true_all[,1,4])
lines(midpoints,N_true_all[,1,5])
title('Global sampling rate')
plot(midpoints,N_true_all[,2,1],type="o",xlim=rev(range(midpoints)))
lines(midpoints,N_true_all[,2,4])
lines(midpoints,N_true_all[,2,5])
title('Period sampling rate')
plot(midpoints,N_true_all[,3,1],type="o",xlim=rev(range(midpoints)),ylim=c(0,2500))
lines(midpoints,N_true_all[,3,4])
lines(midpoints,N_true_all[,3,5])
title('Interval sampling rate')



ymax = max(N_true[,1])
  #max(Nospec)
ymin = 0
ytx  = ymin - 0.125*(ymax-ymin)
yplmin = ymin- .2*(ymax-ymin)
par(mar=c(5,4,4,2)+ 0.1,ps=6)
plot(1,1,xlim=c(max(Bins[,1]),min(Bins[,2])),ylim=c(yplmin,ymax),yaxt='n',xlab="Geological time (Ma)",ylab="Species richness")

lines(midpoints,Nospec,xlab="",ylab="",yaxt='n')
lines(midpoints,N_true[,1],lty=2,yaxt='n')
lines(midpoints,N_true[,2],lty=2,yaxt='n',col="lightgray")
lines(midpoints,N_true[,3],lty=2,yaxt='n',col="lightgray")

sp1 <- seq(0,round(ymax/10)*10,10)
axis(side=2,sp1)
#      ylab = "Species richness")
# abline(h=ytx)
abline(h=0)
text(midpoints,ytx,interval.names,srt=90)
for (ii in 1:27){
  abline(v=Bins[ii,1],col="lightgray")
  #ep(Bins[ii,1],2),c(yplmin,0))
}
abline(v=Bins[1,2],col="lightgray")



## Plotting sampling rates

ymax = 1
#max(Nospec)
ymin = 0
ytx  = ymin - 0.125*(ymax-ymin)
yplmin = ymin- .2*(ymax-ymin)
par(mar=c(5,4,4,2)+ 0.1,ps=6)
plot(1,1,xlim=c(max(Bins[,1]),min(Bins[,2])),ylim=c(yplmin,ymax),yaxt='n',xlab="Geological time (Ma)",ylab="Sampling rate")

lines(midpoints,p_glob[,1],lty=1,yaxt='n')
lines(midpoints,p_glob[,2],lty=2,yaxt='n')
lines(midpoints,p_glob[,3],lty=2,yaxt='n')
sp1 <- seq(0,1,.1)
axis(side=2,sp1)
#      ylab = "Species richness")
# abline(h=ytx)
abline(h=0)
text(midpoints,ytx,interval.names,srt=90)
for (ii in 1:27){
  abline(v=Bins[ii,1],col="lightgray")
  #ep(Bins[ii,1],2),c(yplmin,0))
}
abline(v=Bins[1,2],col="lightgray")






