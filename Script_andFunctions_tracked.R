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

Dinogroups <-c("All dinosaurs","Dinosaurs w/o birds","Ornithischia","Sauropoda","Theropoda","Theropoda w/o birds")
# SHould probably use Sauropodomorpha and not Sauropoda...
dinos  <- pbdb_occurrences(limit="all", base_name="Dinosauria", show=c("phylo", "ident", "time"))
ornits <- pbdb_occurrences(limit="all", base_name="Ornithischia", show=c("phylo", "ident", "time"))
sauros <- pbdb_occurrences(limit="all", base_name="Sauropodomorpha", show=c("phylo", "ident", "time"))
theros <- pbdb_occurrences(limit="all", base_name="Theropoda", show=c("phylo", "ident", "time"))
theroswobird <- pbdb_occurrences(limit="all", base_name="Theropoda",exclude_id=36616, show=c("phylo", "ident", "time"))
dinoswobird  <- pbdb_occurrences(limit="all", base_name="Dinosauria",exclude_id=36616, show=c("phylo", "ident", "time"))

# We could group them nestedly, as they are in a phylogeny;

# Now made new crateDataArrs_v2 which can remove doubles, and returns an array of whether or no
# the entries are in class Aves.
M <- createDataArrs_v2(dinoswobird)

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
period.names = c("Triassic","Jurassic","Cretaceaous")
epoch.names  = c("Triassic","Early Jurassic","Mid Jurassic","Late Jurassic","Early Cretaceaous","Late Cretaceous")

# Periods Triassic, Jurassic, Cretaceous
Bins[,5] = c(3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1)
# Epoch, Triassic, Early Jurassic, Mid Jurassic, Late Jurassic, Early Cretaceaous, Late Cretaceous
Bins[,6] = c(6,6,6,6,6,6,5,5,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,1,1,1,1)
midpoints=(Bins[,1]+Bins[,2])/2
dimnames(Bins)<-list(interval.names,c("Start","End","Duration","PBDB index","Period","Epoch"))


J1 = createDataArrs(dinos)
J2 = createDataArrs(dinoswobird)
J3 = createDataArrs(ornits)
J4 = createDataArrs(sauros)
J5 = createDataArrs(theros)
J6 = createDataArrs(theroswobird)
# Now Data and Times are made with the function above.

# Extracting number of species.
Nospecies = array(NA,c(27,6)); #interval by dinosaurgroup
dimnames(Nospecies)<-list(interval.names,Dinogroups)
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

## Below is for re-generating the data arrays to get a median estimate of observed richness
# This is needed because each occurrence that spans more than one interval is placed into 
# a specified interval based on the duration of the intervals. [i.e an occurrences spanning
# two intervals of duration 5 and 10 Myr are placed with p=1/3 in the first and 2/3 in the
# second]. To account for this when getting species counts, we resample 100 times and use the
# median value for each interval for each group.

Nosp_big = array(NA,c(27,6,100)); #interval by dinosaurgroup by resampled

for (ii in 1:100){
  
  # Regenerating arrays. THis is not optimized code, but works.
  J1 = createDataArrs(dinos)
  J2 = createDataArrs(dinoswobird)
  J3 = createDataArrs(ornits)
  J4 = createDataArrs(sauros)
  J5 = createDataArrs(theros)
  J6 = createDataArrs(theroswobird)
# Extracting number of species.
J = J1;
jj =  1
Nosp_big[,jj,ii] <- getNospec(J$Data)

J = J2;
jj =  2
Nosp_big[,jj,ii] <- getNospec(J$Data)

J = J3;
jj =  3
Nosp_big[,jj,ii] <- getNospec(J$Data)

J = J4;
jj =  4
Nosp_big[,jj,ii] <- getNospec(J$Data)

J = J5;
jj =  5
Nosp_big[,jj,ii] <- getNospec(J$Data)

J = J6;
jj =  6
Nosp_big[,jj,ii] <- getNospec(J$Data)

}


MeanSpec = array(NA,c(27,6,3)); # interval by group by mean/median/mode
# Extracting means and medians
for (ii in 1:27){
  for (jj in 1:6){
    MeanSpec[ii,jj,1] = mean(Nosp_big[ii,jj,])
    MeanSpec[ii,jj,2] = median(Nosp_big[ii,jj,])
    MeanSpec[ii,jj,3] = Mode(Nosp_big[ii,jj,])
  }
}
# writeClipboard(as.character(MeanSpec[,6,]))

dimnames(MeanSpec)<-list(interval.names,Dinogroups,c("Mean","Median","Mode"))
statnames <- c("MLE","ci1","ci2")

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
stat.names = c("MLE","CI1","CI2")
dimnames(p_glob)<-list(interval.names,stat.names,Dinogroups)
for (jj in 1:6){
for (ii in 1:27){
  p_glob[ii,1,jj] <- 1-exp(-poisrates[jj,1]*Bins[ii,3])
  p_glob[ii,2,jj] <- 1-exp(-poisrates[jj,2]*Bins[ii,3])
  p_glob[ii,3,jj] <- 1-exp(-poisrates[jj,3]*Bins[ii,3])
#     1-exp(-coef(fit1)*Bins[ii,3])
  
}
}

# 
# plot(midpoints,p_glob[,1,1],type = "o",lty=1,xlim = rev(range(midpoints)))
# for (jj in 2:6){
#   lines(midpoints,p_glob[,1,jj],type = "o",lty=jj)
# }


# Estimating period specific sampling rates for each of the datasets
poisrates_period = array(NA,c(3,3,6)) # storing them here, period by [mle,ci1,ci2] by dinogroup
dimnames(poisrates_period)<-list(period.names,statnames,Dinogroups)
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
dimnames(p_period)<-list(interval.names,statnames,Dinogroups)
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
dimnames(poisrates_epoch)<-list(epoch.names,stat.names,Dinogroups)
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
dimnames(p_epoch)<-list(interval.names,stat.names,Dinogroups)
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
dimnames(poisrates_interval)<-list(interval.names,stat.names,Dinogroups)
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
dimnames(p_interval)<-list(interval.names,stat.names,Dinogroups)
for (jj in 1:6){
  for (ii in 1:27){
    p_interval[ii,1,jj] <- 1-exp(-poisrates_interval[ii,1,jj]*Bins[ii,3])
    p_interval[ii,2,jj] <- 1-exp(-poisrates_interval[ii,2,jj]*Bins[ii,3])
    p_interval[ii,3,jj] <- 1-exp(-poisrates_interval[ii,3,jj]*Bins[ii,3])
    #     1-exp(-coef(fit1)*Bins[ii,3])
    
  }
}

# Generating large arrays for poisson rates to be plotted
pois_glob = array(NA,c(27,3,6))
pois_period = array(NA,c(27,3,6))
pois_epoch  = array(NA,c(27,3,6))
pois_interval = poisrates_interval;
for (ii in 1:27){
  for (jj in 1:6){
    pois_glob[ii,,jj]<-poisrates[jj,]
  }
  pois_period[ii,,]<-poisrates_period[Bins[ii,5],,]
  pois_epoch[ii,,] <- poisrates_epoch[Bins[ii,6],,]
  
}
