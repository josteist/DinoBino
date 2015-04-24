## loading palobiodb
rm(list=ls())

library(paleobioDB)
library(stats4)

# source("C:/Users/josteist/Documents/DinoBino/R_code/dinbinfun.R")

source("http://www.mn.uio.no/cees/english/people/researcher-postdoc/josteist/r_code/dinbinfun.r")

## loading occurrences of theropods
## We also want to have early_int_no and late_int_no [ein and lin] to place occs inside intervals
## achieved by adding "time" to the show entry.
dinos <- pbdb_occurrences(limit="all", base_name="Sauropodomorpha", show=c("phylo", "ident", "time"))
stegos <- pbdb_occurrences(limit="all", base_name="Stegosauria", show=c("phylo", "ident", "time"))
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
Bins = matrix(data=NA,nrow=27,ncol=5)

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
Bins[,5] = c(3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1)

## Counting occurrences inside each bin for each uniqe species
# A uniqe indexlist of occurrences matched to species rank.
uniqspec <- unique(dinos$mid[dinos$mra==3])
# Data has [species by interval] with number of occurrences per species in each interval.
Data = matrix(data=0,nrow=length(uniqspec),ncol=27)
# Times are the durations in a matrix of same size for ease of computation.
Times = matrix(data=0,nrow=length(uniqspec),ncol=27)
# Inputting the durations in the Times matrix
for (ii in 1:27){
  Times[,ii] <- Bins[ii,3]
  
}
## species by interval matrix
tix = 1 # loop counter
countdoubles = 0; # counting how many occurrences span more than 1 bin
for (ii in uniqspec) {
  ##	dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
  j1<- dinos$ein[dinos$mid==ii]
  j2<- dinos$lin[dinos$mid==ii]
  for (jj in (1:length(j1))) {
    if (j2[jj]>111){
      if (j1[jj]>j2[jj]) {
        countdoubles=countdoubles+1
        bix = seq(j1[jj],j2[jj])
        x = Bins[bix-111,3]			
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
# Number of species per interval
Nospec <- matrix(0,27,1)
for (ii in 1:27){
  Nospec[ii]<-sum(Data[,ii]>0)
  
}
## Testing, estimating the sampling probability using the functions defined below for each 
## interval

# Estimating global sampling rate, i.e. assuming a fixed fossil/sampling intensity for the 
# entire span of time.
Occs <- Data[Data>0] # Occurrences 
dTs  <- Times[Data>0] # List of durations for each of these occurrences
poisrates <- estimatePoiss(dTs,Occs)

p_glob = matrix(NA,27,3)
for (ii in 1:27){
  p_glob[ii,1] <- 1-exp(-poisrates[1]*Bins[ii,3])
  p_glob[ii,2] <- 1-exp(-poisrates[2]*Bins[ii,3])
  p_glob[ii,3] <- 1-exp(-poisrates[3]*Bins[ii,3])
#     1-exp(-coef(fit1)*Bins[ii,3])
  
}



N_true = matrix(NA,27,3)
for (ii in 1:27){
N_true[ii,] <- estimatetrue(Nospec[ii],p_glob[ii])
}
# par(mfrow=c(2,2))
midpoints=(Bins[,1]+Bins[,2])/2


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





