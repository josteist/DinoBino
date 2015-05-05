

# Some testing of the removedoubles alternative

## Some quick loop to estimate the variability in number of species and poisson rates estimated
# due to variability in the emprand placing of multiple inverval spanners


M <- createDataArrs_v2(ornits);

Nosp_big = matrix(NA,100,27)
Pois_big = matrix(NA,100,3);
for (jj in 1:100){
  T <- createDataArrs_v2(ornits);
  Nosp_big[jj,] <- colSums(T$Data>0);
  
  Occs <- T$Data[T$Data>0] # Occurrences 
  dTs  <- T$Times[T$Data>0] # List of durations for each of these occurrences
  Pois_big[jj,] <- estimatePoiss(dTs,Occs)
}



pois_mle_big = array(NA,c(100,6,3));


# Testing for period/epoch
for (jj in 1:100){
  T <- createDataArrs_v2(ornits);
  Nosp_big[jj,] <- colSums(T$Data>0);
  for (ii in 1:6){
    
    tmp = T$Data[,Bins[,6]==ii]
    tmp2 = T$Time[,Bins[,6]==ii]
    Occs = tmp[tmp>0]
    dTs = tmp2[tmp>0]
    if (sum(Occs>1)>1){
      pois_mle_big[jj,ii,]<-estimatePoiss(dTs,Occs)
    }
    
  }
  
  
  
}


hist(pois_mle_big[,3,1])


Td <- createDataArrs_v2(ornits,removedoubles=TRUE)


sum(Md$Data)
sum(M$Data)

colSums(Md$Data>0)
colSums(M$Data>0)
par(mfrow=c(1,1))
plot(midpoints,colSums(Md$Data>0),type="o")
lines(midpoints,colSums(M$Data>0))

Occs <- M$Data[M$Data>0] # Occurrences 
dTs  <- M$Times[M$Data>0] # List of durations for each of these occurrences
poisrates_all <- estimatePoiss(dTs,Occs)
Occs <- Md$Data[Md$Data>0] # Occurrences 
dTs  <- Md$Times[Md$Data>0] # List of durations for each of these occurrences
poisrates_sings <- estimatePoiss(dTs,Occs)

# THese are pretty similar!, Slightly higher for not including the random placements. 
# But the first will depend on the random draws in the Occs.

## COmparing periods
poiscomp_per<- array(NA,c(3,2,3)) # period by dataset by mle,ci1,ci2

for (ii in 1:3){
  tmp = M$Data[,Bins[,5]==ii]
  tmp2 = M$Times[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poiscomp_per[ii,1,] <- estimatePoiss(dTs,Occs)
  tmp = Md$Data[,Bins[,5]==ii]
  tmp2 = Md$Times[,Bins[,5]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poiscomp_per[ii,2,] <- estimatePoiss(dTs,Occs)
}

# These are slightly more different. Also inconsistently different over time



## COmparing periods
poiscomp_epoch<- array(NA,c(6,2,3)) # period by dataset by mle,ci1,ci2

for (ii in 1:6){
  tmp = M$Data[,Bins[,6]==ii]
  tmp2 = M$Times[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poiscomp_epoch[ii,1,] <- estimatePoiss(dTs,Occs)
  tmp = Md$Data[,Bins[,6]==ii]
  tmp2 = Md$Times[,Bins[,6]==ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  poiscomp_epoch[ii,2,] <- estimatePoiss(dTs,Occs)
}
# And, again, using epochs make them even more different;, particularly for Epoch 3 [late Jurassic]

# Comparing intervals
poiscomp_int<- array(NA,c(27,2,3)) # period by dataset by mle,ci1,ci2

for (ii in 1:27){
  tmp = M$Data[,ii]
  tmp2 = M$Times[,ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  if (sum(Occs>1)>1){
  poiscomp_int[ii,1,] <- estimatePoiss(dTs,Occs)}
  tmp = Md$Data[,ii]
  tmp2 = Md$Times[,ii]
  Occs = tmp[tmp>0]
  dTs = tmp2[tmp>0]
  if (sum(Occs>1)>1){
  poiscomp_int[ii,2,] <- estimatePoiss(dTs,Occs)}
}
# And, again, using intervals make them even more different. But they do still correlate well
cor(poiscomp_int[,,1],use="pairwise.complete.obs")



