# Some dirty resampling the some of the dinosaurdata to see if we can get estimates of sampling rates in
# Santonian. We can, but it depends on the 'random' placing of individuals in one of several bins the occurrence
# span. 

tmp = array(0,c(100,1))
tmprate = matrix(0,100,3)
for (ii in 1:100){
  M<- createDataArrs_v2(dinoswobird)
  s<-M$Data[M$Data[,3]>0,3]
  tmp[ii] = max(s)
#   tmprate[ii,]<-estimatePoiss(Bins[3,])
  
  Occs <- M$Data[M$Data[,3]>0,3] # Occurrences 
  dTs  <- M$Times[M$Data[,3]>0,3] # List of durations for each of these occurrences
  if (sum(Occs>1)){
  tmprate[ii,] <- estimatePoiss(dTs,Occs)
  }
}