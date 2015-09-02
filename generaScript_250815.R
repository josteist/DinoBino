
## -- Generating Occurrence count matrices and species count matrices. This is done 100 times and 
# the mode of the occ.count sets are used, due to random placement of occurrence IF spanning more than
# one interval. 

# rm(list=ls())
# load("~/R/DinoBinoGit/DinoBino/PBDB_Download_only.RData")

# This has been run 250815. load Analysis_Genera.Rdata for plotting.

## -- Loading needed libraries.
library(paleobioDB)
library(stats4)
# The library below will be made available online after acceptance of main manuscript.
source("functions_DinoBino.R")

dogenera = TRUE 
# If true analysis is done on genus level
# If false analysis done on species level.

for (dd in 1:4){   # For all downloaded datasets
  # Switchin between datasets
  if (dd==1){
    usenow = dinos;
  } else if (dd==3){
    usenow = sauros;
  } else if (dd==2){
    usenow = ornits;
  } else if (dd==4){
    usenow = theros
  }
  
  # Occurrence counts. Temporary array to get number of unique species/genera
  if (dogener==TRUE){
    tmp = createDataArrs_v4(usenow);
  } else {
    tmp = createDataArrs_v3(usenow);
  }
  
  tmpOccs = array(data=NA,dim=c(100,nrow(tmp$Data),27)); # large array of all 100 drawn occurrence matrices
  tmpSpec = matrix(data=NA,nrow=100,ncol=27);  # large array of 100 drawn occ.matrices' species richness counts.
  for (rr in 1:100){
    # Replicating all collection of data.
    if (dogener==TRUE){
    J = createDataArrs_v4(usenow);
    } else {
      J = createDataArrs_v3(usenow);
    }
    tmpSpec[rr,]=colSums(J$Data>0); # Collecting species number
    tmpOccs[rr,,]=J$Data; # storing the occurrence count matrix
  }
  # Finding the modal occurrence count for each species.
  Occs = array(data=NA,dim=c(nrow(tmp$Data),27));
  for (ss in 1:nrow(tmp$Data)){
    uniqsets = unique(tmpOccs[,ss,])
    repunq = array(0,c(100,1))
    for (iq in 1:length(uniqsets[,1])){
      # For each unique set of occurrence counts
      for (rr in 1:100){
        # For each replicate
        if (sum(abs(tmpOccs[rr,ss,]-uniqsets[iq,]))==0){
          # if this replicate has this combination of occ.counts.
          repunq[rr] = iq;
        }
      }
    }
    # Most common uniqe set of occ.counts
    tiq = Mode(repunq);
    Occs[ss,] = uniqsets[tiq,];
  }
  
  # Replicating sampling rate analysis for interval specific rates.
  p_interval_now = array(data=NA,dim=c(100,27,3));
  Spec_count_now = array(data=NA,dim=c(100,27)); 
  for (rr in 1:100){
    Occs = tmpOccs[rr,,];
    Spec_count_now[rr,] = colSums(Occs>0);
    for (ss in 1:27){
      if (sum((Occs[Occs[,ss]>0,ss])>1)){
        p_interval_now[rr,ss,]  <- estimatePoiss(array(Bins[ss,3],length(Occs[Occs[,ss]>0,ss])),Occs[Occs[,ss]>0,ss])
      }
      
    }
  }
  
  # Tallying range-through diversity.
  tmp_RT = array(NA,dim=c(100,27))
  for (rr in 1:100){
    Occs = tmpOccs[rr,,];
    Range_through = array(NA,dim=dim(Occs))
    for (ss in 1:(dim(Occs)[1])){
      Range_through[ss,seq(min(which(Occs[ss,]>0)),max(which(Occs[ss,]>0)))]=1
    }
    tmp_RT[rr,] = colSums(Range_through,na.rm=TRUE); # Storing range-through richness count.
  }
  
  
  
  colnames(Occs)<-colnames(J$Data);
  rownames(Occs)<-rownames(J$Data);
  
  # Switchin between datasets, store to workspace
  if (dd==1){
    dinos_occs = Occs;
    dinos_occs_tmp = tmpOccs;
    dinos_times = J$Times;
    spec_count_dinos =  Spec_count_now;
    p_interval_dinos = p_interval_now;
    spec_count_rt_dinos = tmp_RT;
  } else if (dd==2){
    ornits_occs=Occs;
    ornits_occs_tmp = tmpOccs;
    ornits_times = J$Times;
    spec_count_ornits =  Spec_count_now;
    p_interval_ornits = p_interval_now;
    spec_count_rt_ornits = tmp_RT;
  } else if (dd==3){
    sauros_occs = Occs;
    sauros_occs_tmp = tmpOccs;
    sauros_times = J$Times;
    spec_count_sauros =  Spec_count_now;
    p_interval_sauros = p_interval_now;
    spec_count_rt_sauros = tmp_RT;
  } else if (dd==4){
    theros_occs=Occs;
    theros_occs_tmp = tmpOccs;
    theros_times = J$Times;
    spec_count_theros =  Spec_count_now;
    p_interval_theros = p_interval_now;
    spec_count_rt_theros = tmp_RT;
  }
  
}

# Number of genera

dim(dinos_occs_tmp)
dim(ornits_occs_tmp)
dim(sauros_occs_tmp)
dim(theros_occs_tmp)


## Getting 'median' rates. 
p_int_median = array(NA,dim=c(27,4,3))
for (ss in 1:27){
  # Dinosauria
  tix = which(p_interval_dinos[,ss,1]==sort(p_interval_dinos[,ss,1])[floor(length(sort(p_interval_dinos[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,1,]=p_interval_dinos[tix[1],ss,]
  # Ornithischia
  tix = which(p_interval_ornits[,ss,1]==sort(p_interval_ornits[,ss,1])[floor(length(sort(p_interval_ornits[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,2,]=p_interval_ornits[tix[1],ss,]
  
  # Sauropodomorpha
  tix = which(p_interval_sauros[,ss,1]==sort(p_interval_sauros[,ss,1])[floor(length(sort(p_interval_sauros[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,3,]=p_interval_sauros[tix[1],ss,]
  
  
  # Theropoda
  tix = which(p_interval_theros[,ss,1]==sort(p_interval_theros[,ss,1])[floor(length(sort(p_interval_theros[,ss,1]))/2)])
  # if more than one hit, just use the first since they are identical
  p_int_median[ss,4,]=p_interval_theros[tix[1],ss,]
  
}

# generating binomial sampling probability array
dimnames(p_int_median)<-list(interval.names,Dinogroups,statnames)
p_binos = array(NA,c(27,4,3));
dimnames(p_binos)<-list(interval.names,Dinogroups,statnames)
for (ss in 1:27){
  for (dd in 1:4){
    p_binos[ss,dd,]=1-exp(-p_int_median[ss,dd,]*Bins[ss,3])
  }
}

# Collecting observed species richness.
SpecRich_test = array(NA,c(27,4));
SpecRich_RT   = array(NA,c(27,4)); #Range through
for (ss in 1:27){
  SpecRich_test[ss,1]=round(median(rowSums(dinos_occs_tmp[,,ss]>0)))
  SpecRich_test[ss,2]=round(median(rowSums(ornits_occs_tmp[,,ss]>0)))
  SpecRich_test[ss,3]=round(median(rowSums(sauros_occs_tmp[,,ss]>0)))
  SpecRich_test[ss,4]=round(median(rowSums(theros_occs_tmp[,,ss]>0)))
  SpecRich_RT[ss,1] = max(spec_count_rt_dinos[,ss])
  SpecRich_RT[ss,2] = max(spec_count_rt_ornits[,ss])
  SpecRich_RT[ss,3] = max(spec_count_rt_sauros[,ss])
  SpecRich_RT[ss,4] = max(spec_count_rt_theros[,ss])
}
dimnames(SpecRich_test)<-list(interval.names,Dinogroups)
dimnames(SpecRich_RT)<-list(interval.names,Dinogroups)


# Estimating true richnesses
N_est_int = array(NA,dim=c(27,4,3))
for (dd in 1:4){
  for (ss in 1:27){
    N_est_int[ss,dd,1]=floor(SpecRich_test[ss,dd]/p_binos[ss,dd,1]);
    N_est_int[ss,dd,2]=max(estimatetrue(SpecRich_test[ss,dd],p_binos[ss,dd,2]));
    N_est_int[ss,dd,3]=min(estimatetrue(SpecRich_test[ss,dd],p_binos[ss,dd,3]));
    
  }
}

# Generating summary table. Here the upper and lower have been confused, but they are obvious.
Dino_out = array(NA,dim=c(27,8))
dimnames(Dino_out)<-list(interval.names,c("Species richness","Species richness - range through","Estimated true richness - MLE","Estimated true richness - upper ci","Estimated true richness - lower ci","Binomial sampling probability - MLE","Binomial sampling probability - lower ci","Binomial sampling probability - upper ci"))
Dino_out[,1] = SpecRich_test[,1];
Dino_out[,2] = SpecRich_RT[,1];
Dino_out[,3:5] = N_est_int[,1,];
Dino_out[,6:8] = p_binos[,1,];
write.table(Dino_out,"Dinosaria_out_genera.txt",sep="\t")

Dino_out = array(NA,dim=c(27,8))
dimnames(Dino_out)<-list(interval.names,c("Species richness","Species richness - range through","Estimated true richness - MLE","Estimated true richness - upper ci","Estimated true richness - lower ci","Binomial sampling probability - MLE","Binomial sampling probability - lower ci","Binomial sampling probability - upper ci"))
Dino_out[,1] = SpecRich_test[,2];
Dino_out[,2] = SpecRich_RT[,2];
Dino_out[,3:5] = N_est_int[,2,];
Dino_out[,6:8] = p_binos[,2,];
write.table(Dino_out,"Ornithischia_out_genera.txt",sep="\t")

Dino_out = array(NA,dim=c(27,8))
dimnames(Dino_out)<-list(interval.names,c("Species richness","Species richness - range through","Estimated true richness - MLE","Estimated true richness - upper ci","Estimated true richness - lower ci","Binomial sampling probability - MLE","Binomial sampling probability - lower ci","Binomial sampling probability - upper ci"))
Dino_out[,1] = SpecRich_test[,3];
Dino_out[,2] = SpecRich_RT[,3];
Dino_out[,3:5] = N_est_int[,3,];
Dino_out[,6:8] = p_binos[,3,];
write.table(Dino_out,"Sauropodomorpha_out_genera.txt",sep="\t")


Dino_out = array(NA,dim=c(27,8))
dimnames(Dino_out)<-list(interval.names,c("Species richness","Species richness - range through","Estimated true richness - MLE","Estimated true richness - upper ci","Estimated true richness - lower ci","Binomial sampling probability - MLE","Binomial sampling probability - lower ci","Binomial sampling probability - upper ci"))
Dino_out[,1] = SpecRich_test[,4];
Dino_out[,2] = SpecRich_RT[,4];
Dino_out[,3:5] = N_est_int[,4,];
Dino_out[,6:8] = p_binos[,4,];
write.table(Dino_out,"Theropoda_out_genera.txt",sep="\t")

midpoints = (Bins[,1]+Bins[,2])/2
# Printing figures
source('SamplingFigures_fordist.R')
# source('RichnessFigure_for_distribution.R')
source('RichnessFigure_withRT.R')

# save.image(file="Analysis_Genera.Rdata")