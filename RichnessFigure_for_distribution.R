# This script prints the estimated richness plots from our paper.

# This file plots the variables 
# N_est_int - 27by4by3 array of estimated richnesses.
# SpecRich_test - the observed richness from PBDB 
# both needs to already be defined.

# Uncomment the call below to print to file.
# png("Richness_Intnam_II.png",width=18,height=12,units='cm',res=2400)
par(mfrow=c(2,2),cex=0.4) 
# Loading a nice color library.
library(RColorBrewer)
color.list <- brewer.pal(4,"Set1")
color.list2 <- color.list
for (ii in 1:4){
  color.list2[ii] = paste(color.list[ii],'44',sep=""); 
  # Making one slightly transparent.
}
par(mar=c(3.5,3.5,3.5,3.5))
labs = c("(a)","(b)","(c)","(d)")
doints = TRUE
dolog =TRUE
plbins = 1:27;

# N_est_int = log10(N_est_int)
# SpecRich_plot = log10(SpecRich_test)
SpecRich_plot = SpecRich_test

for (ii in 1:4){
  if (dolog==TRUE){
    plot(midpoints[plbins],N_est_int[plbins,ii,1], type='l',col=color.list[1],ylim=c(.65,max(c(N_est_int[plbins,ii,]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="",ylab="",log="y",yaxt='n')
    axis(2,at=c(1,10,100,1000))
  } else {
    plot(midpoints[plbins],N_est_int[plbins,ii,1], type='l',col=color.list[1],ylim=c(.65,max(c(N_est_int[plbins,ii,]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="",ylab="",yaxt='n')
    axis(2,at=c(10,100,200,500,800,1000))
    
  }
  # Interval lines
  for (jj in 1:27){
    abline(v=max(Bins[jj,1]),lwd=0.95,col='gray87')
  }
  # Epochs/periods
  for (jj in 1:6){
    abline(v=max(Bins[Bins[,6]==jj,1]),lwd=0.5)
  }
  abline(v=min(Bins[Bins[,6]==6,2]),lwd=0.5)
  lines(midpoints[plbins],SpecRich_plot[plbins,ii],type='o')
  
  for (jj in 1:(length(plbins)-1)){
    polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(N_est_int[plbins[jj:(jj+1)],ii,2]),N_est_int[plbins[jj:(jj+1)],ii,3]),col=color.list2[4],border=NA)
    lines(midpoints[plbins[jj:(jj+1)]],N_est_int[plbins[jj:(jj+1)],ii,1],type='o',col=color.list[4],lwd=1)
    
  }
  title(Dinogroups[ii])
  # 
  if (ii==1){
    mtext("Species richness",side=2,padj=-3,cex=0.5)
  }
  if (ii==3){
    mtext("Species richness",side=2,padj=-3,cex=0.5);
    mtext("Myr",side=1,padj=2.8,cex=.5);
  }
  if (ii==4){
    mtext("Myr",side=1,padj=2.8,cex=.5)
  }
  mtext(labs[ii],side=3,adj=0,padj=-0.3,cex=0.5)
  # trying to draw abbreviations for all intervals in the richness plots.
  
  figlims = par("usr"); # getting figure limits. Try to box the abbreviations from miny to 0
  
  for (ii in 1:27){
    # rect(Bins[ii,1],1,Bins[ii,2],2)
    rect(Bins[ii,1],exp(-3),Bins[ii,2],1,lwd=0.3)
  }
  # rect(xleft, ybottom, xright, ytop,
  for (ii in 1:27){
    if (Bins[ii,3]>4){
      text(midpoints[ii],exp(figlims[3]),substr(interval.names[ii],1,2),cex=.65)
    }
  }
}

# dev.off()