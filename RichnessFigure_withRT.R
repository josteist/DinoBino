# This script prints the estimated richness plots from our paper.

# This file plots the variables 
# N_est_int - 27by4by3 array of estimated richnesses.
# SpecRich_test - the observed richness from PBDB 
# both needs to already be defined.

# Uncomment the call below to print to file.
# jpeg(printfilename,width=8,height=16,units='cm',res=2400)
par(mfrow=c(4,1),cex=0.4) 
# Loading a nice color library.
library(RColorBrewer)
color.list <- brewer.pal(4,"Set1")
color.list2 <- color.list
for (ii in 1:4){
  color.list2[ii] = paste(color.list[ii],'44',sep=""); 
  # Making one slightly transparent.
}
# par(mar=c(3.5,3.5,3.5,3.5))
labs = c("(a)","(b)","(c)","(d)")
doints = TRUE
dolog =TRUE
plbins = 1:27;

# N_est_int = log10(N_est_int)
# SpecRich_plot = log10(SpecRich_test)
SpecRich_plot = SpecRich_test

for (ii in 1:4){
  if (dolog==TRUE){
    par(mar=c(0+(ii>3)*4,4,4-(ii>1)*4,2))
    # par(mar=c(1,4,1,4))
    #4,4,0,2 is last
    #0,4,4,2 is first
    #0,4,0,2 is 2nd and 3rd
    plot(midpoints[plbins],N_est_int[plbins,ii,1], type='l',col=color.list[ii],ylim=c(1-(0.9*(ii==4)),max(c(N_est_int[plbins,ii,]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="",ylab="",log="y",yaxt='n',xaxt="n")
    axis(2,at=c(1,10,100,1000))
  } else {
    plot(midpoints[plbins],N_est_int[plbins,ii,1], type='l',col=color.list[ii],ylim=c(.65,max(c(N_est_int[plbins,ii,]),na.rm=TRUE)*1.1),xlim = rev(range(midpoints[plbins])),xlab="",ylab="",yaxt='n')
    axis(2,at=c(10,100,200,500,800,1000))
    
  }
#   # Interval lines
#   for (jj in 1:27){
#     abline(v=max(Bins[jj,1]),lwd=1,col='gray87')
#   }
#   # Epochs/periods
#   for (jj in 1:6){
#     abline(v=max(Bins[Bins[,6]==jj,1]),lwd=1)
#     
#   }
#   abline(v=min(Bins[Bins[,6]==6,2]),lwd=1)
#   
  for (jj in plbins){
    lines(rep(Bins[jj,1],2),c(.3,1e4),col="lightgray",lwd=1)
  }
  lines(rep(Bins[1,2],2),c(.3,1e4),col="lightgray",lwd=1)
  for (jj in 1:6){
    lines(rep(max(Bins[Bins[,6]==jj,1]),2),c(0.3,1e4),col="grey10",lwd=1)
  }
  lines(rep(min(Bins[Bins[,6]==6,2]),2),c(0.3,1e4),col="grey10",lwd=1)
        
  lines(midpoints[plbins],SpecRich_plot[plbins,ii],type='o')
  lines(midpoints[plbins],SpecRich_RT[plbins,ii],type='o',pch=6,lty=3)
  
  for (jj in 1:(length(plbins)-1)){
    polygon(c(rev(midpoints[plbins[jj:(jj+1)]]),midpoints[plbins[jj:(jj+1)]]),c(rev(N_est_int[plbins[jj:(jj+1)],ii,2]),N_est_int[plbins[jj:(jj+1)],ii,3]),col=color.list2[ii],border=NA)
    lines(midpoints[plbins[jj:(jj+1)]],N_est_int[plbins[jj:(jj+1)],ii,1],type='o',col=color.list[ii],lwd=1)
    
  }
  # title(Dinogroups[ii])
  # 
  # if (ii==1){
  if (dospec==1){
  mtext("Species richness",side=2,padj=-3,cex=0.5)
  } else {
    mtext("Genus richness",side=2,padj=-3,cex=0.5)
  }
  # }
  # if (ii==3){
  # mtext("Species richness",side=2,padj=-3,cex=0.5);
  # mtext("Myr",side=1,padj=2.8,cex=.5);
  
  if (ii==4){
    mtext("Ma",side=1,padj=2.8,cex=.5)
  } else {
    axis(side=1,at=0)
  }
  # mtext(labs[ii],side=2,adj=.3,padj=1.3,cex=0.5)
  mtext(Dinogroups[ii],side=2,adj=.6,padj=1.6,cex=0.5)
  # trying to draw abbreviations for all intervals in the richness plots.
  
  figlims = par("usr"); # getting figure limits. Try to box the abbreviations from miny to 0
  
}
# for (ii in 1:27){
#   # rect(Bins[ii,1],1,Bins[ii,2],2)
#   rect(Bins[ii,1],exp(-3),Bins[ii,2],1,lwd=0.3)
# }
# xtx = array(NA,c(7,1))
# for (ii in 1:6){
#   xtx[ii] = max(Bins[Bins[,6]==ii,1])
#   
# }
# xtx[7] = min(Bins[Bins[,6]==6],1)
axis(1,at=c(240,200,150,100,66))

color.list.eps = c("grey70","grey90")
for (ii in 1:27){
  rect(Bins[ii,1],.3,Bins[ii,2],1,lwd=1,col = color.list.eps[Bins[ii,6] %% 2])
}

for (ii in 1:6){
  rect(max(Bins[Bins[,6]==ii,1]),.01,min(Bins[Bins[,6]==ii,2]),.3,lwd=1,col = color.list.eps[ii %% 2])
}
# for (ii in 1:27){
#   text(midpoints[ii],-.05,substr(interval.names[ii],1,3),cex=.5,srt=90)
# }
for (ii in 1:6){
  text((min(Bins[Bins[,6]==ii,2])+max(Bins[Bins[,6]==ii,1]))/2,.15,epoch.names[ii],cex=1)
}

# 
for (ii in 1:27){
    text(midpoints[ii],.5,substr(interval.names[ii],1,3),cex=.7,srt=90)
}


# dev.off()