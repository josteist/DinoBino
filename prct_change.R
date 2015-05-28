# Making a quick and dirty % change plot over time using both true and MLE richnesses.
whd = 2; #which dinosaurgroup
tmp <- MeanSpec[,whd,1]
tmp2 <- N_est_int[,1,whd]
par(mfrow=c(1,1))
pcdif <- 100*(tmp[1:26]/tmp[2:27]-1)
pcdif2 <- 100*(tmp2[1:26]/tmp2[2:27]-1)

plot(midpoints[2:21],pcdif[1:20],type="o")
lines(midpoints[2:21],pcdif2[1:20],type="o",col="grey")