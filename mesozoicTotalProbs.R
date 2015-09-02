
sum(p_binos[,1,3]*N_est_int[,1,1])/sum(N_est_int[,1,1])
# These are 'estimated richness weighted' average binomial probabilities. They
genbinps <-c(0.804570578,0.742674076,0.860997285)
specbinps<-c(0.586808395,	0.494517311,	0.692982841)

# True genus richness
c(estimatetrue(1272,genbinps[1])[1],min(estimatetrue(1272,genbinps[3])),max(estimatetrue(1272,genbinps[2])))
# True species richness

c(estimatetrue(1751,specbinps[1])[1],min(estimatetrue(1751,specbinps[3])),max(estimatetrue(1751,specbinps[2])))
