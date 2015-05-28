# Strange, because we get an estimate of the sampling rates using theropods only, but not for dinosaurs w/o bird (we even
# get for theropods w/o birds


#)


j = regexpr("Aves",dinos$cll)
find(j==1)
# cll has the class, 
# which here are either reptilia OR Aves. Here we can remove the birds.

j = regexpr("Reptilia",dinos$cll)

use = dinos[dinos$cll=="Reptilia",];
use$odl
# This is slightly strange, there are 5 orders of non-avian dinosaurs listed here,
# Neornithischia (part of Ornithischia)
# <NA>
# Dinosauria
# Saurischia
# Ornithischia

# So, are all theropods in Dinosauria order? Nope, it seems like they are in the <NA> category
# Dinosauria seems to be indeterminate ones
# [4 5]
use[4:5,]
use[c(1,7,9),]

sauros2 = use[use$odl=="Saurischia",];



# So we can group Neornithischia and Ornithischia --> ornits
# Saurichia is itself
# <NA> is what...
# There is some shitty R-stuff here
# Order name
dinos[1:10,]$odl
# Order number
dinos[1:10,]$odn
# Similarly with class
dinos[1:10,]$cll
dinos[1:10,]$cln

use = dinos[dinos$cln==36322,]

## We might as well use these means of generating the different groups, as perhaps the data are more commensurate
# then. Only problem is perhaps to get the theropods....

# Theropod orders are often not given
unique(theros$odl)
# Well, actually here some are classified as Saurischia, which they are though...
theros[51:70,]$odl
# Saurischia as order is order n 38505
theros[theros$odn==38505,]


##
# Thinking abit on why we cannot estimate a rate for the Santonian
# histogramming occurrences for the last three stages
for (ii in 1:3){
Occs = J$Data[J$Data[,ii]>0,ii]
hist(Occs)
}


