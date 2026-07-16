###working with the genoest.txt file from Sargon. The goal is to make it into an input file for a kinship estimator
read.csv("genoest.txt")->genoest

#### each locus for each individual is called from 0 - 2, where 1 is likely to be a heterzygous individual
### need to write a script to make these exactly 0, 1 or 2 for something like plink

### let's look at the distributions first
plot(genoest[,3309])

which(genoest[,2:5001] == 0, arr.ind = TRUE)
###looks like all of the zeros are in col 3308, all in one locus. Can't be estimated there?
plot(genoest[,3309])

### having looked at a bunch, it seems like there's a tight distribution around 1, but lots of variation between 1.3 and 2.
### Let's just do a hard cut off now for everything, and then we can stress test how much this matters later

genoest_calls<-as.data.frame(ifelse(genoest[,2:5001]>1.8, 2, ifelse(genoest[,2:5001]<1.2, 0, 1)))

genoest_calls_matrix<-as.matrix(genoest_calls)
###this looks like it does what we want to some extent, worth working through the vignette

##https://github.com/StoreyLab/popkin
install.packages("popkin")
library(popkin)

kinship<-popkin(genoest_calls_matrix, loci_on_cols = TRUE)
inbreeding <- inbr(kinship)
plot_popkin(kinship)