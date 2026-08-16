###working with the genoest_21k_snps.txt file from Sargon. The goal is to make it into an input file for a kinship estimator
read.csv("genoest_21k_snps.txt")->genoest

#### each locus for each individual is called from 0 - 2, where 1 is likely to be a heterzygous individual
### need to write a script to make these exactly 0, 1 or 2 for something like plink

### let's look at the distributions first
plot(genoest[,2])

which(genoest[,2:21446] == 0, arr.ind = TRUE)
###looks like most of the zeros are in col 13489 and 14085, all in one locus. Can't be estimated there?
plot(genoest[,123]) ##one weird individual
plot(genoest[,8613]) ##very weird distribution overall
plot(genoest[,13490])
plot(genoest[,14086])

###there are so few markers here, I'm just going to get rid of all of them. 

genoest_fixed<-genoest[,c(2:122, 124:8612, 8614:13489, 13491:14085, 14087:21446)]
### having looked at a bunch, it seems like there's a tight distribution around 1, but lots of variation between 1.3 and 2.
### Let's just do a hard cut off now for everything, and then we can stress test how much this matters later

genoest_calls<-as.data.frame(ifelse(genoest_fixed[,2:21441]>1.8, 2, ifelse(genoest_fixed[,2:21441]<1.2, 0, 1)))

genoest_calls_matrix<-as.matrix(genoest_calls)
###this looks like it does what we want to some extent, worth working through the vignette

##https://github.com/StoreyLab/popkin
install.packages("popkin")
library(popkin)

kinship<-popkin(genoest_calls_matrix, loci_on_cols = TRUE)
inbreeding <- inbr(kinship)
plot_popkin(kinship)

pairwise_fst <- pwfst(kinship)


### Let's get a plink and map file so that other analyses can be run.
genoest_ped_sample_info<-cbind.data.frame(c(1:length(genoest_calls[,1])), genoest[,1], NA, NA, NA, NA)
names(genoest_ped_sample_info)<-c("fam", "id", "dad", "mom", "sex", "pheno")
### need to write something here so that I get two columns per SNP, or get a different .ped style out. 
genoest_ped<- cbind.data.frame(genoest_ped_sample_info, genoest_calls_matrix)

write.table(genoest_ped, file="genoest_21k_genocalls.ped", row.names = FALSE, quote = FALSE, sep="\t")

###waiting for metadata from Sargon, but shouldn't matter for ped purposes, so faux .map file
chromosomes<-sample(c(1:24), 21440, replace=T)
names(genoest_calls)
basepair.coordinate<-sample(1:10000000, 21440, replace=T)
dummy_map<-cbind.data.frame(chromosomes, names(genoest_calls), 0, basepair.coordinate )
names(dummy_map)<-c()

write.table(dummy_map, file="genoest_21k_genocalls.map", row.names=FALSE, quote=FALSE, sep="\t")
