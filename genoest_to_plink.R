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
write.table(genoest_calls_matrix, file="genoest_calls_matrix.txt", quote = F, row.names = F, sep="\t")
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


###going to try running gemma to get the kinship from there###

###gemma runs locally as "~/gemma"
SNPs<-t(genoest_calls_matrix)
read.table("peromyscus_mac3_Q30_miss0.25_dp3_ind999_maf001_admixture_ready.2.Q")->Qscores
Qscores[,1]->q_leucopus
Qscores_bimbam<-noquote(c("q_leucopus", mean(q_leucopus), mean(q_leucopus), q_leucopus))
SNPs_bimbam<- cbind(1:length(rowMeans(SNPs)), rowMeans(SNPs), rowMeans(SNPs), SNPs) ### this is what makes it into a bimbam file we need

SNPs_q_bimbam<-noquote(rbind(SNPs_bimbam,Qscores_bimbam))
write.table(SNPs_bimbam, row.names = FALSE, col.names = FALSE, file="SNPs_q.txt")

###needs a phenotype file, but I can give it a dummy one?

y_quant<-rnorm(length(genoest_calls_matrix[,1]), 0, 1) ###obs, this is simulated data, replace with vector when available
hist(y_quant)
write.table(y_quant, row.names=FALSE, col.names=FALSE, file="y_quant.txt")
### better than this would be to add q from admixture so that it's included in the analysis as a non-snp snp. 
###have asked Sargon for this

system("~/gemma -g SNPs.txt -p y_quant.txt -gk 1 -o relmatrix", wait=TRUE) 
system("~/gemma -g SNPs_q.txt -notsnp -p y_quant.txt -gk 1 -o relmatrix_q", wait=TRUE) ###have now added q as a covariate here
system("~/gemma -g SNPs.txt -p y_quant.txt -gk 2 -o relmatrix", wait=TRUE) ###gk 2 gives us standardized matrix

read.table("./output/relmatrix.cXX.txt")->relmatrix
read.table("./output/relmatrix_q.cXX.txt")->relmatrix_q

read.table("./output/relmatrix.sXX.txt")->relmatrix_std


cor(c(as.matrix(relmatrix)), c(as.matrix(relmatrix_q))) ###these two are extremely highly correlated!

summary(c(as.matrix(relmatrix)))
summary(c(as.matrix(relmatrix_q)))
summary(c(as.matrix(relmatrix_std)))
