### want to see how correlated relatednesses are from different variant calling protocols

### miss_0.6 and miss_0.8 have the same number of individuals, but different numbers of markers
read.table("../miss06.relatedness2", header=T)->miss_0.6
read.table("../miss08.relatedness2", header=T)->miss_0.8
paste(miss_0.6$INDV1, miss_0.6$INDV2)->miss_0.6$pair
paste(miss_0.8$INDV1, miss_0.8$INDV2)->miss_0.8$pair

merge(miss_0.6, miss_0.8, by="pair")->merged_miss_0.6_0.8

##plot(merged_miss_0.6_0.8$RELATEDNESS_PHI.x, merged_miss_0.6_0.8$RELATEDNESS_PHI.y)
cor.test(merged_miss_0.6_0.8$RELATEDNESS_PHI.x, merged_miss_0.6_0.8$RELATEDNESS_PHI.y)

### these are only 80% correlated. 
### but there are only 372 SNPs for 0.6, compared to 72K for 0.8

read.table("../miss0.35_ind99.relatedness2", header=T)->miss_0.35_ind99
paste(miss_0.35_ind99$INDV1, miss_0.35_ind99$INDV2)->miss_0.35_ind99$pair
merge(miss_0.35_ind99,miss_0.8, by="pair")->merged_miss_0.8_0.35
cor.test(merged_miss_0.8_0.35$RELATEDNESS_PHI.x, merged_miss_0.8_0.35$RELATEDNESS_PHI.y)
### these are 86% correlated. Still only 1747 SNPs.

merge(miss_0.35_ind99,miss_0.6, by="pair")->merged_miss_0.6_0.35
cor.test(merged_miss_0.6_0.35$RELATEDNESS_PHI.x, merged_miss_0.6_0.35$RELATEDNESS_PHI.y)
### these are 96.7% correlated, as these both have relatively few SNPs. 

##Question - how is relatedness being biased by more SNPs? I would think it's just more precise, but if there's an upward/downward bias, I want to know that. 


### new vcfs from Sargon - all three have ~32K SNPs
read.table("../miss0.25ind999.relatedness2", header=T)->miss_0.25ind999
paste(miss_0.25ind999$INDV1, miss_0.25ind999$INDV2)->miss_0.25ind999$pair
read.table("../miss0.25ind995.relatedness2", header=T)->miss_0.25ind995
paste(miss_0.25ind995$INDV1, miss_0.25ind995$INDV2)->miss_0.25ind995$pair
read.table("../miss0.25indall.relatedness2", header=T)->miss_0.25indall
paste(miss_0.25indall$INDV1, miss_0.25indall$INDV2)->miss_0.25indall$pair


merge(miss_0.25ind995,miss_0.25indall, by="pair")->merged_miss_995_all
merge(miss_0.25ind995,miss_0.25ind999, by="pair")->merged_miss_995_999
merge(miss_0.25ind999,miss_0.25indall, by="pair")->merged_miss_999_all

cor.test(merged_miss_995_all$RELATEDNESS_PHI.x, merged_miss_995_all$RELATEDNESS_PHI.y)
###perfectly correlated

cor.test(merged_miss_999_all$RELATEDNESS_PHI.x, merged_miss_999_all$RELATEDNESS_PHI.y)
### perfectly correlated

cor.test(merged_miss_995_999$RELATEDNESS_PHI.x, merged_miss_995_999$RELATEDNESS_PHI.y)
### perfectly correlated


### So, from this, I'm going to use the 'all'. Because the 'good' individuals give the same relatedness, and the bad individuals can be accounted for statistically in a model?

### interpretation from this tutorial: https://www.kingrelatedness.com/manual.shtml
#Close relatives can be inferred fairly reliably based on the estimated kinship coefficients as shown in the following simple algorithm: 
#an estimated kinship coefficient range >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively


### add some interpretation to the data
length(unique(miss_0.25indall$INDV1)) ## all 1617 individuals

miss_0.25indall$interp<-ifelse(miss_0.25indall$RELATEDNESS_PHI>=0.354, "duplicate/twin", 
                               ifelse(0.177<=miss_0.25indall$RELATEDNESS_PHI & miss_0.25indall$RELATEDNESS_PHI <0.354, "first_degree", 
                                      ifelse(0.0884<=miss_0.25indall$RELATEDNESS_PHI & miss_0.25indall$RELATEDNESS_PHI<0.177, "second_degree", 
                                             ifelse(0.0442<=miss_0.25indall$RELATEDNESS_PHI & miss_0.25indall$RELATEDNESS_PHI<0.0884, "third_degree", "unrelated"))))


###with the cleaner data

miss_0.25ind995$interp<-ifelse(miss_0.25ind995$RELATEDNESS_PHI>=0.354, "duplicate/twin", 
                               ifelse(0.177<=miss_0.25ind995$RELATEDNESS_PHI & miss_0.25ind995$RELATEDNESS_PHI <0.354, "first_degree", 
                                      ifelse(0.0884<=miss_0.25ind995$RELATEDNESS_PHI & miss_0.25ind995$RELATEDNESS_PHI<0.177, "second_degree", 
                                             ifelse(0.0442<=miss_0.25ind995$RELATEDNESS_PHI & miss_0.25ind995$RELATEDNESS_PHI<0.0884, "third_degree", "unrelated"))))
### way, way fewer NAs in this case (only 49) because the data here are much, much better. 

### For an A matrix, can double RELATEDNESS_PHI. 

### let's look at who's represented in the NAs

unique(miss_0.25ind995[which(as.factor(miss_0.25ind995$RELATEDNESS_PHI)== "NaN"),1]) ### there are 7 individuals here
unique(miss_0.25indall[which(as.factor(miss_0.25indall$RELATEDNESS_PHI)== "NaN"),1]) ### there are 426 individuals here
unique(miss_0.25ind999[which(as.factor(miss_0.25ind999$RELATEDNESS_PHI)== "NaN"),1]) ### there are 77 individuals here
length(unique(miss_0.25ind995$INDV1)) ### 1078 individuals, relatedness for 1071
length(unique(miss_0.25indall$INDV1)) ### 1617 individuals, relatedness for 1191
length(unique(miss_0.25ind999$INDV1)) ### 1243 individuals, relatedness for 1166 ### this is the best data set for us to use, too few snps are being used to assume very close relatedness in the all data set. 

write.csv(miss_0.25ind999, file="pero_relatedness.csv", row.names = FALSE)

###look at duplicates
miss_0.25indall[which(miss_0.25indall$INDV1 != miss_0.25indall$INDV2 & as.factor(miss_0.25indall$interp) =="duplicate/twin") , ]
### lots of these probably have too few snps to be able to say this!

#### when I do it with better data, there are no cases like this, with too few SNPs leading to duplicates. 
miss_0.25ind999[which(miss_0.25ind999$INDV1 != miss_0.25ind999$INDV2 & as.factor(miss_0.25ind999$interp) =="duplicate/twin") , ]

miss_0.25ind995[which(miss_0.25ind995$INDV1 != miss_0.25ind995$INDV2 & as.factor(miss_0.25ind995$interp) =="duplicate/twin") , ]
### there's only one dodgy case in this data set. 
