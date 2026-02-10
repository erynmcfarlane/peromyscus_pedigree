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
