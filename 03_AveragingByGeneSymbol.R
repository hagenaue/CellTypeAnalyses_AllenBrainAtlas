#These code documents include all of the code used for the cell type analyses performed on the Allen Brain Atlas 160-region Agilent microarray dataset in the manuscript
##Megan Hagenauer, Alek Pankonin, and Isabelle Birt
##October 27 2016

#*****************************************

#Averaging by Gene Symbol:

# 'data.frame':	30000 obs. of  12 variables:
# str(ProbeInfo)
 # $ X             : int  2 3 4 5 6 7 8 9 10 11 ...
 # $ id            : int  1038360 1053962 1053210 1031532 1070322 1032754 1065870 1053653 1027786 1053652 ...
 # $ name          : Factor w/ 30000 levels "A_23_P100001",..: 6470 4308 27699 6733 18473 17522 1497 23129 20825 3800 ...
 # $ gene.id       : int  57211 4576 5141 122048 1011337 90407 1008239 4821 5701 4821 ...
 # $ gene.symbol   : Factor w/ 18787 levels "A_23_P100522",..: 7715 12772 14069 10627 4227 17913 114 13201 14781 13201 ...
 # $ gene.name     : Factor w/ 18769 levels "-","1-acylglycerol-3-phosphate O-acyltransferase 1 (lysophosphatidic acid acyltransferase, alpha)",..: 7857 17898 13847 13669 4592 17871 479 12381 13947 12381 ...
 # $ entrez.id     : int  81501 4602 5173 386618 0 257313 0 4852 5733 4852 ...
 # $ chromosome    : Factor w/ 29 levels "1","10","11",..: 23 21 14 5 25 17 25 22 1 22 ...
 # $ start.position: Factor w/ 1 level "n/a": 1 1 1 1 1 1 1 1 1 1 ...
 # $ end.position  : Factor w/ 1 level "n/a": 1 1 1 1 1 1 1 1 1 1 ...
 # $ p             : num  3.25e-14 7.57e-40 8.95e-22 2.10e-38 1.82e-26 ...
 # $ fold.change   : num  12.3 12 11.6 11.3 10.7 ...

#Note: Gene Symbol is column 5

sum(is.na(ProbeInfo[,5]))
length(unique(ProbeInfo[,5]))
# [1] 18787
length(unique(ProbeInfo[,4]))
# [1] 18790

#Ok, In the ABA data, there can be multiple probes representing the same gene, which I think (combined with the Join function) is making it so that some of the cell type specific genes are getting more "vote" than others.  So I'm going to try averaging by gene symbol first and then re-running the code. 

AllRegions_OneMatrix_SameOrder_ZscoresbyGene<-matrix(0, length(names(table(ProbeInfo[,5]))), length(AllRegions_OneMatrix_SameOrder_Zscores[1,]))
row.names(AllRegions_OneMatrix_SameOrder_ZscoresbyGene)<-names(table(ProbeInfo[,5]))

dim(AllRegions_OneMatrix_SameOrder_ZscoresbyGene)
# [1] 18787  1014

for(i in c(1:length(AllRegions_OneMatrix_SameOrder_Zscores[1,]))){
AllRegions_OneMatrix_SameOrder_ZscoresbyGene[,i]<-tapply(AllRegions_OneMatrix_SameOrder_Zscores[,i], ProbeInfo[,5], mean)
}


#After averaging by gene symbol, the data needs to be re-scaled so that no particular gene is getting more "vote" in our later analyses than any other gene: 
temp<-t(scale(t(AllRegions_OneMatrix_SameOrder_ZscoresbyGene), center=T, scale=T))
dim(temp)
head(temp)
plot.new()
plot(sort(apply(temp, 1, function(y) mean(y, na.rm=T))))
plot(sort(apply(temp, 1, function(y) sd(y, na.rm=T))))
#good - looks like the re-scaling worked.
AllRegions_OneMatrix_SameOrder_ZscoresbyGene<-temp

head(AllRegions_OneMatrix_SameOrder_ZscoresbyGene)
tail(AllRegions_OneMatrix_SameOrder_ZscoresbyGene)

AllRegions_OneMatrix_SameOrder_Zscores<-AllRegions_OneMatrix_SameOrder_ZscoresbyGene

write.csv(AllRegions_OneMatrix_SameOrder_Zscores, "AllRegions_OneMatrix_SameOrder_Zscores.csv")
