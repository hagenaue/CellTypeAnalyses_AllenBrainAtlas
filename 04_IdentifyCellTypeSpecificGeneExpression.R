#These code documents include all of the code used for the cell type analyses performed on the Allen Brain Atlas 160-region Agilent microarray dataset in the manuscript
##Megan Hagenauer, Alek Pankonin, and Isabelle Birt
##October 27 2016

#*****************************************

CellTypeSpecificGenes_Master3<-read.csv("CellTypeSpecificGenes_Master3.csv", header=T)

colnames(CellTypeSpecificGenes_Master3)

 [1] "Umbrella.Cell.Type"    "Specific.Cell.Type"    "Brain.Region"          "Gene.Symbol..Human."  
 [5] "Gene.Symbol..Mouse."   "Species"               "Age"                   "Statistical.Criterion"
 [9] "Specificity"           "Comparison"            "Platform"              "Citation"             
[13] "Tag"                   "CellType_Primary" 

table(CellTypeSpecificGenes_Master3$CellType_Primary)

colnames(CellTypeSpecificGenes_Master3)[4]<-"GeneSymbol_Human"

CellTypeSpecificGenes_Master3[,4]<-as.character(CellTypeSpecificGenes_Master3[,4])


sum(is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human))
#[1] 364

CellTypeSpecificGenes_Master3NoNA<-CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human)==F,]

#Old version:
# temp<-data.frame(ProbeInfo[,5], AllRegions_OneMatrix_SameOrder_Zscores, stringsAsFactors=F)

#New version (after averaging by gene symbol):
temp<-data.frame(row.names(AllRegions_OneMatrix_SameOrder_Zscores), AllRegions_OneMatrix_SameOrder_Zscores, stringsAsFactors=F)

colnames(temp)[1]<-"GeneSymbol_Human"

sum(is.na(temp[,1]))
#[1] 0

#Unnecessary for this application:
# temp<-temp[is.na(temp[,1])==F,]

sum(temp[,1] %in% CellTypeSpecificGenes_Master3[,4])
# Old version:
# [1] 3137

#New version: (after averaging by gene symbol):
# [1] 1608

sum(CellTypeSpecificGenes_Master3[,4]  %in%  temp[,1])
# [1] 1984


#Note: NAs were causing a serious problem with this join function.  Fixed now. :)
library(plyr)
AllRegions_OneMatrix_SameOrder_Expression_CellType<-join(CellTypeSpecificGenes_Master3, temp, by="GeneSymbol_Human", type="inner")
dim(AllRegions_OneMatrix_SameOrder_Expression_CellType)
# [1] 4013 1028

#In the old version: It was making all possible combinations - so on average, the gene symbols from CellTypeSpecificGenes_Master3 align with an average of 3 probes apiece in the ABA Agilent dataset, and some of the cell type specific probes are found in more than one index.
#from glancing at the data using "head" this does seem to be true.

#Old code - relevant for previous dataset, but not this one.
# colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType)
# cbind(colnames(temp)[-1], colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType)[-c(1:14)] )

#After averaging the z-scores by gene symbol this became:
# [1] 1984 1028

