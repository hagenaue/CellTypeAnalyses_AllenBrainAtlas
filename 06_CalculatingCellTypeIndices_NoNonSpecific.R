#These code documents include all of the code used for the cell type analyses performed on the Allen Brain Atlas 160-region Agilent microarray dataset in the manuscript
##Megan Hagenauer, Alek Pankonin, and Isabelle Birt
##October 27 2016

#*****************************************

#Calculating cell type indices using cell type specific gene lists from the primary cell types after first removing the genes that are listed as "specifically-expressed" in gene lists from different primary cell types (i.e., non-specific genes!)


#I recycled some of the old Affymetrix code quantifying overlap between indices in previous analysis:

AllRegions_OneMatrix_SameOrder_Expression_CellType
AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary

#Note: Some of these indices overlap.

#Making a storage matrix to store information about overlap between primary indices:
 CellTypeSpecificGenes_Master3_Overlap<-matrix(0, length(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)), length(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)) )
 
 colnames(CellTypeSpecificGenes_Master3_Overlap)<-names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary))
row.names(CellTypeSpecificGenes_Master3_Overlap)<-names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary))

#Quantifying overlap between primary cell type indices:
for(i in 1: length(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary))){
	for(j in 1: length(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary))){

CellTypeSpecificGenes_Master3_Overlap[i,j]<-sum(AllRegions_OneMatrix_SameOrder_Expression_CellType[AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary==names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)[i]), 4]%in%AllRegions_OneMatrix_SameOrder_Expression_CellType[AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary==names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)[j]), 4])/length(AllRegions_OneMatrix_SameOrder_Expression_CellType[AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary==names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)[i]), 4])

 }
}

write.csv(CellTypeSpecificGenes_Master3_Overlap, "CellTypeSpecificGenes_Master3_Overlap.csv")



#What happens if we eliminate overlap between primary categories and then make master indices:

dim(AllRegions_OneMatrix_SameOrder_Expression_CellType)
# [1] 4013 1028

#Making an empty first row for the storage matrix:
AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap<-matrix(0, 1, (length(AllRegions_OneMatrix_SameOrder_Expression_CellType[1,])))
colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap)<-colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType)

for(i in 1: length(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary))){

#Choosing all data for a particular primary cell type:
TempCurrentIndexAllInfo<-AllRegions_OneMatrix_SameOrder_Expression_CellType[AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary==names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)[i]), ] 

#All of the gene symbols within the current primary cell type:
TempCurrentIndex<-AllRegions_OneMatrix_SameOrder_Expression_CellType[AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary==names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)[i]), 4] 

#All of the gene symbols within all other primary cell types:
TempAllOtherIndices<-AllRegions_OneMatrix_SameOrder_Expression_CellType[AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary%in%names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary)[-i]), 4]

#Grabs only rows of data with gene symbols not found in other primary cell type indices:
AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap<-rbind(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap, TempCurrentIndexAllInfo[(TempCurrentIndex%in%TempAllOtherIndices)==F,])

}

dim(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap)
# [1] 3317 1028

#removing that one dummy row:
AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap<-AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[-1,]

dim(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap)
# [1] 3316 1028
#After averaging by gene symbol: [1] 1658 1028


write.csv(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap, "AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap.csv")

CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap<-table(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$CellType_Primary)

write.csv(CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap, "CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap.csv")


AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$CellType_Primary)), ncol=(length(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean)<-names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$CellType_Primary))

colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean)<-colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

#Old version of code:
# for(i in c(15:length(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[1,]))){
# AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,i-14]<-tapply(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[,i], AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$CellType_Primary, mean)
# }


#I went back and changed this so that it averaged by tag first, then by primary cell category, to match what we did in the Affymetrix data. In general, I find this produces slightly more valid predictions (because cell type specific gene lists from publications with less strict statistical cut-offs don't dominate the results)

temp<-data.frame(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$Tag, AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$CellType_Primary) 
dim(temp)
# [1] 3316    2

CellTypePrimaryVsTag<-unique(temp)
dim(CellTypePrimaryVsTag)
# [1] 38  2

colnames(CellTypePrimaryVsTag)<-c("Tag","CellType_Primary")
head(CellTypePrimaryVsTag)
                                   # Tag CellType_Primary
# 1     Astrocyte_All_Darmanis_PNAS_2015        Astrocyte
# 50     Astrocyte_All_Cahoy_JNeuro_2008        Astrocyte
# 144    Astrocyte_All_Zhang_JNeuro_2014        Astrocyte
# 187      Astrocyte_All_Doyle_Cell_2008        Astrocyte
# 211  Astrocyte_All_Zeisel_Science_2015        Astrocyte
# 442 Endothelial_All_Darmanis_PNAS_2015      Endothelial

AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$Tag))

colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

for(i in c(15:length(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[1,]))){
AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[,i], AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$Tag, mean)
}

head(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_MeanTag)

temp2<-tapply(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap[,i], AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap$Tag, mean)

Tag<-names(temp2)
CellTypePrimaryVsTag2<-join(as.data.frame(Tag), as.data.frame(CellTypePrimaryVsTag), by="Tag")


for(i in c(1:length(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_MeanTag[1,]))){
AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,i]<-tapply(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_MeanTag[,i], CellTypePrimaryVsTag2[,2], mean)
}

head(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean)


write.csv(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean, "AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean.csv")


is.numeric(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean)


#Outputting a hierarchically-clustered heatmap illustrating the correlations between cell type indices & a standard correlation matrix:

png("Heatmap_CorMatrixPrimaryCellsNoOverlap.png")
heatmap(cor(t(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F])))
dev.off()

CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix<-cor(t(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]))

head(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix)

write.csv(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix, "CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv")



#New violin plots by region:

for(i in c(1:length(names(table(Sample_Region))))){

png(paste(names(table(Sample_Region))[i], "_vsPrimaryCellType.png", sep=""), width=1300, height=750)

Temp<-AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,Sample_Region==names(table(Sample_Region))[i]]

RelativeExpression<-as.numeric(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,Sample_Region==names(table(Sample_Region))[i]])

CellType<-rep(row.names(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean), ncol(Temp))
	myData<-data.frame(RelativeExpression,CellType)

p <- ggplot(myData, aes(x=CellType, y=RelativeExpression, fill=CellType)) + 
  geom_violin()
  
print(p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"), axis.text.x = element_text(angle = 30, hjust = 1)))
  
dev.off()

}

#Making higher resolution plots for figures:
#Regions: Choroid plexus [22], Corpus Callosum [29], Central Glial Substance [17], Dentate Gyrus [35], Globus Pallidus[49, 50]

i<-50

pdf(paste(names(table(Sample_Region))[i], "_vsPrimaryCellType.pdf", sep=""), width=8, height=5)

Temp<-AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,Sample_Region==names(table(Sample_Region))[i]]

RelativeExpression<-as.numeric(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,Sample_Region==names(table(Sample_Region))[i]])

CellType<-rep(row.names(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean), ncol(Temp))
myData<-data.frame(RelativeExpression,CellType)

p <- ggplot(myData, aes(x=CellType, y=RelativeExpression, fill=CellType)) + 
  geom_violin()

print(p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)+theme(axis.text=element_text(size=16),
                                                                        axis.title=element_text(size=16,face="bold"), axis.text.x = element_text(angle = 30, hjust = 1)))

dev.off()



#Outputting a storage matrix of average primary cell type index and standard error per region:

MeanPrimaryCellTypeNoOverlap_ByRegion<-matrix(0, length(AVE_Expression_CellType_Primary_bySampleNoNA[,1]), length(names(table(Region))))
SEPrimaryCellTypeNoOverlap_ByRegion<-matrix(0, length(AVE_Expression_CellType_Primary_bySampleNoNA[,1]), length(names(table(Region))))

row.names(MeanPrimaryCellTypeNoOverlap_ByRegion)<-row.names(AVE_Expression_CellType_Primary_bySampleNoNA)
row.names(SEPrimaryCellTypeNoOverlap_ByRegion)<-row.names(AVE_Expression_CellType_Primary_bySampleNoNA)

colnames(MeanPrimaryCellTypeNoOverlap_ByRegion)<-names(table(Region))
colnames(SEPrimaryCellTypeNoOverlap_ByRegion)<-names(table(Region))

for(i in c(1:length(AVE_Expression_CellType_Primary_bySampleNoNA[,1]))){

MeanPrimaryCellTypeNoOverlap_ByRegion[i,]<-tapply(AVE_Expression_CellType_Primary_bySampleNoNA[i,], Region, function(y) mean(y, na.rm=T))
SEPrimaryCellTypeNoOverlap_ByRegion[i,]<-tapply(AVE_Expression_CellType_Primary_bySampleNoNA[i,], Region, function(y) sd(y, na.rm=T)/sqrt(length(AVE_Expression_CellType_Primary_bySampleNoNA[i,])))

}

write.csv(MeanPrimaryCellTypeNoOverlap_ByRegion, "MeanPrimaryCellTypeNoOverlap_ByRegion.csv")
write.csv(SEPrimaryCellTypeNoOverlap_ByRegion, "SEPrimaryCellTypeNoOverlap_ByRegion.csv")



row.names(AVE_Expression_CellType_Primary_bySampleNoNA)
 [1] "Astrocyte"                "Endothelial"              "Microglia"                "Mural"                   
 [5] "Neuron_All"               "Neuron_Interneuron"       "Neuron_Projection"        "Oligodendrocyte"         
 [9] "Oligodendrocyte_Immature" "RBC"






