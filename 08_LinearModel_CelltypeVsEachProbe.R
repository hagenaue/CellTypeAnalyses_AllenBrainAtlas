#These code documents include all of the code used for the cell type analyses performed on the Allen Brain Atlas 160-region Agilent microarray dataset in the manuscript
##Megan Hagenauer, Alek Pankonin, and Isabelle Birt
##October 27 2016

#*****************************************

#Examining the relationship between each probe and the cell type indices

#First I examined the multicollinearity between the different cell type indices using a variance inflation factor calculation:
library(usdm)

vif(data.frame(AVE_Expression_CellType_Primary_bySampleNoNA[1,], AVE_Expression_CellType_Primary_bySampleNoNA[2,], AVE_Expression_CellType_Primary_bySampleNoNA[3,], AVE_Expression_CellType_Primary_bySampleNoNA[4,], AVE_Expression_CellType_Primary_bySampleNoNA[5,], AVE_Expression_CellType_Primary_bySampleNoNA[6,], AVE_Expression_CellType_Primary_bySampleNoNA[7,], AVE_Expression_CellType_Primary_bySampleNoNA[8,], AVE_Expression_CellType_Primary_bySampleNoNA[9,], AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

# #After averaging by gene symbol, re-scaling, then index, no overlap:
                                            # Variables      VIF
# 1   AVE_Expression_CellType_Primary_bySampleNoNA.1... 3.413506
# 2   AVE_Expression_CellType_Primary_bySampleNoNA.2... 7.300151
# 3   AVE_Expression_CellType_Primary_bySampleNoNA.3... 5.364135
# 4   AVE_Expression_CellType_Primary_bySampleNoNA.4... 9.723203
# 5   AVE_Expression_CellType_Primary_bySampleNoNA.5... 4.970093
# 6   AVE_Expression_CellType_Primary_bySampleNoNA.6... 7.547842
# 7   AVE_Expression_CellType_Primary_bySampleNoNA.7... 2.745026
# 8   AVE_Expression_CellType_Primary_bySampleNoNA.8... 8.054152
# 9   AVE_Expression_CellType_Primary_bySampleNoNA.9... 6.582484
# 10 AVE_Expression_CellType_Primary_bySampleNoNA.10... 1.706381


#Outputting the relationship of each probe with cell type index (untrimmed model):

temp<-data.frame(ProbeInfo[,5], AllRegions_OneMatrix_SameOrder_Zscores, stringsAsFactors=F)

colnames(temp)[1]<-"GeneSymbol_Human"
dim(temp)
temp[c(1:3), c(1:3)]

AllRegions_OneMatrix_SameOrder_ZscoresNoNA<-as.matrix(AllRegions_OneMatrix_SameOrder_Zscores[,is.na(AllRegions_OneMatrix_SameOrder_Zscores[1,])==F])

LM_AllPrimaryCellTypes_Betas<-matrix(0, length(AllRegions_OneMatrix_SameOrder_ZscoresNoNA[,1]), 11)
LM_AllPrimaryCellTypes_Pval<-matrix(0, length(AllRegions_OneMatrix_SameOrder_ZscoresNoNA[,1]), 11)
LM_AllPrimaryCellTypes_SE<-matrix(0, length(AllRegions_OneMatrix_SameOrder_ZscoresNoNA[,1]), 11)
LM_AllPrimaryCellTypes_T<-matrix(0, length(AllRegions_OneMatrix_SameOrder_ZscoresNoNA[,1]), 11)
LM_AllPrimaryCellTypes_Rsquared<-matrix(0, length(AllRegions_OneMatrix_SameOrder_ZscoresNoNA[,1]), 2)

colnames(LM_AllPrimaryCellTypes_Betas)<-c("Intercept", row.names(AVE_Expression_CellType_Primary_bySampleNoNA))
colnames(LM_AllPrimaryCellTypes_Pval)<-c("Intercept", row.names(AVE_Expression_CellType_Primary_bySampleNoNA))
colnames(LM_AllPrimaryCellTypes_SE)<-c("Intercept", row.names(AVE_Expression_CellType_Primary_bySampleNoNA))
colnames(LM_AllPrimaryCellTypes_T)<-c("Intercept", row.names(AVE_Expression_CellType_Primary_bySampleNoNA))
colnames(LM_AllPrimaryCellTypes_Rsquared)<-c("Rsquared", "AdjRsquared")


for(i in 1:length(AllRegions_OneMatrix_SameOrder_ZscoresNoNA[,1])){

temp<-summary.lm(lm(AllRegions_OneMatrix_SameOrder_ZscoresNoNA[i,]~AVE_Expression_CellType_Primary_bySampleNoNA[1,]+AVE_Expression_CellType_Primary_bySampleNoNA[2,]+AVE_Expression_CellType_Primary_bySampleNoNA[3,]+AVE_Expression_CellType_Primary_bySampleNoNA[4,]+AVE_Expression_CellType_Primary_bySampleNoNA[5,]+AVE_Expression_CellType_Primary_bySampleNoNA[6,]+AVE_Expression_CellType_Primary_bySampleNoNA[7,]+AVE_Expression_CellType_Primary_bySampleNoNA[8,]+AVE_Expression_CellType_Primary_bySampleNoNA[9,]+AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

LM_AllPrimaryCellTypes_Betas[i,]<-temp$coefficients[,1]
LM_AllPrimaryCellTypes_SE[i,]<-temp$coefficients[,2]
LM_AllPrimaryCellTypes_T[i,]<-temp$coefficients[,3]
LM_AllPrimaryCellTypes_Pval[i,]<-temp$coefficients[,4]
LM_AllPrimaryCellTypes_Rsquared[i,1]<-temp$r.squared
LM_AllPrimaryCellTypes_Rsquared[i,2]<-temp$adj.r.squared
}

#Version before averaging by gene symbol:
# write.csv(data.frame(ProbeInfo[,5],LM_AllPrimaryCellTypes_Betas), "LM_AllPrimaryCellTypes_Betas.csv")
# write.csv(data.frame(ProbeInfo[,5],LM_AllPrimaryCellTypes_Pval), "LM_AllPrimaryCellTypes_Pval.csv")
# write.csv(data.frame(ProbeInfo[,5],LM_AllPrimaryCellTypes_SE), "LM_AllPrimaryCellTypes_SE.csv")
# write.csv(data.frame(ProbeInfo[,5],LM_AllPrimaryCellTypes_T), "LM_AllPrimaryCellTypes_T.csv")
# write.csv(data.frame(ProbeInfo[,5],LM_AllPrimaryCellTypes_Rsquared), "LM_AllPrimaryCellTypes_Rsquared.csv")
# write.csv(data.frame(ProbeInfo[,5],LM_AllPrimaryCellTypes_Betas,LM_AllPrimaryCellTypes_Pval), "LM_AllPrimaryCellTypes_DF.csv")

#Version after averaging by gene symbol:

write.csv(data.frame(row.names(AllRegions_OneMatrix_SameOrder_ZscoresNoNA),LM_AllPrimaryCellTypes_Betas), "LM_AllPrimaryCellTypes_Betas.csv")
write.csv(data.frame(row.names(AllRegions_OneMatrix_SameOrder_ZscoresNoNA),LM_AllPrimaryCellTypes_Pval), "LM_AllPrimaryCellTypes_Pval.csv")
write.csv(data.frame(row.names(AllRegions_OneMatrix_SameOrder_ZscoresNoNA),LM_AllPrimaryCellTypes_SE), "LM_AllPrimaryCellTypes_SE.csv")
write.csv(data.frame(row.names(AllRegions_OneMatrix_SameOrder_ZscoresNoNA),LM_AllPrimaryCellTypes_T), "LM_AllPrimaryCellTypes_T.csv")
write.csv(data.frame(row.names(AllRegions_OneMatrix_SameOrder_ZscoresNoNA),LM_AllPrimaryCellTypes_Rsquared), "LM_AllPrimaryCellTypes_Rsquared.csv")
write.csv(data.frame(row.names(AllRegions_OneMatrix_SameOrder_ZscoresNoNA),LM_AllPrimaryCellTypes_Betas,LM_AllPrimaryCellTypes_Pval), "LM_AllPrimaryCellTypes_DF.csv")


LM_AllPrimaryCellTypes_DF<-data.frame(row.names(AllRegions_OneMatrix_SameOrder_ZscoresNoNA),LM_AllPrimaryCellTypes_Betas,LM_AllPrimaryCellTypes_Pval)


colnames(LM_AllPrimaryCellTypes_DF)
 [1] "row.names.AllRegions_OneMatrix_SameOrder_ZscoresNoNA." "Intercept"                                            
 [3] "Astrocyte"                                             "Endothelial"                                          
 [5] "Microglia"                                             "Mural"                                                
 [7] "Neuron_All"                                            "Neuron_Interneuron"                                   
 [9] "Neuron_Projection"                                     "Oligodendrocyte"                                      
[11] "Oligodendrocyte_Immature"                              "RBC"                                                  
[13] "Intercept.1"                                           "Astrocyte.1"                                          
[15] "Endothelial.1"                                         "Microglia.1"                                          
[17] "Mural.1"                                               "Neuron_All.1"                                         
[19] "Neuron_Interneuron.1"                                  "Neuron_Projection.1"                                  
[21] "Oligodendrocyte.1"                                     "Oligodendrocyte_Immature.1"                           
[23] "RBC.1" 

Top100GenesAssociatedWEachCelltype<-LM_AllPrimaryCellTypes_DF[1,]

for(i in c(3:12)){
temp<-LM_AllPrimaryCellTypes_DF[LM_AllPrimaryCellTypes_DF[,i]>0,]
temp2<-temp[order(temp[,(i+11)]),]
Top100GenesAssociatedWEachCelltype<-rbind(Top100GenesAssociatedWEachCelltype, temp2[c(1:100),])
}

Top100GenesAssociatedWEachCelltype<-Top100GenesAssociatedWEachCelltype[-1,]
head(Top100GenesAssociatedWEachCelltype)

sum(Top100GenesAssociatedWEachCelltype[,1]%in%CellTypeSpecificGenes_Master3[,4])
#[1] 222

InCellTypeSpecific_Master3<-as.numeric(Top100GenesAssociatedWEachCelltype[,1]%in%CellTypeSpecificGenes_Master3[,4])

write.csv(data.frame(Top100GenesAssociatedWEachCelltype, InCellTypeSpecific_Master3), "Top100GenesAssociatedWEachCelltype.csv")

#Note: I'm pretty sure that my column labels are out of order and the interneurons and projection neurons are reversed.
#This does not seem to be the case for other output. E.g. By Gene Symbol, By Gene Symbol no overlap.
#... but suddenly appears in the output for By GeneSymbolNoOverlapAveIndex
#....this is true if I re-run the analysis using again just gene symbol no overlap without averaging by index, so there is something happening when averaging by index that is making the results funky.
#Ok - taken care of now, and all of the code has been re-run.

