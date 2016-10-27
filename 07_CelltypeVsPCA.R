#These code documents include all of the code used for the cell type analyses performed on the Allen Brain Atlas 160-region Agilent microarray dataset in the manuscript
##Megan Hagenauer, Alek Pankonin, and Isabelle Birt
##October 27 2016

#*****************************************

##PCA: Create a new object that identifies the principal components within the normalized signal data:

pca96<-prcomp(t(AllRegions_OneMatrix_SameOrder_Zscores[,is.na(AllRegions_OneMatrix_SameOrder_Zscores[1,])==F]))
PCA_NormExpression<-pca96$x[,c(1:10)]
rownames(PCA_NormExpression)<-colnames(AllRegions_OneMatrix_SameOrder_Zscores[,is.na(AllRegions_OneMatrix_SameOrder_Zscores[1,])==F]) 

tmp<-cbind(AllRegions_OneMatrix_SameOrder_SampleInfo[is.na(AllRegions_OneMatrix_SameOrder_Zscores[1,])==F,], PCA_NormExpression)
print("PCA results outputted as PC12.txt, dimensions:")
print(dim(tmp))
write.table(tmp,"PC12.txt",sep="\t")


##Looking at relationships between cell type indices and PCA:

for(i in 1:10){
png(paste("PC", i, "_vsRegion.png", sep=""), width=6000, height=750)

PC<-pca96$x[,i]

Region<-Sample_Region[is.na(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]

myData<-data.frame(PC,Region)

p <- ggplot(myData, aes(x=Region, y=PC, fill=3)) + 
  geom_violin()
  
print(p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+stat_summary(fun.y=median, geom="point", size=2, color="red")+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
dev.off()
}




#This first correlation matrix is for primary cell type indices without removing overlapping (non-specific) genes:

PCAvsCellType_corMatrix<-matrix(NA, length(AVE_Expression_CellType_Primary_bySample[,1]), 20)
row.names(PCAvsCellType_corMatrix)<-row.names(AVE_Expression_CellType_Primary_bySample)


for(i in c(1:length(AVE_Expression_CellType_Primary_bySample[,1]))){	
	temp<-AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]
	for (j in c(1:20)){
		PCAvsCellType_corMatrix[i,j]<-cor(pca96$x[,j], temp[i,])
	}
}

write.csv(PCAvsCellType_corMatrix, "PCAvsCellType_corMatrix.csv")
#write.csv(PCAvsCellType_corMatrix, "ExpressionPCAvsCellType_corMatrix.csv")



PCAvsCellTypeTag_corMatrix<-matrix(NA, length(AVE_Expression_CellType_Tag_bySample[,1]), 20)
row.names(PCAvsCellTypeTag_corMatrix)<-row.names(AVE_Expression_CellType_Tag_bySample)


for(i in c(1:length(AVE_Expression_CellType_Tag_bySample[,1]))){
	temp<-AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]	
	for (j in c(1:20)){
		PCAvsCellTypeTag_corMatrix[i,j]<-cor(pca96$x[,j], temp[i,])
	}
}

write.csv(PCAvsCellTypeTag_corMatrix, "PCAvsCellTypeTag_corMatrix.csv")
#write.csv(PCAvsCellTypeTag_corMatrix, "ExpressionPCAvsCellTypeTag_corMatrix.csv")


#This second correlation matrix is for primary cell type indices *after* removing overlapping (non-specific) genes:

PCAvsCellType_corMatrix<-matrix(NA, length(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,1]), 20)
row.names(PCAvsCellType_corMatrix)<-row.names(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean)


for(i in c(1:length(AVE_Expression_CellType_Primary_bySample[,1]))){	
	temp<-AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]
	for (j in c(1:20)){
		PCAvsCellType_corMatrix[i,j]<-cor(pca96$x[,j], temp[i,])
	}
}

write.csv(PCAvsCellType_corMatrix, "PCAvsCellType_corMatrixNoOverlap.csv")
#write.csv(PCAvsCellType_corMatrix, "ExpressionPCAvsCellType_corMatrix.csv")




row.names(AVE_Expression_CellType_Primary_bySample)
[1] "Astrocyte"                "Endothelial"              "Microglia"                "Mural"                   
 [5] "Neuron_All"               "Neuron_Interneuron"       "Neuron_Projection"        "Oligodendrocyte"         
 [9] "Oligodendrocyte_Immature" "RBC"




#I ran two versions of this:
#the first version of this output includes cell type indices without removing overlapping genes:
# AVE_Expression_CellType_Primary_bySampleNoNA<-AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]
#the second version of this output includes cell type indices *after* removing overlapping genes:
	AVE_Expression_CellType_Primary_bySampleNoNA<-AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]



#sink('PCAvsSubjectVarCellTypeNoOverlap.txt')
sink('PCAvsSubjectVarCellType.txt')

print(row.names(AVE_Expression_CellType_Primary_bySampleNoNA))

print("PC1")

summary.lm(lm(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[1,]+AVE_Expression_CellType_Primary_bySampleNoNA[2,]+AVE_Expression_CellType_Primary_bySampleNoNA[3,]+AVE_Expression_CellType_Primary_bySampleNoNA[4,]+AVE_Expression_CellType_Primary_bySampleNoNA[5,]+AVE_Expression_CellType_Primary_bySampleNoNA[6,]+AVE_Expression_CellType_Primary_bySampleNoNA[7,]+AVE_Expression_CellType_Primary_bySampleNoNA[8,]+AVE_Expression_CellType_Primary_bySampleNoNA[9,]+AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

print("PC2")

summary.lm(lm(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[1,]+AVE_Expression_CellType_Primary_bySampleNoNA[2,]+AVE_Expression_CellType_Primary_bySampleNoNA[3,]+AVE_Expression_CellType_Primary_bySampleNoNA[4,]+AVE_Expression_CellType_Primary_bySampleNoNA[5,]+AVE_Expression_CellType_Primary_bySampleNoNA[6,]+AVE_Expression_CellType_Primary_bySampleNoNA[7,]+AVE_Expression_CellType_Primary_bySampleNoNA[8,]+AVE_Expression_CellType_Primary_bySampleNoNA[9,]+AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

print("PC3")

summary.lm(lm(pca96$x[,3]~AVE_Expression_CellType_Primary_bySampleNoNA[1,]+AVE_Expression_CellType_Primary_bySampleNoNA[2,]+AVE_Expression_CellType_Primary_bySampleNoNA[3,]+AVE_Expression_CellType_Primary_bySampleNoNA[4,]+AVE_Expression_CellType_Primary_bySampleNoNA[5,]+AVE_Expression_CellType_Primary_bySampleNoNA[6,]+AVE_Expression_CellType_Primary_bySampleNoNA[7,]+AVE_Expression_CellType_Primary_bySampleNoNA[8,]+AVE_Expression_CellType_Primary_bySampleNoNA[9,]+AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

print("PC4")

summary.lm(lm(pca96$x[,4]~AVE_Expression_CellType_Primary_bySampleNoNA[1,]+AVE_Expression_CellType_Primary_bySampleNoNA[2,]+AVE_Expression_CellType_Primary_bySampleNoNA[3,]+AVE_Expression_CellType_Primary_bySampleNoNA[4,]+AVE_Expression_CellType_Primary_bySampleNoNA[5,]+AVE_Expression_CellType_Primary_bySampleNoNA[6,]+AVE_Expression_CellType_Primary_bySampleNoNA[7,]+AVE_Expression_CellType_Primary_bySampleNoNA[8,]+AVE_Expression_CellType_Primary_bySampleNoNA[9,]+AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

print("PC5")

summary.lm(lm(pca96$x[,5]~AVE_Expression_CellType_Primary_bySampleNoNA[1,]+AVE_Expression_CellType_Primary_bySampleNoNA[2,]+AVE_Expression_CellType_Primary_bySampleNoNA[3,]+AVE_Expression_CellType_Primary_bySampleNoNA[4,]+AVE_Expression_CellType_Primary_bySampleNoNA[5,]+AVE_Expression_CellType_Primary_bySampleNoNA[6,]+AVE_Expression_CellType_Primary_bySampleNoNA[7,]+AVE_Expression_CellType_Primary_bySampleNoNA[8,]+AVE_Expression_CellType_Primary_bySampleNoNA[9,]+AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

print("PC6")

summary.lm(lm(pca96$x[,6]~AVE_Expression_CellType_Primary_bySampleNoNA[1,]+AVE_Expression_CellType_Primary_bySampleNoNA[2,]+AVE_Expression_CellType_Primary_bySampleNoNA[3,]+AVE_Expression_CellType_Primary_bySampleNoNA[4,]+AVE_Expression_CellType_Primary_bySampleNoNA[5,]+AVE_Expression_CellType_Primary_bySampleNoNA[6,]+AVE_Expression_CellType_Primary_bySampleNoNA[7,]+AVE_Expression_CellType_Primary_bySampleNoNA[8,]+AVE_Expression_CellType_Primary_bySampleNoNA[9,]+AVE_Expression_CellType_Primary_bySampleNoNA[10,]))

sink()


png("PC1VsEndothelial.png")
plot(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[2,], ylab="PC1", xlab="Endothelial Index")
BestFitLine<-lm(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[2,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()

png("PC1VsOligodendrocyte.png")
plot(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[8,], ylab="PC1", xlab="Oligodendrocyte Index")
BestFitLine<-lm(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[8,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()

png("PC1VsMural.png")
plot(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[4,], ylab="PC1", xlab="Mural Index")
BestFitLine<-lm(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[4,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()

png("PC1VsRBC.png")
plot(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[10,], ylab="PC1", xlab="RBC Index")
BestFitLine<-lm(pca96$x[,1]~AVE_Expression_CellType_Primary_bySampleNoNA[10,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()


png("PC2VsNeuron.png")
plot(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[5,], ylab="PC2", xlab="Neuron_All Index")
BestFitLine<-lm(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[5,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()

#Wow - there are some serious outliers mixed in there for the Neuron_All index. Weird.


png("PC2VsProjectionNeuron.png")
plot(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[7,], ylab="PC2", xlab="Neuron_Projection Index")
BestFitLine<-lm(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[7,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()


png("PC2VsOPC.png")
plot(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[9,], ylab="PC2", xlab="OPC Index")
BestFitLine<-lm(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[9,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()

png("PC2VsEndothelial.png")
plot(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[2,], ylab="PC2", xlab="Endothelial Index")
BestFitLine<-lm(pca96$x[,2]~AVE_Expression_CellType_Primary_bySampleNoNA[2,])
abline(BestFitLine, col=2)
mtext(paste("R-Squared: ", round(summary.lm(BestFitLine)$r.squared, digits=2), sep=""), line=0)
mtext(paste("P-Value< ", summary.lm(BestFitLine)$coefficients[2,4], sep=""), line=1)
dev.off()



Region<-Sample_Region[is.na(AllRegions_OneMatrix_SameOrder_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]

png("PC1vsPC2_ByRegion.png")
plot(pca96$x[,2]~pca96$x[,1], col=Region, pch=as.numeric(Region)/8)
dev.off()

PCA_MeanByRegion<-matrix(0, length(names(table(Region))), 20)
row.names(PCA_MeanByRegion)<-names(table(Region))

for(i in 1:20){
	PCA_MeanByRegion[,i]<-tapply(X=pca96$x[,i], Region, FUN=function(y) mean(y, na.rm=T))
}

write.csv(PCA_MeanByRegion, "PCA_MeanByRegion.csv")





