#These code documents include all of the code used for the cell type analyses performed on the Allen Brain Atlas 160-region Agilent microarray dataset in the manuscript
##Megan Hagenauer, Alek Pankonin, and Isabelle Birt
##October 27 2016

#*****************************************

#Calculating cell type indices using cell type specific gene lists from the primary cell types: (note: this version is before removing non-specific genes)

AVE_Expression_CellType_Primary_bySample<-matrix(NA, nrow=length(names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary))), ncol=(1028-14))

row.names(AVE_Expression_CellType_Primary_bySample)<-names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary))
# colnames(AVE_Expression_CellType_Primary_bySample)<-colnames(temp)[-1]


for(i in c(15:1028)){
AVE_Expression_CellType_Primary_bySample[,(i-14)]<-tapply(AllRegions_OneMatrix_SameOrder_Expression_CellType[,i], AllRegions_OneMatrix_SameOrder_Expression_CellType$CellType_Primary, mean)
}
#There were warning messages for all columns that contained NA

head(AVE_Expression_CellType_Primary_bySample)

#Outputting a hierarchically-clustered heatmap illustrating the correlation betweeen different cell type indices based on primary cell type:
png("CorrMatrixCellTypeVsCellType_HeatMap.png")
heatmap(cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F])))
dev.off()

CorrelationMatrixCellTypeVsCellType<-cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]))

write.csv(CorrelationMatrixCellTypeVsCellType, "CorrelationMatrixCellTypeVsCellType.csv")
#Huh - all correlations are positive. Perhaps because some samples simply have less signal?  (since ABA doesn't use quantile normalization across samples)

#Calculating cell type indices using cell type specific gene lists from different publications:

AVE_Expression_CellType_Tag_bySample<-matrix(NA, nrow=length(names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$Tag))), ncol=1028-14)
row.names(AVE_Expression_CellType_Tag_bySample)<-names(table(AllRegions_OneMatrix_SameOrder_Expression_CellType$Tag))
# colnames(AVE_Expression_CellType_Tag_bySample)<-colnames(temp)[-1]

for(i in c(15:1028)){
AVE_Expression_CellType_Tag_bySample[,(i-14)]<-tapply(AllRegions_OneMatrix_SameOrder_Expression_CellType[,i], AllRegions_OneMatrix_SameOrder_Expression_CellType$Tag, mean)
}

head(AVE_Expression_CellType_Tag_bySample)

#Outputting a hierarchically-clustered heatmap illustrating the correlation betweeen different cell type indices based on different publications:
png("CorrMatrixCellIndexVsCellIndex_HeatMap.png", width=1000, height=1000)
heatmap(cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F])))
dev.off()

CorrelationMatrixCellIndexVsCellIndex<-cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]))

write.csv(CorrelationMatrixCellIndexVsCellIndex, "CorrelationMatrixCellIndexVsCellIndex.csv")

#Outputting a boxplot depicting the average of each cell type index per brain region:
Sample_Region<-AllRegions_OneMatrix_SameOrder_SampleInfo[,11]

for(i in c(1:length(AVE_Expression_CellType_Primary_bySample[,1]))){
png(paste(row.names(AVE_Expression_CellType_Primary_bySample)[i], "_vsRegion.png"), width=3000, height=480)
par(mar=c(18,3,1,1))
boxplot(AVE_Expression_CellType_Primary_bySample[i,]~Sample_Region, col=2, las=2)
dev.off()
}

for(i in c(1:length(AVE_Expression_CellType_Tag_bySample[,1]))){
png(paste(row.names(AVE_Expression_CellType_Tag_bySample)[i], "_vsRegion.png"), width=3000, height=480)
par(mar=c(18,3,1,1))
boxplot(AVE_Expression_CellType_Tag_bySample[i,]~Sample_Region, col=3, las=2)
dev.off()
}


library(ggplot2)

getwd()

setwd(dir="/Users/mhh/Desktop/UROP Amygdala/ABA Cell type/RegionalCellTypes")


for(i in c(1:length(names(table(Sample_Region))))){

png(paste(names(table(Sample_Region))[i], "_vsPrimaryCellType.png", sep=""), width=1000, height=480)

TempMeans<-apply(AVE_Expression_CellType_Primary_bySample[,Sample_Region==names(table(Sample_Region))[i]], 1, function(y) mean(y, na.rm=T))
TempSE<-apply(AVE_Expression_CellType_Primary_bySample[,Sample_Region==names(table(Sample_Region))[i]], 1, function(y) sd(y, na.rm=T)/sqrt(length(y)))
PrimaryCellTypes<-names(TempMeans)

myData<-data.frame(PrimaryCellTypes, TempMeans,TempSE)

dodge<-position_dodge(width = 0.9)

limits<-aes(ymax=myData$TempMeans+myData$TempSE, ymin=myData$TempMeans-myData$TempSE)

p <- ggplot(data = myData, aes(x=PrimaryCellTypes, y=TempMeans, fill=PrimaryCellTypes))

print(p + geom_bar(stat = "identity", position = dodge) + geom_errorbar(limits, position = dodge, width = 0.25)+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank()))

dev.off()

}


#Very nice - what if we used violin plots instead to illustrate the point?


for(i in c(1:length(names(table(Sample_Region))))){

png(paste(names(table(Sample_Region))[i], "_vsPrimaryCellType.png", sep=""), width=1300, height=750)

Temp<-AVE_Expression_CellType_Primary_bySample[,Sample_Region==names(table(Sample_Region))[i]]
	RelativeExpression<-as.numeric(AVE_Expression_CellType_Primary_bySample[,Sample_Region==names(table(Sample_Region))[i]])
	CellType<-rep(row.names(AVE_Expression_CellType_Primary_bySample), ncol(Temp))
	myData<-data.frame(RelativeExpression,CellType)

p <- ggplot(myData, aes(x=CellType, y=RelativeExpression, fill=CellType)) + 
  geom_violin()
  
print(p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+stat_summary(fun.y=median, geom="point", size=2, color="red"))
  
dev.off()

}


#How about violin plots using cell type tag instead of primary cell type?

for(i in c(1:length(names(table(Sample_Region))))){

png(paste(names(table(Sample_Region))[i], "_vsPrimaryCellType.png", sep=""), width=2000, height=750)

Temp<-AVE_Expression_CellType_Tag_bySample[,Sample_Region==names(table(Sample_Region))[i]]
	RelativeExpression<-as.numeric(AVE_Expression_CellType_Tag_bySample[,Sample_Region==names(table(Sample_Region))[i]])
	CellType<-rep(row.names(AVE_Expression_CellType_Tag_bySample), ncol(Temp))
	myData<-data.frame(RelativeExpression,CellType)

p <- ggplot(myData, aes(x=CellType, y=RelativeExpression, fill=CellType)) + 
  geom_violin()
  
print(p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+stat_summary(fun.y=median, geom="point", size=2, color="red")+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
dev.off()

}




