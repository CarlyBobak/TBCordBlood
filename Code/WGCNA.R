#This file contains computational R code for weighted gene co-expression network analysis (WGCNA).
library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)
allowWGCNAThreads()

nGenes = ncol(expr)
(nSamples = nrow(expr))

expr<-apply(expr,2,as.numeric)

powers=c(1:30) # in practice this should include powers up to 30.
sft0=pickSoftThreshold(expr,powerVector=powers, networkType="signed")

pdf("Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

adjacencyPre = adjacency(expr,power=8 ,type="signed") # Our umbilical cord blood study required a beta power of 8.
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = hclust(as.dist(dissTOMPre), method="average")

#Plot gene tree
pdf("GeneTree.pdf",height=6,width=12)
par(mfrow=c(1,1))
plot(geneTreePre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity", labels=FALSE,hang=0.04);
dev.off()

#MODULE ASSIGNMENTS with cutreeHybrid algorithm
mColorh=NULL
for (ds in 0:4){
  tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
                      minClusterSize = (30), cutHeight = 0.99, 
                      deepSplit = ds, distM = dissTOMPre)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}

#Plot dendrogram and module assignments with deep-split options
pdf("DeepSplit_Choices.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Co-Expression Network",dendroLabels=FALSE);
dev.off()

#deepsplit value 2:

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTreePre, distM = dissTOMPre,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
tiff("clusteringDendrogram.tiff",height=8,width=10,units="in",res=600)
plotDendroAndColors(geneTreePre, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEList = moduleEigengenes(expr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result

tiff("clusteringofeigengenes.tiff",height=4,width=6,units="in",res=300)
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()


#pick a height cut of 0.25 -> correlation of 0.75 to merge. 

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(expr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

plotDendroAndColors(geneTreePre, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

MEs0 = moduleEigengenes(expr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

as.numeric.factor <- function(x) {as.factor(as.numeric(x))}

ctsPheno<-pheno
ctsPheno$pid_child<-NULL
ctsPheno$RIN<-NULL
ctsPheno$Chip<-NULL
ctsPheno$Batch<-NULL
ctsPheno$sex<-as.numeric.factor(ctsPheno$sex)
ctsPheno$conversion_10_at36months<-ifelse(ctsPheno$conversion_10_at36months=="Converter",1,0)
ctsPheno$site<-as.numeric.factor(ctsPheno$site)
ctsPheno$mat_hiv_status<-ifelse(ctsPheno$mat_hiv_status=="Negative",0,1)
ctsPheno$tbdisease_merged<-as.numeric(ctsPheno$tbdisease_merged)-1

#CALCULATE PC FOR VISUALIZATION
PCsPD    = moduleEigengenes((expr),  colors=moduleColors) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(moduleColors))
names = row.names((expr))

#Plot diagnostic module plots
pdf("Module_Visualizationq.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

for (which.module in names(table(moduleColors)))
{
  par(mfrow=c(2,1), mar=c(4, 4.1, 4.1, 2))
  plotMat(t(scale(expr[,moduleColors==which.module])),
          ,cex.axis=2,nrgcols=100,rlabels=F,tck=0, rcols=which.module,main=paste("Heatmap",which.module,"Module"))
  
  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, cex.main=1, ylab="Eigengene Expression",xlab="")
  axis(1,at=n, labels=row.names(expr), las=2, cex.axis=0.4, font=2)
};
dev.off();

#DETERMINE GENE-TRAIT RELATIONSHIPS ACROSS ALL GENES AND MODULES
Conversion = as.data.frame(as.numeric(convPheno$conversion_10_at36months))
names(Conversion)="Conversion"
GS.Conversion=as.numeric(cor(convExpr,Conversion,use="p"))
GS.ConversionColor=numbers2colors(GS.Conversion,signed=T)

TB = as.data.frame(as.numeric(convPheno$tbdisease_merged))
names(TB)="TB"
GS.TB=as.numeric(cor(convExpr,TB,use="p"))
GS.TBColor=numbers2colors(GS.TB,signed=T)

Conversion.TB = as.data.frame(as.numeric(convPheno$conversion_10_at36months)*as.numeric(convPheno$tbdisease_merged))
names(Conversion.TB)="Conversion.TB"
GS.Conversion.TB=as.numeric(cor(convExpr,Conversion.TB,use="p"))
GS.Conversion.TBColor=numbers2colors(GS.Conversion.TB,signed=T)

maternalTB = as.data.frame(as.numeric(convPheno$tb_everdx))
names(maternalTB)="Maternal TB"
GS.Maternal.TB=as.numeric(cor(convExpr,maternalTB,use="p"))
GS.Maternal.TBColor=numbers2colors(GS.Maternal.TB,signed=T)

datColors0=data.frame(moduleColors, GS.ConversionColor,GS.TBColor,GS.Conversion.TBColor, GS.Maternal.TBColor)

#Plot gene tree and color bars
pdf("Drakenstein_ColoredModules.pdf",height=8,width=14)
plotDendroAndColors(geneTreePre, colors=datColors0, main="Gene Dendrogram and Module Colors (Beta=8)", groupLabels=c("Module colors", "Conversion", "Tuberculosis","ConversionxTB", "Maternal TB"), dendroLabels=FALSE, hang=0.03, addGuide=FALSE, guideHang=0.05) 
dev.off()

moduleTraitCor = cor(MEs, ctsPheno, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(convExpr)) 

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");


dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(10, 9, 3, 3))
# Display the correlation values within a heatmap plot
tiff("moduleHeatmap.tiff",height=6,width=12,units="in",res=900)
par(mar = c(10, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(ctsPheno),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


# Define variable of interest
conversion = as.data.frame(ctsPheno$conversion_10_at36months);
names(conversion) = "conversion"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(expr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(expr, conversion, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(conversion), sep="");
names(GSPvalue) = paste("p.GS.", names(conversion), sep="");

#variable with names of modules of interest
ConversionModules<-rownames(moduleTraitPvalue)[which(moduleTraitPvalue[,"conversion_10_at36months"]<=0.1)]

midnightBlueGenes<-colnames(expr)[moduleColors=="midnightblue"]
length(midnightBlueGenes) #187

redGenes<-colnames(expr)[moduleColors=="red"]
length(redGenes) #744

greyGenes<-colnames(expr)[moduleColors=="grey60"]
length(greyGenes) #111

lightGreenGenes<-colnames(expr)[moduleColors=="lightgreen"]
length(lightGreenGenes) #99

salmonGenes<-colnames(expr)[moduleColors=="salmon"]
length(salmonGenes) #1731

turquoiseGenes<-colnames(expr)[moduleColors=="turquoise"]
length(turquoiseGenes) #1617

lightcyanGenes<-colnames(expr)[moduleColors=="lightcyan"]
length(lightcyanGenes) #120

magentaGenes<-colnames(expr)[moduleColors=="magenta"]
length(magentaGenes) #120

brownGenes<-colnames(expr)[moduleColors=="brown"]
length(brownGenes) #984

greenYellowGenes<-colnames(expr)[moduleColors=="greenyellow"]
length(greenYellowGenes) #1275

royalBlueGenes<-colnames(expr)[moduleColors=="royalblue"]
length(royalBlueGenes) #37

blackGenes<-colnames(expr)[moduleColors=="black"]
length(blackGenes) #903

greenGenes<-colnames(expr)[moduleColors=="green"]
length(greenGenes) #1248

pinkGenes<-colnames(expr)[moduleColors=="pink"]
length(pinkGenes) #1073




write.csv(colnames(expr),"backgroundGenes.csv",row.names = F)
write.csv(royalBlueGenes,"royalBlueGenes.csv",row.names=F)
write.csv(greenYellowGenes,"greenYellowGenes.csv",row.names=F)
write.csv(lightGreenGenes,"lightGreenGenes.csv",row.names=F)
write.csv(brownGenes,"brownGenes.csv",row.names=F)
write.csv(pinkGenes,"pinkGenes.csv",row.names=F)
write.csv(turquoiseGenes,"turquoiseGenes.csv",row.names=F)


write.table(colnames(expr),"backgroundGenes.txt",row.names = F,sep="\t")
write.table(royalBlueGenes,"royalBlueGenes.txt",row.names=F,sep="\t")
write.table(greenYellowGenes,"greenYellowGenes.txt",row.names=F,sep="\t")
write.table(lightGreenGenes,"lightGreenGenes.txt",row.names=F,sep="\t")
write.table(brownGenes,"brownGenes.txt",row.names=F,sep="\t")
write.table(pinkGenes,"pinkGenes.txt",row.names=F,sep="\t")
write.table(turquoiseGenes,"turquoiseGenes.txt",row.names=F,sep="\t")

library(ggpubr)

testPheno<-c("conversion_10_at36months","tbdisease_merged","tb_everdx","pneumonia_any","smoke_current","prem","lbw")
phenoMod<-pheno[,testPheno]
MEdf<-data.frame(MEs,phenoMod)
MEdf<-MEdf[-which(is.na(pheno$conversion_10_at36months)),]

testMEs<-c(3,4,6,11,14)

p1<-ggboxplot(MEdf,x="conversion_10_at36months",y=colnames(MEs)[11],color="conversion_10_at36months",palette = c('goldenrod1','blue3'),
          add = "jitter", xlab=F)
p1 <- p1 + stat_compare_means()

p1


p2<-ggboxplot(MEdf,x="tbdisease_merged",y=colnames(MEs)[14],color="tbdisease_merged",palette = c('goldenrod1','blue3'),
              add = "jitter", xlab=F)
p2 <- p2 + stat_compare_means()

p2

subMEDf<-MEdf[MEdf$conversion_10_at36months=="Converter",]

p3<-ggboxplot(subMEDf,x="tbdisease_merged",y=colnames(MEs)[14],color="tbdisease_merged",palette = c('goldenrod1','blue3'),
              add = "jitter", xlab=F)
p3 <- p3 + stat_compare_means()

p3

converterTB<-as.character(MEdf$tbdisease_merged)
converterTB[MEdf$tbdisease_merged=="No TB Diagnosis" & MEdf$conversion_10_at36months=="Converter"]<- "Converter: No TB Disease"
converterTB[MEdf$conversion_10_at36months!="Converter"]<-"Nonconverter"

MEdf$converterTB<-as.factor(converterTB)

comps<-list(c("TB Disease","Converter: No TB Disease"),c("TB Disease","Nonconverter"),c("Converter: No TB Disease","Nonconverter"))

testModules<-c(4,6,9,10,11,14) #0.05 cutoff any stat

moduleKey<-paste("M",1:14,sep="")
names(moduleKey)<-colnames(MEs)

pdf("interestingModules.pdf")

for(i in testModules){

  p4<-ggboxplot(MEdf,x="converterTB",y=colnames(MEs)[i],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA, ylab=moduleKey[i])
  p4 <- p4 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none") 

  print(p4)
}
dev.off()

for(i in testModules){
  
  p4<-ggboxplot(MEdf,x="converterTB",y=colnames(MEs)[i],color="converterTB",
                palette = c('red','goldenrod1','blue3'),
                order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
                add = "jitter", xlab=F, outlier.shape=NA, ylab=moduleKey[i],title = paste("M",i,sep=""))
  p4 <- p4 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none")

  png(paste("M",i,"byTB.png",sep=""),height=4,width=6,units="in",res=200)
  print(p4)
  dev.off()
}



testModules<-c(4,6,9,10,11,14) 

MEdf$tb_everdx<-ifelse(MEdf$tb_everdx==1,"Prior TB","No Prior TB Diagnosis")


p4<-ggboxplot(MEdf,x="converterTB",y=colnames(MEs)[10],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA)
p4 <- p4 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none",text = element_text(size = 15)) 

MEdf$tb_everdx<-as.factor(MEdf$tb_everdx)

p5<-ggboxplot(MEdf,x="tb_everdx",y=colnames(MEs)[10],color="tb_everdx",
              palette = c("turquoise","salmon"),
              order=c("Prior TB","No Prior TB Diagnosis"),
              add = "jitter", xlab="Maternal TB Diagnosis", outlier.shape=NA)
p5 <- p5 +stat_compare_means()+ theme(legend.position = "none") 


tiff("MEgreenyellowdist.tiff",height=6,width=6, units="in",res=600)
ggarrange(p4,p5,nrow=2)
dev.off()


p4<-ggboxplot(MEdf,x="converterTB",y=colnames(MEs)[11],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Disease","Nonconverter"),
              add = "jitter", xlab=F, ylab="M11 Eigengene", outlier.shape=NA)
p4 <- p4 +stat_compare_means(method="wilcox",comparisons = comps)+ theme(legend.position = "none") 

png("cleanM11byTB.png",height=3.5,width=7,units="in",res=500)
p4
dev.off()

p5<-ggboxplot(MEdf,x="tb_everdx",y=colnames(MEs)[11],color="tb_everdx",
              palette = c("salmon","turquoise"),
              order=c("Prior TB","No Prior TB Diagnosis"),
              add = "jitter", xlab="Maternal TB Diagnosis", ylab="M11 Eigengene")
p5 <- p5 +stat_compare_means()+ theme(legend.position = "none",text = element_text(size = 15)) 

png("M11EverTb.png",height=4,width=5,units="in",res=500)
p5
dev.off()

MEdf$lbw<-ifelse(MEdf$lbw==0,"Normal","Low")

p6<-ggboxplot(MEdf,x="lbw",y=colnames(MEs)[11],color="lbw",
              palette = c("dark orange","dark green"),
              order=c("Low","Normal"),
              add = "jitter", xlab="Birthweight", ylab="M11 Eigengene", outlier.shape=NA)
p6 <- p6 +stat_compare_means()+ theme(legend.position = "none") 


tiff("MERoyalBluedist.tiff",height=8,width=8, units="in",res=600)
ggarrange(p4,ggarrange(p5,p6,ncol=2),nrow=2,heights=c(1.5,1))
dev.off()

tiff("RoyalBlueBoxplot.tiff",height=4,width=8, units="in",res=600)
p4
dev.off()

tiff("additionalRoyalBlueBoxplots.tiff",height=4,width=8, units="in",res=600)
ggarrange(p5,p6,ncol=2)
dev.off()

intModulesNames<-colnames(MEs[testModules])

p4<-ggboxplot(MEdf,x="tb_everdx",y=colnames(MEs)[14],color="tb_everdx",
              palette = c("salmon","turquiose"),
              order=c("Prior TB","No Prior TB Diagnosis"),
              add = "jitter", xlab="Maternal TB Diagnosis")
p4 <- p4 +stat_compare_means()+ theme(legend.position = "none") 


#checked others -> no association

#low birthweight associated with royal blue

p4<-ggboxplot(MEdf,x="lbw",y=colnames(MEs)[14],color="lbw",
                       palette = c("turquoise","salmon"),
                       #order=c("Ever TB","Never TB"),
                       add = "jitter", xlab="low birth weight")
p4 <- p4 +stat_compare_means()+ theme(legend.position = "none") 

print(p4)

converterTB<-as.factor(converterTB)

Name<-colnames(convExpr)
Description<-Name

writeExpr<-data.frame(Name,Description,t(convExpr))
write.table(writeExpr,file="convergenceExprData.txt",sep="\t", col.names = T, row.names = F, append = F, quote = F)

ciberDat<-writeExpr
ciberDat$Description<-NULL
colnames(ciberDat)[1]<-"Gene"

write.table(ciberDat,file="ciberExprData.txt",sep="\t", col.names = T, row.names = F, append = F, quote = F)

cat(c(as.character(c(length(converterTB),2,1)),"\n"),file="converterTBLabels.cls",sep="\t")
cat(c("#",levels(convPheno$conversion_10_at36months),"\n"),file="converterTBLabels.cls",sep="\t",append=T)
cat(as.character(convPheno$conversion_10_at36months),file="converterTBLabels.cls",sep="\t",append=T)

#DETERMINE MODULE SIG. USING -LOG PVALUE FOR EACH INDIVIDUAL GENE
datSummary <- read.csv("ResultsMaster.csv", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE) #THIS FILE CONTAINS A VECTOR OF DGE P-VALUES FROM DIFFERENTIAL GENE EXPRESSION ANALYSES, GENES NEED TO BE SORTED IN THE EXACT SAME ORDER AS MATRIX INPUT
GS_Conversion=-log10(datSummary$Converter) #Conversion DGE P-VALUES
GS_TB=-log10(datSummary$TB.Disease) #TB P-VALUES
GS_convTB=-log10(datSummary$ConverterxTB)#TB within Converters DGE P-VALUES


ColModuleKey<-moduleKey
names(ColModuleKey)<-gsub("ME", "", names(moduleKey))

orderedModuleKey<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n")
names(orderedModuleKey)<-names(ColModuleKey)
names(orderedModuleKey)[14]<-"purple2"
moduleColors[moduleColors=="pink"]<-"purple2"

barPlotDF<-as.data.frame(cbind(GS_Conversion,ColModuleKey[moduleColors],moduleColors,orderedModuleKey[moduleColors]))

#PLOT MODULE SIGNIFICANCE
pdf("ModuleSignificance.pdf")
par(mfrow=c(3,1))
verboseBarplot(GS_Conversion,orderedModuleKey[moduleColors],color=names(orderedModuleKey), border="black",main="Conversion" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9, ylim=c(0,1.2),names.arg=moduleKey)
abline(h=0.6,col="black",lty = 2)   
verboseBarplot(GS_TB,orderedModuleKey[moduleColors],color=names(orderedModuleKey),border="black",main="TB Disease" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = T,ylim=c(0,1.2),names.arg=moduleKey)
abline(h=0.6,col="black",lty = 2)   
verboseBarplot(GS_convTB,orderedModuleKey[moduleColors],color=names(orderedModuleKey),border="black",main="TB Disease within Converters" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = T,ylim=c(0,1.2),names.arg=moduleKey)
abline(h=0.6,col="black",lty = 2)    
dev.off()

library(latex2exp)

png("justCoversionBarPlot.png",height=4,width=8,units="in",res=600)
par(mfrow=c(1,1))
verboseBarplot(GS_Conversion,orderedModuleKey[moduleColors],color=names(orderedModuleKey), border="black",main="Conversion" ,xlab="",ylab=TeX("$-log_{10}(p-value)$") ,las=2, cex.axis=1,cex.lab=1.4, cex=1,cex.main=1.5, ylim=c(0,1.2),names.arg=moduleKey)
abline(h=0.6,col="black",lty = 2)
dev.off()


##CELL-TYPE AND GENE-SET ENRICHMENT
Gene  = colnames(expr)
enrichments = userListEnrichment(Gene, moduleColors,
                                 fnIn = NULL,
                                 catNmIn = NULL,
                                 #fnIn = c("GeneList","ModuleColors"),
                                 #catNmIn =  c("Genes","Modules"),
                                 nameOut = "ModuleEnrichment2.csv", 
                                 useBrainLists = TRUE,
                                 useBloodAtlases = TRUE,
                                 omitCategories = "grey", 
                                 outputCorrectedPvalues = TRUE,
                                 useStemCellLists = TRUE, 
                                 outputGenes = TRUE, 
                                 minGenesInCategory = 2, 
                                 useBrainRegionMarkers = TRUE,
                                 useImmunePathwayLists = TRUE,
                                 usePalazzoloWang = TRUE)

#PLOT CONNECTIVITY VS GS.TRAIT MEASURES
CM_Pre=intramodularConnectivity(adjacencyPre,colors=moduleColors)
names(CM_Pre)
intModulesNames[6]<-"MEpurple2"
pdf("ConversionConnectivity.pdf")
for(i in intModulesNames){
  whichmodule=i;
  restrict1=paste("ME",moduleColors,sep="")==whichmodule
  verboseScatterplot (CM_Pre$kWithin[restrict1],GS.Conversion[restrict1],col=moduleColors[restrict1],
                      xlab=paste(i,"Connectivity (k) ",sep=""),ylab="G.S. TST Conversion(p)")
}
dev.off()

# CALCULATE MODULE MEMBERSHIP VALUES (aka. module eigengene based connectivity kME)
datKME=signedKME(expr, MEs)
colorOfColumn=substring(names(datKME),4)
colorOfColumn[14]<-"purple2"

pdf("Regress_Modules.pdf",w=11,h=8.5)
par(mfrow = c(3,3))
moduleNo=1
names(moduleKey2)[14]<-"purple2"
for (module in names(moduleKey2[testModules]))
{
  column = match(module,colorOfColumn) 
  restModule=moduleColors==module
  
  verboseScatterplot(datKME[restModule,column],GS.Conversion[restModule],
                     xlab=paste("M",moduleNo,"Gene Connectivity (kME)"),ylab=TeX("Correlation with TST Conversion"),
                     main=paste("kME vs Conversion Gene Significance"),col=module,abline=T,las=1, cex.axis=0.9, pch=16, cex.main = 0.75, cex.lab = 0.75)
  
  verboseScatterplot(datKME[restModule,column],GS.TB[restModule],
                     xlab=paste("M",moduleNo,"Gene Connectivity (kME)"),ylab=TeX("Correlation with TB Diagnosis"),
                     main=paste("kME vs TB Gene Significance"),col=module,abline=T,las=1, cex.axis=0.9, pch=16, cex.main = 0.75, cex.lab = 0.75)

  verboseScatterplot(datKME[restModule,column],GS.Maternal.TB[restModule],
                     xlab=paste("M",moduleNo,"Gene Connectivity (kME)"),ylab=TeX("Correlation with Maternal TB History"),
                     main=paste("kME vs Maternal TB History Gene Significance"),col=module, abline=T,las=1, cex.axis=0.9, pch=16, cex.main = 0.75, cex.lab = 0.75)
  moduleNo=moduleNo+1
}
dev.off()

column = match("royalblue",colorOfColumn) 
module="royalblue"
restModule=moduleColors==module

png("royalBlueRegression.png",height=6,width=6, units="in",res=900)
par(mar=c(5.1, 5.1, 4.1, 2.1))
verboseScatterplot(datKME[restModule,column],GS.Conversion[restModule],
                   xlab=paste("M11 Gene Connectivity (kME) "),ylab=TeX("Correlation with TST Conversion"),
                   col=module,abline=T,las=1, pch=16,cex.axis = 0.9,cex.lab = 1.5)
dev.off()

#MODULE MEMBERSHIP (kME) KME is defined as correlation between expression and modules
geneModuleMembership1 = signedKME((expr), MEs)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(expr)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(expr))
kMEtable1  = cbind(Gene,Gene,moduleColors)
for (i in 1:length(colorsPD))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"Drakenstein_kMEtable.csv",row.names=FALSE)

topGenesKME = NULL
for (c in 1:length(colorsPD)){
  kMErank1    = rank(-geneModuleMembership1[,c])
  maxKMErank  = rank(apply(cbind(kMErank1+.00001),1,max))
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=20])
}; 

colnames(topGenesKME) = colorsPD
topGenesKME #OUTPUT TOP HUB GENES

moduleKey2<-moduleKey
names(moduleKey2)<-modNames

newColNames<-moduleKey2[colnames(topGenesKME)]

colnames(topGenesKME)<-as.character(newColNames)
topGenesKME<-topGenesKME[,as.character(ColModuleKey)]

write.csv(topGenesKME,"TBConversionHubGenes.csv")

#WRITE OUTPUT
Pre<-as.data.frame(as.numeric(pheno$conversion_10_at36months))
names(Pre)<-"Conversion"
modNames = substring(names(MEs), 3)
nGenes = ncol(expr) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(expr) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(expr,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(expr, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

Colors=moduleColors
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="Network_WGCNA_OUTPUT_q.csv")