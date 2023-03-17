#library(DESeq)
library(RColorBrewer)
library(gplots)
library(ggplot2)
#library(RUVSeq)
#library(EDASeq)
library(limma)
library(tidyverse)
library(ggpubr)

expr <- read.delim("ExprsMat.txt", header=FALSE)
dim(expr)
expr<-t(expr)
colnames(expr)<-expr[1,]
colnames(expr)[1]<-"pid_child"

pheno <- read.csv("cordPheno.csv")
dim(pheno)

cordDataSet<-merge(pheno,expr)
dim(cordDataSet)

###Check breast feeding ever.

meta<-read.delim("/Volumes/GoogleDrive/My Drive/PhD QBS/Collaborations/Leo Cord Gene Expression/Umbilical-cord-blood-transcriptome-analysis/MetaData.txt")
colnames(meta)[4]<-"pid_child"

#datacleaning
cordDataSet<-merge(meta,cordDataSet)
cordDataSet<-cordDataSet[,-c(2:4,6:8,12:19)]
colnames(cordDataSet)[2]<-"sex"

empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

library(magrittr)
library(tidyverse)


## transform all columns
cordDataSet %<>% mutate_each(funs(empty_as_na))
cordDataSet$time_conversion_new_years<-as.numeric(cordDataSet$time_conversion_new_years)
cordDataSet$timeconversion_year_cutoff<-as.numeric(cordDataSet$timeconversion_year_cutoff)
cordDataSet$conversion_10_at36months<-as.factor(cordDataSet$conversion_10_at36months)
cordDataSet$tbdisease_merged<-as.factor(cordDataSet$tbdisease_merged)
cordDataSet$time_tb_new_years<-as.numeric(cordDataSet$time_tb_new_years)
cordDataSet$NEW_tb_momonrx<-as.factor(cordDataSet$NEW_tb_momonrx)
cordDataSet$tb_everdx<-as.factor(cordDataSet$tb_everdx)
cordDataSet$pneumonia_any<-as.factor(cordDataSet$pneumonia_any)
cordDataSet$site<-as.factor(cordDataSet$site)
cordDataSet$mat_hiv_status<-as.factor(cordDataSet$mat_hiv_status)
cordDataSet$smoke_current<-as.factor(cordDataSet$smoke_current)
cordDataSet$anyPriorTB<-as.factor(ifelse(cordDataSet$tb_priorepisodes>0,1,0))
cordDataSet$tb_household<-as.factor(cordDataSet$tb_household)
cordDataSet$prem<-as.factor(cordDataSet$prem)
cordDataSet$lbw<-as.factor(cordDataSet$lbw)
cordDataSet$sex<-as.factor(cordDataSet$sex)
cordDataSet$Batch<-as.factor(cordDataSet$Batch)

#saveRDS(cordDataSet,"mergedTBDat.RData")

#######Load cleaned data
cordDataSet<-readRDS("mergedTBDat.RData")
expr<-cordDataSet[,-c(1:24,ncol(cordDataSet))]
pheno<-cordDataSet[,c(1:24,ncol(cordDataSet))]
summary(pheno)


write.csv(pheno,"pheno.csv",row.names = F)

#Generate Table One
library(table1)

#table1 by site

pheno$Outcome<-as.character(pheno$conversion_10_at36months)
pheno$Outcome[pheno$tbdisease_merged=="TB Disease"]<-"TB Disease"
pheno$Outcome<-as.factor(pheno$Outcome)

table1(~NEW_tb_momonrx+tb_everdx+tb_household+mat_hiv_status+smoke_current+site+sex+birthweight_gms|Outcome, data=pheno[is.na(pheno$conversion_10_at36months)==F,])



table1(~mat_hiv_status+birthweight_gms+Breast_feeding_months+smoke_current+tb_household+NEW_tb_momonrx+tb_everdx+sex+site+tbdisease_merged|conversion_10_at36months,data=pheno[is.na(pheno$conversion_10_at36months)==F,])
chisq.test(pheno$site,pheno$conversion_10_at36months)
chisq.test(pheno$mat_hiv_status,pheno$conversion_10_at36months) #p value 0.1238
t.test(pheno$birthweight_gms~pheno$conversion_10_at36months)
chisq.test(pheno$smoke_current,pheno$conversion_10_at36months)
t.test(pheno$Breast_feeding_months~pheno$conversion_10_at36months) #pvalue is 0.1377
chisq.test(pheno$tb_household,pheno$conversion_10_at36months)
chisq.test(pheno$tb_everdx,pheno$conversion_10_at36months)
chisq.test(pheno$sex, pheno$conversion_10_at36months)
chisq.test(pheno$site,pheno$conversion_10_at36months)
chisq.test(pheno$tbdisease_merged,pheno$conversion_10_at36months)

#make nicer factor
pheno$tbdisease_merged<-as.factor(ifelse(pheno$tbdisease_merged=="1","TB Disease","No TB Diagnosis"))

table1(~site+mat_hiv_status+gestation_delivery+birthweight_gms+smoke_current+prem+lbw+Breast_feeding_months+tb_household+NEW_tb_momonrx+tb_everdx+pneumonia_any+pneum_count+sex+conversion_10_at36months|tbdisease_merged,data=pheno[is.na(pheno$tbdisease_merged)==F,])
chisq.test(pheno$site,pheno$tbdisease_merged)
chisq.test(pheno$mat_hiv_status,pheno$tbdisease_merged)
chisq.test(pheno$sex,pheno$tbdisease_merged)
chisq.test(pheno$NEW_tb_momonrx,pheno$tbdisease_merged)
chisq.test(pheno$tb_household,pheno$tbdisease_merged)
chisq.test(pheno$tb_everdx,pheno$tbdisease_merged)
t.test(birthweight_gms~tbdisease_merged,data=pheno) #pvalue 0.1366
chisq.test(pheno$smoke_current,pheno$tbdisease_merged) #pvalue is 0.06966
t.test(pneum_count~tbdisease_merged,data=pheno) #pvalue is 0.06282
chisq.test(pheno$pneumonia_any,pheno$tbdisease_merged) #pvalue is 0.09833
chisq.test(pheno$conversion_10_at36months,pheno$tbdisease_merged) #pvalue is 0.09833
prop.test(table(pheno$lbw,pheno$tbdisease_merged)) #pvalue is 0.1274
t.test(pheno$Breast_feeding_months~pheno$tbdisease_merged) #pvalue is 0.09865



#for within converters
withinConvertersPheno<-pheno[pheno$conversion_10_at36months=="Converter",]
table1(~site+mat_hiv_status+gestation_delivery+birthweight_gms+smoke_current+prem+lbw+Breast_feeding_months+tb_household+NEW_tb_momonrx+tb_everdx+pneumonia_any+pneum_count+sex|tbdisease_merged,data=withinConvertersPheno[is.na(withinConvertersPheno$tbdisease_merged)==F,])
chisq.test(withinConvertersPheno$site,withinConvertersPheno$tbdisease_merged)
chisq.test(withinConvertersPheno$mat_hiv_status,withinConvertersPheno$tbdisease_merged)
t.test(birthweight_gms~tbdisease_merged,data=withinConvertersPheno) 
chisq.test(withinConvertersPheno$smoke_current,withinConvertersPheno$tbdisease_merged) 
t.test(pneum_count~tbdisease_merged,data=withinConvertersPheno)
chisq.test(withinConvertersPheno$pneumonia_any,withinConvertersPheno$tbdisease_merged) 
prop.test(table(withinConvertersPheno$lbw,withinConvertersPheno$tbdisease_merged)) 
t.test(withinConvertersPheno$Breast_feeding_months~withinConvertersPheno$tbdisease_merged) #pvalue is 0.05612

#for pneumonia
#for within converters
pheno$pneumonia_any<-as.factor(ifelse(pheno$pneumonia_any==0,"No Pneumonia","Pneumonia"))
table1(~site+mat_hiv_status+gestation_delivery+birthweight_gms+smoke_current+prem+lbw+Breast_feeding_months+tb_household+NEW_tb_momonrx+tb_everdx+pneumonia_any+sex|pneumonia_any,data=pheno)
chisq.test(pheno$site,pheno$pneumonia_any)
chisq.test(pheno$mat_hiv_status,pheno$pneumonia_any)
t.test(birthweight_gms~pneumonia_any,data=pheno) 
chisq.test(pheno$smoke_current,pheno$pneumonia_any) 
prop.test(table(pheno$tb_household,pheno$pneumonia_any)) #0.0592
prop.test(table(pheno$lbw,pheno$pneumonia_any)) 
prop.test(table(pheno$sex,pheno$pneumonia_any)) #0.01934
t.test(pheno$Breast_feeding_months~pheno$pneumonia_any) 

#start with converters vs not converters
ix<-which(is.na(pheno$conversion_10_at36months))
convPheno<-pheno[-ix,]
convExpr<-expr[-ix,]

convExpr<-apply(convExpr,2,as.numeric)

# Check covariates
Group <- factor(convPheno$conversion_10_at36months, levels=c("Nonconverter","Converter"))
RIN <- as.numeric(convPheno$RIN)
Batch <- factor(convPheno$Batch)
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Group+RIN+Batch)
colnames(design)

fit <- lmFit(t(convExpr), design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
Diff<- topTable(fit2, coef="GroupConverter", adjust="BH", n=25000)
#write.csv(Diff, file="CoverterDEResults.csv", quote=FALSE, row.names=TRUE)

#create volcano plots of DGE signatures
plotDiff<-Diff
plotDiff$sig<-ifelse(Diff$P.Value<0.005,"Significant","Non-significant")
plotDiff$lab<-ifelse(Diff$P.Value<0.001,"Yes","No")
plotDiff$dir<-ifelse(Diff$logFC<0,"Down","Up")
plotDiff$gene<-rownames(Diff)
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$colr<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
levels(plotDiff$colr)
plotDiff$colr<-paste(plotDiff$sig,plotDiff$dir,sep=" ")
plotDiff$colr[plotDiff$sig=="Non-significant"]<-"Non-significant"
plotDiff$colr<-as.factor(plotDiff$colr)
levels(plotDiff$colr)

library(latex2exp)
library(ggrepel)

png("converterVolcanoLabels.png",width=6,height=6,units="in",res=600)
  volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
    geom_point(aes(col=colr),size=1) + #add points colored by significance
    scale_color_manual(values=c("red3","cornflower blue", "black"),labels=c("Up Regulated","Down Regulated","Non-significant"),breaks = rev(levels(plotDiff$colr)))+
    geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
    annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
    geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
    annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
    xlab(TeX("$log_2 Fold Change$"))+ylab(TeX("$-log_{10}(p-value)$"))+
    theme_classic()+theme(legend.title=element_blank(), text = element_text(size = 15))
  volc+geom_text_repel(data=plotDiff[plotDiff$lab=="Yes",], aes(label=gene),size=3)
dev.off()

tiff("converterVolcano.tiff",width=6,height=6,units="in",res=600)
  volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
    geom_point(aes(col=colr),size=1) + #add points colored by significance
    scale_color_manual(values=c("black","cornflower blue", "red3"),labels=c("Up Regulated","Down Regulated","Insignificant"),breaks = rev(levels(plotDiff$colr)))+
    geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
    annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
    geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
    annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
    xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
    theme_classic()+theme(legend.title=element_blank())
  volc
dev.off()

#check the PC components for correlations with variables approaching significants from Table 1
testFactors<-c("site","mat_hiv_status","smoke_current","prem","lbw","pneumonia_any","sex","conversion_10_at36months","tb_household","NEW_tb_momonrx","tb_everdx","tbdisease_merged")
testContinuous<-c("gestation_delivery","birthweight_gms","Breast_feeding_months","pneum_count")


DGEs<-row.names(Diff[Diff$P.Value<0.005,])
DGEMatrix<-convExpr[,DGEs]

pcs<-prcomp(DGEMatrix)

library(ggbiplot)

colors = c('goldenrod1','blue3')
names(colors) = c("Converter","Nonconverter")

g <- ggbiplot(pcs, choices=c(1,2),
              groups = as.factor(convPheno$conversion_10_at36months),
              circle = T, var.axes = F,
              ellipse=T)
g <- g + geom_point(aes(col=as.factor(convPheno$conversion_10_at36months)), size=1)
g<- g + scale_colour_manual(values = colors, name="TB Category")
g <- g + theme(legend.direction = 'vertical',
               legend.position = 'right')
g <- g + theme_classic()+theme(text = element_text(size = 15))

g

png("ConverterPCA_ellipse.png",height=6,width=6,units="in",res=600)
g
dev.off()


library(ggpubr)

PCdata<-data.frame(pcs$x[,1:3])
PCdata$Conversion<-convPheno$conversion_10_at36months
PCdata$TBDisease<-convPheno$tbdisease_merged
PCdata$MaternalPriorTB<-as.factor(ifelse(convPheno$tb_everdx=="1","Prior TB","No Prior TB Diagnosis"))
PCdata$BreastFeedingDuration<-convPheno$Breast_feeding_months

p1<-ggboxplot(PCdata,x="Conversion",y="PC1", color="Conversion",
              palette=c("goldenrod1","blue3"),
              order=c("Converter","Nonconverter"),
              add="jitter",xlab=F,outlier.shape=NA)


p1<-p1+stat_compare_means()+theme(legend.position = "none",text = element_text(size = 15))

p2<-ggboxplot(PCdata,x="TBDisease",y="PC1", color="TBDisease",
              palette=c("red","dodger blue"),
              order=c("TB Disease","No TB Diagnosis"),
              add="jitter",xlab=F,outlier.shape=NA)


p2<-p2+stat_compare_means()+theme(legend.position = "none")

p3<-ggboxplot(PCdata,x="MaternalPriorTB",y="PC2", color="MaternalPriorTB",
              palette=c("salmon","turquoise"),
              order=c("Prior TB","No Prior TB Diagnosis"),
              add="jitter",xlab=F,outlier.shape=NA)


p3<-p3+stat_compare_means()+theme(legend.position = "none",text = element_text(size = 15))

png("conversionPCBox.png",height=4.5, width=4.5, units="in", res=600)
p1
dev.off()

tiff("TBPCBox.tiff",height=4, width=4, units="in", res=600)
p2
dev.off()

png("everTBPCBox.png",height=4.5, width=4.5, units="in", res=600)
p3
dev.off()


tiff("ConversionPCAwithBoxPlots.tiff",height=6, width=10, units="in", res=600)
ggarrange(g,p1,p2,p3,ncol=2,common.legend=T,legend="top")
dev.off()

#df<-data.frame(convPheno$conversion_10_at36months,convPheno$site,convPheno$smoke_current,convPheno$pneumonia_any,convPheno$tbdisease_merged,convPheno$tb_everdx)
#colnames(df)<-c("Conversion","Site","Smoking Status","Pneumonia","Tuberculosis Disease","Maternal Tuberculosis")

df<-data.frame(convPheno$conversion_10_at36months)#,convPheno$tbdisease_merged)
colnames(df)<-c("Conversion")#,"Tuberculosis Disease")


mycols <- colorRamp2(breaks = c(-1,0,1), 
                     colors = c("navy","white","red"))

dists<-c("euclidean","manhattan","maximum","canberra","minkowski","kendall","pearson","spearman","binary")
clusts<-c("complete","average","median","single","ward.D","ward.D2","mcquitty","centroid")

library(ComplexHeatmap)
library(circlize)

setwd("/Volumes/GoogleDrive/My Drive/PhD QBS/Collaborations/Leo Cord Gene Expression/HeatmapConverge")

row.names(convExpr)<-convPheno$pid_child

dists<-c("euclidean","manhattan","maximum","canberra","minkowski","kendall","pearson","spearman")
clusts<-c("complete","average","median","single","ward.D","ward.D2","mcquitty","centroid")

top_annotation = HeatmapAnnotation(name="TB Cord Blood",df=df,
                                   col=list("Conversion"=colors),
                                   show_legend = T,
                                   annotation_legend_param = list(title_gp = gpar(fontsize = 12), 
                                                                  labels_gp = gpar(fontsize = 12)))

convExprCenter<-apply(DGEMatrix,2,function(x){x<-x-median(x)})

for(i in dists){
  for(j in clusts){
    for(k in clusts){
      for(l in dists){

    
    HM<-Heatmap(t(convExprCenter), 
                name = " ", #title of legend
                column_title = NULL, row_title = NULL,
                show_column_names=T,
                clustering_distance_columns=i,
                clustering_method_columns = j,
                clustering_method_rows = k,
                clustering_distance_rows = l,
                row_names_gp = gpar(fontsize = 8),
                top_annotation = top_annotation,
                show_heatmap_legend = T,
                show_row_dend = F,
                heatmap_legend_param = list(color_bar = "continuous",labels_gp=gpar(fontsize=12)),
                #cluster_rows = color_branches(row_dend, k = 3),
                col=mycols
    )
    
    tiff(paste(i,j,k,".tiff",sep=""),height=12,width=12,units="in",res=600)
    draw(HM,  heatmap_legend_side = "left", annotation_legend_side="bottom")
    dev.off()
  }}}
}



HM<-Heatmap(t(convExprCenter), 
            name = " ", #title of legend
            column_title = NULL, row_title = NULL,
            show_column_names=T,
            clustering_distance_columns="canberra",
            clustering_method_columns = "ward.D2",
            clustering_method_rows = "ward.D2",
            clustering_distance_rows = "canberra",
            row_names_gp = gpar(fontsize = 8),
            top_annotation = top_annotation,
            show_heatmap_legend = T,
            show_row_dend = F,
            heatmap_legend_param = list(color_bar = "continuous",labels_gp=gpar(fontsize=12)),
            #cluster_rows = color_branches(row_dend, k = 3),
            col=mycols
)

tiff("conversionHeatmap.tiff",height=16,width=12,units="in",res=900)
draw(HM,  heatmap_legend_side = "bottom", annotation_legend_side="bottom")
dev.off()



tiff("PCDistConv.tiff",height=5,width=7,units="in", res=300,)
plot(pcs)
dev.off()
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]

library(ggpubr)
library(ggbiplot)

tiff("converterPCA.tiff",height=5,width=5,units="in",res=300)
g <- ggbiplot(pcs, choices=c(1,2),
              groups = convPheno$conversion_10_at36months, 
              circle = T, var.axes = F)
g <- g + geom_point(aes(col=convPheno$conversion_10_at36months), size=3)
g<- g + scale_color_manual(values = c("Converter"="red","Nonconverter"="navy"))
g <- g + theme(legend.direction = 'vertical', 
               legend.position = 'right')+theme_classic()
print(g)
dev.off()

PCdata<-data.frame(pcs$x[,1:3])
PCdata$Conversion<-convPheno$conversion_10_at36months


p1 <- ggboxplot(PCdata, x = "Conversion", y = "PC1",
                color = "Conversion", palette = c('red','navy'),
                add = "jitter", xlab=F)
#  Add p-value
p1 <- p1 + stat_compare_means()

#Prepare for GSEA

#create the rank file
keepCol<-c("logFC","P.Value")
output<-topTable(fit2, coef="GroupConverter", n=Inf, sort.by="P")[,keepCol]
rank<-sign(output$logFC)*(-1)*log10(output$P.Value)
#names(rank)<-row.names(output)

rank<-data.frame(row.names(output),rank)
colnames(rank)<-c("SYMBOL","RANK")
rank$SYMBOL<-as.factor(rank$SYMBOL)
rownames(rank)<-NULL

write.table(rank, "convrank.rnk.txt", sep = "\t", row.names = F, quote=F)
rank_conv<-rank


library(cogena)
library(fgsea)

ourGmtList<-gmt2list("/Users/carlybobak/Library/CloudStorage/GoogleDrive-carly.a.bobak@dartmouth.edu/Shared drives/Hill Research Group/STUDENT AND STAFF FILES_non project files only/Bobak/TB Transcription Integration Project/Final Data/MyCustomizedSet2.gmt")
allPathwayNames<-names(ourGmtList)

set.seed(0)
conversiongsea<-fgsea(ourGmtList,rank,eps=1e-100,minSize=50,maxSize=1000)
View(conversiongsea)
conversiongsea$magNES<-abs(conversiongsea$NES)
conversiongsea_sub<-conversiongsea[conversiongsea$pval<0.005,]
conversiongsea_sub <- conversiongsea_sub %>% arrange(desc(magNES))
#collage leading edge
for(i in 1:nrow(conversiongsea_sub)){
  conversiongsea_sub$leadingEdge[i]<-paste(conversiongsea_sub$leadingEdge[[i]],collapse = ", ")
}
conversiongsea_sub$leadingEdge<-unlist(conversiongsea_sub$leadingEdge)
write.csv(conversiongsea_sub,"conversionGSEAresults.csv",row.names=F)

## quick barplot
conversiongseadf<-conversiongsea_sub[1:5,]

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

conversiongseadf$xlab<-sapply(strsplit(conversiongseadf$pathway,"%"),`[`, 1)
conversiongseadf$xlab<-capwords(tolower(conversiongseadf$xlab))
conversiongseadf$xlab[4]<-"TCR Signaling"

#save for a bigger figure later?
ggbarplot(conversiongseadf,x="xlab",y="NES")


#all this for WGCNA
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
BiocManager::install(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
library(WGCNA)

moduleTraitCor = cor(PC1, convPheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(convExpr)) 

for(i in testFactors){
  p<-summary(aov(PC1~convPheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC1<-t(moduleTraitPvalue)

tiff("PC1Converter.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC1[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC1 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC2, convPheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(convExpr)) 

for(i in testFactors){
  p<-summary(aov(PC2~convPheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC2<-t(moduleTraitPvalue)

tiff("PC2Converter.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC2[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC2 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC3, convPheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(convExpr)) 

for(i in testFactors){
  p<-summary(aov(PC3~convPheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC3<-t(moduleTraitPvalue)

tiff("PC3Converter.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC3[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC3 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

biplot(pcs) #biggest drivers here are DEFA1 and DEFA3

write.csv(pcs$rotation[,1:3],"converterPCALoadings.csv")


################### TB disease vs not

#start with converters vs not converters
ix<-which(is.na(pheno$tbdisease_merged))
TBPheno<-pheno[-ix,]
TBExpr<-expr[-ix,]

TBExpr<-apply(TBExpr,2,as.numeric)

Group <- factor(TBPheno$tbdisease_merged, levels=c("No TB Diagnosis","TB Disease"))
Group
design <- model.matrix(~Group)
head(design)
colnames(design) <- c("No TB Diagnosis","TB Disease")
fit <- lmFit(t(TBExpr), design)


set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
head(colnames(expr),50)
tail(colnames(expr),50)


#Add covariates
RIN <- as.numeric(TBPheno$RIN)
Batch <- factor(TBPheno$Batch)
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Group+RIN+Batch)
colnames(design)

fit <- lmFit(t(TBExpr), design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
Diff<- topTable(fit2, coef="GroupTB Disease", adjust="BH", n=25000)
write.csv(Diff, file="TBDEResults.csv", quote=FALSE, row.names=TRUE)

#create volcano plots of DGE signatures

plotDiff<-Diff
plotDiff$sig<-ifelse(Diff$P.Value<0.005,"Significant","Non-significant")
plotDiff$lab<-ifelse(Diff$P.Value<0.001,"Yes","No")
plotDiff$dir<-ifelse(Diff$logFC<0,"Down","Up")
plotDiff$gene<-rownames(Diff)
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$colr<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
levels(plotDiff$colr)
plotDiff$colr<-paste(plotDiff$sig,plotDiff$dir,sep=" ")
plotDiff$colr[plotDiff$sig=="Non-significant"]<-"Non-significant"
plotDiff$colr<-as.factor(plotDiff$colr)
levels(plotDiff$colr)

library(latex2exp)
library(ggrepel)

tiff("TBVolcanoLabels.tiff",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("red3","cornflower blue","black"),labels=c("Up Regulated","Down Regulated","Non-significant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
  theme_classic()+theme(legend.title=element_blank())
volc+geom_text_repel(data=plotDiff[plotDiff$lab=="Yes",], aes(label=gene),size=3)
dev.off()

plotDiff$lab<-ifelse(Diff$P.Value<0.001,"Yes","No")

tiff("TBVolcanoLabels2.tiff",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("red3","cornflower blue", "black"),labels=c("Up Regulated","Down Regulated","Insignificant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
  theme_classic()+theme(legend.title=element_blank())
volc+geom_text_repel(data=plotDiff[plotDiff$lab=="Yes",], aes(label=gene),size=3)
dev.off()

tiff("TBVolcano.tiff",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("black","cornflower blue", "red3"),labels=c("Up Regulated","Down Regulated","Insignificant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
  theme_classic()+theme(legend.title=element_blank())
volc
dev.off()

DGEs<-row.names(Diff[Diff$P.Value<0.005,])
DGEMatrix<-TBExpr[,DGEs]

pcs<-prcomp(DGEMatrix)
summary(pcs)
tiff("PCDistTB.tiff",height=5,width=7,units="in", res=300,)
plot(pcs)
dev.off()
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]

moduleTraitCor = cor(PC1, TBPheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(TBExpr)) 

for(i in testFactors){
  p<-summary(aov(PC1~TBPheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC1<-t(moduleTraitPvalue)

tiff("PC1TB.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC1[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC1 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC2, TBPheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(TBExpr)) 

for(i in testFactors){
  p<-summary(aov(PC2~TBPheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC2<-t(moduleTraitPvalue)

tiff("PC2TB.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC2[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC2 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC3, TBPheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(TBExpr)) 

for(i in testFactors){
  p<-summary(aov(PC3~TBPheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC3<-t(moduleTraitPvalue)

tiff("PC3TB.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC3[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC3 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

write.csv(pcs$rotation[,1:3],"diseasePCALoadings.csv")

library(ggbiplot)

colors = c("red","dodger blue")
names(colors) = c("TB Disease","No TB Diagnosis")

g <- ggbiplot(pcs, choices=c(1,2),
              groups = as.factor(TBPheno$tbdisease_merged),
              circle = T, var.axes = F, ellipse=T)
g <- g + geom_point(aes(col=as.factor(TBPheno$tbdisease_merged)), size=1)
g<- g + scale_colour_manual(values = colors, name="TB Disease")
g <- g + theme(legend.direction = 'vertical',
               legend.position = 'right')
g <- g + theme_classic()

g

tiff("TBPCA_ellipse.tiff",height=6,width=6,units="in",res=600)
g
dev.off()

library(ggpubr)

PCdata<-data.frame(pcs$x[,1:3])
PCdata$Conversion<-TBPheno$conversion_10_at36months
PCdata$TBDisease<-TBPheno$tbdisease_merged
PCdata$MaternalPriorTB<-as.factor(ifelse(TBPheno$tb_everdx=="1","Prior TB","No Prior TB Diagnosis"))
PCdata$Smoking<-as.factor(ifelse(TBPheno$smoke_current=="1","Current Smoker","Current Nonsmoker"))

p1<-ggboxplot(PCdata,x="TBDisease",y="PC1", color="TBDisease",
              palette=c("red","dodger blue"),
              order=c("TB Disease","No TB Diagnosis"),
              add="jitter",xlab=F,outlier.shape=NA)


p1<-p1+stat_compare_means()+theme(legend.position = "none")

p2<-ggboxplot(PCdata[is.na(PCdata$Smoking)==F,],x="Smoking",y="PC2", color="Smoking",
              palette=c("maroon","darkgrey"),
              order=c("Current Smoker","Current Nonsmoker"),
              add="jitter",xlab=F,outlier.shape=NA)


p2<-p2+stat_compare_means()+theme(legend.position = "none")


tiff("TBDGEPCBox.tiff",height=4, width=4, units="in", res=600)
p1
dev.off()

tiff("TBDGEPCSmokingBox.tiff",height=4, width=4, units="in", res=600)
p2
dev.off()

#GSEA time


#create the rank file
output<-topTable(fit2, coef="GroupTB Disease", n=Inf, sort.by="P")[,keepCol]
rank<-sign(output$logFC)*(-1)*log10(output$P.Value)
#names(rank)<-row.names(output)

rank<-data.frame(row.names(output),rank)
colnames(rank)<-c("SYMBOL","RANK")
rank$SYMBOL<-as.factor(rank$SYMBOL)
rownames(rank)<-NULL

write.table(rank, "TBrank.rnk.txt", sep = "\t", row.names = F, quote=F)
rank_TB<-rank


set.seed(0)
TBgsea<-fgsea(ourGmtList,rank,eps=1e-100,minSize=50,maxSize=1000)
View(TBgsea)
TBgsea$magNES<-abs(TBgsea$NES)
TBgsea_sub<-TBgsea[TBgsea$pval<0.005,]
TBgsea_sub <- TBgsea_sub %>% arrange(desc(magNES))
#collage leading edge
for(i in 1:nrow(TBgsea_sub)){
  TBgsea_sub$leadingEdge[i]<-paste(TBgsea_sub$leadingEdge[[i]],collapse = ", ")
}
TBgsea_sub$leadingEdge<-unlist(TBgsea_sub$leadingEdge)
write.csv(TBgsea_sub,"TBGSEAresults.csv",row.names=F)

################### TB disease vs not within converters

#start with converters vs not converters
ix<-which(is.na(pheno$tbdisease_merged)|pheno$conversion_10_at36months=="Nonconverter")
TBPheno2<-pheno[-ix,]
TBExpr2<-expr[-ix,]

TBExpr2<-apply(TBExpr2,2,as.numeric)

Group <- factor(TBPheno2$tbdisease_merged, levels=c("No TB Diagnosis","TB Disease"))
Group
design <- model.matrix(~Group)
head(design)
colnames(design) <- c("No TB Diagnosis","TB Disease")
fit <- lmFit(t(TBExpr2), design)


set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")



#Add covariates
RIN <- as.numeric(TBPheno2$RIN)
Batch <- factor(TBPheno2$Batch)
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Group+RIN+Batch)
colnames(design)

fit <- lmFit(t(TBExpr2), design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
Diff<- topTable(fit2, coef="GroupTB Disease", adjust="BH", n=25000)
write.csv(Diff, file="TBDEResultsJustConverters.csv", quote=FALSE, row.names=TRUE)

#create volcano plots of DGE signatures

plotDiff<-Diff
plotDiff$sig<-ifelse(Diff$P.Value<0.005,"Significant","Non-significant")
plotDiff$lab<-ifelse(Diff$P.Value<0.005,"Yes","No")
plotDiff$dir<-ifelse(Diff$logFC<0,"Down","Up")
plotDiff$gene<-rownames(Diff)
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$colr<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
levels(plotDiff$colr)
plotDiff$colr<-paste(plotDiff$sig,plotDiff$dir,sep=" ")
plotDiff$colr[plotDiff$sig=="Non-significant"]<-"Non-significant"
plotDiff$colr<-as.factor(plotDiff$colr)
levels(plotDiff$colr)

library(latex2exp)
library(ggrepel)

png("TBVolcanoLabelsJustConverters.png",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("red3","cornflower blue", "black"),labels=c("Up Regulated","Down Regulated","Non-significant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("$log_2 Fold Change$"))+ylab(TeX("$-log_{10}(p-value)$"))+
  theme_classic()+theme(legend.title=element_blank(),text = element_text(size = 15))
volc+geom_text_repel(data=plotDiff[plotDiff$lab=="Yes",], aes(label=gene),size=3)
dev.off()

plotDiff$lab<-ifelse(Diff$P.Value<0.005,"Yes","No")

tiff("TBVolcanoLabelsJustConverters2.tiff",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("red3","cornflower blue", "black"),labels=c("Up Regulated","Down Regulated","Non-significant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
  theme_classic()+theme(legend.title=element_blank())
volc+geom_text_repel(data=plotDiff[plotDiff$lab=="Yes",], aes(label=gene),size=3)
dev.off()

tiff("TBVolcanoJustConverters.tiff",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("black","cornflower blue", "red3"),labels=c("Up Regulated","Down Regulated","Non-significant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
  theme_classic()+theme(legend.title=element_blank())
volc
dev.off()

DGEs<-row.names(Diff[Diff$P.Value<0.005,])
DGEMatrix<-TBExpr2[,DGEs]

pcs<-prcomp(DGEMatrix)
summary(pcs)
tiff("PCDistTB2.tiff",height=5,width=7,units="in", res=300,)
plot(pcs)
dev.off()
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]

moduleTraitCor = cor(PC1, TBPheno2[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(TBExpr2)) 

for(i in testFactors){
  if(i=="conversion_10_at36months"|i=="tb_everdx"){next} #no variability here!
  p<-summary(aov(PC1~TBPheno2[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC1<-t(moduleTraitPvalue)

tiff("PC1TBJustConverters.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC1[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC1 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC2, TBPheno2[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(TBExpr2)) 

for(i in testFactors){
  if(i=="conversion_10_at36months"|i=="tb_everdx"){next} #no variability here!
  p<-summary(aov(PC1~TBPheno2[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC2<-t(moduleTraitPvalue)

tiff("PC2TBJustConverters.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC2[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC2 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC3, TBPheno2[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(TBExpr2)) 

for(i in testFactors){
  if(i=="conversion_10_at36months"|i=="tb_everdx"){next} #no variability here!
  p<-summary(aov(PC1~TBPheno2[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}
moduleTraitPvaluePC3<-t(moduleTraitPvalue)

tiff("PC3TBJustConverters.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC3[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC3 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

write.csv(pcs$rotation[,1:3],"diseaseAmongConvertersPCALoadings.csv")

DGEs<-row.names(Diff[Diff$P.Value<0.005,])
DGEMatrix<-convExpr[,DGEs]

#pcs<-prcomp(DGEMatrix)

library(ggbiplot)

colors = c('red','goldenrod1')
names(colors) = c("TB Disease","TST Converter")

TBPheno2$Outcome<-as.character(TBPheno2$tbdisease_merged)
TBPheno2$Outcome[TBPheno2$Outcome=="No TB Diagnosis"]<-"TST Converter"
TBPheno2$Outcome<-as.factor(TBPheno2$Outcome)
  
  
g <- ggbiplot(pcs, choices=c(1,2),
              groups = as.factor(TBPheno2$Outcome),
              circle = T, var.axes = F, ellipse=T)
g <- g + geom_point(aes(col=as.factor(TBPheno2$Outcome)), size=1)
g<- g + scale_colour_manual(values = colors, name="TB Category")
g <- g + theme(legend.direction = 'vertical',
               legend.position = 'right')
g <- g + theme_classic() +theme(text = element_text(size = 15))

g

png("TBinConvertersPCA_ellipse.png",height=6,width=6,units="in",res=600)
g
dev.off()

PCdata<-data.frame(pcs$x[,1:3])
PCdata$TbCategory<-TBPheno2$Outcome
PCdata$MaternalPriorTB<-as.factor(ifelse(TBPheno2$tb_everdx=="1","Prior TB","No Prior TB Diagnosis"))
PCdata$Smoking<-as.factor(ifelse(TBPheno2$smoke_current=="1","Current Smoker","Current Nonsmoker"))

p1<-ggboxplot(PCdata,x="TbCategory",y="PC1", color="TbCategory",
              palette=c("red","goldenrod1"),
              order=c("TB Disease","TST Converter"),
              add="jitter",xlab=F,outlier.shape=NA)


p1<-p1+stat_compare_means()+theme(legend.position = "none",text = element_text(size = 15))

png("TbinConvPCAOutcome.png",height=6, width=6, units="in",res=600)
p1
dev.off()

p2<-ggboxplot(PCdata[is.na(PCdata$Smoking)==F,],x="Smoking",y="PC1", color="Smoking",
              palette=c("maroon","darkgrey"),
              order=c("Current Smoker","Current Nonsmoker"),
              add="jitter",xlab=F,outlier.shape=NA)


p2<-p2+stat_compare_means()+theme(legend.position = "none",text = element_text(size = 15))

png("TbinConvPCASmoking.png",height=6, width=6, units="in",res=600)
p2
dev.off()

df<-data.frame(TBPheno2$Outcome, PCdata$Smoking)#,convPheno$tbdisease_merged)
colnames(df)<-c("TB Category","Smoking Category")


tbCatCol<-c("red","goldenrod1")
names(tbCatCol)<-c("TB Disease","Converter")

smokeCol<-c("maroon","darkgrey")
names(smokeCol)<-c("Current Smoker","Current Nonsmoker")

top_annotation = HeatmapAnnotation(name="TB Cord Blood",df=df,
                                   col=list("TB Category"=c("TB Disease"="red","TST Converter"="goldenrod1"),
                                            "Smoking Category"=c("Current Smoker"="maroon","Current Nonsmoker"="darkgrey")),
                                   show_legend = T,
                                   annotation_legend_param = list(title_gp = gpar(fontsize = 12), 
                                                                  labels_gp = gpar(fontsize = 12)))

convExprCenter<-apply(DGEMatrix,2,function(x){x<-x-median(x)})

HM<-Heatmap(t(convExprCenter), 
            name = " ", #title of legend
            column_title = NULL, row_title = NULL,
            show_column_names=T,
            clustering_distance_columns="canberra",
            clustering_method_columns = "ward.D2",
            clustering_method_rows = "ward.D2",
            clustering_distance_rows = "canberra",
            row_names_gp = gpar(fontsize = 8),
            top_annotation = top_annotation,
            show_heatmap_legend = T,
            show_row_dend = F,
            heatmap_legend_param = list(color_bar = "continuous",labels_gp=gpar(fontsize=12)),
            #cluster_rows = color_branches(row_dend, k = 3),
            col=mycols
)

tiff("TbinConversionHeatmap.tiff",height=10,width=10,units="in",res=900)
draw(HM,  heatmap_legend_side = "bottom", annotation_legend_side="bottom")
dev.off()

#create the rank file
output<-topTable(fit2, coef="GroupTB Disease", n=Inf, sort.by="P")[,keepCol]
rank<-sign(output$logFC)*(-1)*log10(output$P.Value)
#names(rank)<-row.names(output)

rank<-data.frame(row.names(output),rank)
colnames(rank)<-c("SYMBOL","RANK")
rank$SYMBOL<-as.factor(rank$SYMBOL)
rownames(rank)<-NULL

write.table(rank, "TBinconvrank.rnk.txt", sep = "\t", row.names = F, quote=F)
rank_TBinconv<-rank


set.seed(0)
TBinConvgsea<-fgsea(ourGmtList,rank,eps=1e-100,minSize=50,maxSize=1000)
TBinConvgsea$magNES<-abs(TBinConvgsea$NES)
TBinConvgsea_sub<-TBinConvgsea[TBinConvgsea$pval<0.005,]
TBinConvgsea_sub <- TBinConvgsea_sub %>% arrange(desc(magNES))
#collage leading edge
for(i in 1:nrow(TBinConvgsea_sub)){
  TBinConvgsea_sub$leadingEdge[i]<-paste(TBinConvgsea_sub$leadingEdge[[i]],collapse = ", ")
}
TBinConvgsea_sub$leadingEdge<-unlist(TBinConvgsea_sub$leadingEdge)
write.csv(TBinConvgsea_sub,"TBinConvGSEAresults.csv",row.names=F)








###################################pneumonia

expr<-apply(expr,2,as.numeric)

Group <- factor(pheno$pneumonia_any, levels=c("No Pneumonia","Pneumonia"))
Group
design <- model.matrix(~Group)
head(design)
colnames(design) <- c("Nonconverter","Converter")
fit <- lmFit(t(expr), design)


set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
head(colnames(expr),50)
tail(colnames(expr),50)


#Add covariates
RIN <- as.numeric(pheno$RIN)
Batch <- factor(pheno$Batch)
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Group+RIN+Batch)
colnames(design)

fit <- lmFit(t(expr), design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
Diff<- topTable(fit2, coef="GroupPneumonia", adjust="BH", n=25000)
write.csv(Diff, file="PneumoniaDEResults.csv", quote=FALSE, row.names=TRUE)

#create volcano plots of DGE signatures

plotDiff<-Diff
plotDiff$sig<-ifelse(Diff$P.Value<0.005,"Significant","Insignificant")
plotDiff$lab<-ifelse(Diff$P.Value<0.001,"Yes","No")
plotDiff$dir<-ifelse(Diff$logFC<0,"Down","Up")
plotDiff$gene<-rownames(Diff)
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$col<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
plotDiff$colr<-as.factor(paste(plotDiff$sig,plotDiff$dif,sep=" "))
levels(plotDiff$colr)
plotDiff$colr<-paste(plotDiff$sig,plotDiff$dir,sep=" ")
plotDiff$colr[plotDiff$sig=="Insignificant"]<-"Insignificant"
plotDiff$colr<-as.factor(plotDiff$colr)
levels(plotDiff$colr)

library(latex2exp)
library(ggrepel)

tiff("pneumoniaVolcanoLabels.tiff",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("black","cornflower blue", "red3"),labels=c("Up Regulated","Down Regulated","Insignificant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
  theme_classic()+theme(legend.title=element_blank())
volc+geom_text_repel(data=plotDiff[plotDiff$lab=="Yes",], aes(label=gene),size=3)
dev.off()

tiff("pneumoniaVolcano.tiff",width=6,height=6,units="in",res=600)
volc = ggplot(plotDiff, aes(logFC, -log10(P.Value))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=colr),size=1) + #add points colored by significance
  scale_color_manual(values=c("black","cornflower blue", "red3"),labels=c("Up Regulated","Down Regulated","Insignificant"),breaks = rev(levels(plotDiff$colr)))+
  geom_hline(yintercept = -log10(0.005),colour="red",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.005)+0.075,label="p-value=0.005",colour="red",size=2)+
  geom_hline(yintercept = -log10(0.001),colour="orange",linetype="dashed")+
  annotate("text",x=-0.78,y=-log10(0.001)+0.075,label="p-value=0.001",colour="orange",size=2)+
  xlab(TeX("log_2 Fold Change"))+ylab(TeX("-log_{10}(p-value)"))+
  theme_classic()+theme(legend.title=element_blank())
volc
dev.off()

#check the PC components for correlations with variables approaching significants from Table 1
testFactors<-c("site","mat_hiv_status","smoke_current","prem","lbw","pneumonia_any","sex","conversion_10_at36months","tb_household","NEW_tb_momonrx","tb_everdx","tbdisease_merged")
testContinuous<-c("gestation_delivery","birthweight_gms","Breast_feeding_months","pneum_count")


DGEs<-row.names(Diff[Diff$P.Value<0.005,])
DGEMatrix<-expr[,DGEs]

pcs<-prcomp(DGEMatrix)
tiff("PCDistPneumonia.tiff",height=5,width=7,units="in", res=300,)
plot(pcs)
dev.off()
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]

moduleTraitCor = cor(PC1, pheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(expr)) 

for(i in testFactors){
  p<-summary(aov(PC1~pheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC1<-t(moduleTraitPvalue)

tiff("PC1Pneumonia.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC1[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC1 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC2, pheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(expr)) 

for(i in testFactors){
  p<-summary(aov(PC2~pheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC2<-t(moduleTraitPvalue)

tiff("PC2Pneumonia.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC2[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC2 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()

moduleTraitCor = cor(PC3, pheno[,testContinuous], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(expr)) 

for(i in testFactors){
  p<-summary(aov(PC3~pheno[,i]))
  nm<-colnames(moduleTraitPvalue)
  nm<-append(nm,i)
  moduleTraitPvalue<-cbind(moduleTraitPvalue,p[][[1]]$`Pr(>F)`[1])
  colnames(moduleTraitPvalue)<-nm
}

moduleTraitPvaluePC3<-t(moduleTraitPvalue)

tiff("PC3Pneumonia.tiff",height=8,width=12, units="in", res=600)
par(mar=c(5, 8.5, 4, 2),xpd=F)
barplot(-log10(moduleTraitPvaluePC3[,1]), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(moduleTraitPvaluePC1),
        xlab="DGE PC3 association P-value (-log10) ")
abline(v=-log10(0.05), col="red")
dev.off()
