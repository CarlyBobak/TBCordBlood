#analysing ciberdata results

ciberResults <- read.csv("~/Dropbox (Dartmouth College)/Drakenstein, Carly:Leo/cibersortResults/CIBERSORTx_Job3_Results.csv")

dim(ciberResults)

ciberResults<-data.frame(ciberResults,phenoMod[-which(is.na(phenoMod$conversion_10_at36months)),])
ciberResults$converterTB<-MEdf$converterTB
ciberResults$tb_everdx<-MEdf$tb_everdx

comps<-list(c("TB Disease","Converter: No TB Disease"),c("TB Disease","Nonconverter"),c("Converter: No TB Disease","Nonconverter"))


pdf("ciberSortResultsAbsolute_anova.pdf")

for(i in 2:23){
  
  p4<-ggboxplot(ciberResults,x="converterTB",y=colnames(ciberResults)[i],color="converterTB",
                palette = c('red','goldenrod1','blue3'),
                order=c("TB Disease","Converter: No TB Disease","Nonconverter"),
                add = "jitter", xlab=F, outlier.shape=NA)
  p4 <- p4 +stat_compare_means(method="wilcox",comparisons = comps)+ theme(legend.position = "none") 
  
  print(p4)
}

dev.off()

pdf("ciberSortResultsPriorTB.pdf")

for(i in 2:23){
  
  p4<-ggboxplot(ciberResults,x="tb_everdx",y=colnames(ciberResults)[i],color="tb_everdx",
                palette = c("turquoise","salmon"),
                order=c("Prior TB","No Prior TB Diagnosis"),
                add = "jitter", xlab="Maternal Prior TB Diagnosis", outlier.shape=NA)
  p4 <- p4 +stat_compare_means(method="t.test")+ theme(legend.position = "none") 
  
  print(p4)
}

dev.off()

i=23
t.test(ciberResults[ciberResults$converterTB=="TB Disease",i],ciberResults[ciberResults$converterTB=="Converter: No TB Disease",i])
t.test(ciberResults[ciberResults$converterTB=="TB Disease",i],ciberResults[ciberResults$converterTB=="Nonconverter",i])
t.test(ciberResults[ciberResults$converterTB=="Converter: No TB Disease",i],ciberResults[ciberResults$converterTB=="Nonconverter",i])

#try combining groups

ciberResults$BCells<-ciberResults$B.cells.naive+ciberResults$B.cells.memory
ciberResults$CD4<-ciberResults$T.cells.CD4.naive+ciberResults$T.cells.CD4.memory.resting+ciberResults$T.cells.CD4.memory.activated+ciberResults$T.cells.follicular.helper+ciberResults$T.cells.regulatory..Tregs.
ciberResults$Dendritic<-ciberResults$Dendritic.cells.activated+ciberResults$Dendritic.cells.resting
ciberResults$Macrophages<-ciberResults$Macrophages.M0+ciberResults$Macrophages.M1+ciberResults$Macrophages.M2
ciberResults$NK<-ciberResults$NK.cells.activated+ciberResults$NK.cells.resting
ciberResults$Mast<-ciberResults$Mast.cells.activated+ciberResults$Mast.cells.resting
ciberResults$Other<-ciberResults$Plasma.cells+ciberResults$Eosinophils

#KEEP MONOCYTES, NEUTROPHILS, CD8, AND GAMMA DELTA

t.test(ciberResults[ciberResults$converterTB=="TB Disease","Other"],ciberResults[ciberResults$converterTB=="Converter: No TB Disease","Other"])
t.test(ciberResults[ciberResults$converterTB=="TB Disease","Other"],ciberResults[ciberResults$converterTB=="Nonconverter","Other"])
t.test(ciberResults[ciberResults$converterTB=="Converter: No TB Disease","Other"],ciberResults[ciberResults$converterTB=="Nonconverter","Other"])


## TRY TOAST

BiocManager::install("TOAST")
library(TOAST)

CellType <- list(Bcells = 1,
                 CD8T = 2,
                 CD4T = 3,
                 NK = 4,
                 Monocytes = 5,
                 Neutrophils = 6)
## choose (up to 20) significant markers 
## per cell type
myMarker <- ChooseMarker(LM_22, 
                         CellType, 
                         nMarkCT = 20,
                         chooseSig = TRUE,
                         verbose = FALSE)
lapply(myMarker, head, 3)
  
