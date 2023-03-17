# Validation with Anderson data

library(WGCNA)
options(stringsAsFactors = FALSE);

metaTBFolder<-"/Users/carlybobak/Library/CloudStorage/GoogleDrive-carly.a.bobak@dartmouth.edu/Shared drives/Hill Research Group/STUDENT AND STAFF FILES_non project files only/Bobak/TB Transcription Integration Project/Final Data"
allTBData<-readRDS(paste(metaTBFolder,"tidiedProcessedExpression.RData",sep="/"))
allTBPhenoData<-readRDS(paste(metaTBFolder,"completeIntegrationObject.RData",sep="/"))
allTBPhenoData<-allTBPhenoData$pheno

#Kenya
GSE39939<-t(allTBData$GSE39939)

#Malawi
GSE39940<-t(allTBData$GSE39940)

#previously loaded don't rerun won't work now
GSE39939_pheno<-allTBPhenoData[allTBPhenoData$dataset=="GSE39939",]
GSE39940_pheno<-allTBPhenoData[allTBPhenoData$dataset=="GSE39940",]

#validate conversion,TB, TB in conversion, and blue module genes in both Anderson sets
convGenes<-c("UBAP2L",	"CD37",	"LOC100130552",	"PEX6",	"FLJ37786",	"LOC644214",	"LOC100129466",	"PRB3",	"RPS6KB2",	"KIF6",	"RTEL1",	"KIAA0495",	"CXorf48",	"LOC100133678",	"JARID2",	"SURF6",	"IWS1",	"FLJ23584",	"INTS7",	"C19orf28",	"ZNF622",	"ARHGEF1",	"LOC284422",	"XRCC5",	"LOC644037",	"SIRPD",	"FAM3A",	"CD74",	"LOC345041",	"MYL5",	"ZNF362",	"CLEC12A",	"CUL4B",	"HLA-DOA",	"SOCS7",	"LOC100129795",	"MED14",	"SUGT1L1",	"ARHGAP15",	"FLJ37307",	"MGC72080",	"HLA-DMA",	"HLA-DQA1",	"ADPRH",	"MLST8",	"FAM65A",	"LOC653080",	"GHDC",	"STRN",	"LOC647450",	"FBXL11",	"SPNS3",	"C14orf121",	"PABPN1",	"ROGDI",	"T-SP1",	"SUSD5",	"LRPPRC",	"DEFA1",	"ZNF653",	"SP110",	"BAZ2A",	"TTLL1",	"LYL1",	"LOC399491",	"LOC642113",	"TMF1",	"PALMD",	"PPIG",	"PLUNC",	"C17orf103",	"LOC642367",	"TWSG1",	"LOC389293",	"LOC391692",	"NSUN5B",	"API5",	"CD38",	"HLA-DPA1",	"LOC388248",	"C3orf10",	"UBE2F",	"LGALS9",	"SPSB3",	"MGC52498",	"LOC648749",	"USPL1",	"LPAR3",	"PSMB2",	"NIPBL",	"LOC100131381",	"IRF2BP2",	"TEX11",	"CCDC6",	"C7orf70",	"VPREB3",	"HLA-DRB4",	"MIR106A",	"LOC100129141",	"ZBTB44",	"DEFA3",	"LOC441642",	"LOC387753",	"BTN3A2",	"IFITM3",	"CD8B",	"JARID1A",	"C6orf176",	"CDYL2",	"ATPBD3",	"LOC643873",	"TRAPPC6B",	"SCD",	"LOC727732")
#71 out of 114

TBGenes<-c("SULT1A3",	"LOC643357",	"HMBS",	"LOC100130561",	"NCOA3",	"CYBA",	"CSTB",	"ASXL1",	"HIST1H2AC",	"ACTG1",	"CD37",	"YY1",	"PSMC1",	"LOC440589",	"LOC390466",	"FAM21D",	"NOP14",	"PRKAB1",	"PARP1",	"LOC643319",	"SNURF",	"RPL36AL",	"ARPC4",	"LOC643176",	"ICOS",	"GGA2",	"LOC389342",	"VNN2",	"C1orf41",	"MAPK1",	"STAG3",	"PPP1CA",	"LOC651143",	"TACC3",	"LOC100134537",	"KIAA1949",	"HERC6",	"TMEM101",	"ARPC1B",	"ZBTB44",	"BANF1",	"ARRB2",	"ARHGEF11",	"BRD9",	"VPS11",	"EWSR1",	"LGALS9",	"RPLP0",	"ACRC",	"GRWD1",	"HOXA6",	"LOC644153",	"FBXO7",	"KLRD1",	"LOC641750",	"EFS",	"HIST1H2BD",	"ZNF695",	"LOC646791",	"PLXDC1")
#45 out of 60

TBinConvGenes<-c("PARP1",	"WDR4",	"KLRD1",	"YY1",	"BCL11A",	"N4BP2L1",	"RPL36AL",	"NCOA3",	"LOC100131989",	"FCRLA",	"KLF12",	"SEC61A1",	"MYO9B",	"TUBB1",	"MAPK1",	"TMEM77",	"LOC100132717",	"BIRC3",	"EWSR1",	"MGEA5",	"C18orf10",	"GGA2",	"LOC92755",	"C2orf25",	"ARHGEF7",	"NARS",	"LOC728014",	"SH3BGRL3",	"TRIM38",	"DIS3L")
#23 out of 30 in Anderson

royalBlueGenes #exists from WGCNA code
#32 out of 37 genes available in Anderson

AllGenesForValidation<-unique(c(convGenes,TBGenes,TBinConvGenes,royalBlueGenes)) #228 genes
AllGenesForValidation<-AllGenesForValidation[AllGenesForValidation%in%colnames(GSE39939)] #158 genes

#validate on LTBI vs OD, TB vs LTBI, TB vs OD

GSE39939_val<-t(GSE39939[,AllGenesForValidation])
GSE39940_val<-t(GSE39940[,AllGenesForValidation])

#LTBI vs OD
GSE39939_val_LTBI<-GSE39939_val[,GSE39939_pheno$disease%in%c("LTBI","Other Disease")]
GSE39940_val_LTBI<-GSE39940_val[,GSE39940_pheno$disease%in%c("LTBI","Other Disease")]

Outcome <- factor(GSE39939_pheno$disease[GSE39939_pheno$disease%in%c("LTBI","Other Disease")], levels=c("Other Disease","LTBI"))
HIV <- as.factor(GSE39939_pheno$HIV[GSE39939_pheno$disease%in%c("LTBI","Other Disease")])
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Outcome)
colnames(design)

fit <- lmFit(GSE39939_val_LTBI, design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
GSE39939_Diff_LTBI<- topTable(fit2, coef="OutcomeLTBI", adjust="BH", n=Inf)
write.csv(GSE39939_Diff_LTBI, file="GSE39939_LTBIvsOD.csv", quote=FALSE, row.names=TRUE)

cat(paste(rownames(GSE39939_Diff_LTBI[GSE39939_Diff_LTBI$adj.P.Val<0.05,]),collapse="\n"))

Outcome <- factor(GSE39940_pheno$disease[GSE39940_pheno$disease%in%c("LTBI","Other Disease")], levels=c("Other Disease","LTBI"))
HIV <- as.factor(GSE39940_pheno$HIV[GSE39940_pheno$disease%in%c("LTBI","Other Disease")])
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Outcome)
colnames(design)

fit <- lmFit(GSE39940_val_LTBI, design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
GSE39940_Diff_LTBI<- topTable(fit2, coef="OutcomeLTBI", adjust="BH", n=Inf)
write.csv(GSE39940_Diff_LTBI, file="GSE39940_LTBIvsOD.csv", quote=FALSE, row.names=TRUE)

cat(paste(rownames(GSE39940_Diff_LTBI[GSE39939_Diff_LTBI$adj.P.Val<0.05,]),collapse="\n"))

#conversion result data frame

conv_val<-convResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
conv_val$GSE39939FC<-GSE39939_Diff_LTBI[rownames(conv_val),"logFC"]
conv_val$GSE39939Pval<-GSE39939_Diff_LTBI[rownames(conv_val),"P.Value"]
conv_val$GSE39939AdjPval<-GSE39940_Diff_LTBI[rownames(conv_val),"adj.P.Val"]
conv_val$GSE39940FC<-GSE39940_Diff_LTBI[rownames(conv_val),"logFC"]
conv_val$GSE39940Pval<-GSE39940_Diff_LTBI[rownames(conv_val),"P.Value"]
conv_val$GSE39940AdjPval<-GSE39940_Diff_LTBI[rownames(conv_val),"adj.P.Val"]
View(conv_val)

write.csv(conv_val,"TSTConversionValidation.csv,",row.names = T)

conv_val_plot<-convResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
conv_val_plot2<-convResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
conv_val_plot$val_FC<-GSE39939_Diff_LTBI[rownames(conv_val),"logFC"]
conv_val_plot$val_P<-GSE39939_Diff_LTBI[rownames(conv_val),"P.Value"]
conv_val_plot$val_AdjP<-GSE39940_Diff_LTBI[rownames(conv_val),"adj.P.Val"]
conv_val_plot2$val_FC<-GSE39940_Diff_LTBI[rownames(conv_val),"logFC"]
conv_val_plot2$val_P<-GSE39940_Diff_LTBI[rownames(conv_val),"P.Value"]
conv_val_plot2$val_AdjP<-GSE39940_Diff_LTBI[rownames(conv_val),"adj.P.Val"]

conv_val_plot<-rbind(conv_val_plot,conv_val_plot2)
conv_val_plot$val_country<-c(rep("Kenya",114),rep("Malawi",114))
conv_val_plot$label_gene<-ifelse(abs(conv_val_plot$logFC)>0.5|abs(conv_val_plot$val_FC)>0.5,1,0)
conv_val_plot$sig<-ifelse(conv_val_plot$val_AdjP<0.05,"Significant","Not Significant")
conv_val_plot$symbol<-rep(row.names(conv_val),2)

p = ggplot(na.omit(conv_val_plot), aes(logFC, val_FC)) + 
  geom_point(aes(col=val_country,shape=sig),size=1) + 
  scale_color_manual(values=c("gold","turquoise"),labels=c("Kenya","Malawi"))+
  geom_hline(yintercept = 0,colour="grey",linetype="dashed")+
  geom_vline(xintercept=0,colour="grey",linetype="dashed")+
  xlab("Cord Blood TST Conversion logFC")+ylab("Anderson LTBI vs Other Disease logFC")+
  theme_classic()+theme(legend.title=element_blank(), text = element_text(size = 10))

tiff("convergence_validation.tiff",height=5,width=5,units="in",res=400)
p+geom_text_repel(data=conv_val_plot[conv_val_plot$label_gene==1,], aes(label=symbol),size=3)
dev.off()

validated_LTBI_OD<-unique(conv_val_plot$symbol[conv_val_plot$val_AdjP<0.05])

# TB vs LTBI

GSE39939_val_TBLTBI<-GSE39939_val[,GSE39939_pheno$disease%in%c("LTBI","TB")]
GSE39940_val_TBLTBI<-GSE39940_val[,GSE39940_pheno$disease%in%c("LTBI","TB")]

Outcome <- factor(GSE39939_pheno$disease[GSE39939_pheno$disease%in%c("LTBI","TB")], levels=c("LTBI","TB"))
HIV <- as.factor(GSE39939_pheno$HIV[GSE39939_pheno$disease%in%c("LTBI","TB")])
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Outcome)
colnames(design)

fit <- lmFit(GSE39939_val_TBLTBI, design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
GSE39939_Diff_TBLTBI<- topTable(fit2, coef="OutcomeTB", adjust="BH", n=Inf)
write.csv(GSE39939_Diff_TBLTBI, file="GSE39939_LTBIvsTB.csv", quote=FALSE, row.names=TRUE)

cat(paste(rownames(GSE39939_Diff_TBLTBI[GSE39939_Diff_TBLTBI$adj.P.Val<0.05,]),collapse="\n"))

Outcome <- factor(GSE39940_pheno$disease[GSE39940_pheno$disease%in%c("LTBI","TB")], levels=c("LTBI","TB"))
HIV <- as.factor(GSE39940_pheno$HIV[GSE39940_pheno$disease%in%c("LTBI","TB")])
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Outcome)
colnames(design)

fit <- lmFit(GSE39940_val_TBLTBI, design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
GSE39940_Diff_TBLTBI<- topTable(fit2, coef="OutcomeTB", adjust="BH", n=Inf)
write.csv(GSE39940_Diff_TBLTBI, file="GSE39940_LTBIvsTB.csv", quote=FALSE, row.names=TRUE)

cat(paste(rownames(GSE39940_Diff_TBLTBI[GSE39939_Diff_TBLTBI$adj.P.Val<0.05,]),collapse="\n"))

#TB result data frame

#read in TB results

TBResults<-read.csv(file="TBDEResults.csv",row.names=1)

TB_val<-TBResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val$GSE39939FC<-GSE39939_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val$GSE39939Pval<-GSE39939_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val$GSE39939AdjPval<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]
TB_val$GSE39940FC<-GSE39940_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val$GSE39940Pval<-GSE39940_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val$GSE39940AdjPval<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]
View(TB_val)

write.csv(TB_val,"TBValidation_TBvsLTBI.csv,",row.names = T)

TB_val_plot<-TBResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot2<-TBResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot$val_FC<-GSE39939_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val_plot$val_P<-GSE39939_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val_plot$val_AdjP<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]
TB_val_plot2$val_FC<-GSE39940_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val_plot2$val_P<-GSE39940_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val_plot2$val_AdjP<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]

TB_val_plot<-rbind(TB_val_plot,TB_val_plot2)
TB_val_plot$val_country<-c(rep("Kenya",60),rep("Malawi",60))
TB_val_plot$label_gene<-ifelse(abs(TB_val_plot$logFC)>0.5|abs(TB_val_plot$val_FC)>0.5,1,0)
TB_val_plot$sig<-ifelse(TB_val_plot$val_AdjP<0.05,"Significant","Not Significant")
TB_val_plot$symbol<-rep(row.names(TB_val),2)

p = ggplot(na.omit(TB_val_plot), aes(logFC, val_FC)) + 
  geom_point(aes(col=val_country,shape=sig),size=1) + 
  scale_color_manual(values=c("gold","turquoise"),labels=c("Kenya","Malawi"))+
  geom_hline(yintercept = 0,colour="grey",linetype="dashed")+
  geom_vline(xintercept=0,colour="grey",linetype="dashed")+
  xlab("Cord Blood TB Disease logFC")+ylab("Anderson TB vs LTBI logFC")+
  theme_classic()+theme(legend.title=element_blank(), text = element_text(size = 10))

tiff("TBDis_TBLTBI_validation.tiff",height=5,width=5,units="in",res=400)
p+geom_text_repel(data=TB_val_plot[TB_val_plot$label_gene==1,], aes(label=symbol),size=3)
dev.off()

validated_TB_LTBI<-unique(TB_val_plot$symbol[TB_val_plot$val_AdjP<0.05])

#read in TB in Conv results

TBinConvResults<-read.csv(file="TBDEResultsJustConverters.csv",row.names=1)

TB_val<-TBinConvResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val$GSE39939FC<-GSE39939_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val$GSE39939Pval<-GSE39939_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val$GSE39939AdjPval<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]
TB_val$GSE39940FC<-GSE39940_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val$GSE39940Pval<-GSE39940_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val$GSE39940AdjPval<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]
View(TB_val)

write.csv(TB_val,"TBinConversionValidation_TBvsLTBI.csv,",row.names = T)

TB_val_plot<-TBinConvResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot2<-TBinConvResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot$val_FC<-GSE39939_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val_plot$val_P<-GSE39939_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val_plot$val_AdjP<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]
TB_val_plot2$val_FC<-GSE39940_Diff_TBLTBI[rownames(TB_val),"logFC"]
TB_val_plot2$val_P<-GSE39940_Diff_TBLTBI[rownames(TB_val),"P.Value"]
TB_val_plot2$val_AdjP<-GSE39940_Diff_TBLTBI[rownames(TB_val),"adj.P.Val"]

TB_val_plot<-rbind(TB_val_plot,TB_val_plot2)
TB_val_plot$val_country<-c(rep("Kenya",30),rep("Malawi",30))
TB_val_plot$label_gene<-ifelse(abs(TB_val_plot$logFC)>0.5|abs(TB_val_plot$val_FC)>0.5,1,0)
TB_val_plot$sig<-ifelse(TB_val_plot$val_AdjP<0.05,"Significant","Not Significant")
TB_val_plot$symbol<-rep(row.names(TB_val),2)

p = ggplot(na.omit(TB_val_plot), aes(logFC, val_FC)) + 
  geom_point(aes(col=val_country,shape=sig),size=1) + 
  scale_color_manual(values=c("gold","turquoise"),labels=c("Kenya","Malawi"))+
  geom_hline(yintercept = 0,colour="grey",linetype="dashed")+
  geom_vline(xintercept=0,colour="grey",linetype="dashed")+
  xlab("Cord Blood TB Disease logFC")+ylab("Anderson TB vs LTBI logFC")+
  theme_classic()+theme(legend.title=element_blank(), text = element_text(size = 10))

tiff("TBinConv_TBLTBI_validation.tiff",height=5,width=5,units="in",res=400)
p+geom_text_repel(data=TB_val_plot[TB_val_plot$label_gene==1,], aes(label=symbol),size=3)
dev.off()

validated_TB_LTBI<-unique(append(validated_TB_LTBI,unique(TB_val_plot$symbol[TB_val_plot$val_AdjP<0.05])))

# TB vs OD

GSE39939_val_TBOD<-GSE39939_val[,GSE39939_pheno$disease%in%c("Other Disease","TB")]
GSE39940_val_TBOD<-GSE39940_val[,GSE39940_pheno$disease%in%c("Other Disease","TB")]

Outcome <- factor(GSE39939_pheno$disease[GSE39939_pheno$disease%in%c("Other Disease","TB")], levels=c("Other Disease","TB"))
HIV <- as.factor(GSE39939_pheno$HIV[GSE39939_pheno$disease%in%c("Other Disease","TB")])
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Outcome)
colnames(design)

fit <- lmFit(GSE39939_val_TBOD, design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
GSE39939_Diff_TBOD<- topTable(fit2, coef="OutcomeTB", adjust="BH", n=Inf)
write.csv(GSE39939_Diff_TBOD, file="GSE39939_ODvsTB.csv", quote=FALSE, row.names=TRUE)

cat(paste(rownames(GSE39939_Diff_TBOD[GSE39939_Diff_TBOD$adj.P.Val<0.05,]),collapse="\n"))

Outcome <- factor(GSE39940_pheno$disease[GSE39940_pheno$disease%in%c("Other Disease","TB")], levels=c("Other Disease","TB"))
HIV <- as.factor(GSE39940_pheno$HIV[GSE39940_pheno$disease%in%c("Other Disease","TB")])
#check PC components later: smoke, sex, site, HIV, pneumonia

design <- model.matrix(~Outcome)
colnames(design)

fit <- lmFit(GSE39940_val_TBOD, design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
GSE39940_Diff_TBOD<- topTable(fit2, coef="OutcomeTB", adjust="BH", n=Inf)
write.csv(GSE39940_Diff_TBOD, file="GSE39940_ODvsTB.csv", quote=FALSE, row.names=TRUE)

cat(paste(rownames(GSE39940_Diff_TBOD[GSE39939_Diff_TBOD$adj.P.Val<0.05,]),collapse="\n"))

#TB result data frame

TB_val<-TBResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val$GSE39939FC<-GSE39939_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val$GSE39939Pval<-GSE39939_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val$GSE39939AdjPval<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]
TB_val$GSE39940FC<-GSE39940_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val$GSE39940Pval<-GSE39940_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val$GSE39940AdjPval<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]
View(TB_val)

write.csv(TB_val,"TBValidation_TBvsOD.csv,",row.names = T)

TB_val_plot<-TBResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot2<-TBResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot$val_FC<-GSE39939_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val_plot$val_P<-GSE39939_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val_plot$val_AdjP<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]
TB_val_plot2$val_FC<-GSE39940_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val_plot2$val_P<-GSE39940_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val_plot2$val_AdjP<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]

TB_val_plot<-rbind(TB_val_plot,TB_val_plot2)
TB_val_plot$val_country<-c(rep("Kenya",60),rep("Malawi",60))
TB_val_plot$label_gene<-ifelse(abs(TB_val_plot$logFC)>0.15&abs(TB_val_plot$val_FC)>0.15,1,0)
TB_val_plot$sig<-ifelse(TB_val_plot$val_AdjP<0.05,"Significant","Not Significant")
TB_val_plot$symbol<-rep(row.names(TB_val),2)

p = ggplot(na.omit(TB_val_plot), aes(logFC, val_FC)) + 
  geom_point(aes(col=val_country,shape=sig),size=1) + 
  scale_color_manual(values=c("gold","turquoise"),labels=c("Kenya","Malawi"))+
  geom_hline(yintercept = 0,colour="grey",linetype="dashed")+
  geom_vline(xintercept=0,colour="grey",linetype="dashed")+
  xlab("Cord Blood TB Disease logFC")+ylab("Anderson TB vs Other Disease logFC")+
  theme_classic()+theme(legend.title=element_blank(), text = element_text(size = 10))

tiff("TBDis_TBOD_validation.tiff",height=5,width=5,units="in",res=400)
p+geom_text_repel(data=TB_val_plot[TB_val_plot$label_gene==1,], aes(label=symbol),size=3)
dev.off()

validated_TB_OD<-unique(TB_val_plot$symbol[TB_val_plot$val_AdjP<0.05])

#read in TB in Conv results

TB_val<-TBinConvResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val$GSE39939FC<-GSE39939_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val$GSE39939Pval<-GSE39939_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val$GSE39939AdjPval<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]
TB_val$GSE39940FC<-GSE39940_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val$GSE39940Pval<-GSE39940_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val$GSE39940AdjPval<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]
View(TB_val)

write.csv(TB_val,"TBinConversionValidation_TBvsOD.csv,",row.names = T)

TB_val_plot<-TBinConvResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot2<-TBinConvResults %>% select(logFC,P.Value,adj.P.Val) %>% filter(P.Value<0.005)
TB_val_plot$val_FC<-GSE39939_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val_plot$val_P<-GSE39939_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val_plot$val_AdjP<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]
TB_val_plot2$val_FC<-GSE39940_Diff_TBOD[rownames(TB_val),"logFC"]
TB_val_plot2$val_P<-GSE39940_Diff_TBOD[rownames(TB_val),"P.Value"]
TB_val_plot2$val_AdjP<-GSE39940_Diff_TBOD[rownames(TB_val),"adj.P.Val"]

TB_val_plot<-rbind(TB_val_plot,TB_val_plot2)
TB_val_plot$val_country<-c(rep("Kenya",30),rep("Malawi",30))
TB_val_plot$label_gene<-ifelse(abs(TB_val_plot$logFC)>0.5|abs(TB_val_plot$val_FC)>0.5|TB_val_plot$val_AdjP<0.05,1,0)
TB_val_plot$sig<-ifelse(TB_val_plot$val_AdjP<0.05,"Significant","Not Significant")
TB_val_plot$symbol<-rep(row.names(TB_val),2)

p = ggplot(na.omit(TB_val_plot), aes(logFC, val_FC)) + 
  geom_point(aes(col=val_country,shape=sig),size=1) + 
  scale_color_manual(values=c("gold","turquoise"),labels=c("Kenya","Malawi"))+
  geom_hline(yintercept = 0,colour="grey",linetype="dashed")+
  geom_vline(xintercept=0,colour="grey",linetype="dashed")+
  xlab("Cord Blood TB Disease logFC")+ylab("Anderson TB vs LTBI logFC")+
  theme_classic()+theme(legend.title=element_blank(), text = element_text(size = 10))

tiff("TBinConv_TBOD_validation.tiff",height=5,width=5,units="in",res=400)
p+geom_text_repel(data=TB_val_plot[TB_val_plot$label_gene==1,], aes(label=symbol),size=3)
dev.off()

validated_TB_OD<-unique(append(validated_TB_OD,unique(TB_val_plot$symbol[TB_val_plot$val_AdjP<0.05])))

allTBGenes<-unique(c(convGenes,TBGenes,TBinConvGenes))
allTBGenes<-allTBGenes[allTBGenes%in%colnames(GSE39939)]
length(allTBGenes)

allValidatedGene<-unique(c(validated_LTBI_OD,validated_TB_OD,validated_TB_LTBI))
intersect(allTBGenes,allValidatedGene)


#list expressions, samples X genes
cordBlood<-list()
cordBlood$data<-expr
Anderson1<-list()
Anderson1$data<-GSE39939
Anderson2<-list()
Anderson2$data<-GSE39940

multiExpr<-list(cordBlood=cordBlood,Anderson1=Anderson1,Anderson2=Anderson2)

#reference colour, list of each gene and its corresponding colour
names(moduleColors)<-Name

multiColor<-list(cordBlood=moduleColors)

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
 # Save the result

saveRDS(mp,"modulePreservation.Rdata")

# check kME connectivity vs significance in M11 and TST vs LTBI
# Both GSE39939 and GSE39940

whichmodule="MEroyalblue";

#kWithin find connectivity of royal blue genes that are in validation
#Match to correlation in conversion

GS.GSE39939.LTBI=as.numeric(cor(GSE39939,as.numeric(factor(GSE39939_pheno$disease,levels=c("LTBI","TB","Other Disease"))),use="p"))
names(GS.GSE39939.LTBI)<-colnames(GSE39939)
connectivity<-CM_Pre$kWithin
names(connectivity)<-colnames(expr)

M11Validate<-royalBlueGenes[royalBlueGenes%in%colnames(GSE39939)]

tiff("KenyaConnectivity.tiff",height=5,width=5,units="in",res=300)
verboseScatterplot(connectivity[M11Validate],GS.GSE39939.LTBI[M11Validate],col="royalblue",
                      xlab="M11 Gene Connectivity (kME)",ylab="Anderson Kenya G.S.",pch=20)
abline(lm(GS.GSE39939.LTBI[M11Validate]~connectivity[M11Validate]))
dev.off()


GS.GSE39940.LTBI=as.numeric(cor(GSE39940,as.numeric(factor(GSE39940_pheno$disease,levels=c("LTBI","Other Disease","TB"))),use="p"))
names(GS.GSE39940.LTBI)<-colnames(GSE39940)

tiff("MalawiConnecivity.tiff",height=5,width=5,units="in",res=300)
verboseScatterplot(connectivity[M11Validate],GS.GSE39940.LTBI[M11Validate],col="royalblue",
                   xlab="M11 Gene Connectivity (kME)",ylab="Anderson Malawi G.S.",pch=20)
abline(lm(GS.GSE39940.LTBI[M11Validate]~connectivity[M11Validate]))
dev.off()

#check the eigengene for M11 across the diseases
M11Eigen_GSE39939<-moduleEigengenes(GSE39939[,M11Validate],rep("royalblue",length(M11Validate)))$eigengenes[,1]

Val_Eigen_GSE39939<-data.frame(M11Eigen_GSE39939,GSE39939_pheno$disease)
colnames(Val_Eigen_GSE39939)<-c("M11 Eigengene","Disease Label")

validation_colors<-c('red','goldenrod1','darkgrey')
names(validation_colors)<-c("TB","LTBI","Other Disease")

tiff("EigengeneKenya.tiff",height=5,width=5,units="in",res=300)
ggboxplot(Val_Eigen_GSE39939,x='Disease Label',y='M11 Eigengene',
          color="Disease Label",add="jitter",outlier.shape=NA, xlab="Disease Label (Kenya)",
          palette = validation_colors,order=c("LTBI","TB","Other Disease"))+
  stat_compare_means(comparisons = list(c("LTBI","TB"),c("TB","Other Disease"),c("LTBI","Other Disease")))+theme(legend.position = "none")
dev.off()

M11Eigen_GSE39940<-moduleEigengenes(GSE39940[,M11Validate],rep("royalblue",length(M11Validate)))$eigengenes[,1]

Val_Eigen_GSE39940<-data.frame(M11Eigen_GSE39940,GSE39940_pheno$disease)
colnames(Val_Eigen_GSE39940)<-c("M11 Eigengene","Disease Label")

validation_colors<-c('red','goldenrod1','darkgrey')
names(validation_colors)<-c("TB","LTBI","Other Disease")

tiff("EigengeneMalwai.tiff",height=5,width=5,units="in",res=300)
ggboxplot(Val_Eigen_GSE39940,x='Disease Label',y='M11 Eigengene',
          color="Disease Label",add="jitter",outlier.shape=NA, xlab="Disease Label (Malawi)",
          palette = validation_colors,order=c("LTBI","TB","Other Disease"))+
  stat_compare_means(comparisons = list(c("LTBI","TB"),c("TB","Other Disease"),c("LTBI","Other Disease")))+theme(legend.position = "none")
dev.off()



