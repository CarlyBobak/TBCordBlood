#fgsea files for enrichment map

head(conversiongsea_sub)

# GO.ID      {tab} Description                     {tab} p.Val {tab} FDR  {tab} Phenotype

#name, as in GMT
#description: before the first '%'
#pvalue
#fdr pvalue
#pheno (1 for up, -1 for down)

conversiongsea_out<-data.frame(conversiongsea_sub$pathway)
colnames(conversiongsea_out)<-"GO.ID"
conversiongsea_out$Description<-sapply(str_split(conversiongsea_sub$pathway,'%'),`[`,1)
conversiongsea_out$p.Val<-conversiongsea_sub$pval
conversiongsea_out$FDR<-conversiongsea_sub$padj
conversiongsea_out$Phenotype<-ifelse(conversiongsea_sub$NES>0,'+1','-1')
View(conversiongsea_out)

write_tsv(conversiongsea_out,"conversionGSEAforCytoscape.tsv")

TBgsea_out<-data.frame(TBgsea_sub$pathway)
colnames(TBgsea_out)<-"GO.ID"
TBgsea_out$Description<-sapply(str_split(TBgsea_sub$pathway,'%'),`[`,1)
TBgsea_out$p.Val<-TBgsea_sub$pval
TBgsea_out$FDR<-TBgsea_sub$padj
TBgsea_out$Phenotype<-ifelse(TBgsea_sub$NES>0,'+1','-1')

write_tsv(TBgsea_out,"TBGSEAforCytoscape.tsv")

TBinConvgsea_out<-data.frame(TBinConvgsea_sub$pathway)
colnames(TBinConvgsea_out)<-"GO.ID"
TBinConvgsea_out$Description<-sapply(str_split(TBinConvgsea_sub$pathway,'%'),`[`,1)
TBinConvgsea_out$p.Val<-TBinConvgsea_sub$pval
TBinConvgsea_out$FDR<-TBinConvgsea_sub$padj
TBinConvgsea_out$Phenotype<-ifelse(TBinConvgsea_sub$NES>0,'+1','-1')

write_tsv(TBinConvgsea_out,"TBinConvGSEAforCytoscape.tsv")

Smokegsea_out<-data.frame(Smokegsea_sub$pathway)
colnames(Smokegsea_out)<-"GO.ID"
Smokegsea_out$Description<-sapply(str_split(Smokegsea_sub$pathway,'%'),`[`,1)
Smokegsea_out$p.Val<-Smokegsea_sub$pval
Smokegsea_out$FDR<-Smokegsea_sub$padj
Smokegsea_out$Phenotype<-ifelse(Smokegsea_sub$NES>0,'+1','-1')

write_tsv(Smokegsea_out,"SmokingGSEAforCytoscape.tsv")









