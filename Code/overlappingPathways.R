library(ggplot2)

#get overlapping pathways

overlapping_paths<-intersect(conversiongsea_sub$pathway,TBgsea_sub$pathway)
overlapping_paths<-intersect(overlapping_paths,TBinConvgsea_sub$pathway)
overlapping_paths<-intersect(overlapping_paths,Smokegsea_sub$pathway)

overlapping_gsea_df<-data.frame(c(rep("TST Conversion",5),rep("TB Disease",5),rep("TB Disease in\nTST Converters",5),rep("Maternal Smoking",5)))
colnames(overlapping_gsea_df)<-"Hypothesis"
overlapping_gsea_df$pathway<-rep(overlapping_paths,4)

library(tidyverse)
overlapping_gsea_df <- overlapping_gsea_df %>% arrange(pathway)

magNES<-c()
neglog10pval<-c()

for(path in overlapping_paths){
  magNES<-append(magNES,conversiongsea_sub$magNES[conversiongsea_sub$pathway==path])
  magNES<-append(magNES,TBgsea_sub$magNES[TBgsea_sub$pathway==path])
  magNES<-append(magNES,TBinConv_sub$magNES[TBinConv_sub$pathway==path])
  magNES<-append(magNES,Smokegsea_sub$magNES[Smokegsea_sub$pathway==path])
  
  neglog10pval<-append(neglog10pval,(-1)*log10(conversiongsea_sub$pval[conversiongsea_sub$pathway==path]))
  neglog10pval<-append(neglog10pval,(-1)*log10(TBgsea_sub$pval[TBgsea_sub$pathway==path]))
  neglog10pval<-append(neglog10pval,(-1)*log10(TBinConv_sub$pval[TBinConv_sub$pathway==path]))
  neglog10pval<-append(neglog10pval,(-1)*log10(Smokegsea_sub$pval[Smokegsea_sub$pathway==path]))
}

overlapping_gsea_df$magNes<-magNES
overlapping_gsea_df$neglog10pval<-neglog10pval

overlapping_gsea_df$pathway<-tolower(sapply(str_split(overlapping_gsea_df$pathway,"%"),`[`,1))


set.seed(42)

library(latex2exp)

tiff("gseabubbleplot.tiff",height=4,width=8,units="in",res=600)
ggplot(overlapping_gsea_df, aes(Hypothesis, forcats::fct_rev(pathway), fill = neglog10pval, size = magNes)) +
  geom_point(shape = 21, stroke = 0) +
  geom_hline(yintercept = seq(.5, 4.5, 1), size = .2) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(4, 10)) +
  scale_fill_gradient(low = "orange", high = "blue", breaks = c(0, 10, 20), labels = c(0, 10, 20), limits = c(0, 20)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust=0)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "top", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Magnitude NES", fill = TeX("$-log_{10}$p-value"), x = NULL, y = NULL)
dev.off()

#pathway z score for defense response across TB categories, with point type by maternal smoking

DefenseGenes<-ourGmtList["DEFENSE RESPONSE%GOBP%GO:0006952"][[1]]
ImmuneResponseGenes<-ourGmtList["REGULATION OF IMMUNE RESPONSE%GOBP%GO:0050776"][[1]]

#We'll need our DE results
convResults<-read.csv(file="CoverterDEResults.csv", row.names=1)

#And our Zscore Function
SigZScore<-function(M,up,down){
  if(length(up)>0){
    upScore<-apply(M[up,],2,mean)
  }
  else(upScore=0)
  if(length(down)>0){
    downScore<-apply(M[down,],2,mean)
  }
  else(downScore=0)
  
  Score<-(upScore-downScore)
  return(Score)
}

ConvDefenseGeneResults<-na.omit(convResults[DefenseGenes,])
DefenseDown<-rownames(ConvDefenseGeneResults[which(ConvDefenseGeneResults$logFC<0),])
DefenseUp<-rownames(ConvDefenseGeneResults[which(ConvDefenseGeneResults$logFC>0),])

DefenseUp<-DefenseUp[DefenseUp%in%colnames(convExpr)]
DefenseDown<-DefenseDown[DefenseDown%in%colnames(convExpr)]


DefenseSigScore<-SigZScore(t(convExpr),DefenseUp,DefenseDown)
(centerOn<-median(DefenseSigScore[convPheno$conversion_10_at36months=="Nonconverter"&convPheno$tbdisease_merged=="No TB Diagnosis"]))
DefenseSigScore = DefenseSigScore - centerOn

library(ggpubr)

converterTB<-as.character(convPheno$tbdisease_merged)
converterTB[convPheno$tbdisease_merged=="No TB Diagnosis" & convPheno$conversion_10_at36months=="Converter"]<- "Converter: No TB Disease"
converterTB[convPheno$conversion_10_at36months!="Converter"]<-"Nonconverter"
converterTB<-as.factor(converterTB)

ZScoreDF<-data.frame(converterTB,DefenseSigScore,convPheno$smoke_current)
colnames(ZScoreDF)<-c("Status","Score","Smoker")
ZScoreDF$Smoker<-as.factor(ifelse(ZScoreDF$Smoker==1, "Current Smoker","Non-smoker"))
ZScoreDF$Smoker<-relevel(ZScoreDF$Smoker,ref = "Current Smoker")

comps<-list(c("TB Disease","Converter: No TB Disease"),c("TB Disease","Nonconverter"),c("Converter: No TB Disease","Nonconverter"))

colors = c('red','goldenrod1','blue3')
names(colors) = c("TB Disease","Converter: No TB Disease","Nonconverter")

xlabs<-c("TB Disease","Converter:\nNo TB Disease","Nonconverter")


ZScoreDF<-na.omit(ZScoreDF)

library(ggpubr)

p1 <- ggboxplot(ZScoreDF, x = "Status", y = "Score",
                color = "Status", palette = colors,
                order=c("TB Disease","Converter: No TB Disease","Nonconverter"),
                outlier.shape=NA, xlab=F,ylab="Defense Response Pathway Z-Score") + guides(fill="none",color="none") +
  geom_jitter(aes(shape = Smoker,color=Status),alpha=0.5) + 
  scale_x_discrete(labels= xlabs) +labs(shape="Maternal Smoking\nStatus")
  
p1 <- p1 +stat_compare_means(method="wilcox.test",comparisons = comps) + theme(legend.position = "right",text = element_text(size = 10),axis.text.x=element_text(color = "black",angle=30, vjust=.8, hjust=0.8))
p1

tiff("DefensePathwayZscore",height=6,width=6,units="in",res=400)
p1
dev.off()

ConvDefenseGeneResults<-na.omit(convResults[ImmuneResponseGenes,])
DefenseDown<-rownames(ConvDefenseGeneResults[which(ConvDefenseGeneResults$logFC<0),])
DefenseUp<-rownames(ConvDefenseGeneResults[which(ConvDefenseGeneResults$logFC>0),])

DefenseUp<-DefenseUp[DefenseUp%in%colnames(convExpr)]
DefenseDown<-DefenseDown[DefenseDown%in%colnames(convExpr)]


DefenseSigScore<-SigZScore(t(convExpr),DefenseUp,DefenseDown)
(centerOn<-median(DefenseSigScore[convPheno$conversion_10_at36months=="Nonconverter"&convPheno$tbdisease_merged=="No TB Diagnosis"]))
DefenseSigScore = DefenseSigScore - centerOn


ZScoreDF<-data.frame(converterTB,DefenseSigScore,convPheno$smoke_current)
colnames(ZScoreDF)<-c("Status","Score","Smoker")
ZScoreDF$Smoker<-as.factor(ifelse(ZScoreDF$Smoker==1, "Current Smoker","Non-smoker"))
ZScoreDF$Smoker<-relevel(ZScoreDF$Smoker,ref = "Current Smoker")

ZScoreDF<-na.omit(ZScoreDF)

library(ggpubr)

p1 <- ggboxplot(ZScoreDF, x = "Status", y = "Score",
                color = "Status", palette = colors,
                order=c("TB Disease","Converter: No TB Disease","Nonconverter"),
                outlier.shape=NA, xlab=F,ylab="Regulation of Immune Response Pathway Z-Score") + guides(fill="none",color="none") +
  geom_jitter(aes(shape = Smoker,color=Status),alpha=0.5) + 
  scale_x_discrete(labels= xlabs) +labs(shape="Maternal Smoking\nStatus")

p1 <- p1 +stat_compare_means(method="wilcox.test",comparisons = comps) + theme(legend.position = "right",text = element_text(size = 10),axis.text.x=element_text(color = "black",angle=30, vjust=.8, hjust=0.8))
p1

tiff("ImmuneResponsePathwayZscore.tiff",height=6,width=6,units="in",res=400)
p1
dev.off()

select_pathways<-read.csv("select_pathways.csv")
select_pathways$Pathway<-toupper(select_pathways$Pathway)

tiff("WGCNA_pathway_viz.tiff",height=8,width=12,units="in",res=400)
ggbarplot(select_pathways,x='Pathway',y='log10pval',fill = 'Pathway',palette = "npg")+ylab(TeX("$-log_{10}$p-value"))+
  geom_hline(yintercept = -log10(0.05),color="darkgrey",linetype="dashed")+
  annotate("text",x=0.65,y=-log10(0.05)+0.25,label="p-value=0.05",colour="navy",size=3, fontface =2)+
  theme(text = element_text(size = 11),axis.text.x=element_text(color = "black",angle=30, vjust=.8, hjust=0.8),legend.position="none")
dev.off()

