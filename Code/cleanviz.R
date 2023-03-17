DGEs<-row.names(Diff[Diff$P.Value<0.005,])
DGEMatrix<-convExpr[,DGEs]

pcs<-prcomp(DGEMatrix)

library(ggbiplot)

colors = c('red','goldenrod1','blue3')
names(colors) = c("TB Disease","Converter","Nonconverter")

g <- ggbiplot(pcs, choices=c(1,2),
              groups = as.factor(convPheno$conversion_10_at36months),
              circle = T, var.axes = F)
g <- g + geom_point(aes(col=as.factor(convPheno$conversion_10_at36months)), size=1)
g<- g + scale_colour_manual(values = colors, name="TB Category")
g <- g + theme(legend.direction = 'vertical',
               legend.position = 'right')
g <- g + theme_classic()

g

tiff("ConverterPCA.tiff",height=6,width=6,units="in",res=600)
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


p1<-p1+stat_compare_means()+theme(legend.position = "none")

p2<-ggboxplot(PCdata,x="TBDisease",y="PC1", color="TBDisease",
              palette=c("red","dodger blue"),
              order=c("TB Disease","No TB Diagnosis"),
              add="jitter",xlab=F,outlier.shape=NA)


p2<-p2+stat_compare_means()+theme(legend.position = "none")

p3<-ggboxplot(PCdata,x="MaternalPriorTB",y="PC2", color="MaternalPriorTB",
              palette=c("salmon","turquoise"),
              order=c("Prior TB","No Prior TB Diagnosis"),
              add="jitter",xlab=F,outlier.shape=NA)


p3<-p3+stat_compare_means()+theme(legend.position = "none")

tiff("conversionPCBox.tiff",height=4, width=4, units="in", res=600)
p1
dev.off()

tiff("TBPCBox.tiff",height=4, width=4, units="in", res=600)
p2
dev.off()

tiff("everTBPCBox.tiff",height=4, width=4, units="in", res=600)
p3
dev.off()


df<-data.frame(convPheno$conversion_10_at36months)#,convPheno$tbdisease_merged)
colnames(df)<-c("Conversion")#,"Tuberculosis Disease")

top_annotation = HeatmapAnnotation(name="TB Cord Blood",df=df,
                                   col=list("Conversion"=colors),
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

tiff("conversionHeatmap.tiff",height=16,width=12,units="in",res=900)
draw(HM,  heatmap_legend_side = "bottom", annotation_legend_side="bottom")
dev.off()

