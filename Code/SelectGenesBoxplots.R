DGEdf<-data.frame(convExpr[,c("VNN2","PNKD","EIF3A","PARP8","DCTN4","PGAM1")])
DGEdf$converterTB<-MEdf$converterTB
DGEdf$tb_everdx<-MEdf$tb_everdx

p1<-ggboxplot(DGEdf,x="converterTB",y=colnames(DGEdf)[1],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA, ylab=colnames(DGEdf)[1])
p1 <- p1 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none") 
p1

p2<-ggboxplot(DGEdf,x="converterTB",y=colnames(DGEdf)[2],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA, ylab=colnames(DGEdf)[2])
p2 <- p2 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none") 
p2

p3<-ggboxplot(DGEdf,x="converterTB",y=colnames(DGEdf)[3],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA, ylab=colnames(DGEdf)[3])
p3 <- p3 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none") 
p3


p4<-ggboxplot(DGEdf,x="converterTB",y=colnames(DGEdf)[4],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA, ylab=colnames(DGEdf)[4])
p4 <- p4 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none") 
p4

p5<-ggboxplot(DGEdf,x="converterTB",y=colnames(DGEdf)[5],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA, ylab=colnames(DGEdf)[5])
p5 <- p5 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none") 
p5

p6<-ggboxplot(DGEdf,x="converterTB",y=colnames(DGEdf)[6],color="converterTB",
              palette = c('red','goldenrod1','blue3'),
              order=c("TB Disease","Converter: No TB Diagnosis","Nonconverter"),
              add = "jitter", xlab=F, outlier.shape=NA, ylab=colnames(DGEdf)[6])
p6 <- p6 +stat_compare_means(method="wilcox.test",comparisons = comps)+ theme(legend.position = "none") 
p6

