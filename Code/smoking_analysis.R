# Smoking Analysis

################### TB disease vs not

#start with converters vs not converters
ix<-which(is.na(pheno$smoke_current))
SmokePheno<-pheno[-ix,]
SmokeExpr<-expr[-ix,]

SmokeExpr<-apply(SmokeExpr,2,as.numeric)
SmokeOutcome<-ifelse(SmokePheno$smoke_current==1,"Current Smoker","Not Current Smoker")

Group <- factor(SmokeOutcome,levels = c("Not Current Smoker","Current Smoker"))
Group

#Add covariates
RIN <- as.numeric(SmokePheno$RIN)
Batch <- factor(SmokePheno$Batch)
HIV <- factor(SmokePheno$mat_hiv_status)


design <- model.matrix(~Group+RIN+Batch+HIV)
colnames(design)

fit <- lmFit(t(SmokeExpr), design)
set.seed(324)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")
Diff<- topTable(fit2, coef="GroupCurrent Smoker", adjust="BH", n=25000)
write.csv(Diff, file="SmokingDEResults.csv", quote=FALSE, row.names=TRUE)

#create the rank file
output<-topTable(fit2, coef="GroupCurrent Smoker", n=Inf, sort.by="P")[,keepCol]
rank<-sign(output$logFC)*(-1)*log10(output$P.Value)
#names(rank)<-row.names(output)

rank<-data.frame(row.names(output),rank)
colnames(rank)<-c("SYMBOL","RANK")
rank$SYMBOL<-as.factor(rank$SYMBOL)
rownames(rank)<-NULL

write.table(rank, "Smokerank.rnk.txt", sep = "\t", row.names = F, quote=F)
rank_smoke<-rank



set.seed(0)
Smokegsea<-fgsea(ourGmtList,rank,eps=1e-100,minSize=50,maxSize=1000)
View(Smokegsea)
Smokegsea$magNES<-abs(Smokegsea$NES)
Smokegsea_sub<-Smokegsea[Smokegsea$pval<0.005,]
Smokegsea_sub <- Smokegsea_sub %>% arrange(desc(magNES))
#collage leading edge
for(i in 1:nrow(Smokegsea_sub)){
  Smokegsea_sub$leadingEdge[i]<-paste(Smokegsea_sub$leadingEdge[[i]],collapse = ", ")
}
Smokegsea_sub$leadingEdge<-unlist(Smokegsea_sub$leadingEdge)
write.csv(Smokegsea_sub,"SmokeGSEAresults.csv",row.names=F)



