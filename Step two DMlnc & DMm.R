##Organize expression spectrum
load("cescmyNorm.RData")
keep<-colnames(myNorm)
cescFPKM1<-cescFPKM
###According to the methylation level, the samples were divided into high and low methylation groups to screen for differentially expressed genes. The sample is a shared sample of methylation profile and transcriptome expression profile.
###Only use cancer samples for subsequent analysis
colnames(cescFPKM1)<-gsub("\\.","-",substring(colnames(cescFPKM1),1,15))
cescFPKM1<-cescFPKM1[,keep]
###Only retain cancer samples
cescFPKM1<-cescFPKM1[,-c(298,299,300)]
myNorm<-myNorm[,-c(298,299,300)]
####Extract lnc and mRNA expression profiles
cesc_lnc<-cescFPKM1[which(rownames(cescFPKM1)%in%lnc_id[,1]),]
cesc_mrna<-cescFPKM1[which(rownames(cescFPKM1)%in%mrna_id[,1]),]

####Data preprocessing
####Remove genes with gene expression values exceeding 30% from 0。
index<-which(unlist(apply(cesc_lnc,1,function(x){
  sum(x == 0)/length(x) > 0.3
})))
cesc_lnc2<-cesc_lnc[-index,]

index<-which(unlist(apply(cesc_mrna,1,function(x){
  sum(x == 0)/length(x) > 0.3
})))
cesc_mrna2<-cesc_mrna[-index,]
###K-means
library(impute)
beta<-t(apply(cesc_lnc2,1,function(x){ifelse(x == 0,NA,x)}))
beta<-impute.knn(beta)
betaData=beta$data
betaData<-betaData+0.00001
cesc_lnc2<-betaData
cesc_lnc2<-as.data.frame(cesc_lnc2)

beta<-t(apply(cesc_mrna2,1,function(x){ifelse(x == 0,NA,x)}))
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
cesc_mrna2=betaData
cesc_mrna2<-as.data.frame(cesc_mrna2)
###Find differentially methylated lnc: DMlnc (upregulated) 
num<-ceiling(dim(cescFPKM1)[2]*0.4)
DMlnc<-data.frame(gene=NA,fc=NA,p.value=NA)
lnc_low_mean<-data.frame(gene=NA,mean=NA)
lnc_top_mean<-data.frame(gene=NA,mean=NA)
lnc_low_me<-data.frame(gene=NA,mean=NA)
lnc_top_me<-data.frame(gene=NA,mean=NA)
for (i in 1:dim(lnc_me_exp)[1]) {
  index<-order(lnc_me_exp[i,],decreasing = T)
  top<-index[1:num]
  low<-index[179:297]
  gene<-rownames(lnc_me_exp)[i]
  lnc_low_mean[i,]<-c(gene,rowMeans(cesc_lnc2[gene,low]))
  lnc_top_mean[i,]<-c(gene,rowMeans(cesc_lnc2[gene,top]))
  lnc_low_me[i,]<-c(gene,rowMeans(lnc_me_exp[gene,low]))
  lnc_top_me[i,]<-c(gene,rowMeans(lnc_me_exp[gene,top]))
  fc<-mean(as.numeric(cesc_lnc2[gene,low]))/mean(as.numeric(cesc_lnc2[gene,top]))
  cesc_lnc2[gene,]<-log2(cesc_lnc2[gene,])
  p<-t.test(as.numeric(cesc_lnc2[gene,top]),as.numeric(cesc_lnc2[gene,low]))$p.value
  DMlnc[i,]<-c(gene,fc,p)
}
DMlnc$logfc<-log2(as.numeric(DMlnc$fc))
DMlnc$fdr<-p.adjust(as.numeric(DMlnc$p.value),method = "fdr")
index<-which(DMlnc$logfc>0)
index2<-which(DMlnc$fdr<0.05)
DMlnc<-DMlnc[intersect(index,index2),]
write.table(lnc_low_mean,"lnc_low_mean.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(lnc_top_mean,"lnc_top_mean.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(lnc_low_me,"lnc_low_me.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(lnc_top_me,"lnc_top_me.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(DMlnc,"DMlnc.txt",quote = F,row.names = F,col.names = T,sep = "\t")
###Find differentially methylated mrna: DMm (upregulated) 

DMmrna<-data.frame(gene=NA,fc=NA,p.value=NA)
mrna_low_me<-data.frame(gene=NA,mean=NA)
mrna_top_me<-data.frame(gene=NA,mean=NA)
for (i in 1:dim(mrna_me_exp)[1]) {
  index<-order(mrna_me_exp[i,],decreasing = T)
  top<-index[1:num]
  low<-index[179:297]
  gene<-rownames(mrna_me_exp)[i]
  mrna_low_me[i,]<-c(gene,rowMeans(mrna_me_exp[gene,low]))
  mrna_top_me[i,]<-c(gene,rowMeans(mrna_me_exp[gene,top]))
  fc<-mean(as.numeric(cesc_mrna2[gene,low]))/mean(as.numeric(cesc_mrna2[gene,top]))
  cesc_mrna2[gene,]<-log2(cesc_mrna2[gene,])
  p<-t.test(as.numeric(cesc_mrna2[gene,top]),as.numeric(cesc_mrna2[gene,low]))$p.value
  DMmrna[i,]<-c(gene,fc,p)
}
DMmrna$logfc<-log2(as.numeric(DMmrna$fc))
DMmrna$fdr<-p.adjust(as.numeric(DMmrna$p.value),method = "fdr")
index<-which(DMmrna$logfc>0)
index2<-which(DMmrna$fdr<0.05)
DMmrna<-DMmrna[intersect(index,index2),]
write.table(DMmrna,"DMmrna.txt",quote = F,row.names = F,col.names = T,sep = "\t")