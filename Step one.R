library(TCGAbiolinks)
library(SummarizedExperiment)
library(clusterProfiler)
library(ChAMP)
library(rtracklayer)
library(stringr)
####Download clinical data
TCGAbiolinks::getGDCprojects()$project_id
cancer_type<-"TCGA-CESC"
clinical<-GDCquery_clinic(project=cancer_type,type="clinical")
######Transcriptome expression profile
query <- GDCquery(project = "TCGA-CESC", 
                  legacy = FALSE, 
                  experimental.strategy = "RNA-Seq", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")
GDCdownload(query)#Download expression spectrum data
###Directly process as expression spectrum
data<-GDCprepare(query = query)
cescFPKM<-assay(data)
#####Download methylation spectrum
query_me <- GDCquery(project = "TCGA-CESC", 
                     legacy = TRUE, 
                     data.category = "DNA methylation")
GDCdownload(query_me)#Download expression spectrum data
data<-GDCprepare(query = query_me)
cescMe<-assay(data)
#### Extracted the geneid of lncrna
gtf <- rtracklayer::import("./gencode.v36.long_noncoding_RNAs.gtf") 
gtf_df=as.data.frame(gtf)
lnc_id<-gtf_df$gene_id
###Remove duplicates geneid
lnc_id<-as.data.frame(lnc_id[!duplicated(lnc_id)])
####Remove version number
lnc_id<-apply(lnc_id,1,function(x){
  strsplit(x,"\\.")[[1]][1]
})
lnc_id<-as.data.frame(lnc_id)
lnc_id<-as.data.frame(lnc_id[!duplicated(lnc_id),])
###Distinguish between disease and normal samples, with cancer samples in front (306) and normal samples behind (3)
barcode<-substring(colnames(cescFPKM),14,15)
table(barcode)
tumor_samp1<-cescFPKM[,which(barcode=="01")]
tumor_samp2<-cescFPKM[,which(barcode=="06")]
normal_samp<-cescFPKM[,which(barcode=="11")]
cescFPKM<-cbind(tumor_samp1,tumor_samp2)
cescFPKM<-cbind(cescFPKM,normal_samp)
####Extracted the geneid of mrna
library(rtracklayer)
gtf <- rtracklayer::import("./gencode.v36.basic.annotation.gtf") 
gtf_df=as.data.frame(gtf)
mrna_data<-gtf_df[which(gtf_df$gene_type=="protein_coding"),]
mrna_id<-mrna_data$gene_id
###Remove duplicates geneid
mrna_id<-as.data.frame(mrna_id[!duplicated(mrna_id)])
mrna_id<-apply(mrna_id,1,function(x){
  strsplit(x,"\\.")[[1]][1]
})
mrna_id<-as.data.frame(mrna_id)
mrna_id<-as.data.frame(mrna_id[!duplicated(mrna_id),])
####Extract lnc and mRNA expression profiles
cesc_lnc<-cescFPKM1[lnc_id[,1],]
cesc_lnc<-na.omit(cesc_lnc)
cesc_mrna<-cescFPKM1[mrna_id[,1],]
cesc_mrna<-na.omit(cesc_mrna)
####Processing methylation expression profiles
####450k
####Remove samples that do not exist in the expression profile data
sam_mrna<-gsub("\\.","-",substring(colnames(cesc_mrna),1,15))
sam_me<-substring(colnames(cescMe),1,15)
sam<-intersect(sam_me,sam_mrna)
cescMe<-cescMe[,match(sam,sam_me)]
temp<-substring(colnames(cescMe),1,15)
cescMe1<-cescMe
for (i in 1:length(temp)) {
  index<-which(temp[i]==sam_mrna)
  cescMe1[,index]<-cescMe[,i]
}
colnames(cescMe1)<-sam_mrna
cescMe<-cescMe1
######Processing clinical data
######Corresponding the order of clinical data with the order of gene expression profile samples, where the last five samples are duplicate samples, including two 06 and three normal samples
clinical2<-as.data.frame(matrix(nrow = 309,ncol = 73))
for(i in 1:length(colnames(cescMe))){
  index<-which(substring(colnames(cescMe)[i],1,12)==clinical[,1])
  if(index>0){
    clinical2[i,]<-clinical[index,]
  }
  
}
clinical<-clinical2
clinical2[,1]<-paste0(clinical2[,1],rep(c("-01","-06","-11"),c(304,2,3)))
clinical2[,2]<-clinical[,1]
clinical2<-clinical2[,c(1,2)]
clinical2$group_class<-rep(c("Tumor","Normal"),c(306,3))
colnames(clinical2)<-c("sampleID","patient","group_class")

#####Organize as ChAMP object
beta<-as.matrix(cescMe)
beta=impute.knn(beta) 
sum(is.na(beta))
beta=beta$data
beta=beta+0.00001
myLoad=champ.filter(beta = beta,pd=clinical2) #This step has already been automatically filtered
####Quality Control
QC = champ.QC(beta = myLoad$beta,pheno = myLoad$pd$group_class)
####normalization
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=4)
###Normalization will generate a lot of NA. Let's see which samples contain NA
num.na<-apply(myNorm, 2, function(x){sum(is.na(x))})
table(num.na)
#####There are 300 samples with no NA, and 9 samples with a lot of NA generated. Remove these samples.
names(num.na)<-colnames(myNorm)
dt=names(num.na[num.na>0])
keep<-setdiff(colnames(myNorm),dt)
myNorm<-myNorm[,keep]
pd=myLoad$pd
pd<-pd[pd$sampleID %in% colnames(myNorm),]
##Export processed methylation information for subsequent analysis
save(myNorm,file = './cescmyNorm.RData')