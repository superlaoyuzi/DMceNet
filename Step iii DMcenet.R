###Calculate the correlation between DMmRNA and DMlnc, and extract significant positive correlations
mrna_hash<-hash()
for (i in 1:dim(cesc_mrna2)[1]) {
  mrna_hash[[rownames(cesc_mrna2)[i]]]<-cesc_mrna2[i,]
}
lnc_hash<-hash()
for (i in 1:dim(cesc_lnc2)[1]) {
  lnc_hash[[rownames(cesc_lnc2)[i]]]<-cesc_lnc2[i,]
}
cesc_cerna<-hash()
for(i in 1:dim(DMmrna)[1]){
  for (j in 1:dim(DMlnc)[1]) {
    temp<-cor.test(as.numeric(mrna_hash[[DMmrna[i,1]]]),as.numeric(lnc_hash[[DMlnc[j,1]]]))
    index<-dim(DMlnc)[1]*(i-1)+j
    if(temp$estimate>0.3 & temp$p.value<0.05)
      cesc_cerna[[paste(DMmrna[i,1],DMlnc[j,1])]]<-paste(temp$estimate,temp$p.value)
  }
  print(i)
}
cesc_cerna2<-data.frame(mrna=NA,lnc=NA,cor=NA,p.value=NA)
keys<-keys(cesc_cerna)
values<-values(cesc_cerna)
for (i in 1:length(keys)) {
  cesc_cerna2[i,]<-c(strsplit(keys[i]," ")[[1]][1],strsplit(keys[i]," ")[[1]][2],strsplit(values[i]," ")[[1]][1],strsplit(values[i]," ")[[1]][2])
  print(i)
}
length(table(cesc_cerna2[,1]))
length(table(cesc_cerna2[,2]))
lnc_exp<-data.frame(lnc=rownames(cesc_lnc2),exp=rowMeans(cesc_lnc2))

###Extract the average methylation levels of DMlnc and DMmra from the ceRNA network
rownames(lnc_top_me)<-lnc_top_me[,1]
rownames(lnc_low_me)<-lnc_low_me[,1]
rownames(mrna_top_me)<-mrna_top_me[,1]
rownames(mrna_low_me)<-mrna_low_me[,1]
DMlnc_top_me<-lnc_top_me[DMlnc[,1],]
DMlnc_low_me<-lnc_low_me[DMlnc[,1],]
DMmrna_top_me<-mrna_top_me[DMmrna[,1],]
DMmrna_low_me<-mrna_low_me[DMmrna[,1],]

##Read candidate ceRNAs collected from the database
load("candidate_cerna2.RData")
mikeys2<-keys(candidate_cerna2)
for (i in mikeys2) {
  if(candidate_cerna2[[i]]<1){
    del(i,candidate_cerna2)
  }
}
mikeys2<-keys(candidate_cerna2)
##Methylation driven ceRNA recognition
cesc<-apply(cesc_cerna2, 1, function(x){
  paste(x[1],x[2])
})
cesc_cerna2<-do.call(rbind,strsplit(intersect(mikeys2,cesc)," "))
length(table(cesc_cerna2[,1]))
length(table(cesc_cerna2[,2]))

