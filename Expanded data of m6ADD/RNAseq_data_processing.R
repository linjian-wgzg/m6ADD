rm(list=ls())
##Read the expression profile data of rectal cancer
shuju=read.table("C:\\Users\\Administrator\\Desktop\\TCGA-READ.htseq_fpkm.tsv",sep="\t",header=T,stringsAsFactors = F)
#Remove the version number of the gene name
library(stringr)
x=as.data.frame(str_split_fixed(shuju$Ensembl_ID, "[.]", 2))
shuju[,1]=x[,1]

colnames(shuju)[1]="name"
rownames(shuju)=shuju$name
shuju=shuju[,-1]

######Delete all 0 rows
aa=apply(shuju,1,sum)
loc=which(aa==0)
shuju=shuju[-loc,]
#######Delete the lines that are greater than 70% as 0, and round to the nearest integer.
loc1=NULL
for(i in 1:dim(shuju)[1]){
  loc=which(shuju[i,]==0)
  if(length(loc)>round(dim(shuju)[2]*0.70)){
    loc1=c(loc1,i)
  }
}
shuju=shuju[-loc1,]
#####Replace the 0 value of the remaining sites with the minimum value
shuju1=as.vector(as.matrix(shuju))
num=sort(unique(shuju1))
for(i in 1:dim(shuju)[1]){
  loc=which(shuju[i,]==0)
  shuju[i,loc]=num[2]
}
write.csv(shuju,"C:\\Users\\Administrator\\Desktop\\Pretreatment.csv")




####Set the control group and case groupSet the control group and case group
t=grep(".11",colnames(shuju))
data_control=as.data.frame(shuju[,t])

t1=grep(".01",colnames(shuju))
data_case=shuju[,t1]

###t2=grep(".06",colnames(shuju))
###data_case=shuju[,c(t1,t2)]
###Manually delete columns that do not meet the criteria
colnames(data_control)
data_control=data_control[,-c(2,3,4,10,13,16)]
colnames(data_case)
#data_case=data_case[,-c(13)]


mean_case=apply(data_case,1,mean)
mean_control=apply(data_control,1,mean)
fold_change=mean_case/mean_control
loc=which(fold_change>=2|fold_change<=1/2)
data=cbind(data_case,data_control)
data$fc=fold_change
rownames(data)=rownames(shuju)
#############t-test fdr
p=0
sig_fold_change=which(fold_change>=2|fold_change<=1/2)#
for(i in 1:dim(data)[1]){
  p_value=t.test(data_case[i,],data_control[i,])$p.value
  p=c(p,p_value)
}
p=p[-1]

head(p)
data$p=p
data$padj=p.adjust(data$p,method="fdr",length(data$p))
sig_padj=which(data$padj<0.05)#
sig_m=intersect(sig_fold_change,sig_padj)###
fc.padj=data.frame(id=rownames(data),fc=data$fc,p=data$p,padj=data$padj)
sig.data=data[sig_m,]#
sig.shuju=rownames(data)[sig_m]
###Delete the last three columns
sig.data=sig.data[,-c(177,178,179)]


write.table(sig.shuju,file="C:\\Users\\Administrator\\Desktop\\sig_shuju.txt",sep="\t",quote=F,row.names = F)#####
write.csv(sig.data,file="C:\\Users\\Administrator\\Desktop\\sig_data.csv",sep=",")





#########Screen the expression data of differentially expressed m6A regulatory proteins and differentially expressed mRNA
m6A_id=read.csv("C:\\Users\\Administrator\\Desktop\\m6A_id2.csv",header=T)
colnames(m6A_id)=c("X")
sig.data=read.csv("C:\\Users\\Administrator\\Desktop\\sig_data.csv")
library(dplyr)
sig_m6A=semi_join(sig.data,m6A_id,by="X")
new_data=anti_join(sig.data,m6A_id,by="X")
write.csv(sig_m6A,"C:\\Users\\Administrator\\Desktop\\sigm6A.csv")
write.csv(new_data,"C:\\Users\\Administrator\\Desktop\\sig.csv")
