expr=read.table("C:\\Users\\Administrator\\Desktop\\TCGA-READ.htseq_fpkm.tsv",header=T)
library(stringr)
x=as.data.frame(str_split_fixed(expr$Ensembl_ID, "[.]", 2))
expr[,1]=x[,1]
rownames(expr)=expr[,1]
expr=expr[,-1]
##Screen cancer samples
t1=grep(".01",colnames(expr))
expr=expr[,t1]
colnames(expr)
##Convert gene name from ENSEMBL ID to SYMBOL ID
library(org.Hs.eg.db)
detach("package:dplyr")
writers <- rownames(expr)

cols <- c("SYMBOL")

expr_id=select(org.Hs.eg.db, keys=writers, columns=cols, keytype="ENSEMBL")
expr_id =na.omit(expr_id)
expr_id=expr_id[!duplicated(expr_id$ENSEMBL), ] 
expr_id=expr_id[!duplicated(expr_id$SYMBOL), ] 

library(dplyr)
expr[,167]=rownames(expr)
colnames(expr)[167]="ENSEMBL"
expr_symbol=semi_join(expr,expr_id,by="ENSEMBL")
rownames(expr_symbol)=expr_id[,2]
expr_symbol=expr_symbol[,-167]
write.csv(expr_symbol,"C:\\Users\\Administrator\\Desktop\\Full_site_expression(symbol).csv")

##Convert the ensembl id of the difference site to symbol
sig=read.csv("C:\\Users\\Administrator\\Desktop\\sig.csv",header=T)
rownames(sig)=sig[,2]
sig=sig[,-c(1,2)]
sig=sig[,-c(167:176)]

ENSG=as.data.frame(rownames(sig))
colnames(ENSG)=c("id")
sig[,167]=ENSG[,1]
colnames(sig)[167]="id"
detach("package:dplyr")
library(org.Hs.eg.db)

writers <- c(as.character(ENSG$id))

cols <- c("SYMBOL")

symbol=select(org.Hs.eg.db, keys=writers, columns=cols, keytype="ENSEMBL")

library(dplyr)
symbol=na.omit(symbol)
colnames(symbol)=c("id","id2")

symbol=symbol[!duplicated(symbol$id2), ] 
symbol=symbol[!duplicated(symbol$id), ] 

sig_symbol=semi_join(sig,symbol,by="id")

rownames(sig_symbol)=symbol[,2]
sig_symbol=sig_symbol[,-167]

write.csv(sig_symbol,"C:\\Users\\Administrator\\Desktop\\sig_symbol.csv")

##PPI network (ppi network comes from the protein interaction network integrated by STRING and HPRD database)
ppi=read.csv("C:\\Users\\Administrator\\Desktop\\ppi.csv",header=T)
m6A_id=read.csv("C:\\Users\\Administrator\\Desktop\\m6A_id.csv",header=T)


signame=as.data.frame(rownames(sig_symbol))
m6Aname=as.data.frame(rownames(m6A_id))
colnames(signame)="name"
colnames(m6Aname)="name"
signame=rbind(signame,m6Aname)
signame=as.data.frame(signame[!duplicated(signame$name), ])
signame[,2]=signame[,1]
colnames(signame)=c("name1","name2")
##Screen the sites that are one-step neighbors of the differential gene in the PPI network, and construct a one-step neighbor sub-network
onestep1=semi_join(ppi,signame,by="name1")
onestep2=semi_join(ppi,signame,by="name2")
Onestep=rbind(onestep1,onestep2)
Onestep=base::unique(Onestep)

onestep1=as.data.frame(Onestep[,1])
onestep2=as.data.frame(Onestep[,2])
colnames(onestep1)="name"
colnames(onestep2)="name"
onestep_node=rbind(onestep1,onestep2)
onestep_node=base::unique(onestep_node)


onestep_expr=expr_symbol
onestep_expr[,167]=rownames(onestep_expr)

colnames(onestep_expr)[167]="name"
expr_symbol=semi_join(onestep_expr,onestep_node,by="name")


colnames(m6A_id)[1]="name"
m6A_id=semi_join(expr_symbol,m6A_id,by="name")
expr_symbol=anti_join(expr_symbol,m6A_id,by="name")

rownames(expr_symbol)=expr_symbol[,167]
sig=expr_symbol[,-167]

rownames(m6A_id)=m6A_id[,167]
m6A_id=m6A_id[,-167]

######Delete the rows where sig is all 0
aa=apply(sig,1,sum)
loc=which(aa==0)
sig=sig[-loc,]
#######Delete the lines that are greater than 70% as 0, and round to the nearest integer.
loc1=NULL
for(i in 1:dim(sig)[1]){
  loc=which(sig[i,]==0)
  if(length(loc)>round(dim(sig)[2]*0.70)){
    loc1=c(loc1,i)
  }
}
sig=sig[-loc1,]

#####Replace the 0 value of the remaining sites with the minimum value
sig1=as.vector(as.matrix(sig))
num=sort(unique(sig1))
for(i in 1:dim(sig)[1]){
  loc=which(sig[i,]==0)
  sig[i,loc]=num[2]
}
##Calculate the correlation between m6A regulatory protein and differentially expressed genes


m6A_id=as.data.frame(t(m6A_id))
sig=as.data.frame(t(sig))


test_m6A=sample_n(m6A_id,166)
test_sig=sample_n(sig,166)
random_cor=cor(test_sig,test_m6A,method = c("pearson"))
random_corr=as.data.frame(random_cor)
random_cor_a=random_corr[c(1:8342),c(1:6)]

random_cor_a=as.matrix(random_cor_a)
random_cor_list=data.frame(name1=rep(colnames(random_cor_a),each=dim(random_cor_a)[1]),name2=rep(rownames(random_cor_a),dim(random_cor_a)[2]),value=matrix(random_cor_a)) 
random_cor_list[,3]=abs(random_cor_list[,3])
random_cornum=mean(random_cor_list[,3])

cor=cor(sig,m6A_id,method = c("pearson"))
###cor=corr.test(sig,m6A,method = c("pearson"))
###corr=cor$r
###corr=as.data.frame(corr)
corr=as.data.frame(cor)
cor_a=corr[c(1:8342),c(1:6)]

cor_a=as.matrix(cor_a)
cor_list=data.frame(name1=rep(colnames(cor_a),each=dim(cor_a)[1]),name2=rep(rownames(cor_a),dim(cor_a)[2]),value=matrix(cor_a)) 
cor_list[,3]=abs(cor_list[,3])

cor_list3=cor_list

cor_list=subset(cor_list, value > 0.4) 



##Calculate the distance of the protein in the network
edge=Onestep


library(igraph)

node=as.matrix(edge)
dim(node)=c(300669*2,1)
node=as.data.frame(node)
node=base::unique(node)

node=as.character(node[,1])


sigid=colnames(sig)
num1=which(sigid %in% node==TRUE)
sigid=sigid[num1]

m6Aid=colnames(m6A_id)
num2=which(m6Aid %in% node==TRUE)
m6Aid=m6Aid[num2]


m6Aid2=as.data.frame(m6Aid)
colnames(m6Aid2)="name1"
sigid2=as.data.frame(sigid)
colnames(sigid2)="name2"

corlist1=semi_join(cor_list,m6Aid2,by="name1")
corlist1=semi_join(corlist1,sigid2,bu="name2")


##Replace the protein name with a number
corlist2=corlist1[,c(1:3)]
name=node
for(i in 1:nrow(corlist2)){
  corlist2[i,4]=which(name==as.character(corlist1[i,1]))
  corlist2[i,5]=which(name==as.character(corlist1[i,2]))
}
corlist2[,1]=corlist2[,4]
corlist2[,2]=corlist2[,5]
corlist2=corlist2[,-c(4,5)]



m6Anum=as.data.frame(corlist2[,1])
m6Anum=base::unique(m6Anum)
m6Anum=as.character(m6Anum[,1])

signum=as.data.frame(corlist2[,2])
signum=base::unique(signum)
signum=as.character(signum[,1])

###Perturb the network a thousand times to calculate the weight and take the average as the threshold
test1=vector()
test2=vector()
test3=vector()
test4=vector()
test5=vector()
test6=vector()

for(i in 1:1000){
  set.seed(i)
  tubiao=erdos.renyi.game(nrow(onestep_node),nrow(edge), type="gnm")
  x=diameter(tubiao)
  ceshi <- shortest.paths(tubiao,m6Anum,signum)
  ceshi=as.data.frame(ceshi)
  rownames(ceshi)=m6Anum
  colnames(ceshi)=signum
  ceshi=ceshi[c(1:6),c(1:2759)]
  ceshi=as.matrix(ceshi)
  juli=data.frame(name1=rep(rownames(ceshi),each=dim(ceshi)[2]),name2=rep(colnames(ceshi),dim(ceshi)[1]),value=matrix(ceshi)) 
  juli[,1]=as.integer(as.character(juli[,1]))
  juli[,2]=as.integer(as.character(juli[,2]))
  juli=semi_join(juli,corlist2,by=c("name1","name2"))
  corlist2[,4]=juli[,3]
  corlist2[,5]=1-(1-corlist2[,3])*((corlist2[,4]-1)/(x-1))
  test1[i]=min(corlist2[c(1:196),5])
  test2[i]=min(corlist2[c(197:1755),5])
  test3[i]=min(corlist2[c(1756:2899),5])
  test4[i]=min(corlist2[c(2900:4362),5])
 test5[i]=min(corlist2[c(4363:6323),5])
 test6[i]=min(corlist2[c(6324:6366),5])

}


mean1=mean(test1)
mean2=mean(test2)
mean3=mean(test3)
mean4=mean(test4)
mean5=mean(test5)
mean6=mean(test6)

#Calculate the original network weight
pic=graph_from_data_frame(edge, directed = F, vertices = node)
n=diameter(pic)###the Longest shortest distance
ceshi <- shortest.paths(pic,m6Aid,sigid)
ceshi=as.data.frame(ceshi)
ceshi=ceshi[c(1:6),c(1:8342)]
ceshi=as.matrix(ceshi)
juli=data.frame(name1=rep(rownames(ceshi),each=dim(ceshi)[2]),name2=rep(colnames(ceshi),dim(ceshi)[1]),value=matrix(ceshi)) 


cor_list3[,4]=juli[,3]
cor_list3[,5]=1-(1-cor_list3[,3])*((cor_list3[,4]-1)/(n-1))

m6A1=corlist2[c(1:196),c(1,2,5)]
m6A2=corlist2[c(197:1755),c(1,2,5)]
m6A3=corlist2[c(1756:2899),c(1,2,5)]
m6A4=corlist2[c(2900:4362),c(1,2,5)]
m6A5=corlist2[c(4363:6323),c(1,2,5)]
m6A6=corlist2[c(6324:6366),c(1,2,5)]



m6A1=subset(m6A1, V5 > mean1) 
m6A2=subset(m6A2, V5 > mean2) 
m6A3=subset(m6A3, V5 > mean3) 
m6A4=subset(m6A4, V5 > mean4) 
m6A5=subset(m6A5, V5 > mean5) 
m6A6=subset(m6A6, V5 > mean6) 



m6A_total=rbind(m6A1,m6A2,m6A3,m6A4,m6A5,m6A6)


m6A_total2=m6A_total[,c(1:3)]
name=node
for(i in 1:nrow(m6A_total2)){
  m6A_total2[i,4]=name[as.numeric(m6A_total2[i,1])]
  m6A_total2[i,5]=name[as.numeric(m6A_total2[i,2])]
}
m6A_total2[,1]=m6A_total2[,4]
m6A_total2[,2]=m6A_total2[,5]
m6A_total2=m6A_total2[,-c(4,5)]



write.csv(m6A_total2,"C:\\Users\\Administrator\\Desktop\\m6A_net.csv")
write.table(m6Aid,"C:\\Users\\Administrator\\Desktop\\m6A protein mapped to PPI.txt")

###m6A_net+One-step neighbor network is a network mined by modules in cytoscape, which is used for the subsequent module mining process