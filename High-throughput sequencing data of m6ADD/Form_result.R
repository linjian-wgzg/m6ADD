anno<-read.table("\\path of ensg_symbol file\\\ensg_symbol.txt",header = F,sep="\t")
MTC <-read.table("\\path of diff_site file\\MetDiff_diff_site.xls",header = T,sep="\t")
for(i in 1:length(MTC$name)){
  MTC$Symbol_id[i] = as.character(anno$V2[which(as.character(anno$V1)==as.character(MTC$name[i]))])
}
write.table(MTC,"\\path of diff_site file\\MetDiff_diff_site.xls",sep ="\t", col.names = T,row.names = F,quote = F)

RTC <-read.table("\\path of diff_site file\\RADAR_diff_site.xls",header = T,sep="\t")
for(i in 1:length(RTC$name)){
  RTC$Symbol_id[i] = as.character(anno$V2[which(as.character(anno$V1)==as.character(RTC$name[i]))])
}

write.table(RTC,"\\path of diff_site file\\RADAR_diff_site.xls",sep ="\t", col.names = T,row.names = F,quote = F)

met <-read.table("\\path of diff_site file\\MetDiff_diff_site.xls",header = T,sep="\t")
Rd <-read.table("\\path of diff_site file\\RADAR_diff_site.xls",header = T,sep="\t")

Jj <- intersect(x=as.character(met$Symbol_id), y = as.character(Rd$name))
for(i in 1:length(met$Symbol_id)){
  if(all(is.element(as.character(met$Symbol_id[i]), Jj))){
    met$share[i] =="MetDiff & RADAR"
  }else { met$share[i]="MetDiff"}
  
}
for(i in 1:length(Rd$Symbol_id)){
  if(all(is.element(as.character(Rd$Symbol_id[i]), Jj))){
    Rd$share[i] = "MetDiff & RADAR"
  }else { Rd$share[i]="RADAR"}
  
}

##MeTDiff Adjustment
met<-met[,-which(names(met)%in%c("score","blockCount","blockSizes","blockStarts","lg.p","lg.fdr","fold_enrchment","diff.lg.fdr","thickStart","thickEnd","itemRgb","name"))]
names(met)<-c("chr","start","end","strand","p_value","logFC","Symbol_id","share")
met <- met[,c(7,1,4,2,3,5,6,8)]
met$p_value<-pow(2,met$p_value)
GSEm <-rep("GSE119168",length(met$Symbol_id))
met$GSE <-GSEm
Tm <-rep("NA",length(met$Symbol_id))
met$Treat_method <-Tm
Dm <-rep("omental tumor",length(met$Symbol_id))
met$Disease <-Dm
met <- met[,c(11,1,2,3,4,5,6,7,9,10,8)]
names(met)<-c("Disease","Symbol_id","chr","strand","start","end","p_value","logFC","GSE","Treat_method","share")

##RADAR  Adjustment
Rd<-Rd[,-which(names(Rd)%in%c("score","blockCount","blockSizes","blockStarts","thickStart","thickEnd","itemRgb","name"))]
Rd <- Rd[,c(7,1,4,2,3,6,5,8)]
GSER <-rep("GSE119168",length(Rd$Symbol_id))
Rd$GSE <-GSER
TR <-rep("NA",length(Rd$Symbol_id))
Rd$Treat_method <-TR
DR <-rep("GSE119168",length(Rd$Symbol_id))
Rd$Disease <-DR
Rd <- Rd[,c(11,1,2,3,4,5,6,7,9,10,8)]
ALL <-rbind(met,Rd)
write.table(met,"\\path of diff_site file\\GSE119168_MetDiffdiffsite.xls",sep ="\t", col.names = T,row.names = F,quote = F)
write.table(Rd,"\\path of diff_site file\\GSE119168_RADARdiff.xls",sep ="\t", col.names = T,row.names = F,quote = F)

##SCORE
##MeTDiff
met<-met[which(as.character(met$share)=="MetDiff & RADAR"),]
met<-met[,-which(names(InsertM)%in%c("share"))]
M1<-length(which(met$logFC>=0))
M2<-length(which(met$logFC<0))
Mean_m1 <-sum(met[which(met$logFC>=0),8])/M1
Mean_m2 <-sum(met[which(met$logFC<0),8])/M2
##RADAR
Rd<-Rd[which(as.character(Rd$share)=="MetDiff & RADAR"),]
Rd<-InsertR[,-which(names(InsertR)%in%c("share"))]
R1<-length(which(Rd$logFC>=0))
R2<-length(which(Rd$logFC<0))
Mean_r1 <-sum(Rd[which(Rdt$logFC>=0),8])/R1
Mean_r2 <-sum(Rd[which(Rd$logFC<0),8])/R2
ALL_mean1<-(sum(Rd[which(met$logFC>=0),8])+sum(met[which(met$logFC>=0),8]))/(R1+M1)
ALL_mean2<-(sum(Rd[which(met$logFC<0),8])+sum(met[which(met$logFC<0),8]))/(R2+M2)
MeT1<-met[which(met$logFC>=0),]
MeT2<-met[which(met$logFC<0),]
for(i in 1:length(which(met$logFC>=0)){
  score_m1[i]<-(MeT1[i,8]-Mean_m1)/ALL_mean1
}
MeT1$score<-score_m1
for(i in 1:length(which(met$logFC<0)){
  score_m2[i]<-(MeT2[i,8]-Mean_m2)/ALL_mean2
}
MeT2$score<-score_m2
met <-rbind(MeT1,MeT2)
Rd1<-met[which(Rd$logFC>=0),]
Rd2<-met[which(Rd$logFC<0),]
for(i in 1:length(which(Rd$logFC>=0)){
  score_r1[i]<-(Rd1[i,8]-Mean_r1)/ALL_mean1
}
Rd1$score<-score_r1
for(i in 1:length(which(Rd$logFC<0)){
  score_r2[i]<-(Rd2[i,8]-Mean_r2)/ALL_mean2
}
Rd2$score<-score_r2
Rd <-rbind(Rd1,Rd2)
Score<-rbind(Rd,met)
write.table(Score,"\\path of diff_site file\\GSE119168_Score.xls",sep ="\t", col.names = T,row.names = F,quote = F)