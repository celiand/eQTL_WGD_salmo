#### Import data
eQTL_data<-read.table(file="Data/rna_LL_all_august_version_01.txt",header=TRUE)
eQTL_data_nona<-eQTL_data[!is.na(eQTL_data$duplicatestatus),]

##filter only Intra-chromosomal non syntenic (ICNS) connections
ICNS_table<-eQTL_data_nona[eQTL_data_nona$duplicatestatus=="Intra-chromosomal non syntenic" & eQTL_data_nona$kind=="trans" ,c(1,2,7,8,9,11,12,16)]
### compute raw distance between eQTL and gene
ICNS_table$rawdist<-abs(ICNS_table$snppos-ICNS_table$start)



### import information about duplicated regions
duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")

###format table for easier manipulation
duplication_tableXinfo<-duplication_table[,2:4]
colnames(duplication_tableXinfo)<-c("chrom","start","end")

duplication_tableYinfo<-duplication_table[,5:7]
colnames(duplication_tableYinfo)<-c("chrom","start","end")

formattedtable<-rbind.data.frame(duplication_tableXinfo,duplication_tableYinfo,stringsAsFactors = FALSE)


### add region information

for(i in 1:nrow(ICNS_table)){
  genechr<-ICNS_table[i,"chr"]
  genepos<-ICNS_table[i,"start"]
  
  generegioninfo<-formattedtable[formattedtable$chrom==genechr & formattedtable$start<=genepos & formattedtable$end>=genepos,c("start","end")]
  
  
  eQTLchr<-ICNS_table[i,"snpchrom"]
  snppos<-ICNS_table[i,"snppos"]
  
  snpregioninfo<-formattedtable[formattedtable$chrom==eQTLchr & formattedtable$start<=snppos & formattedtable$end>=snppos,c("start","end")]
  
  
  if(nrow(generegioninfo)!=0 & nrow(snpregioninfo)!=0){
    ### the maximum distance possible is then the maximum of the 4 values - the minimum of the 4 values (since we are in the same chromosome)
    allvalues<-c(generegioninfo[1,1],generegioninfo[1,2],snpregioninfo[1,1],snpregioninfo[1,2])
    maxdistance<-max(allvalues)-min(allvalues)
  }else{
    maxdistance<-NA
  }
  
  ICNS_table[i,"maxidistance"]<-maxdistance
}

ICNS_table<-ICNS_table[!is.na(ICNS_table$maxidistance),]

summary(ICNS_table$maxidistance)


### we can do the same for Intra-chromosomal syntenic (ICS) region


ICS_table<-eQTL_data_nona[eQTL_data_nona$duplicatestatus=="Intra-chromosomal syntenic" & eQTL_data_nona$kind=="trans",c(1,2,7,8,9,11,12,16)]
ICS_table$rawdist<-abs(ICS_table$snppos-ICS_table$start)


for(i in 1:nrow(ICS_table)){
  genechr<-ICS_table[i,"chr"]
  genepos<-ICS_table[i,"start"]
  
  generegioninfo<-formattedtable[formattedtable$chrom==genechr & formattedtable$start<=genepos & formattedtable$end>=genepos,c("start","end") ]
  
  
  eQTLchr<-ICS_table[i,"snpchrom"]
  snppos<-ICS_table[i,"snppos"]
  
  snpregioninfo<-formattedtable[formattedtable$chrom==eQTLchr & formattedtable$start<=snppos & formattedtable$end>=snppos,c("start","end")]
  
  
  if(nrow(generegioninfo)!=0 & nrow(snpregioninfo)!=0){
    ### the maximum distance possible is then the maximum of the 4 values - the minimum of the 4 values (since we are in the same chromosome)
    allvalues<-c(generegioninfo[1,1],generegioninfo[1,2],snpregioninfo[1,1],snpregioninfo[1,2])
    maxdistance<-max(allvalues)-min(allvalues)
  }else{
    maxdistance<-NA
  }
  
  ICS_table[i,"maxidistance"]<-maxdistance
}

ICS_table<-ICS_table[!is.na(ICS_table$maxidistance),]

summary(ICS_table$maxidistance)


### we can scale the distance

ICNS_table$scaled_dist<-ICNS_table$rawdist/ICNS_table$maxidistance
ICS_table$scaled_dist<-ICS_table$rawdist/ICS_table$maxidistance

mergedtable<-rbind.data.frame(ICNS_table,ICS_table,stringsAsFactors = FALSE)

library(ggplot2)

ggplot(ICNS_table,aes(x=scaled_dist,fill=duplicatestatus))+geom_histogram(color="black")+xlab(c(0,1))+scale_fill_manual(values=c("#B5AFDE"))

ggplot(ICS_table,aes(x=scaled_dist,fill=duplicatestatus))+geom_histogram(color="black")+xlab(c(0,1))+scale_fill_manual(values=c("#B6D8B6"))


ggplot(mergedtable, aes(scaled_dist, color=duplicatestatus))+stat_ecdf(linewidth=2)+scale_color_manual(values=c("#B5AFDE","#B6D8B6"))


ggplot(mergedtable, aes(x=scaled_dist, color=duplicatestatus))+geom_density(linewidth=2)+scale_color_manual(values=c("#B5AFDE","#B6D8B6"))+theme_bw(15)

