## import data
real_data<-read.table(file="Data/rna_LL_all_august_version_01.txt",header=TRUE)
real_data_filter<-real_data[!is.na(real_data$duplicatestatus),]
only_homeologous<-real_data_filter[real_data_filter$duplicatestatus=="Inter-chromosomal ohnolog",]

## get duplication data
duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")
duplication_table$typeregion<-ifelse(duplication_table$aore_count>=duplication_table$lore_count,"AORe","LORe")


## define is region is Aore or Lore
for(i in 1:nrow(only_homeologous)){
  possnp<-only_homeologous[i,"snppos"]
  chromsnp<-only_homeologous[i,"snpchrom"]
  
  posgene<-only_homeologous[i,"start"]
  chromgene<-only_homeologous[i,"chr"]
  
  typeregion<-duplication_table[(duplication_table$chromosome_x==chromsnp & duplication_table$begin_x<=possnp & duplication_table$end_x>=possnp & duplication_table$chromosome_y==chromgene & duplication_table$begin_y<=posgene & duplication_table$end_y>=posgene) | (duplication_table$chromosome_y==chromsnp & duplication_table$begin_y<=possnp & duplication_table$end_y>=possnp & duplication_table$chromosome_x==chromgene & duplication_table$begin_x<=posgene & duplication_table$end_x>=posgene),"typeregion"]
  only_homeologous[i,"typeregion"]<-typeregion
}

data_circos_homeologous<-data.frame(chromgene=only_homeologous$chr,startgene=only_homeologous$start,endgene=only_homeologous$end,chromsnp=only_homeologous$snpchrom,startsnp=only_homeologous$snppos,endsnp=only_homeologous$snppos,typeregion=only_homeologous$typeregion)

### save data
write.table(data_circos_homeologous,file="Data/data_circos_homeologous.txt",row.names = FALSE)


## test enrichment

## first get number of bp in both regions
duplication_table$nbbp<-(duplication_table$end_x-duplication_table$begin_x)+(duplication_table$end_y-duplication_table$begin_y)
aorebp<-sum(duplication_table[duplication_table$typeregion=="AORe","nbbp"])
lorebp<-sum(duplication_table[duplication_table$typeregion=="LORe","nbbp"])

## then get the data
tablecircos<-read.table(file="Data/data_circos_homeologous.txt",header=TRUE)
nbconnectionAORE<-nrow(tablecircos[tablecircos$typeregion=="AORe",])
nbconnectionLORE<-nrow(tablecircos[tablecircos$typeregion=="LORe",])


## get the probability of a connection to fall in AORe region
prob_AORE<-as.numeric(as.character(aorebp))/(as.numeric(as.character(aorebp))+as.numeric(as.character(lorebp)))



vector_AORE_random<-c()
vector_LORE_random<-c()

## for 1000 time, sample randomly according to compited probability
for(i in 1:1000){
  samples <- rbinom(nrow(tablecircos), 1, prob_AORE)
  result <- ifelse(samples == 1, "AORe", "LORe")
  vector_AORE_random[i]<-length(result[result=="AORe"])
  vector_LORE_random[i]<-length(result[result=="LORe"])
}

observed <- c(nbconnectionAORE, nbconnectionLORE)  # Observed frequencies 
expected <- c(mean(vector_AORE_random), mean(vector_LORE_random))  # Expected frequencies 

chisq.test(matrix(c(observed, expected), nrow = 2))