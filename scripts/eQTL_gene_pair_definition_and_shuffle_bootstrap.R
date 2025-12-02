## Script use to define duplicate region of Figure1 and shuffling code
newdata<-read.csv(file="Data/cistrans_gene_August.csv",header=TRUE)

#### -- Define regions -- ####

## change chromosome name
newdata$snpchrom<-as.numeric(as.character(ifelse(substr(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 1),4,4)==0,substr(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 1),5,5),substr(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 1),4,5))))
newdata$snppos<-as.numeric(as.character(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 2)))

## import duplicated region information
duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")


## convert position to number

### these are the length of the 29 chromosomes in Atlantic salmon
lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector)[1:28])

## convert position for duplication table
chromname<-c(1:29)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)

duplication_table$startxnb<-duplication_table$begin_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]

duplication_table$endxnb<-duplication_table$end_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]

duplication_table$startynb<-duplication_table$begin_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]

duplication_table$endynb<-duplication_table$end_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]

##and for the data

newdata$snpposnb<-newdata$snppos+listval[match(newdata$snpchrom,listval$chromname),3]

newdata$genestartnb<-newdata$start+listval[match(newdata$chr,listval$chromname),3]

newdata$geneendnb<-newdata$end+listval[match(newdata$chr,listval$chromname),3]


####end of transformation


newdata2<-newdata
lengthdata<-length(newdata2)

### for each pair of Regions (each duplicated region), we run the following test
for(i in 1:length(duplication_table$id)){
  duplication_tablerow1<-duplication_table[i,]
  ## this test is defining duplicated region (see decision_tree and explications below for more information)
  test<-ifelse(newdata$snpposnb>=duplication_tablerow1$startxnb & newdata$snpposnb<=duplication_tablerow1$endxnb,ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(newdata$snpposnb>=duplication_tablerow1$startynb & newdata$snpposnb<=duplication_tablerow1$endynb,ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Gene present",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Gene present",ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Gene present",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Gene present","NA"))))))
  newdata2[,lengthdata+i]<-test
}

## notice that this part does not define Inter-chromosomal non ohnolog, which are defined below when gene and eQTL are present but does not have any of the other labels

## here is a decomposition of the test, in order of their apperance in the test above

## newdata$snpposnb>=duplication_tablerow1$startxnb & newdata$snpposnb<=duplication_tablerow1$endxnb // this part check if the eQTL is located in the Region 1 of the current homeologous pair 
## newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog" // This part check if the start of the gene is located in the region 2 of the current homeologous pair
## newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog" // Same than before but with the end of the gene
## newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1] // This test if the gene is on the same chromosome as the eQTL
## newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic" // This check if the start of the gene is in Region 1
## newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic" // this does the same than before but with the end of the gene


## newdata$snpposnb>=duplication_tablerow1$startynb & newdata$snpposnb<=duplication_tablerow1$endynb // this part check if the eQTL is located in the Region 2 of the current homeologous pair
## newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog" // This part check if the start of the gene is in region 1
## newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog" // same as before but with the end of the gene
## newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1] // this part check if the gene is on the same chromosome as the eQTL
## newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic" // this part check if the start of the gene is in Region 2
## newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic" // same as before but with the end of the gene


## newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Gene present" // This part just check is the start of the gene is present in Region 2
## newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Gene present" // same as before but for the end of the gene
## newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Gene present" // same as before but in Region 1
## newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb // same as before but with the start of the gene



            

## then finish process to get one value output (indeed, the part above is just a part of the process)

list_of_values <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic")

# Define a function to apply row-wise
process_row <- function(row) {
  if (any(row %in% list_of_values)) {
    return(row[row %in% list_of_values][1])  # If any values in the list is found, get the first one
  } else if ("Gene present" %in% row && "eQTL present" %in% row) {  ## else if the eQTl and the gene are found, define as Inter-chromosomal non ohnolog
    return("Inter-chromosomal non ohnolog")
  } else { 
    return(NA)
  }
}

# Apply the function to each eQTL - gene pair
newdata2$result <- apply(newdata2, 1, process_row)

## extract the defined status for the eQTL-gene pairs
newdata$duplicatestatus<-newdata2[,103]

## save the results
write.table(newdata,file="Data/rna_LL_all_august_version_01.txt",row.names = FALSE)


#### -- Shuffle code -- ####

shuffle_data<- function(table){
  cols_to_shuffle <- c("chr","start","end")
  cols_to_shuffle2<-c("snpchrom","snppos")
  
  # combine the columns to shuffle into one matrix
  shuffle_matrix <- table[, cols_to_shuffle]
  shuffle_matrix2 <- table[, cols_to_shuffle2]
  
  
  # shuffle the rows of the matrix
  shuffle_matrix <- shuffle_matrix[sample(nrow(shuffle_matrix)),]
  shuffle_matrix2 <- shuffle_matrix2[sample(nrow(shuffle_matrix2)),]
  
  # split the shuffled matrix back into the original columns
  newtab<- cbind.data.frame(shuffle_matrix, shuffle_matrix2)
  return(newtab)
  
}

## shuffle the data (one time)
shuffleddata<-shuffle_data(newdata)


## Bootstrap function (same process as before, but repeated after shuffling)

Bootstrap_trans<-function(nb,table_init){
  
  library(dplyr)
  library(tidyr)
  
  duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")
  
  
  
  lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                  114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                  55994222,45305548,41468476,43051128)
  chromadd<-c(0,cumsum(lengthvector)[1:28])
  
  chromname<-c(1:29)
  chrommid<-chromadd+lengthvector/2
  listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
  
  # Define a function to apply row-wise
  process_row <- function(row) {
    if (any(row %in% list_of_values)) {
      return(row[row %in% list_of_values][1])  # Output the first found value
    } else if ("Gene present" %in% row && "eQTL present" %in% row) {
      return("Inter-chromosomal non ohnolog")
    } else {
      return(NA)
    }
  }
  
  dupstatus_list <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic","Inter-chromosomal non ohnolog")
  df_bootstrap<-c()
  for (y in 1:nb){
    
    shuffleddata<-shuffle_data(table_init)
    
    ## convert position to number
    
    duplication_table$startxnb<-duplication_table$begin_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]
    
    duplication_table$endxnb<-duplication_table$end_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]
    
    duplication_table$startynb<-duplication_table$begin_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]
    
    duplication_table$endynb<-duplication_table$end_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]
    
    shuffleddata$snpposnb<-shuffleddata$snppos+listval[match(shuffleddata$snpchrom,listval$chromname),3]
    
    shuffleddata$genestartnb<-shuffleddata$start+listval[match(shuffleddata$chr,listval$chromname),3]
    
    shuffleddata$geneendnb<-shuffleddata$end+listval[match(shuffleddata$chr,listval$chromname),3]
    
    
    
    shuffleddata2<-shuffleddata
    lengthdata<-length(shuffleddata2)
    for(i in 1:length(duplication_table$id)){
      duplication_tablerow1<-duplication_table[i,]
      
      test<-ifelse(shuffleddata$snpposnb>=duplication_tablerow1$startxnb & shuffleddata$snpposnb<=duplication_tablerow1$endxnb,ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$genestartnb>=chromadd[shuffleddata$snpchrom] & shuffleddata$genestartnb<=chromadd[shuffleddata$snpchrom+1],ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(shuffleddata$snpposnb>=duplication_tablerow1$startynb & shuffleddata$snpposnb<=duplication_tablerow1$endynb,ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$genestartnb>=chromadd[shuffleddata$snpchrom] & shuffleddata$genestartnb<=chromadd[shuffleddata$snpchrom+1],ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Gene present",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Gene present",ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Gene present",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Gene present","NA"))))))
      shuffleddata2[,lengthdata+i]<-test
    }
    
    
    ## then process to get one value output
    
    list_of_values <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic")
    
    # Apply the function to each row 
    shuffleddata2$result <- apply(shuffleddata2, 1, process_row)
    
    shuffleddata$duplicatestatus<-shuffleddata2[,96]
    
    shuffleddata$data<-"shuffled"
    shuffleddata<-cbind.data.frame(shuffleddata,table_init[,c(1:6,10)])
    
    fulltab<-rbind.data.frame(table_init,shuffleddata,stringsAsFactors = FALSE)
    
    
    dfcount <- fulltab%>% filter(!is.na(duplicatestatus))  %>% filter(kind!="cis")  %>% group_by(duplicatestatus,data) %>% 
      summarise(total_count=n(),.groups = 'drop') %>%
      as.data.frame() %>% filter(duplicatestatus %in% dupstatus_list) %>% group_by(data) %>% 
      mutate(prop = total_count / sum(total_count))
    
    df_wide <- pivot_wider(dfcount, id_cols = duplicatestatus, names_from = data, values_from = total_count)
    df_wide$enrichment<-df_wide$real/df_wide$shuffled
    df_wide$bootstrapnb<-y
    
    df_bootstrap<-rbind.data.frame(df_bootstrap,df_wide,stringsAsFactors = FALSE)
    
  }
  return(df_bootstrap)
}

Bootstrap_data<-Bootstrap_trans(nb=1000,table_init = newdata)


## save results
write.table(Bootstrap_data,file="Data/bootstrap_rnaLL_1000_01.txt",row.names = FALSE)





#### -- Plots and proportion -- ####
newdata<-read.table(file="Data/rna_LL_all_august_version_01.txt",header=TRUE)
shuffleddata<-read.table(file="Data/bootstrap_rnaLL_1000_01.txt",header=TRUE)

newdata$data<-"real"
shuffleddata$data<-"shuffled"
shuffleddata<-cbind.data.frame(shuffleddata,newdata[,c(1:6,10)])

fulltab<-rbind.data.frame(newdata,shuffleddata,stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(tidyr)

## Note: Inter-chromosomal ohnolog is synonymous for Inter-chromosomal homeolog

dupstatus_list <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic","Inter-chromosomal non ohnolog")

dfcount <- fulltab%>% filter(!is.na(duplicatestatus))  %>% filter(kind!="cis")  %>% group_by(duplicatestatus,data) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame() %>% filter(duplicatestatus %in% dupstatus_list) %>% group_by(data) %>% 
  mutate(prop = total_count / sum(total_count))



ggplot(dfcount,aes(x=data,y=prop,fill=duplicatestatus))+geom_bar(position = "fill", stat = "identity")+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+theme_bw(base_size = 30)+ylab("Proportion") +xlab("Data type")+guides(fill=guide_legend("Interaction type"))

dfcountreal<-dfcount[dfcount$data=="real",]
ggplot(dfcountreal,aes(x="",y=prop,fill=duplicatestatus))+geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+
  theme(legend.text = element_text(size=15))

##prop :
## ICNO: 61
## ICO : 2
## ICNS : 20
## ICS : 17

df_wide <- pivot_wider(dfcount, id_cols = duplicatestatus, names_from = data, values_from = total_count)
df_wide$enrichment<-df_wide$real/df_wide$shuffled







# Compute mean and standard error
result <- shuffleddata %>%
  group_by(duplicatestatus) %>%
  summarize(
    mean_enrichment = mean(enrichment),
    se_enrichment = sd(enrichment) / sqrt(n())
  )

ggplot(result,aes(x=duplicatestatus,y=mean_enrichment,fill=duplicatestatus))+geom_bar(stat="identity")+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+theme_bw(30)+geom_hline(yintercept=1,color="darkgrey",linetype = "dashed",size=4)+ scale_y_continuous(trans = "log2")+
  xlab("Interaction type") +ylab("Log2 of enrichment")+geom_errorbar(aes(x=duplicatestatus, ymin = mean_enrichment - 1.96*se_enrichment, ymax = mean_enrichment + 1.96*se_enrichment), color = "black")+theme(legend.position="none")#+theme(axis.text.x = element_text(angle = 45, hjust=1)) 
