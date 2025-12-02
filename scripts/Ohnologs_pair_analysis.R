#### Celian Diblasi
#### 15/03/2024
#### Investigate ohnlogs gene expression



## -- function to find pair -- ##

## import fine mapped data
data<-read.table(file="Data/synchroLL_finemap.txt",header=TRUE)

cisdata<-read.table(file="Data/synchroLL_cis_finemap.txt",header=TRUE)

data$kind<-"trans"
cisdata$kind<-"cis"


cistransdata<-rbind.data.frame(cisdata,data,stringsAsFactors = FALSE)



##function to find pair

## import ohnolog pair data

genetab<-read.table(file="Data/ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
genetab<-genetab[genetab$type=="ss4r",]
genetab<-genetab[,-c(1,2,6)]


## this function find the status of the ohnolog pair
getPairstatus2<-function(threshold){
  
  ## filter data based on a threshold a significance
  data_subset<- cistransdata[cistransdata$posterior_prob >= threshold, ]
  
  ## loop on all pairs
  for(y in 1:length(genetab$ID)){
    
    ## get gene name and eQTl regulating these genes
    gene1row1<-genetab[y,"gene1"]
    tabpair1<-data_subset[data_subset$gene==gene1row1,]
    gene2row1<-genetab[y,"gene2"]
    tabpair2<-data_subset[data_subset$gene==gene2row1,]  
    
    ltab1<-length(tabpair1$SNP)
    ltab2<-length(tabpair2$SNP)
    
    ## define base value (if nothing is found)
    reg<-"No eQTL"
    bothdiffreg<-"nobothdiffreg"
    bothsamereg<-"nobothsamereg"
    
    if(ltab1==0 && ltab2==0){
    }
    if(ltab1==0 && ltab2!=0){ ##if a regulator is find for only one gene of the pair
      reg<-"Single copy eQTL"
    }
    if(ltab1!=0 && ltab2==0){ ##if a regulator is find for only one gene of the pair
      reg<-"Single copy eQTL"
    }
    if(ltab1!=0 && ltab2!=0){ ##if a regulator is find for both genes
      issamesnp<-ifelse(tabpair1$SNP%in%tabpair2$SNP,"yes","no") ##ifa similar SNP is regulating both genes
      if("no"%in%issamesnp){
        bothdiffreg<-"Distinct eQTL"
      }
      if("yes"%in%issamesnp){
        bothsamereg<-"Shared eQTL"
      }
    }
    
    
    
    ##say that if regulated by same snp then we keep only this result
    if(bothdiffreg=="Distinct eQTL" && bothsamereg=="Shared eQTL"){
      reg<-"Shared eQTL"
    }else if(bothdiffreg=="nobothdiffreg" && bothsamereg=="Shared eQTL"){
      reg<-"Shared eQTL"
    }else if(bothdiffreg=="Distinct eQTL" && bothsamereg=="nobothsamereg"){
      reg<-"Distinct eQTL"
    }
    
    ##select a random common SNP if shared
    if(reg=="Shared eQTL"){
      common_SNP<-sample(tabpair1[issamesnp=="yes","SNP"],1)
    }else{
      common_SNP<-NA
    }
    
    
    genetab[y,"regulation"]<-reg
    genetab[y,"Common_SNP"]<-common_SNP
  }
  
  return(genetab)
}


### this part is used to make a table to see the effect of the threshold (see Supplementary figure 3)
for(y in seq(from=0,to=1,by=0.05)){
  result<-getPairstatus2(y)
  result_reg<-result$regulation
  result_commonSNP<-result$Common_SNP
  namecol_reg<-paste("regulation",y,sep="_")
  namecol_SNP<-paste("Common_SNP",y,sep="_")
  genetab[,namecol_reg]<-result_reg
  genetab[,namecol_SNP]<-result_commonSNP
}

## save the data
write.table(genetab,file="Data/genetabgenepair05.txt",row.names=FALSE)

## make a smaller table to use with only regulation and common SNP 0

##genetabsub<-genetab[,c(1:5)]

##write.table(genetabsub,file="Data/genetabgenepair05_0cutoff.txt",row.names=FALSE)



### Get gene expression for both copies -- ###

library(tidyverse)

###load expression, genotype and covariate data

Expression <- read_table2("Data/Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
Genotype <- read_table2("Data/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
Covariate <-  read_table2("Data/covariate_synchro.txt") 
## fix the name in genotype
names(Genotype)[2:ncol(Genotype)] <- str_remove(names(Genotype)[2:ncol(Genotype)], "0_")
# scale the expression
Expression[, 2:ncol(Expression)] <- t(apply(Expression[, 2:ncol(Expression)], 1, scale))


## import the pair status data
genepair<-read.table(file="Data/genetabgenepair05.txt",header=TRUE)


genepair1<-Expression[Expression$pid%in%genepair$gene1,]
genepair2<-Expression[Expression$pid%in%genepair$gene2,]



# Pivot genotype table to long format
#geno_table_long_p1 <- genepair1 %>%
#pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")

#Pivot expression table to long format
expr_table_long1 <- genepair1 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

expr_table_long2 <- genepair2 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

covar_kept<-Covariate[c(3,5),]
id<-colnames(covar_kept)
tcovar<-t(covar_kept)
long_covar<-cbind.data.frame(id,tcovar)
long_covar<-long_covar[-1,]
colnames(long_covar)<-c("id","sex","LL")

long_covar_LL<-long_covar[,c(1,3)]

#longcovarll<-long_covar_LL[long_covar_LL$LL=="LL",] ## 907 individuals



mergedpair1 <- merge(expr_table_long1 , long_covar, by = "id")
mergedpair2 <- merge(expr_table_long2 , long_covar, by = "id")

##filter

mergedpair1_filter<- mergedpair1[mergedpair1$LL=="LL",]
mergedpair2_filter<- mergedpair2[mergedpair2$LL=="LL",]

##can keep more column if wanted
mergedpair1_filter<-mergedpair1_filter[,1:3]
mergedpair2_filter<-mergedpair2_filter[,1:3]


library(dplyr)

# Left join table 1 with table 2 based on gene1 and pid
joined_table1_2 <- left_join(genepair, mergedpair1_filter, by = c("gene1" = "pid"))
joined_table1_2<-joined_table1_2[!is.na(joined_table1_2$expression),]

# Left join result with table 3 based on gene2 and pid
joined_table1_2_3 <- left_join(joined_table1_2, mergedpair2_filter, by = c("id"="id","gene2" = "pid"))
joined_table1_2_3 <- joined_table1_2_3 [!is.na(joined_table1_2_3$expression.y),]


## save data
colnames(joined_table1_2_3)[47:48]<-c("Expression_gene1","Expression_gene2")
write.table(joined_table1_2_3,file="Data/expression_gene_pair02.txt",row.names=FALSE)


### -- Analysis and plots -- ###

library(ggplot2)
library(dplyr)
library(tidyr)

## load data

gene_pair_tab<-read.table(file="Data/expression_gene_pair02.txt",header=TRUE)

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))

genepair<-read.table(file="Data/genetabgenepair05.txt",header=TRUE)


results2<-left_join(result, genepair, by = c("ID"="ID"))

### extract the pairs that share eQTL to check for GO enrichment

sharedpairs<-results2[results2$regulation_0=="Shared eQTL",c("gene1","gene2")]

sharedpairs_export <- c(sharedpairs$gene1, sharedpairs$gene2)

writeLines(sharedpairs_export, "Data/sharedpairs_export.txt")

# same for background

allgenes <- results2[,c("gene1","gene2")]

allgenes_export <- c(allgenes$gene1, allgenes$gene2)

writeLines(allgenes_export, "Data/allgenes_export.txt")

## Check how count change according to the threshold
genetab<-results2

library(dplyr)

# Extract regulation column names
regulation_columns <- grep("^regulation_", colnames(genetab), value = TRUE)

# Extract unique eQTL types
eqtl_types <- unique(genetab$regulation_0)

# Initialize an empty data frame to store results
summary_table <- data.frame(regulation_number = numeric(0))

# Add columns for each eQTL type
for (eqtl_type in eqtl_types) {
  summary_table[[eqtl_type]] <- integer(0)
}

# Loop through regulation columns
for (reg_col in regulation_columns) {
  regulation_number <- as.numeric(gsub("^regulation_", "", reg_col))
  row_data <- list(regulation_number)
  
  for (eqtl_type in eqtl_types) {
    nb_eqtl_type <- sum(genetab[[reg_col]] == eqtl_type, na.rm = TRUE)
    row_data[[eqtl_type]] <- nb_eqtl_type
  }
  
  # Add a row to the summary table
  summary_table <- rbind(summary_table, do.call(cbind, row_data))
}


# Reshape the data
summary_table_long <- summary_table %>%
  pivot_longer(cols = -V1, names_to = "eQTL_type", values_to = "count")

# Rename the eQTL types 
eqtl_type_names <- c("No eQTL", "Distinct eQTL", "Single copy eQTL", "Shared eQTL")
summary_table_long$eQTL_type <- factor(summary_table_long$eQTL_type, levels = eqtl_type_names)
colnames(summary_table_long)[1]<-"threshold"

summary_table_long$eQTL_type<-as.factor(summary_table_long$eQTL_type)
summary_table_long$eQTL_type<-factor(summary_table_long$eQTL_type,levels=c("Distinct eQTL","No eQTL","Shared eQTL", "Single copy eQTL"))


library(ggplot2)

ggplot(summary_table_long,aes(x=threshold,y=count,color=eQTL_type))+geom_line(size=2)+geom_point(size=4)+scale_color_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+theme_bw(20)




### now make analysis and plot using only p=0

# Define the desired order of levels
desired_order <- c("Shared eQTL", "Distinct eQTL", "Single copy eQTL", "No eQTL")

# Convert regulation_0 to a factor with the desired order
results2$regulation_0 <- factor(results2$regulation_0, levels = desired_order)



### get the average coefficient for each category
df2 <- results2 %>%
  group_by(regulation_0) %>%
  summarise(Mean = mean(spearman_coeff))

df2$regulation_0<- factor(df2$regulation_0, levels = desired_order)


ggplot(results2,aes(x=spearman_coeff,fill=regulation_0))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation_0,scales = 'free_y',ncol=4)+xlab("Spearman coefficient")+
  scale_fill_manual(values=c("#8CE86D","#F9B27C","#D186D8","#99B8EA"))+geom_vline(data = df2, mapping = aes(xintercept = Mean),color="firebrick3",size=1.5)+coord_flip()


results2<-left_join(result, genepair, by = c("ID"="ID"))

dfcount2 <- results2%>% group_by(regulation_0) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()  %>% 
  mutate(prop = total_count / sum(total_count))

ggplot(dfcount2,aes(x="",y=prop,fill=regulation_0))+geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+
  theme(legend.text = element_text(size=25))


### prop:

#Distinct: 29%
#No eQTL: 24%
## Single eQTL: 42%
## Shared eQTL: 5%


##make a slightly different plot with top10% and top 25% information

percent10<-round(0.1*length(results2$gene2))
percent25<-round(0.25*length(results2$gene2))

df_ordered <- results2[order(results2$spearman_coeff, decreasing = TRUE), ]

df_extracted <- head(df_ordered, n = percent10)

dfcounttop10 <- df_extracted %>% group_by(regulation_0) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()  %>% 
  mutate(prop = total_count / sum(total_count))

df_extracted2 <- head(df_ordered, n = percent25)

dfcounttop25 <- df_extracted2 %>% group_by(regulation_0) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()  %>% 
  mutate(prop = total_count / sum(total_count))

dfcount2$data<-"all"
dfcounttop10$data<-"Top 10%"
dfcounttop25$data<-"Top 25%"

catdfcount<-rbind.data.frame(dfcount2,dfcounttop10,dfcounttop25)


# Define the desired order of levels
desired_order2 <- c("Top 10%", "Top 25%", "all")

# Convert regulation_0 to a factor with the desired order
catdfcount$data <- factor(catdfcount$data, levels = desired_order2)


ggplot(catdfcount,aes(x = data, y = prop, fill = regulation_0))+geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~regulation_0)+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+
  theme_bw(25)+geom_text(aes(label = signif(total_count)), nudge_y = 0.015, size=7)+theme(axis.text.x = element_text(angle = 45, hjust=1))




### make ks test to see differences between distribution

results3<-as.data.frame(results2)

ks.test(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"],results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"],results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"],results3[results3$regulation_0=="No eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"],results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"],results3[results3$regulation_0=="No eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"],results3[results3$regulation_0=="No eQTL","spearman_coeff"])

## same but with absolute coefficient
ks.test(abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]))


x<-abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"])
x<-abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"])
x<-abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"])
x<-abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]) 

mean(x)
sd(x)/sqrt(length(x))



#### null expectations - correlation between random gene pairs
library(tidyverse)
gene_pair_tab<-read.table(file="Data/expression_gene_pair02.txt",header=TRUE)

head(gene_pair_tab)

listgene1<-unique(gene_pair_tab$gene1)
listgene2<-unique(gene_pair_tab$gene2)


expressiongene1<-gene_pair_tab[,c("gene1","id","Expression_gene1")]
expressiongene2<-gene_pair_tab[,c("gene2","id","Expression_gene2")]
colnames(expressiongene1)<-c("gene1","sampleid","Expression_gene1")
colnames(expressiongene2)<-c("gene2","sampleid","Expression_gene2")

### shuffle the genes to make random pairs
shuffled_list1 <- sample(listgene1)
shuffled_list2 <- sample(listgene2)

pairs <- data.frame(gene1 = shuffled_list1, gene2 = shuffled_list2)
pairs$ID<-paste(pairs$gene1,"Ssal",pairs$gene2,"Ssal",sep="_")


all_combinations <- pairs %>%
  tidyr::crossing(sampleid = unique(c(expressiongene1$sampleid, expressiongene2$sampleid)))

# merge tables
shuffled_genepair_expression <- all_combinations %>%
  left_join(expressiongene1, by = c("gene1", "sampleid")) %>%
  left_join(expressiongene2, by = c("gene2", "sampleid"))



# calculate spearman coefficient as for observed data
result_shuffled <- shuffled_genepair_expression %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))

summary(result_shuffled$spearman_coeff)



x<-abs(result_shuffled$spearman_coeff)

mean(x)
sd(x)/sqrt(length(x))

### actually the trend is even stronger with non-absolute correlation

### compared to observed data


genepair<-read.table(file="Data/genetabgenepair05.txt",header=TRUE)

results2<-left_join(result, genepair, by = c("ID"="ID"))


x<-abs(results2$spearman_coeff)

mean(x)
sd(x)/sqrt(length(x))




#### Compute the number of pair in AORe and LOre region for each category


library(ggplot2)
library(dplyr)
library(tidyr)

gene_pair_tab<-read.table(file="expression_gene_pair02.txt",header=TRUE)

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))


genepair<-read.table(file="Data/genetabgenepair05.txt",header=TRUE)

results2<-left_join(result, genepair, by = c("ID"="ID"))
results2<-results2[,c(1:6)]

genetab<-read.table(file="Data/ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
genetab<-genetab[genetab$type=="ss4r",]
genetab<-genetab[,-c(1,2)]


mergeddata<-merge(results2,genetab,by="ID")
mergeddatafilter<-mergeddata[!is.na(mergeddata$redip.class),] ### from 4627 pairs to 3020


aorepairs<-mergeddatafilter[mergeddatafilter$redip.class=="AORe",]
lorepairs<-mergeddatafilter[mergeddatafilter$redip.class=="LORe",]

table(lorepairs$regulation_0)
table(aorepairs$regulation_0)

dfprop<-data.frame(data=c("Distinct eQTL","Distinct eQTL","Single copy eQTL","Single copy eQTL","No eQTL","No eQTL","Shared eQTL","Shared eQTL"),
                   redip.class=c("AORe","LORe","AORe","LORe","AORe","LORe","AORe","LORe"),prop=c(0.71,0.29,0.73,0.27,0.78,0.22,0.62,0.38))


dfprop$data<-as.factor(dfprop$data)

dfprop$data<-factor(dfprop$data,levels=c("No eQTL","Distinct eQTL","Single copy eQTL","Shared eQTL"))

ggplot(dfprop,aes(x=data,y=prop,fill=redip.class))+geom_bar(position = "fill", stat = "identity")+theme_bw(25)+ylab("Proportion")+scale_fill_manual(values=c("#5387b4","#db6d48"))+coord_flip()



ggplot(mergeddatafilter,aes(x=spearman_coeff,fill=redip.class))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~redip.class,scales = 'free_y',ncol=4)+xlab("Spearman coefficient")+
  scale_fill_manual(values=c("#5387b4","#db6d48"))+coord_flip()

chidata<-data.frame(LORe=c(59,241,342,194),AORe=c(97,603,936,548))

chisq.test(chidata)  ## X-squared = 9.7479, df = 3, p-value = 0.02084


#### how many eQTL present for all the shared pairs ?


sharedpairs<-mergeddata[mergeddata$regulation_0=="Shared eQTL",]
length(unique(sharedpairs$Common_SNP_0)) ### 60 unique eQTLs

table(sharedpairs$Common_SNP_0) ## ssa14_51354372 and ssa14_51370830 responsible for 73+91 = 164 pairs


##### check % of shared pairs in LOre + check expression level if we take only AORe pairs from all regions

library(ggplot2)
library(dplyr)
library(tidyr)

gene_pair_tab<-read.table(file="Data/expression_gene_pair02.txt",header=TRUE)

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))


genepair<-read.table(file="Data/genetabgenepair05.txt",header=TRUE)

results2<-left_join(result, genepair, by = c("ID"="ID"))
results2<-results2[,c(1:6)]

genetab<-read.table(file="Data/ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
genetab<-genetab[genetab$type=="ss4r",]
genetab<-genetab[,-c(1,2)]


mergeddata<-merge(results2,genetab,by="ID")
mergeddatafilter<-mergeddata[!is.na(mergeddata$redip.class),] ### from 4627 pairs to 3020


aorepairs<-mergeddatafilter[mergeddatafilter$redip.class=="AORe",]
lorepairs<-mergeddatafilter[mergeddatafilter$redip.class=="LORe",]



x<-abs(aorepairs[aorepairs$regulation_0=="Shared eQTL","spearman_coeff"]) #### 0.75 +-0.01
x<-abs(aorepairs[aorepairs$regulation_0=="Distinct eQTL","spearman_coeff"]) #### 0.49 +-0.01
x<-abs(aorepairs[aorepairs$regulation_0=="Single copy eQTL","spearman_coeff"]) #### 0.48 +- 0.01
x<-abs(aorepairs[aorepairs$regulation_0=="No eQTL","spearman_coeff"]) #### 0.56 +- 0.01

mean(x)
sd(x)/sqrt(length(x))


x<-abs(lorepairs[lorepairs$regulation_0=="Shared eQTL","spearman_coeff"]) #### 0.70 +-0.03
x<-abs(lorepairs[lorepairs$regulation_0=="Distinct eQTL","spearman_coeff"]) #### 0.52 +-0.01
x<-abs(lorepairs[lorepairs$regulation_0=="Single copy eQTL","spearman_coeff"]) #### 0.53 +- 0.01
x<-abs(lorepairs[lorepairs$regulation_0=="No eQTL","spearman_coeff"]) #### 0.56 +- 0.02
