#### Marie Saitou And Celian Diblasi - 20/05/25

library(tidyverse)
library(ggplot2)
#####Make expression bin plot
df <- read_csv("Data/geneLL_class.csv")
#https://www.dropbox.com/s/776gg2bwmv3evly/geneLL_class.csv?dl=0
df %>%
  mutate(bin = ntile(median_all_TPM, n=10)) %>% 
  dplyr::count(cistrans,bin)  -> df2

desired_order <- c("trans", "cis", "both")

# Convert to a factor with the desired order
df2$cistrans<- factor(df2$cistrans, levels = desired_order)



ggplot(df2, aes(x = as.factor(bin), y = n, fill = cistrans)) + 
  geom_bar(position = "fill",stat = "identity") +theme_minimal()+ 
  scale_fill_manual(values = c("#f1b047","#71aaaa","#b890db","#e2e2e2"))+
  labs(x = "Gene expression fraction",y="eQTL presence")

### Make the barplot of proportion of eQTL

### load Lead SNP data
df_eQTL<-read.csv(file="Data/cistrans_gene_August.csv",header=TRUE)


### cluster eQTLs in "cis", "trans" or "both", and add proportion
df_eQTL_summary <- df_eQTL %>%
  group_by(lead_SNP) %>%
  summarise(
    kind = if (n_distinct(kind) > 1) "both" else unique(kind)
  )%>% count(kind) %>%
  mutate(proportion = n / sum(n), 
         ypos = cumsum(proportion) - proportion / 2 )

df_eQTL_summary$kind<- factor(df_eQTL_summary$kind, levels = desired_order)

### plot
ggplot(df_eQTL_summary, aes(x ="", y = proportion, fill = kind)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of types of eQTL",
    x = "All eQTLs",
    y = "Proportion of eQTL classes"
  ) +
  theme_bw(15)+scale_fill_manual(values = c("#f1b047","#71aaaa","#b890db"))+geom_text(aes(y = ypos, label = n), color = "black", size = 5)


#### effect size

mean(abs(df_eQTL[df_eQTL$kind=="cis","slope"]))
mean(abs(df_eQTL[df_eQTL$kind=="trans","slope"]))

###Chromosomal distribution of eQTL



library(dplyr)
library(tidyverse)
library(data.table)


## import data
imputed <-  read_table("Data/SNPinfo.imputed.txt") 

wgs <- read_table("Data/SNPinfo.wgsparents.txt") 

cis <-  read_table("Data/LLsmolt.cispermute_fig1.txt") 

trans <- read_table("Data/LLsmolt.trans.adjust.fig1.txt") 


### vector length of chromosomes

lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)


library(dplyr)
library(ggplot2)



window_size <- 1e7

#remove ssa and ssa0
imputed <- imputed %>%
  dplyr::mutate(CHR=str_remove_all(CHR,"ssa")) %>% 
  dplyr::mutate_at('CHR',as.numeric) 


result_list <- list()


start_pos <- 0

for (chr in unique(imputed$CHR)) {
  sub_data <- imputed[imputed$CHR == chr, ]
  max_window <- ceiling(max(sub_data$POS) / window_size)
  
  # Count by window
  for (window in 1:max_window) {
    window_start <- (window - 1) * window_size
    window_end <- window * window_size
    count <- sum(sub_data$POS > window_start & sub_data$POS <= window_end)
    result_list[[length(result_list) + 1]] <- data.frame(CHR = chr, window = window, counts = count, pos = start_pos + window)
  }
  
  # Update the start position of the next chromosome
  start_pos <- start_pos + max_window
}

# Combine results
result <- bind_rows(result_list)

# Odd/Even Chromosome Flags
result$odd_even <- ifelse(result$CHR %% 2 == 1, "Odd", "Even")



# Calculate the center coordinate for the resulting data frame
chromosome_labels <- result %>%
  group_by(CHR) %>%
  summarise(center_pos = mean(pos)) %>%
  ungroup()

# Plot
p <-  ggplot(result, aes(x = pos, y = counts, fill = odd_even)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Odd" = "black", "Even" = "gray60")) +
  labs(
    title = "imputed SNP counts/1M windows") +
  theme_minimal()

# Add x-axis labels to the center of each chromosome
imputed_plot <- p + scale_x_continuous(breaks = chromosome_labels$center_pos, labels = chromosome_labels$CHR)




##############  wgs


#remove ssa and ssa0
wgs <- wgs %>%
  dplyr::mutate(CHR=str_remove_all(CHR,"ssa")) %>% 
  dplyr::mutate_at('CHR',as.numeric) 

# Prepare an empty list to store the results
result_list <- list()

# Initialize starting position
start_pos <- 0

# Counting by window for each chromosome
for (chr in unique(wgs$CHR)) {
  sub_data <- wgs[wgs$CHR == chr, ]
  max_window <- ceiling(max(sub_data$POS) / window_size)
  
  # Count by window
  for (window in 1:max_window) {
    window_start <- (window - 1) * window_size
    window_end <- window * window_size
    count <- sum(sub_data$POS > window_start & sub_data$POS <= window_end)
    result_list[[length(result_list) + 1]] <- data.frame(CHR = chr, window = window, counts = count, pos = start_pos + window)
  }
  
  # Update the start position of the next chromosome
  start_pos <- start_pos + max_window
}

# Combine results
result <- bind_rows(result_list)

# Odd/Even Chromosome Flags
result$odd_even <- ifelse(result$CHR %% 2 == 1, "Odd", "Even")




# Calculate the center coordinate for the resulting data frame
chromosome_labels <- result %>%
  group_by(CHR) %>%
  summarise(center_pos = mean(pos)) %>%
  ungroup()

# plot
p <-  ggplot(result, aes(x = pos, y = counts, fill = odd_even)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Odd" = "black", "Even" = "gray60")) +
  labs(
    title = "WGS variants/1M windows") +
  theme_minimal()

# Add x-axis labels to the center of each chromosome
 
wgs_plot <- p + scale_x_continuous(breaks = chromosome_labels$center_pos, labels = chromosome_labels$CHR)



library(patchwork)
final_plot <- (imputed_plot / wgs_plot )



###### cis



#remove ssa and ssa0


cis <- cis %>%
  dplyr::mutate(CHR=str_remove_all(CHR,"ssa")) %>% 
  dplyr::mutate_at('CHR',as.numeric) 

# Prepare an empty list to store the results
result_list <- list()

# Initialize starting position
start_pos <- 0

# Counting by window for each chromosome
for (chr in unique(cis$CHR)) {
  sub_data <- cis[cis$CHR == chr, ]
  max_window <- ceiling(max(sub_data$POS) / window_size)
  
  # Count by window
  for (window in 1:max_window) {
    window_start <- (window - 1) * window_size
    window_end <- window * window_size
    count <- sum(sub_data$POS > window_start & sub_data$POS <= window_end)
    result_list[[length(result_list) + 1]] <- data.frame(CHR = chr, window = window, counts = count, pos = start_pos + window)
  }
  
  # Update the start position of the next chromosome
  start_pos <- start_pos + max_window
}

# Combine results
result <- bind_rows(result_list)

# Odd/Even Chromosome Flags
result$odd_even <- ifelse(result$CHR %% 2 == 1, "Odd", "Even")


# Calculate the center coordinate for the resulting data frame
chromosome_labels <- result %>%
  group_by(CHR) %>%
  summarise(center_pos = mean(pos)) %>%
  ungroup()

# plot
p <-  ggplot(result, aes(x = pos, y = counts, fill = odd_even)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Odd" = "black", "Even" = "gray60")) +
  labs(
    title = "cis eQTL/1M windows") +
  theme_minimal()

# Add x-axis labels to the center of each chromosome
 
cis_plot <- p + scale_x_continuous(breaks = chromosome_labels$center_pos, labels = chromosome_labels$CHR)


cis_plot


library(patchwork)
final_plot <- (cis_plot/ imputed_plot / wgs_plot )




###### trans


trans <- trans %>%
  dplyr::mutate(CHR=str_remove_all(CHR,"ssa")) %>% 
  dplyr::mutate_at('CHR',as.numeric) 

# Prepare an empty list to store the results
result_list <- list()

# Initialize starting position
start_pos <- 0

# Counting by window for each chromosome
for (chr in unique(trans$CHR)) {
  sub_data <- trans[trans$CHR == chr, ]
  max_window <- ceiling(max(sub_data$POS) / window_size)
  
  # Count by window
  for (window in 1:max_window) {
    window_start <- (window - 1) * window_size
    window_end <- window * window_size
    count <- sum(sub_data$POS > window_start & sub_data$POS <= window_end)
    result_list[[length(result_list) + 1]] <- data.frame(CHR = chr, window = window, counts = count, pos = start_pos + window)
  }
  
  # Update the start position of the next chromosome
  start_pos <- start_pos + max_window
}

# Combine results
result <- bind_rows(result_list)

# Odd/Even Chromosome Flags
result$odd_even <- ifelse(result$CHR %% 2 == 1, "Odd", "Even")

# Calculate the center coordinate for the resulting data frame
chromosome_labels <- result %>%
  group_by(CHR) %>%
  summarise(center_pos = mean(pos)) %>%
  ungroup()

# plot
p <- ggplot(result, aes(x = pos, y = counts, fill = odd_even)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Odd" = "black", "Even" = "gray60")) +
  labs(
    title = "trans eQTL/1M windows") +
  theme_minimal()

# Add x-axis labels to the center of each chromosome
 
trans_plot <- p + scale_x_continuous(breaks = chromosome_labels$center_pos, labels = chromosome_labels$CHR)



trans_plot + theme_bw(20)

library(patchwork)
final_plot <- (trans_plot / cis_plot) +
  plot_layout(heights = c(3, 1,1))
final_plot





#### TSS analysis 

### important gff file from assembly v3.1
gffile<-read.table("Salmo_salar.Ssal_v3.1.106_filtered.gff.gz",fill=TRUE)

gffilegene<-gffile[gffile$V3=="gene",]

df_eQTL<-read.csv(file="Data/cistrans_gene_August.csv",header=TRUE)
df_eQTL$snppos<-as.numeric(as.character(sapply(strsplit(df_eQTL$lead_SNP,"_"),FUN = `[[`, 2)))

df_eQTL_cis<-df_eQTL[df_eQTL$kind=="cis",]

library(tidyverse)
library(ggplot2)

### use strand information to know which position to keep 

df_eQTL_cis <- df_eQTL_cis %>%
  mutate(
    TSS = pmap_chr(
      list(gene, start, end),
      function(val, st, en) {
        match_row <- gffilegene %>%
          filter(str_detect(V9, fixed(val))) %>%
          slice_head(n = 1)
        
        if (nrow(match_row) == 0) {
          return(NA_character_)
        }
        
        strand <- match_row$V7
        
        if (strand == "+") {
          return(as.character(st))
        } else if (strand == "-") {
          return(as.character(en))
        } else {
          return(NA_character_)
        }
      }
    )
  )

df_eQTL_cis$TSS<-as.numeric(as.character(df_eQTL_cis$TSS))
df_eQTL_cis$distfromTSS<-df_eQTL_cis$TSS-df_eQTL_cis$snppos

ggplot(df_eQTL_cis,aes(x=distfromTSS))+geom_density(adjust=0.5)
