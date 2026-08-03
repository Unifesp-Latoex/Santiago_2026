
########################## Univariate Analysis of metabolomics data
# In this script, we perform statistical analysis on metabolomics dataset that 
# includes: 
# 1) Input and check files
# 2) Getting p-value with Kruskal–Wallis test (non-parametric version of ANOVA)
# 3) Differential Analysis by Dunn's test 
# 4) Boxplots and heatmap of significant analysis





############################################################################
############################################################################
############################################################################
###################################### Updates, install and call packages 

install.packages("ggsci")
install.packages("matrixStats")
install.packages("ggrepel")
install.packages("tidyverse")
if(!require(FSA)) install.packages("FSA")


suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(matrixStats)
  library(ggrepel)
  library(VennDiagram)
  library(magick)
  library(dplyr)
  library(ComplexHeatmap)
  library(dendextend)
  library(NbClust)
  library(factoextra)
  library(FSA)
})

Sys.time()
rm(list=ls())











############################################################################### 
############################################################################### 
############################################################################### 
######################################################### Input and check files:

### In the following, we import the cleaned and imputed feature table for the dataset.
# See Script_DataCleaning for more information.
Directory<- c("")

dirs <- dir(path=paste(getwd(), sep=""), full.names=TRUE, recursive=TRUE)
folders <- unique(dirname(dirs))
files <- list.files(folders, full.names=TRUE)
files_1 <- basename((files))
files_2 <- dirname((files))
# Creating a Result folder
dir.create(path=paste(files_2[[1]], "/Univariate_Results", sep=""), showWarnings = TRUE)
fName <-paste(files_2[[1]], "/Univariate_Results", sep="")


# We need to import the feature table and metadata
# If you performed data cleaning, we import the cleaned and imputed feature table
ft_url <- paste0(Directory, "/DataCleanup_Results/Step5_Normalised_quantile.csv")
# metadata
md_url <-  paste0(Directory, "/DataCleanup_Results/metadata2.csv")


# read data
ft <- read.csv(ft_url, header = T, check.names = F, row.names = 1)
md <- read.csv(md_url, header = T, check.names = F, row.names = 1)

# dimension of the data
dim(ft)
dim(md)

head(ft)
head(md)


# how many files in the metadata are also present in the feature table
table(rownames(md) %in% colnames(ft))
# which file names in the metadata are not in the feature table?
setdiff(rownames(md),colnames(ft))
md <- md[rownames(md) %in% colnames(ft),]
dim(md)


# transpose ft
ft <- t(ft)
head(ft)
dim(ft)





################  Filter and Process data
### 1) Concatenate ft and md and filter the data based on the desired metadata
# (here you could remove outliers)
library(dplyr)
DataI <- cbind.data.frame(md,ft)
dim(DataI)

Datai2<-DataI %>%
  filter(ATTRIBUTE_Condition!="POOL")



### 2) Remove the features that have zero intensities for all the samples in the dataset tested.
# This step increases statistical power by removing the samples that undetected features
nmeta<-ncol(md)+1
b <- Datai2[ ,nmeta:ncol(Datai2)]

# calculate the mean of each feature in all the samples
Mean_Total <-mean(b[,1])
for(i in 2:ncol(b)) { 
  meanvalue<-mean(b[,i])
  Mean_Total <- cbind(Mean_Total, meanvalue)
}

# get the features that have zero absence in all the features
# for quantile normalization, the intensity that represents the zero absence is 
# the value you selected during imputation
# for example, 3000 or 5000
Index_Nonzeros = NULL
for(i in 1:ncol(Mean_Total)) { 
  if (Mean_Total[i]>10000) {
    Index_Nonzeros <- cbind(Index_Nonzeros, i)
  }
}

# Get new feature table and metadata
d<-b[,Index_Nonzeros]
metai2<-Datai2[,1:ncol(md)]


### 3) log-transformation of the data
df_log <- log(d + 1)


### 4) Concatenate final data
Data<-  cbind.data.frame(metai2,df_log)









############################################################################### 
############################################################################### 
############################################################################### 
################## Getting p-value with Kruskal–Wallis test

### Run Kruskal–Wallis test
# We can run a for loop to pass each feature column into the first argument of the
# aov function, while the second argument, time point, is constant.
nmeta<-ncol(md)+1
a <- nmeta:ncol(Data)

# Run Kruskal–Wallis tests
kruskal_out <- vector("list", length(a))

for (i in seq_along(a)) {
  kruskal_out[[i]] <- kruskal.test(
    Data[[a[i]]] ~ as.factor(Data$ATTRIBUTE_concat3)
  )
}

names(kruskal_out) <- a   # label results by metabolite name


##### Convert results to data.frame ----
kruskal_df <- data.frame(
  metabolite = names(kruskal_out),
  statistic  = sapply(kruskal_out, function(x) x$statistic),
  p_value    = sapply(kruskal_out, function(x) x$p.value),
  stringsAsFactors = FALSE
)

kruskal_df
rownames(kruskal_df) <- colnames(Data[nmeta:ncol(Data)])
kruskal_df["p_fdr"] <- p.adjust(kruskal_df$p_value,method="fdr")
kruskal_df["significant"] <- ifelse(kruskal_df$p_fdr<0.05,"Significant","Nonsignificant")
kruskal_df$metabolite<-colnames(Data[nmeta:ncol(Data)])


# Save Kruskal–Wallis Table
kruskal_df_sorted <- kruskal_df %>% arrange (p_fdr)
write.csv(kruskal_df_sorted, file.path(fName,'Output_kruskal_df_sorted.csv'),row.names =TRUE)



### Plot kruskal_df results
plot_kruskal_df <- ggplot(kruskal_df,aes(x=log(statistic,base=10),y=-log(p_value,base=10),color=significant))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = c("pink", "royalblue"))+
  #scale_color_jama()+
  ylab("-log(p)")+
  xlab("log(F)")+
  geom_text_repel(data=kruskal_df %>% arrange(p_value) %>% slice_head(n=10),
                  aes(label=metabolite),size=3,show.legend = FALSE,max.overlaps = 100)+
  theme(legend.title = element_blank())

plot_kruskal_df

# Save ANOVA plot
pdf(file = paste0(fName, "/plot_kruskal_df.pdf"))
plot_kruskal_df
dev.off()








###############################################################################
###############################################################################
###############################################################################
######################### Post-hoc analysis: Dunn's test + Volcano plots


# Identify feature columns
feature_cols <- colnames(Data)[nmeta:ncol(Data)]

# Extract group levels
groups <- unique(Data$ATTRIBUTE_concat3)
group_pairs <- combn(groups, 2, simplify = FALSE)


## Run Dunn's test for each metabolite
dunn_results_list <- list()
for (met in feature_cols) {
  
  dunn_res <- dunnTest(
    x = Data[[met]],
    g = Data$ATTRIBUTE_concat3,
    method = "bh" #Benjamini–Hochberg (most common FDR method)
  )
  
  dunn_table <- dunn_res$res
  dunn_table$metabolite <- met
  
  dunn_results_list[[met]] <- dunn_table
}

# Combine all results
dunn_results <- bind_rows(dunn_results_list)

colnames(dunn_results) <- c("Comparison", "Z", "p_value", "p_adj", "metabolite")

write.csv(dunn_results,
          file.path(fName, "Output_Dunn_full_results.csv"),
          row.names = FALSE)







###############################################################################
###############################################################################
###############################################################################
######################### Volcano plots for C-section

###### Get pair comparison
  pair<-unlist(group_pairs[1])
  g1 <- pair[1]
  g2 <- pair[2]
  comparison_name <- paste(g1, "vs", g2, sep = "_")
  
  
###### Filter Dunn results for this comparison
  pair_res <- dunn_results %>%
    filter(Comparison == paste(g1, "-", g2) |
             Comparison == paste(g2, "-", g1))
  
###### Compute log2 Fold Change (mean difference)
  ColumnfromMetadata <- 9
  DataFirst<-Data %>%
    filter(Data[ColumnfromMetadata]==g1)
  
  DataFirst<-DataFirst[,nmeta:ncol(DataFirst)]
  
  DataSecond<-Data %>%
    filter(Data[ColumnfromMetadata]==g2)
  DataSecond<-DataSecond[,nmeta:ncol(DataSecond)]
  
  FCm=NULL
  mUP=NULL
  mDOWN=NULL
  for (i in 1:ncol(DataFirst)){
    mUP[i]<-mean(DataFirst[,i])
    mDOWN[i]<-mean(DataSecond[,i])
    xx <- mUP[i]/mDOWN[i]
    FCm <- cbind(FCm, xx)
  }
  colnames(FCm)<-colnames(Data[nmeta:ncol(Data)])
  
  
  FC <- sapply(FCm, as.numeric)
  logFC <- log(abs(FC),base=2)

  volcano_df <- data.frame(
    metabolite = pair_res$metabolite,
    log2FC = logFC,
    p_adj = pair_res$p_adj
  )
  volcano_df$negLog10FDR <- -log10(volcano_df$p_adj)
  volcano_df$Significant <- ifelse(volcano_df$p_adj < 0.05 & abs(volcano_df$log2FC) > 0,
                                   "Significant", "Not Significant")
  
  
  

# Add results
output_dunn_SPF<-pair_res  
output_dunn_SPF$log2FC<- volcano_df$log2FC
output_dunn_SPF$Group<- Grouping
output_dunn_SPF$Significant<- volcano_df$Significant
  
write.csv(output_dunn_SPF,
          file.path(fName, "Output_Dunn_C-SECTION.csv"),
          row.names = FALSE)

  
###### Get significance for each metabolite based on p-value ad FC
t<-output_dunn_SPF


# Coloring
vector=NULL
limit<-0
for(i in 1:nrow(t)) { 
  if(t$Significant[i] == "Significant" && t$log2FC[i]< 0){
    vector[i]<-c(g2)} 
  else if (t$Significant[i] == "Not Significant"){
    vector[i]<-c("Not Significant")}
  else if (t$Significant[i] == "Significant" && t$log2FC[i]> 0){
    vector[i]<-c(g1)}
} 

t$stats_volcano <- vector

 # Get values for plot
temp=-log(t$p_adj)
temp[is.infinite(temp)]<-NA
my_max=max(temp, na.rm=TRUE)
stratigh<-seq(0, my_max, by=2)
log2FC<-t$log2FC

limit<--1
f1<-seq(0,0)
# f2<-seq(limit,limit)
# f3<-seq(1,1)


plot_tukey <- ggplot(t,aes(x=log2FC,y=-log(p_value,base=10),color=stats_volcano))+ 
  geom_point()+
  theme_minimal()+
  scale_colour_manual(values = c("darkorange",  "mediumpurple3","gray"))+
  ylab("-log(p)")+
  ggtitle(comparison_name) +
  theme(legend.title = element_blank())+
  geom_text_repel(data=t %>% arrange(p_adj) %>% slice_head(n=0),
                  aes(label=metabolite),size=3,show.legend = FALSE,max.overlaps = 20)+
  geom_vline(xintercept = f1, linetype="dashed", 
             color = "gray",linewidth=1)
  # geom_vline(xintercept = f2, linetype="dashed", 
  #          color = "gray",linewidth=1)+
  # geom_vline(xintercept = f3, linetype="dashed", 
  #          color = "gray",linewidth=1)

plot_tukey


  # Save plot
  ggsave(filename = paste0(fName,"/Volcano_", comparison_name, ".pdf"),
         plot = plot_tukey,
         width = 7,
         height = 6)

  
  
  t_down_CS_late <- output_dunn_SPF %>%
    filter(log2FC<0, Significant=='Significant') %>%
    arrange (p_adj)
  
  t_up_CS_early <- output_dunn_SPF %>%
    filter(log2FC>0, Significant=='Significant') %>%
    arrange (p_adj)


