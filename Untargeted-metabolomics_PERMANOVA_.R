################### Unsupervised multivariate analyses: PCoA and PERMANOVA
# In this notebook we will perform a Principal Coordinate Analysis (PCoA), also 
# known as metric or classical Multidimensional Scaling (metric MDS) to explore 
# and visualize patterns in an untargeted mass spectromtery-based metabolomics 
# dataset. We will then assess statistical significance of the patterns and 
# dispersion of different sample types using permutational multivariate analysis 
# of variance (PERMANOVA).

################### Updates, install and call packages 
install.packages('vegan')
install.packages('ggplot2')
install.packages("RColorBrewer")
install.packages("Polychrome")

library(Polychrome)
library(RColorBrewer)
library(vegan)
library(ggplot2)

# last update
Sys.time()

################### Input and check files:
### In the following, we import the cleaned and imputed feature table for the dataset.
# See Script_DataCleaning for more information.
rm(list=ls())

# setting the current directory as the working directory
Directory<- c("")
setwd(Directory)


dirs <- dir(path=paste(getwd(), sep=""), full.names=TRUE, recursive=TRUE)
folders <- unique(dirname(dirs))
files <- list.files(folders, full.names=TRUE)
files_1 <- basename((files))
files_2 <- dirname((files))
# Creating a Result folder
dir.create(path=paste(files_2[[1]], "/PERMANOVA_Results", sep=""), showWarnings = TRUE)
fName <-paste(files_2[[1]], "/PERMANOVA_Results", sep="")

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




############################################################################
#####################################  Filter and Process data
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

### 4) Z-scoring Data (scaling and centering)
fts <- scale(df_log, center = TRUE, scale = TRUE)

### 5) Concatenate final data
Data_fil<-  cbind.data.frame(metai2,fts)





############################################################################
##################################### PERMANOVA
nmeta<-ncol(md)+1
a<-Data_fil[, nmeta:ncol(Data_fil)]
md2<-Data_fil[,1:ncol(metai2)]
distm2 <- dist(a, method = 'euclidean')

PERMANOVA_results=NULL
adonres <- adonis2(distm2 ~ md2[,5])
R2<-adonres$R2[1]*100
p<-adonres$`Pr(>F)`[1]
temp1<-cbind(R2, p)
PERMANOVA_results<-temp1

for (x in 6:ncol(md2)) {
  adonres <- adonis2(distm2 ~ md2[,x])
  R2<-adonres$R2[1]*100
  p<-adonres$`Pr(>F)`[1]
  temp<-cbind(R2, p)
  PERMANOVA_results<-rbind(PERMANOVA_results,temp)
  
}


metadata<-colnames(md2)
metadata<-metadata[5:length(metadata)]
PERMANOVA_results<-data.frame(PERMANOVA_results)
rownames(PERMANOVA_results)<-metadata
PERMANOVA_results<-cbind(PERMANOVA_results,metadata)
colnames(PERMANOVA_results)<-c("R2","pvalue","metadata")

PERMANOVA_results2 <- PERMANOVA_results %>% 
  filter(pvalue!="NA")


Significance=NULL
for (i in 1:nrow(PERMANOVA_results2)) {
  temp<-PERMANOVA_results2$pvalue[i]
  if (temp>0.05) {
    Significance[i]="Non-Significant"}
  else {
    Significance[i]="Significant"}
}

PERMANOVA_results2$Significance<-Significance

temp<-as.character(paste("p-value = ", as.character(PERMANOVA_results2$pvalue)))
PERMANOVA_results2$pvalue2<-temp


PERMANOVA_results3 <- PERMANOVA_results2 %>% 
  filter(Significance!="Non-Significant")%>%
  arrange(R2)

