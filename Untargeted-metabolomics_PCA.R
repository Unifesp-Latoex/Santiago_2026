################### Unsupervised multivariate analyses: PCoA and PERMANOVA
# In this notebook we will perform a Principal Component Analysis (PCA), also 
# known as metric or classical Multidimensional Scaling (metric MDS) to explore 
# and visualize patterns in an untargeted mass spectromtery-based metabolomics 
# dataset. We will then assess statistical significance of the patterns and 
# dispersion of different sample types using permutational multivariate analysis 
# of variance (PERMANOVA).

################### Updates, install and call packages 


# install.packages("FactoMineR")
# install.packages("ggcorrplot")
# install.packages("corrr")
# install.packages('vegan')
# install.packages('ggplot2')
# install.packages("RColorBrewer")
# install.packages("Polychrome")
# install.packages("factoextra")
# install.packages("patchwork")

library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(Polychrome)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(patchwork)
library(factoextra)
library(dplyr)

# last update
Sys.time()


############################################################################
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
dir.create(path=paste(files_2[[1]], "/PCA_Results", sep=""), showWarnings = TRUE)
fName <-paste(files_2[[1]], "/PCA_Results", sep="")

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
##################################### PCA of all samples



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

### 4) Z-scoring Data (scaling and centering)
fts <- scale(df_log, center = TRUE, scale = TRUE)

### 5) Concatenate final data
Data_fil<-  cbind.data.frame(metai2,fts)

write.csv(fts, file.path(Directory,'DataCleanup_Results/Step6_Filter_Logtransfor_Scaled.csv'),row.names =TRUE)





################ Compute PCA
nmeta<-ncol(md)+1
aa<-Data_fil[,nmeta:ncol(Data_fil)]
res.pca <- prcomp(aa, scale = F)


# Results for Variables
res.var <- get_pca_var(res.pca)
Coordinates_var<-res.var$coord # Coordinates 
Contribution_var<-res.var$contrib  # Contributions to the PCs
quality_var<-res.var$cos2 # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(res.pca)
Coordinates_ind <- res.ind$coord  # Coordinates of the samples in the score
Contribution_ind<-res.ind$contrib # Contributions to the PCs of the samples in the score
quality_ind<-res.ind$cos2 # Quality of representation 

# Get Pc coordinates
PcoA_points <- as.data.frame(Coordinates_ind)
PcoA_points<-PcoA_points[,1:4]
names(PcoA_points)[1:4] <- c('PCoA1', 'PCoA2', 'PCoA3', 'PCoA4')
a<-data.frame(res.pca$sdev)
variance_temp<-get_eig(res.pca)
variance<-as.numeric(variance_temp$variance.percent)

# Get explained variance for PC1, PC2 and PC3
temp1<-variance[1]
temp2<-format.default(temp1, digits = 4)
temp3<-variance[2]
temp4<-format.default(temp3, digits = 4)
temp5<-variance[3]
temp6<-format.default(temp5, digits = 4)

# Save PCA table
write.csv(PcoA_points, file.path(fName,'PcoA_points.csv'),row.names =TRUE)


######## Screeplot
# Check the explained variance in each of the PCs
#Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(res.pca)
# SAVE
pdf(file = paste0(fName, "/PCA_ExplainedVariable.pdf"))
fviz_eig(res.pca)
dev.off()

######## Score Plot  
md_fil<-Data_fil[,1:ncol(md)]
my_colors <- c(
  "#d62728", "#9467bd","#1f77b4", "#ff7f0e", "#2ca02c",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#5254a3", "#9c9ede", "#8ca252", "#bd9e39", "#ad494a",
  "#a55194", "#6b6ecf", "#cedb9c", "#e7cb94", "#e7969c"
)



####################### Score Plot  
md_fil<-Data_fil[,1:ncol(md)]
# Using the ggplot2 library, we can plot our PCoA using the Euclidean distance (=PCA).
plotPCA1<-ggplot(PcoA_points, aes(x = PCoA1, y = PCoA2, colour = md_fil$ATTRIBUTE_TimePoint2)) +
  geom_point(aes(shape=md_fil$ATTRIBUTE_Condition),size=4)+
  #geom_point(size=4)+
  scale_color_manual(values = my_colors) +
  ggtitle("All birth modes, all postnatal days")+
  xlab(paste('PCA1, Explained variance',temp2,'%', sep = ' ')) +
  ylab(paste('PCA2, Explained variance,',temp4,'%', sep = ' ')) +
  theme(legend.title=element_blank())+
  stat_ellipse(aes(group = md_fil$ATTRIBUTE_concat3))

plotPCA1

# SAVE
pdf(file = paste0(Directory, "PCA_Results/PCA_Score_PC1PC2_allsamples2.pdf"))
plotPCA1
dev.off()






