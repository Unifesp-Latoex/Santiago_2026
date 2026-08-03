################### DATA CLEANING
# This script is used for cleaning the feature table, an output of metabolomics 
# experiment, containing all the features with their corresponding intensities.
# The data cleanup steps involved are: 
# 1) Organize data
# 2) Weight correction
# 3) Blank subtraction
# 4) Removal of bad injections
# 5) Imputation 
# 6) Normalization (quantile and total area)





##################################################################################
######################################### Package installation
# installing and calling the necessary packages:
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("IRdisplay")) install.packages("IRdisplay")
if (!require("BiocManager", quietly = TRUE))
{
  install.packages("BiocManager")
}
if(!require("preprocessCore", quietly = TRUE))
{
  BiocManager::install("preprocessCore")
}




library(BiocManager)
library(preprocessCore)
library(ggplot2)
library(dplyr)
library(IRdisplay)
library(tibble)


rm(list=ls())




##################################################################################
######################################### Input files needed:
### 1) Feature table: An output of metabolomics experiment, containing all the 
# features or peaks (LC-MS peaks here) with their corresponding intensities. 
# Ex: The feature table used in the test data is obtained by MZmine or MPP. (Filetype: .csv file)
### 2) Metadata: Created by the user about the files used obtaining the feature 
# table (It can be a csv/txt/tsv file). The columns in a metadata should be created
# with the following format: filename (1st column having all the filenames in 
# the same order as the columns in feature table), all the other columns with 
# column name such as: ATTRIBUTE_yourDesiredAttribute.

# setting the current directory as the working directory
# Here you should put the directory you want to save your data.
rm(list=ls())
Directory<- c("")
setwd(Directory)

# Getting all the files in the folder
dirs <- dir(path=paste(getwd(), sep=""), full.names=TRUE, recursive=TRUE)
folders <- unique(dirname(dirs))
files <- list.files(folders, full.names=TRUE)
files_1 <- basename((files))
files_2 <- dirname((files))
# Creating a Result folder
dir.create(path=paste(files_2[[1]], "/DataCleanup_Results", sep=""), showWarnings = TRUE)
fName <-paste(files_2[[1]], "/DataCleanup_Results", sep="")

# We need to import the feature table and metadata
ft_url <- paste0(Directory, "Filtered_Entities.csv")
# metadata
md_url <-paste0(Directory, "metadata.csv")


# Upload data and check data (lets check if the data has been read correctly)
ft <- read.csv(ft_url, header = T, check.names = F, row.names = 1)
md <- read.csv(md_url, header = T, check.names = F, row.names = 1)

head(ft)
dim(ft)
head(md)
dim(md)





##################################################################################
######################################### Creating Functions:
# Before getting into the Data cleanup steps, we have created a function that can
# be used later for data summarization. By creating functions, we don't have to 
# write these big codes multiple times. Instead, we just use the function name. 
# The following cells in this section will not produce any outputs here. 
# The outputs will be produced when we give input variables to these functions in
# the later sections.

# Using this function, we get an idea of the multiple levels in each of the mentioned
# attributes in the metadata as well as the datatype of each attribute. 
# This function takes metadata table as its input.

#Function: InsideLevels
InsideLevels <- function(metatable){
  lev <- c()
  typ<-c()
  for(i in 1:ncol(metatable)){
    x <- levels(droplevels(as.factor(metatable[,i])))
    if(is.double(metatable[,i])==T){x=round(as.double(x),2)}
    x <-toString(x)
    lev <- rbind(lev,x)
    
    y <- class(metatable[,i])
    typ <- rbind(typ,y)
  }
  out <- data.frame(INDEX=c(1:ncol(metatable)),ATTRIBUTES=colnames(metatable),LEVELS=lev,TYPE=typ,row.names=NULL)
  return(out)
}





##################################################################################
#########################################   Arranging metadata and feature table in the same order:
# Next, we are trying to bring the feature table and metadata in the
# correct format such as the rownames of metadata and column names of feature table
# are the same. They both are the file names and they need to be the same, as from
# now on, we will call the columns in our feature table based on our metadata 
# information. Thus, using the metadata, the user can filter their data easily. 
# You can also directly deal with your feature table without metadata by getting
# your hands dirty with some coding!! But having a metadata improves the 
# user-experience greatly.

# 1) First, let's have a look at the different conditions within each attribute of our metadata.
InsideLevels(md)

new_ft <- ft
new_md <- md

#Removing Peak area extensions
colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft))
rownames(new_md) <- gsub(' Peak area','',rownames(new_md))
new_md <- new_md[,colSums(is.na(new_md))<nrow(new_md)] #Removing if any NA columns present in the md file
rownames(new_md) <- trimws(rownames(new_md), which = c("both")) #remove the (front & tail) spaces, if any present, from the rownames of md

#Removing the spaces (if any) from each column of md and converting them all to UPPERCASE
for(i in 1:ncol(new_md)){
  if(is.factor(new_md[,i]) | is.character(new_md[,i]) == T){
    new_md[,i] <- trimws(new_md[,i], which = c("both")) #First remove spaces in the front and end of each column of md
    new_md[,i] <- gsub(' ','_', new_md[,i]) # Replace the spaces (in the middle) to underscore
    new_md[,i] <- factor(casefold(new_md[,i], upper=T)) #convert all to UPPERCASE
  } else if (is.numeric(new_md[,i]) | is.integer(new_md[,i]) | is.double(new_md[,i]) == T){
    new_md[,i] <- new_md[,i]
  }
}

#Changing the row names of the files
rownames(new_ft) <- paste(new_ft$'row ID',round(new_ft$'Mass',digits = 3),round(new_ft$'Retention Time',digits = 3),new_ft$Compound, sep = '_')

#Picking only the files with column names containing 'mzXML'
new_ft <- new_ft[,grep('mzXML',colnames(new_ft))]
InsideLevels(new_md)
new_ft<- new_ft[,order(colnames(new_ft))] #ordering the ft by its column names
new_md <-new_md[order(rownames(new_md)),] #ordering the md by its row names

#lists the colnames(ft) which are not present in md
unmatched_ft <- colnames(new_ft)[which(is.na(match(colnames(new_ft),rownames(new_md))))] 
cat("These", length(unmatched_ft),"columns of feature table are not present in metadata:")
if((length(unmatched_ft) %% 2) ==0)
{matrix(data=unmatched_ft,nrow=length(unmatched_ft)/2,ncol=2)}else
{matrix(data=unmatched_ft,nrow=(length(unmatched_ft)+1)/2,ncol=2)}

flush.console()
Sys.sleep(0.2)

#lists the rownames of md which are not present in ft
unmatched_md <- rownames(new_md)[which(is.na(match(rownames(new_md),colnames(new_ft))))] 
cat("These", length(unmatched_md),"rows of metadata are not present in feature table:")
if((length(unmatched_md) %% 2) ==0)
{matrix(data=unmatched_md,nrow=length(unmatched_md)/2,ncol=2)}else
{matrix(data=unmatched_md,nrow=(length(unmatched_md)+1)/2,ncol=2)}

#Removing those unmatching columns and rows:
if(length(unmatched_ft)!=0){new_ft <- subset(ft, select = -c(which(is.na(match(colnames(ft),rownames(md))))) )}
if(length(unmatched_md)!=0){new_md <- md[-c(which(is.na(match(rownames(md),colnames(ft))))),]}

#checking the dimensions of our new ft and md:
cat("The number of rows and columns in our original ft is:",dim(ft),"\n")
cat("The number of rows and columns in our new ft is:",dim(new_ft),"\n")
cat("The number of rows and columns in our new md is:",dim(new_md))

# Notice that the number of columns of feature table is same as the number of rows
# in our metadata. Now, we have both our feature table and metadata in the same order.
#checking if they are the same
if(identical(colnames(new_ft),rownames(new_md))==T){
  print("The column names of ft and rownames of md are the same")}else{print("The column names of ft and rownames of md are not the same")}

#Check data once again
head(new_ft)
dim(new_ft)
head(new_md)
dim(new_md)



##################################################################################
######################################### STEP 1) Correct intensities based on tissue weigth

#Getting the weights
weight <- as.numeric(new_md$ATTRIBUTE_Weight)
#Normalizing the weights
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
weights_norm <-range01(weight)
weights_norm[weights_norm==0]<-1

#Dividing each element of a particular column with its column sum
Normalized_data_weights = NULL
for (i in 1:ncol(new_ft)){
  x <- new_ft[,i] / weights_norm[i]
  Normalized_data_weights <- cbind(Normalized_data_weights, x)
}


Normalized_data_weights<-data.frame(Normalized_data_weights)
rownames(Normalized_data_weights)<-rownames(new_ft)
colnames(Normalized_data_weights) <- colnames(new_ft)
dim(Normalized_data_weights)


#SAVE
write.csv(Normalized_data_weights, file.path(fName,'Step1_weightnormalization.csv'),row.names =TRUE)








##################################################################################
######################################### STEP 2) Blank Subtraction 

## Splitting the data into Blanks and Samples using Metadata:
# For the first step: Blank removal, we need to split the data as spectra obtained 
# from blanks and samples respectively using the metadata. More about Blank removal
# in the next section.
InsideLevels(new_md)

#  split the blanks from the sample 
#If subset_data exists, it will take it as "data", else take new_md as "data"
if(exists("subset_data")==T){data <-subset_data}else{data <-new_md}
InsideLevels(data)

#Condition <- as.double(unlist(readline("Enter the index number of the attribute to split sample and blank:")))
Condition <- 1

Levels_Cdtn <- levels(droplevels(as.factor(data[,Condition[1]])))
Levels_Cdtn 

#Among the shown levels of an attribute, select the ones to keep
Blk_id <- 3
####### ATTENTION: HERE YOU MUST WRITE THE NUMBER OF THE INDEX OF THE ATTRIBUTE TO SLIT SAMPLE AND BLANK
paste0('You chosen blank is:',Levels_Cdtn[Blk_id])

#Splitting the data into blanks and samples based on the metadata
md_Blank1 <- data[(data[,Condition] == Levels_Cdtn[Blk_id]),]
Blank1 <- Normalized_data_weights[,which(colnames(Normalized_data_weights)%in%rownames(md_Blank1)),drop=F] 
md_Samples1 <- data[(data[,Condition] != Levels_Cdtn[Blk_id]),]
Samples <- Normalized_data_weights[,which(colnames(Normalized_data_weights)%in%rownames(md_Samples1)),drop=F] 

head(Blank1)
dim(Blank1)
head(Samples)
dim(Samples)

# Check intensity of the blanks (if there is no carryover, it should be much lower than samples)
sum_Blank=NULL
for(i in 1:ncol(Blank1)) { 
  sum_Blank[i] <- sum(Blank1[,i])}
indexB<-seq(1, nrow(Blank1), by=1)

mean_Blank=NULL
for(i in 1:ncol(Blank1)) { 
  mean_Blank[i] <- mean(Blank1[,i])}

log_base<-2
IBlank <- t(Blank1)
IBlank<-data.frame(IBlank)
IBlank$log<-log(sum_Blank, as.numeric(log_base))
IBlank$SampleID<-rownames(IBlank)
IBlank$SumI<-sum_Blank
IBlank$MeanI<-mean_Blank



ggplot_theme <- theme(panel.background = element_blank(),
                      legend.background = element_blank(),
                      legend.key = element_blank(),
                      axis.line = element_line(linewidth = 2, lineend = "square"),
                      axis.text = element_text(size = 12),
                      axis.title = element_text(size = 18))


### Blank plot
tmp_plot <- ggplot(data = IBlank, aes(x = as.character(SampleID), y = MeanI)) +
  geom_boxplot()+
  geom_jitter()+
  ggplot_theme +
  ylab(paste0("log", log_base, "(peak intensity)")) +
  xlab("Samples") +
  theme(axis.text.x = element_blank())

tmp_plot

minI<-min(mean_Blank)
minI
toleranceB<-minI*200

#Remove bad blanks
goodBlanks<-IBlank %>%
  filter(IBlank$MeanI<toleranceB)
namesBlank<-rownames(goodBlanks)
namesBlank

Blank<-Blank1 %>%
  select(all_of(namesBlank))



## Blank Removal
#  LC-MS/MS, we use solvents called Blanks which are usually injected time-to-time
# to prevent carryover of the sample. The features coming from these Blanks would also be detected by LC-MS/MS instrument. Our goal here is to remove these features from our samples. The other blanks that can be removed are: Signals coming from growth media alone in terms of microbial growth experiment, signals from the solvent used for extraction methods and so on. Therefore, it is best practice to measure mass spectra of these blanks as well in addition to your sample spectra.
# How do we remove these blank features?
# Since we have the feature table split into Control blanks and Sample groups now,
# we can compare blanks to the sample to identify the background features coming
# from blanks. A common filtering method is to use a cutoff to remove features 
# that are not present sufficient enough in our biological samples.

# The steps followed in the next few cells are:
## 1) We find an average for all the feature intensities in your blank set and 
# sample set. Therefore, for n no.of features in a blank or sample set, we get 
# n no.of averaged features.
## 2) Next, we get a ratio of this average_blanks vs average_sample. This ratio 
# Blank/sample tells us how much of that particular feature of a sample gets its
# contribution from blanks. If it is more than 30% (or Cutoff as 0.3), we consider
# the feature as noise.
## 3) The resultant information (if ratio > Cutoff or not) is stored in a bin such 
# as 1 = Noise or background signal, 0 = Feature Signal
## 4) We count the no.of features in the bin that satisfies the condition ratio > 
# cutoff, and consider those features as 'noise or background features' and remove them.

#When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
#Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3
Cutoff <-0.2

#Getting mean for every feature in blank and Samples
Avg_blank <- rowMeans(Blank, na.rm= FALSE, dims = 1) # set na.rm = FALSE to check if there are NA values. When set as TRUE, NA values are changed to 0
Avg_samples <- rowMeans(Samples, na.rm= FALSE, dims = 1)

#Getting the ratio of blank vs Sample
Ratio_blank_Sample <- (Avg_blank+1)/(Avg_samples+1)

# Creating a bin with 1s when the ratio>Cutoff, else put 0s
Bg_bin <- ifelse(Ratio_blank_Sample > Cutoff, 1, 0 )
Blank_removal <- cbind(Samples,Bg_bin)

# Checking if there are any NA values present. Having NA values in the 4 variables will affect the final dataset to be created
temp_NA_Count <-cbind(Avg_blank ,Avg_samples,Ratio_blank_Sample,Bg_bin)

print('No of NA values in the following columns:')
print(colSums(is.na(temp_NA_Count)))

#Calculating the number of background features and features present
print(paste("No.of Background or noise features:",sum(Bg_bin ==1,na.rm = TRUE)))
print(paste("No.of features after excluding noise:",(nrow(Samples) - sum(Bg_bin ==1,na.rm = TRUE)))) 

Blank_removal <- Blank_removal %>% filter(Bg_bin == 0) # Taking only the feature signals
Blank_removal <- as.matrix(Blank_removal[,-ncol(Blank_removal)]) # removing the last column Bg_bin 


## You will get an output saying the:
# "No of NA values in the following columns:"
# "No.of Background or noise features:"
# No.of features after excluding noise:"

write.csv(Blank_removal,file.path(fName,'Step2_Blanks_Removed.csv'),row.names =TRUE)
head(Blank_removal)
dim(Blank_removal)











##################################################################################
######################################### STEP 3) Remove bad injections

# Check the intensity sum of all samples
sum_Samples=NULL
for(i in 1:ncol(Blank_removal)) { 
  sum_Samples[i] <- sum(Blank_removal[,i])}
indexS<-seq(1, nrow(Blank_removal), by=1)

mean_Samples=NULL
for(i in 1:ncol(Blank_removal)) { 
  mean_Samples[i] <- mean(Blank_removal[,i])}

log_base<-2
ISamples <- t(Blank_removal)
ISamples<-data.frame(ISamples)
ISamples$log<-log(sum_Samples, as.numeric(log_base))
ISamples$SampleID<-rownames(ISamples)
#ISamples$SampleID<-md_Samples1$ATTRIBUTE_InjOrd
ISamples$SumI<-sum_Samples
ISamples$MeanI<-mean_Samples


ggplot_theme <- theme(panel.background = element_blank(),
                      legend.background = element_blank(),
                      legend.key = element_blank(),
                      axis.line = element_line(linewidth = 2, lineend = "square"),
                      axis.text = element_text(size = 12),
                      axis.title = element_text(size = 18))

#### Boxplot, first version
tmp_plot2 <- ggplot(data = ISamples, aes(x = SampleID, y = MeanI)) +
  geom_boxplot(outlier.size = 2) +
  ggplot_theme +
  ylab(paste0("Mean of peak intensity")) +
  xlab("Samples") +
  theme(axis.text.x = element_blank())
tmp_plot2


#Remove bad samples
badSamples<-ISamples %>%
  filter(ISamples$MeanI<toleranceB*10)
namesbadSamples<-rownames(badSamples)
namesbadSamples


goodSamples<-ISamples %>%
  filter(ISamples$MeanI>toleranceB*10)
namesSamples<-rownames(goodSamples)
namesSamples


tempp2<-as_tibble(Blank_removal)
rownames(tempp2)<- rownames(Blank_removal)

Samples1<-tempp2 %>%
  select(all_of(namesSamples))
rownames(Samples1)<-rownames(tempp2)

md_Samples <- md_Samples1 %>% rownames_to_column(var = "File") %>% 
  as_tibble() %>% 
  filter(File %in% namesSamples)
temp2<-md_Samples$File
md_Samples<-md_Samples[,-1]
rownames(md_Samples)<-temp2



write.csv(Samples1,file.path(fName,'Step3_Badinjections_removed.csv'),row.names =TRUE)
head(Samples1)
dim(Samples1)









##################################################################################
######################################### STEP 4) Imputation
# For several reasons, real world datasets might have some missing values in it,
# in the form of NA, NANs or 0s. Eventhough the gapfilling step of MZmine fills 
# the missing values, we still end up with some missing values or 0s in our 
# feature table. This could be problematic for statistical analysis.
# In order to have a better dataset, we cannot simply discard those rows or 
# columns with missing values as we will lose a chunk of our valuable data. 
# Instead we can try imputing those missing values. Imputation involves replacing
# the missing values in the data with a meaningful, reasonable guess. There are 
# several methods, such as:

## 1) Mean imputation (replacing the missing values in a column with the mean or
# average of the column)
## 2) Replacing it with the most frquent value
## 3) Several other machine learning imputation methods such as k-nearest neighbors
# algorithm(k-NN), Hidden Markov Model(HMM)

# Here, we use ft and see the frequency distribution of its features with the ggplot.
# It shows where the features are present in higher number.
Cutoff=0
#creating bins from -1 to 10^10 using sequence function seq()
bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

#cut function cuts the give table into its appropriate bins
scores_gapfilled <- cut(as.matrix(Samples1),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10')) 

#transform function convert the tables into a column format: easy for visualization 
FreqTable<-transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values

colnames(FreqTable)[1] <- 'Range_Bins'
FreqTable #Uncomment the line if you want to see the FreqTable used for the following ggplot.

## GGPLOT2
ggplot(FreqTable, aes(Range_Bins, Log_Freq))+ 
  geom_bar(stat="identity",position = "dodge", width=0.3) + 
  scale_fill_brewer(palette = "Set1") +
  ggtitle(label="Frequency plot - Gap Filled") +
  xlab("Range") + ylab("(Log)Frequency") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5))
# The plot will show you the minimum value and this minimum value will be used for imputation.

#Create the imputed table
Imputed <- Samples1
Imputed<-data.frame(Imputed)
rownames(Imputed)<-rownames(Samples1)
colnames(Imputed)<-colnames(Samples1)

# Set the imputation value
Cutoff_LOD <- 3000
print(paste0("The minimum value greater than 0 in gap-filled table: ",Cutoff_LOD)) 

# Find the indexes that have values lower than the Cutoff_LOD
index<-which(Imputed <Cutoff_LOD)

#generate one random number between 2990 and 3010
min_cutoff<-Cutoff_LOD-10
max_cutoff<-Cutoff_LOD+10
Cutoff_LOD_random <- runif(n=length(index), min=min_cutoff, max=max_cutoff)

#imput values on matrix
Imputed[Imputed <Cutoff_LOD]<- Cutoff_LOD_random
head(Imputed)
dim(Imputed)

#removing all the rows with only cutoff values:
Imputed<-Imputed[rowMeans(Imputed)!= Cutoff_LOD,]
dim(Imputed)

#save imputed table
write.csv(Imputed,file.path(fName,paste0('Step4_Imputed_QuantTable_filled_with_',Cutoff_LOD,'_CutOff_Used_',Cutoff,'.csv')),row.names =TRUE)










##################################################################################
######################################### STEP 6) Normalization


##### 1) Total area normalization
o<-ncol(md)+1
input_table<-Imputed
saved_colnames<-colnames(Imputed)
saved_rownames<-rownames(Imputed)

#Getting column-wise sums of the input-data
sample_sum <- as.numeric(colSums(input_table, na.rm= TRUE, dims = 1))
sample_sum


#Dividing each element of a particular column with its column sum
Normalized_data_TA=NULL
Normalized_data_TA <- input_table[,1] / sample_sum[1]
for (i in 2:ncol(input_table)){
  x <- input_table[,i] / sample_sum[i]
  Normalized_data_TA <- cbind(Normalized_data_TA, x)
}

Normalized_data_TA<-data.frame(Normalized_data_TA)
colnames(Normalized_data_TA) <- colnames(input_table)
rownames(Normalized_data_TA)<-rownames(input_table)
print(paste('No.of NA values in Normalized data:',sum(is.na(Normalized_data_TA)== TRUE)))
dim(Normalized_data_TA)


# SAVE
write.csv(Normalized_data_TA, file.path(fName,'Step5_Normalised_Quant_table_TA.csv'),row.names =TRUE)




##### 2) Quantile Normalization
input_table<-Imputed
saved_colnames<-colnames(Imputed)
saved_rownames<-rownames(Imputed)

# Normalize input table
Normalized_data_quantile <- data.frame(normalize.quantiles(as.matrix(input_table), copy = TRUE))
colnames(Normalized_data_quantile) <- saved_colnames
rownames(Normalized_data_quantile) <- saved_rownames


write.csv(Normalized_data_quantile, file.path(fName,'Step5_Normalised_quantile.csv'),row.names =TRUE)
write.csv(md_Samples, file.path(fName,'Metadata2.csv'),row.names =TRUE)


