library(tidyverse)
library(dplyr)
library(edgeR)
library(corrr)
library(factoextra)
library(ggfortify)

# Set working directory
setwd('~/counts/')

#######################################
# Merge count tables in a single one #
#####################################
setwd('~/counts/pepper/')
list <- list.files(pattern = '.txt$')

# Create a base table
df <- read.table('counts.11441.txt', header = T) %>%
  dplyr::select(Geneid,
                -Chr, 
                -Start,
                -End,
                -Strand,
                -Length, 
                -c(starts_with('out'))) %>% column_to_rownames('Geneid')

## PASTE ALL THE COUNTS IN A SINGLE COUNT MATRIX FOR PEPPER ##
for (i in list) {
df <-  cbind(df, 
                  read.table(i, header = T) %>%
                                      dplyr::select(Geneid,
                                                    -Chr, 
                                                    -Start,
                                                    -End,
                                                    -Strand,
                                                    -Length, starts_with('out')) %>%
                    column_to_rownames('Geneid'))
}

# Changing the sample names
colnames(df) <- gsub('^out\\.([S|D]RR.*[0-9])\\..*$', '\\1', colnames(df))


df_pepper <- df

dim(df_pepper) #779 samples

## PASTE ALL THE COUNTS IN A SINGLE COUNT MATRIX FOR TOMATO ##
setwd('~/counts/tomato/')
list <- list.files(pattern = '.txt$')

# Create a base table
df <- read.table('counts.213528.txt', header = T) %>%
  dplyr::select(Geneid,
                -Chr, 
                -Start,
                -End,
                -Strand,
                -Length, 
                -c(starts_with('out'))) %>% column_to_rownames('Geneid')

## PASTE ALL THE COUNTS IN A SINGLE COUNT MATRIX FOR PEPPER ##
for (i in list) {
  df <-  cbind(df, 
               read.table(i, header = T) %>%
                 dplyr::select(Geneid,
                               -Chr, 
                               -Start,
                               -End,
                               -Strand,
                               -Length, starts_with('out')) %>%
                 column_to_rownames('Geneid'))
}

# Changing the sample names
colnames(df) <- gsub('^out\\.([S|D]RR.*[0-9])\\..*$', '\\1', colnames(df))

df_tomato <- df

dim(df_tomato) #130 samples

##################################
## Create METADATA table pepper ##
#################################

metadata.list <- list.files('~/counts/pepper/metadata/')
metadata <- data.frame() #empty data frame

for (i in metadata.list) {
  
  temp <- read.csv(paste('~/counts/pepper/metadata/',i, sep = '/'), header = T, fill = T) %>% 
    select(Run, 
           LibraryName,
           BioProject, 
           ScientificName,
           avgLength, 
           LibraryStrategy, 
           Platform, 
           Sample, 
           LibraryLayout,
           LibrarySelection)
  metadata <- rbind(metadata, temp)
}
dim(metadata) #881 samples 
colnames(metadata)

View(metadata)

##################################
## Create METADATA table tomato ##
##################################

metadata.list <- list.files('~/counts/tomato/metadata/')
metadata <- data.frame() #empty data frame

for (i in metadata.list) {
  
  temp <- read.csv(paste('~/counts/tomato/metadata/',i, sep = '/'), header = T, fill = T) %>% 
    select(Run, 
           LibraryName,
           BioProject, 
           ScientificName,
           avgLength, 
           LibraryStrategy, 
           Platform, 
           Sample, 
           LibraryLayout,
           LibrarySelection)
  metadata <- rbind(metadata, temp)
}
dim(metadata) #881 samples 
colnames(metadata)

View(metadata)

# FILTER RNA-Seq data!!
metadata$LibraryStrategy %>% table
metadata <- metadata %>% filter(metadata$LibraryStrategy == 'RNA-Seq')

#EXPLORE THE DATA
metadata$Platform %>% table # ALL is illumina
metadata$BioProject %>% table %>% length() #31 projects
metadata$ScientificName %>% table 
metadata$LibraryLayout %>% table
metadata$LibrarySelection %>% table
metadata$avgLength %>% summary()

########################################################
# Check that all the samples in df and metadata match #
######################################################

df <- df[,!(colnames(df) %>% table >= 2)] 

df[!(colnames(df) %in% metadata$Run)]  %>% colnames()
df <- df[!(grepl('\\.1', colnames(df)))]
dim(df) #619

metadata$Run %in% colnames(df) %>% table #695 TRUE 
metadf <- metadata
metadata <- metadata %>% filter(metadata$Run %in% colnames(df)) 

# ARE THE SAME SAMPLES ???
nrow(metadata[1:561,]) == ncol(df) #MUST BE TRUE <<<<<

################################################################################
# Normalization

########################
## GETTING CPM COUNTS ##
########################

# Filter libraries with less than 1 M reads
table(apply(df, 2, sum) > 1000000) #628
df <- df[,apply(df, 2, sum) > 1000000]


d <- DGEList(counts = df)
norm.counts <- cpm(d, normalized.lib.sizes = T, log = F)

## OMITTING GENES THAT CPM < 1 IN ALL SAMPLES
table(rowSums(norm.counts > 1) == 0) #7394
norm.counts.filt <- norm.counts[rowSums(norm.counts > 1) != 0, ]
variances <- apply(norm.counts.filt, 1, var)
variances<- round(variances, digits = 1)

table(variances == 0) #2066 genes with variance 0 

#################
# Create a PCA #
###############

# Transpose the matrix so that rows = samples and columns = variables
pca_matrix <- t(log2(norm.counts+1))
pca <- prcomp(pca_matrix)

# Look at the first 10 rows and first 5 columns of the matrix
pca_matrix[1:10, 1:5]

pc_eigenvalues <- pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

fviz_eig(pca, addlabels = TRUE, ncp = 10)

pc_scores <- pca$x
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") 


metadata <- read.csv('metadata.csv', header = T) 
dim(metadata)


pc_scores_2 <- merge(pc_scores, metadata, by = 'sample')


colnames(pc_scores_2)
# Print pca
pc_scores_2 %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC4, color=Stage)) +
  geom_point()

#autoplot(pca, data = metadata, colour = 'avgLength')

View(metadata)
write.csv(metadata, 'metadata.csv')

autoplot(pca, data = metadata, colour = 'Stage', loadings = TRUE)
