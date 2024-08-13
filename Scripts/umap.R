library(tidyverse)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(dplyr)
library(corrr)
library(factoextra)
library(ggfortify)
library(patchwork)

###############################
## Load COUNTS table pepper ##
#############################

df <- read.table('/Users/alejandra/pepper_counts/counts/counts_matrix.txt')

#################################
## Load METADATA table pepper ##
###############################

metadata <- read.csv('/Users/alejandra/pepper_counts/metadata/metadata_clean.csv')


########################################################
# Check that all the samples in df and metadata match #
######################################################

# ARE THE SAME SAMPLES ???
nrow(metadata) == ncol(df) #MUST BE TRUE <<<<<

################################################################################
# Normalization

########################
## GETTING CPM COUNTS ##
########################

# Filter libraries with less than 1 M reads
d <- DGEList(counts = df)
d <- calcNormFactors(d, method = "TMM")

norm.counts <- cpm(d, normalized.lib.sizes = T, log = F)


## OMITTING GENES THAT CPM < 1 IN ALL SAMPLES
temp <- norm.counts[!rowSums(norm.counts) == 0,] 
dim(temp)

#### PCA!!!

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
#pc_eigenvalues %>% View

fviz_eig(pca, addlabels = TRUE, ncp = 10)

pc_scores <- pca$x
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") 
colnames(pc_scores)

metadata2 <- metadata %>% dplyr::rename(sample = Run)
pc_scores_2 <- merge(pc_scores, metadata2, by = 'sample')

# Print pca
pc_scores_2 %>% 
  # create the plot
  ggplot(aes(x = PC2, y = PC6, color=Organism)) +
  geom_point() 

#  stat_ellipse()

temp <- pc_scores_2 %>% filter(PC2 > 50) %>% select(sample, PC2) %>% select(sample)
metadata2 %>% filter(sample %in% temp$sample) %>% View




#autoplot(pca, data = metadata, colour = 'avgLength')

View(metadata)
write.csv(metadata, 'metadata.csv')

autoplot(pca, data = metadata, colour = 'Stage', loadings = TRUE)


norm.counts.filt <- norm.counts[rowSums(norm.counts > 1) != 0, ]
variances <- apply(norm.counts.filt, 1, sd)
variances <- round(variances, digits = 2)
summary(variances)

temp <- norm.counts[variances < 1,] 
rowSums(temp) %>% summary


## UMAP
# Calcular UMAP
library(umap)
umap_result <- umap(norm.counts)

# Plot UMAP
plot(umap_result$layout, main = "UMAP Plot")
