# Libraries
library(dplyr)
library(tidyverse)
library(DESeq2)
library(edgeR)

# Set working directory
setwd('~/Documents/Pepper/')

# Load data
counts <- read.table('counts_matrix.txt',header = T)
metadata <- read.csv('metadata.csv', header = T)

counts <- counts %>% select(metadata$Run) 

summary(colSums(counts))
table(colSums(counts) > 10000000)


# Filtrar las columnas cuya suma sea mayor a 100000
filtered_sums <- colSums(counts)
# Crear un dataframe adecuado para el gráfico
data_to_plot <- data.frame(Sum = filtered_sums)

# Agrupar las columnas en rangos de suma
data_to_plot <- data_to_plot %>%
  mutate(Range = cut(Sum, breaks = seq(1000000, max(Sum), by = 1000000)))

# Convertir los valores de Value a millones
data_to_plot

boxplot <- 
  ggplot(data_to_plot, aes(y = Sum)) +
  geom_boxplot() +
  ylab("Library Size") +
  coord_flip() +
  theme_classic() 

density <-
  ggplot(data_to_plot, aes(y = Sum)) +
  geom_density() +
  ylab("Library Size") +
  theme_classic() +
  coord_flip() 


boxplot + density + plot_layout(ncol = 1)

########################
## GETTING CPM COUNTS ##
########################

# Filter libraries with less than 1 M reads
d <- DGEList(counts = counts)
norm.counts <- cpm(d, normalized.lib.sizes = T, log = F)


## OMITTING GENES THAT CPM < 1 IN ALL SAMPLES
table(rowSums(norm.counts == 0) == 0) #4107
norm.counts.filt <- norm.counts[rowSums(norm.counts == 0) != 0, ]
variances <- apply(norm.counts.filt, 1, var)
variances <- round(variances, digits = 2)
hist(log2(variances))
summary(var)
(variances == 0) %>% table

temp <- norm.counts[variances == 0,] 
View(temp)


# Just fruit
metadata_fruit <- metadata %>% filter(Tissue_type == 'Fruit')
counts_fruit <- norm.counts[,metadata_fruit$Run]


#################
# Create a PCA #
###############

# Transpose the matrix so that rows = samples and columns = variables
pca_matrix <- t(log2(counts_fruit+1))
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
min(pc_scores$PC1)
pc_scores %>% filter(PC1 < -500) %>% select(sample, PC1)
#metadata <- read.csv('metadata.csv', header = T) 
#dim(metadata)

metadata_fruit$create_date <- as.factor(metadata_fruit$create_date)
metadata2 <- metadata_fruit %>% 
  dplyr::rename(sample = Run) %>% 
  select(sample, Tissue, Tissue_type, Organism, BioProject,
         LibrarySelection, Fruit_Stage_V3,
         create_date, Variety,Cultivar.2, Fruit_Stage_V2)

pc_scores_2 <- merge(pc_scores, metadata2, by = 'sample')

pc_scores_2 %>% filter(Tissue_type == 'Root') %>% View

# Print pca
pc_scores_2 %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC3, color= Fruit_Stage_V3)) +
  geom_point() 


