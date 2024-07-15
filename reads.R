library(tidyverse)
library(dplyr)
library(edgeR)
library(corrr)
library(factoextra)
library(ggfortify)
library(patchwork)

# Set working directory
setwd('/Users/alejandra/pepper_counts/counts')

#######################################
# Merge count tables in a single one #
#####################################

list <- list.files(pattern = '^counts.*')
length(list)

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

dim(df) # 35845 genes, 1766 samples


# Changing the sample names - Removing the 'out.' from the begging
colnames(df) <- gsub('^out\\.([S|D|E]RR.*[0-9])\\..*$', '\\1', colnames(df))
colnames(df) <- gsub('(^.*)\\.1','\\1',colnames(df))

# Check duplicated columns
colnames(df)[duplicated(colnames(df))] %>% length()

# Eliminar columnas duplicadas, manteniendo la primera aparición
df <- df[ , !duplicated(colnames(df))]
dim(df)

# Save Samples ID
samples <- colnames(df) %>% as.data.frame() 
write.csv(samples,'muestras_descargadas.csv', row.names = F)


# How many samples do not have any read?
head(df)
#df <- as.matrix(df)
summary(colSums(df))
table(colSums(df) > 1000000) # Keep in mind that we have libraries with very few reads

# Filtrar las columnas cuya suma sea mayor a 100000
df[, colSums(df) > 1000000] #81 FALSE
df <- df[, colSums(df) > 1000000] #81 FALSE
dim(df) #1470

summary(colSums(df))

# Filtrar las columnas cuya suma sea mayor a 100000
filtered_sums <- colSums(df)
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
  

# Backup
df_pepper <- df


##################################
## Load METADATA table pepper ##
#################################

## ---> ARREGLAR AQUI!!

metadata <- read.csv('/Users/alejandra/pepper_counts/metadata/metadata_clean.csv')

df[!(colnames(df) %in% metadata$Run)] %>% colnames %>% as.data.frame %>% View




df <- df[colnames(df) %in% metadata$Run]  # Which are not in the metadata? 


########################################################
# Check that all the samples in df and metadata match #
######################################################



colnames(df) %in% metadata$Run %>% table
metadata$Run %in% colnames(df) %>% table #747 TRUE 
metadata_backup <- metadata
metadata <- metadata %>% filter(metadata$Run %in% colnames(df)) 

# ARE THE SAME SAMPLES ???
nrow(metadata) == ncol(df) #MUST BE TRUE <<<<<

################################################################################
# Normalization

########################
## GETTING CPM COUNTS ##
########################

# Filter libraries with less than 1 M reads
d <- DGEList(counts = df)
norm.counts <- cpm(d, normalized.lib.sizes = T, log = F)

## OMITTING GENES THAT CPM < 1 IN ALL SAMPLES
table(rowSums(norm.counts > 1) == 0) #4107
norm.counts.filt <- norm.counts[rowSums(norm.counts > 1) != 0, ]
variances <- apply(norm.counts.filt, 1, var)
variances <- round(variances, digits = 2)

(variances == 0) %>% table

temp <- norm.counts[variances == 0,] 
View(temp)



# Clean this section of the code <--------------

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
min(pc_scores_2$PC1)
pc_scores_2 %>% filter(PC1 < -500) %>% select(sample, PC1)
#metadata <- read.csv('metadata.csv', header = T) 
#dim(metadata)

metadata2 <- metadata %>% rename(sample = Run) %>% select(sample, Tissue, Tissue_overall, Tissue_Type, Organism, BioProject)
pc_scores_2 <- merge(pc_scores, metadata2, by = 'sample')


colnames(pc_scores_2)
# Print pca
pc_scores_2 %>% 
  # create the plot
  ggplot(aes(x = PC3, y = PC5, color=Organism)) +
  geom_point()

#autoplot(pca, data = metadata, colour = 'avgLength')

View(metadata)
write.csv(metadata, 'metadata.csv')

autoplot(pca, data = metadata, colour = 'Stage', loadings = TRUE)

