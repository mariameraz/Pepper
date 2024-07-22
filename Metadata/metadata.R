## This scripts was created to clean the metadata table obtained from NCBI

# Load libraries
library(tidyverse)
library(dplyr)

##################################
## Create METADATA table pepper ##
#################################

# Define the columns we want to keep in metadata

columns_to_select <- read.csv('/Users/alejandra/pepper_counts/counts/colnames_sra_metadata.csv')
columns_to_select <- columns_to_select$x 

# Initialize an empty data frame with the correct structure
metadata <- data.frame(matrix(ncol = length(columns_to_select), nrow = 0))
colnames(metadata) <- columns_to_select

metadata.list <- list.files('/Users/alejandra/pepper_counts/metadata/Projects/', pattern = '.csv')
length(metadata.list) # 111 projects
metadata.list
coln <- list()
for (i in metadata.list) {
  temp <- read.csv(paste('/Users/alejandra/pepper_counts/metadata/Projects/', i, sep = ''), header = TRUE, fill = TRUE)
  
  # Print the first few rows of the temporary data frame
  coln[[i]] <- colnames(temp)
  
  # Check which of the columns are available in the current data frame
  available_columns <- columns_to_select[columns_to_select %in% colnames(temp)]
  
  # Select only the available columns
  temp <- temp %>% select(all_of(available_columns))
  
  # Add missing columns with NA values
  missing_columns <- setdiff(columns_to_select, available_columns)
  if (length(missing_columns) > 0) {
    temp[missing_columns] <- NA
  }
  
  # Ensure the columns are in the correct order
  temp <- temp %>% select(all_of(columns_to_select))
  
  # Print the first few rows of the modified data frame
  #print(head(temp))
  
  # Bind the rows to the metadata data frame
  metadata <- rbind(metadata, temp)
}


# Delete duplicated samples
duplicated(metadata$Run) %>% table()
metadata <- metadata[!(duplicated(metadata$Run)),]
dim(metadata) #2413 x 147

# Back up 
metadata2 <- metadata
# FILTER RNA-Seq data!!
metadata <- metadata %>% filter(metadata$Assay.Type == 'RNA-Seq') # 2327

### Cleaning the data

# Clean ReleaseDate column:
metadata$ReleaseDate <- gsub('(^.*)-.*-.*', '\\1', metadata$ReleaseDate)
metadata$ReleaseDate %>% table(exclude = F) # Remove the 2012 samples, they are too old.
metadata <- metadata %>% filter(!grepl('2012', ReleaseDate))

# Cleaning create_data colun:
# Just keep the year
metadata$create_date <- gsub('(^.*)-.*-.*', '\\1', metadata$create_date)

# Clean Organism column:
metadata <- metadata %>% filter(grepl('Capsicum', Organism))


# Eliminar proyectos que no nos sirven:
metadata <- metadata %>% 
                  filter(!BioProject %in% c('PRJEB26324', 
                                            'PRJNA704710', 
                                            'PRJNA235215', 
                                            'PRJNA509740',
                                            'PRJNA641558',
                                            'PRJNA635538'))

# Keeping only samples that are in the count matrix
counts <- read.table(file = '/Users/alejandra/pepper_counts/counts/counts_matrix.txt')

# Obtener los nombres de las columnas que no están en metadata$Run
setdiff(colnames(counts), metadata$Run)

# Filtering samples in metadata that are in the count matrix
metadata <- metadata %>% filter(Run %in% colnames(counts))
dim(metadata)

#EXPLORE THE DATA
metadata$Assay.Type %>% table
metadata$Platform %>% table(exclude = F) # Not all is illumina
metadata$BioProject %>% table() %>% length() #106 projects
metadata$Organism %>% table(exclude=F)
metadata$LibraryLayout %>% table(exclude = F)
metadata$LibrarySelection %>% table(exclude = F)
metadata$Instrument %>% table(exclude = F)

metadata$Instrument <- gsub('^HiSeq X Ten$', 'Illumina HiSeq X Ten', metadata$Instrument)

metadata$Center.Name %>% table(exclude = F)
metadata$Sample_name %>% table(exclude = F) # Esta columna esta extrana
metadata$tissue_type %>% table(exclude = F)
metadata$Cultivar %>% table(exclude = F)
metadata$variety %>% table(exclude = F)
metadata$replicate %>% table(exclude = F)
metadata$sample_type %>% table(exclude = F)
metadata$tissue %>% table(exclude = F)
metadata$geo_loc_name_country %>% table(exclude = F)
metadata$Ecotype %>% table(exclude = F)
metadata$treatment %>% table(exclude = F)
metadata$Sample.Name %>% table(exclude = F)
(metadata$Bytes %>% as.numeric %>% na.omit() %>% sum())
metadata$Bytes %>% is.na() %>% table()
metadata$Bases %>% as.numeric %>% na.omit() %>% sum()
metadata$Bases %>% is.na() %>% table()
metadata$dev_stage %>% table(exclude = F)
metadata$isolate %>% table(exclude = F)
summary(metadata$AvgSpotLen %>% na.exclude()) 
metadata$Age %>% table(exclude = F)
metadata$Collection_Date %>% table(exclude = F) 
metadata$Consent %>% table(exclude = F)
metadata$LibrarySource %>% table(exclude = F)
metadata$create_date %>% table(exclude = F)
metadata$Host_disease %>% table(exclude = F) 

metadata <- metadata %>% select(-c(SAMPLE, 
                       sample_description,
                       sample_disease_status,
                       sample_health_state,
                       Experimental_Factor._infect..exp.,
                       ENA.FIRST.PUBLIC..run.,
                       ENA.LAST.UPDATE..run.,
                       ENA_last_update,
                       ENA.LAST.UPDATE..run.,
                       ena_first_public, Extension_Project_of_GeneChip,
                       Diet, Biological_Replicate, Barcode, cell_type,
                       Host_disease,
                       host_infra_specific_name,
                       host_sex,
                       host_taxid,
                       Host_Diet,
                       Replica,
                       collected_time_after_inoculation,
                       altitude, Aliquote, BREED, source_material_id,
                       samp_store_temp,
                       samp_mat_process, samp_collect_device,
                       Consent,
                       time, disease_stage, DATASTORE.provider,
                       DATASTORE.region, SRA.Study, version,
                       GEO_Accession..exp., geographic_location_.country_and.or_sea., 
                       geographic_location_.latitude.,
                       geographic_location_.longitude.))


metadata <- metadata %>% 
  mutate(Organism = ifelse(BioProject %in% 
                             c('PRJB85966', 
                               'PRJNA837070',
                               'PRJNA922936', 
                               'PRJNA846120',
                               'PRJNA987024',
                               'PRJNA675380') & 
                             Organism == "Capsicum",
                           "Capsicum annuum", Organism)) %>%
  mutate(Organism = ifelse(BioProject %in% 
                             c('PRJDB11441') &
                             Organism == "Capsicum",
                           "Capsicum chinense", Organism)) %>%
  mutate(Organism = ifelse(Organism %in% 
                             c('Capsicum annuum var. annuum',
                               'Capsicum annuum var. glabriusculum'),
                           "Capsicum annuum", Organism)) %>%
  mutate(Organism = ifelse(BioProject %in%  
                             c('PRJNA770309') | 
                             Organism %in% c('Capsicum chinense x Capsicum baccatum var. baccatum',
                                           'Capsicum baccatum var. baccatum x Capsicum annuum var. annuum',
                                           'Capsicum baccatum var. baccatum x Capsicum chinense', 
                                           'Capsicum annuum var. annuum x Capsicum baccatum var. baccatum',
                                           'Capsicum annuum var. annuum x Capsicum chinense'),
                           "Capsicum Hibrid", Organism)) %>%
  mutate(Organism = ifelse(Organism %in% "Capsicum baccatum var. baccatum",
                           "Capsicum baccatum", Organism)) %>%
  mutate(Organism = ifelse(BioProject %in% 
                             c('PRJNA641558') &
                             Organism == "Capsicum cardenasii",
                           "Capsicum chinense", Organism)) 


metadata$Organism %>% table(exclude = F) 


metadata$Tissue <- metadata$tissue
metadata$Tissue <- stringr::str_to_title(metadata$Tissue)

metadata$Tissue <- gsub('(Root)*[0-9]', 'Root', metadata$Tissue)
metadata$Tissue <- gsub('Root Root', 'Root', metadata$Tissue)
metadata$Tissue <- gsub('Roots', 'Root', metadata$Tissue)
metadata$Tissue <- gsub('FruitRoot', 'Root', metadata$Tissue)
metadata$Tissue <- gsub('RootRoot', 'Root', metadata$Tissue)
metadata$Tissue <- gsub('Fruit Root-Root', 'Root', metadata$Tissue)
metadata$Tissue <- gsub('Root Transcriptome Sequence', 'Root', metadata$Tissue)
metadata$Tissue <- gsub('Pepper leaves', 'Leaf', metadata$Tissue)
metadata$Tissue <- gsub('Pepper Leaf', 'Leaf', metadata$Tissue)
metadata$Tissue <- gsub('Leaves', 'Leaf', metadata$Tissue)
metadata$Tissue <- gsub('(Fruit)*[0-9]', 'Fruit', metadata$Tissue)
metadata$Tissue <- gsub('Steam', 'Stem', metadata$Tissue)
metadata$Tissue <- gsub('Apical Stem Tissue Rootmm From Tip', 'Stem', metadata$Tissue)
metadata$Tissue <- gsub('Stems', 'Stem', metadata$Tissue)
metadata$Tissue <- gsub('Whole PlantRoot', 'Seedling', metadata$Tissue)
metadata$Tissue <- gsub('Whole RootRoot Days Seedlings', 'Seedling', metadata$Tissue)
metadata$Tissue <- gsub('The Exocarp Of Pepper', 'Exocarp', metadata$Tissue)
metadata$Tissue <- gsub('Fruits', 'Fruit', metadata$Tissue)
metadata$Tissue <- gsub('Fruit \\(Pericarp\\)', 'Pericarp', metadata$Tissue)
metadata$Tissue <- gsub('FlowerRoot', 'Flower', metadata$Tissue)


metadata <- metadata %>% 
  mutate(Tissue = ifelse(Experiment %in% 
                           c('SRX10073984',
                             'SRX10073982',
                             'SRX10073980',
                             'SRX10073978'),
                         "Root", Tissue)) %>%
  mutate(Tissue = ifelse(Experiment %in% 
                           c('SRX10073983',
                             'SRX10073981',
                             'SRX10073979',
                             'SRX10073977'),
                         "Leaf", Tissue)) 



metadata$Tissue %>% table()
metadata %>% filter(Tissue %in% 'Leaf\\, Stem\\, Root\\, Flower\\, Immature Fruit\\, Mature Fruit') %>% View


metadata$Collection_Date %>% table(exclude = F) # Some samples were collected in 2013 



# Save the metadata table
write.csv(metadata, '/Users/alejandra/pepper_counts/metadata/metadata_clean.csv')
dim(metadata) #2214 muestras




## Get the Ids of the samples that we need to process
# samples <- read.csv('muestras_descargadas.csv') %>% select('.') %>% setNames('Run')
# metadata$Run %in% samples$Run %>% table
# new_samples <- metadata %>% filter(!(Run %in% samples$Run)) %>% select(Run) 
# 
# write.csv(new_samples,'/Users/alejandra/pepper_counts/metadata/SRA_por_bajar.csv', 
#           row.names = F, quote = F,)
