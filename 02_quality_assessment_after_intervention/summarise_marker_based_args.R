# Aline Cuénod, 2021
# This script summarises the marker based identifications, which have been performed per group
# The following arguments are required: 
# (i) path to the common dir of the different output files (e.g. './Species_identification/PAPMID/')
# (ii) path to the file 'Strains_numbering_first_Shipment_2021-02.csv'
# (iii) output path to the summary csv file

# Load packages
library('tidyverse')
library('dplyr')

# Define arguments
args = commandArgs(trailingOnly=TRUE)

# Import all PAPMID output files
files = list.files(path = args[1],        # directory to search within
                   pattern = ".*export.csv", # regex pattern, some explanation below
                   recursive = T,          # search subdirectories
                   full.names = TRUE          # return the full path
)

PAPMID_export = lapply(files, read.csv)  # read all the matching files

# Convert all columns to character
for (i in 1:length(PAPMID_export)){
  PAPMID_export[[i]][]<-lapply(PAPMID_export[[i]], as.character)
  PAPMID_export[[i]]['file']<-files[i]
}


# Merge all PAPMID outputs
PAPMID_export <- PAPMID_export %>% reduce(full_join, by = colnames(PAPMID_export[[1]]))

# Keep only 1 row per spectra, filter for the one with the highest MatchCount Score
PAPMID_export_filtered<-as.data.frame(PAPMID_export %>% 
                                        group_by(Sample) %>% filter(DataCount == max(as.numeric(as.character(DataCount)))) %>% ungroup())


# Rename Winkia neuii -> Actinomyces neii in order to compare
PAPMID_export_filtered['Genus']<-ifelse(PAPMID_export_filtered$Genus == 'Winkia' & grepl('neuii', PAPMID_export_filtered$Species), 'Actinomyces', PAPMID_export_filtered$Genus)

# In order to evaluate the species identification result, create 'species' column, not including subspecies info
PAPMID_export_filtered['species']<-gsub('\\-*subsp.*', '',PAPMID_export_filtered$Species)


# Count number of species and genus with highest match, add which species have been identified per sample and keep one entry per spectrum
PAPMID_export_filtered <- PAPMID_export_filtered %>%
  group_by(Sample) %>%
  mutate(n_species_identified = n_distinct(paste(Genus, species, sep=' '))) %>% 
  mutate(n_genera_identified = n_distinct(Genus)) %>%
  mutate(species_identified = paste(unique(paste(Genus, species, sep=' ')), collapse='-')) %>%
  mutate(genus_identified = paste(unique(Genus), collapse='-')) %>%
  filter(row_number() == 1)

# Remove empty or redundant columns columns
PAPMID_export_filtered <- PAPMID_export_filtered[!PAPMID_export_filtered$Sample=='',c("Sample","Run","species_identified","genus_identified","DataCount","X.Subunit.Mass.","n_species_identified","n_genera_identified", "file")]

# Some 'qnt' files have been analysed with PAPMID and classifer, keep only classifier results 
PAPMID_export_filtered<- PAPMID_export_filtered[!grepl('BRU-(17|18|39|40)-1-125', PAPMID_export_filtered$Sample),]

# remove the axima ones in the microflex output and vice versa
PAPMID_export_filtered['drop'] <- if(any(grepl('01_Spectratest_Microflex', args[1]))){
  ifelse(grepl('BRU', PAPMID_export_filtered$Sample), 'keep', 'drop')
}else{
  ifelse(grepl('BRU', PAPMID_export_filtered$Sample), 'drop', 'keep')
} 

PAPMID_export_filtered<-PAPMID_export_filtered[PAPMID_export_filtered$drop == 'keep',]
PAPMID_export_filtered$drop <- NULL

# Import all classifier outputs. (Same approach, but subtyping module)
# List all files
files = list.files(path = args[1], 
                   pattern = ".*Identification Report.*csv$", # regex pattern, some explanation below
                   recursive = TRUE,          # search subdirectories
                   full.names = TRUE          # return the full path
)

# Define delimiter
delimiter <- ','

# Read all files
output = lapply(files, read.csv, sep = delimiter)  # read all the matching files

# Convert all columns to chracter and add genus
for (i in 1:length(output)){
  colnames(output[[i]])<-c("X","spectra","sample","match_count","species_identified","profiles_matched")
  output[[i]]['Genus']<-ifelse(strsplit(files[i], '\\/')[[1]][13] != "", strsplit(files[i], '\\/')[[1]][13], strsplit(files[i], '\\/')[[1]][14])
  output[[i]][]<-lapply(output[[i]], as.character)
  output[[i]]['file']<-files[i]
}


# Merge all classifier outputs
output <- output %>% reduce(full_join, by = colnames(output[[1]]))
output$X<-NULL
output$spectra<-NULL

# some 'quantity' spectra have wrongly been analysed with the strep viridans classifier, although they belong to other genera, remove these
output<-output[!(output$Genus == 'S_viridans' & grepl('qty', output$sample)),]

# remove the axima ones in the microflex output and vice versa
output['drop'] <- if(any(grepl('01_Spectratest_Microflex', args[1]))){
  ifelse(grepl('BRU', output$sample), 'keep', 'drop')
}else{
  ifelse(grepl('BRU', output$sample), 'drop', 'keep')
} 

output<-output[output$drop == 'keep',]
output$drop <- NULL

# for the sample where the ID was repeated, keep the newer one
#Find all duplicates
output<-output%>% 
  group_by(sample) %>%
  mutate(repeated = ifelse(n_distinct(file) > 1, TRUE, FALSE))

#Subset to relevant
output['drop'] <- ifelse(output$repeated == TRUE & !grepl('20210219|20210224|Odense', output$file), 1, 0)
output <- output[output$drop == 0,]
output$drop<-NULL
output$repeated <- NULL

# Add 'Genus' column
output['Genus']<-ifelse(output$Genus == 'S_aureus_complex', 'Staphylococcus', 
                        ifelse(output$Genus == 'S_viridans', 'Streptococcus', output$Genus))


#count species identified
output['species_identified_count']<-gsub('-like', '.like', output$species_identified)
# xx-like: two different species
output['species_identified_count']<-gsub('(like\\-)(\\d{1,2})', 'like.\\2', output$species_identified_count)
# xx-cluster1: distinguishing within two clusters of the same species
output['species_identified_count']<-gsub('-cluster\\d{1,2}', '', output$species_identified_count)
output['species_identified_count']<-gsub('-subsp-[[:alnum:]]+', '', output$species_identified_count)
output['species_identified_count']<-gsub(' or ', '-', output$species_identified_count)

n_species_identified<-list()
for (i in 1:nrow(output)){
  n_species_identified<-append(n_species_identified, length(unique(strsplit(output$species_identified_count[[i]], "-")[[1]])))
}

output['n_species_identified']<-unlist(n_species_identified)
output$species_identified_count<-NULL
output['n_genera_identified']<-1

# Merge PAPMID and Classifier Outputs
output<-merge(output, PAPMID_export_filtered, by.x=c('sample','file','match_count', 'species_identified', 'Genus','n_species_identified','n_genera_identified'), by.y = c('Sample', 'file','DataCount', 'species_identified', 'genus_identified','n_species_identified','n_genera_identified'), all=TRUE)

# Add strainnumber
output['strainnumber'] <- gsub('(.*\\.\\d\\.)(\\d{2})(\\.\\d\\..*)', '\\2', output$sample)
output['strainnumber']<-ifelse(nchar(output$strainnumber)==1, as.character(paste0('0', output$strainnumber)), as.character(output$strainnumber))

# Remove entries with weird strainnumber
output<-output[!nchar(output$strainnumber)>2,]

# Remove 'contaminated'
output<-output[!grepl('conta', output$sample),]

# change species name from bronchoseptica to bronchiseptica
output$species_identified<-gsub('bronchoseptica', 'bronchiseptica', output$species_identified)

# In order to evaluate whether the identification is correct, add specied by NGS
numbering<-read.csv2(args[2], sep=',')
numbering['species_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\3', numbering$Strain)
numbering['genus_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\1', numbering$Strain)
numbering['genus_NGS']<-gsub('Winkia', 'Actinomyces', numbering$genus_NGS)
numbering['Numbering_Shipment_1']<-ifelse(nchar(numbering$Numbering_Shipment_1)==1, paste0('0',numbering$Numbering_Shipment_1), numbering$Numbering_Shipment_1)

# Add strainnumber
output<-merge(output, numbering, by.x = 'strainnumber', by.y = 'Numbering_Shipment_1')


# Check if correct genus amongst identified genera
output['correct_genus_identified']<-mapply(grepl, output$genus_NGS, output$Genus)
output['correct_species_identified']<-mapply(grepl, output$species_NGS, output$species_identified)

# some spectra were first analysed with the module of the wrong phylogenetic group, remove these
output['in_correct_input_dir']<-ifelse(output$strainnumber %in% c('01', '02', '03', '04', '05', '06') & grepl('Klebsiella', output$file), TRUE,
                                        ifelse(output$strainnumber %in% c('07', '08') & grepl('Listeria', output$file), TRUE,
                                               ifelse(output$strainnumber %in% c('09', '10', '11', '12', '13', '14', '15', '16') & grepl('Shigella', output$file), TRUE,
                                                      ifelse(output$strainnumber %in% c('17', '18', '19', '20') & grepl('Burkholderia',output$file), TRUE,
                                                             ifelse(output$strainnumber %in% c('21', '22', '23') & grepl('Bordetella', output$file), TRUE,
                                                                    ifelse(output$strainnumber %in% c('24', '25', '29') & grepl('viridans', output$file), TRUE,
                                                                           ifelse(output$strainnumber %in% c('26', '27', '28', '30', '31', '32') & grepl('others', output$file), TRUE,
                                                                                  ifelse(output$strainnumber %in% c('33', '34') & grepl('Bacterioides', output$file), TRUE, 
                                                                                         ifelse(output$strainnumber %in% c('35', '36', '37', '38') & grepl('Enterobacter', output$file), TRUE, 
                                                                                                ifelse(output$strainnumber %in% c('39', '40', '41') & grepl('aureus', output$file), TRUE, 
                                                                                                       ifelse(output$strainnumber %in% c('42', '43') & grepl('Coryne', output$file), TRUE, 
                                                                                                              ifelse(output$strainnumber == '44' & grepl('Gardnarella', output$file), TRUE, 
                                                                                                                     ifelse(output$strainnumber %in% c('45', '46') & grepl('Actino', output$file), TRUE, 
                                                                                                                            ifelse(output$strainnumber == '47' & grepl('Pasteurella', output$file), TRUE, FALSE))))))))))))))

# remove the ones which were in wrong subtyping module
output<-output[output$in_correct_input_dir == TRUE,]
output$in_correct_input_dir<-NULL
output$file<-NULL

write.csv(output, args[3], row.names = F)
