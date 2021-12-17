library('dplyr')
library('purrr')
library('ggplot2')
library('stringr')
library('tidyr')
library('ggpubr')
library('fastDummies')
library('caret')


# plot evaluation from all participating laboratories
# read in spectra eval
dir<-'./Spectra_and_SpeciedID/02_quality_assessment_after_intervention/01_spectra/03_readout/'
eval_dirs <-list.files(path = dir,pattern="*read_out_sum.csv", recursive = T)
eval_list <-lapply(paste0(dir, eval_dirs), read.csv, sep=';', stringsAsFactors = FALSE, header=TRUE)

eval_list<-lapply(eval_list, function(df) mutate_at(df, .vars = c("position","frac_peaks_repr"), as.character))

eval<- eval_list %>% reduce(full_join, by = colnames(eval_list[1]))
eval<-eval[!is.na(eval$spectra),]
colnames(eval)[colnames(eval) == 'strain_number'] <- 'strainnumber'

# add lab
eval['lab'] <- gsub('\\..*', '',eval$spectra)

# check position
# add 'method'
eval['method'] <- gsub('(.*\\.)(\\d{1})(\\.\\d{2}\\.\\d{1}.*)','\\2', eval$spectra)
eval['replicate'] <- gsub('(.*\\.\\d{1}\\.\\d{2}\\.)(\\d{1})(.*)','\\2', eval$spectra)
# remove '02' tags for spectra acquired on device 44
eval$replicate <- ifelse((eval$lab== 'device\\_44' & grepl('\\_02$', eval$spectra)), '2', 
                          ifelse((eval$lab== 'device\\_44' & !grepl('\\_02$', eval$spectra)), '1', eval$replicate))


# add brukercode if available
eval['brukercode']<-ifelse(grepl('[[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12}', eval$spectra), gsub('(.*)([[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12})(.*)', '\\2', eval$spectra), NA)

# make strainnumber 2 digits
eval['strainnumber']<-ifelse(nchar(eval$strainnumber)==1, as.character(paste0('0', eval$strainnumber)), as.character(eval$strainnumber))

# for a few labs, what is assigned as 'run' is not unique and the combination of strainnnumber, position, run and lab results in duplicates. make the 'run' unique by including the method
eval$run <- ifelse(eval$lab == 'device_05', paste0(gsub('(.*\\.)(\\d)(\\.\\d{2}.*)', '\\2',eval$spectra), '20210414'), eval$run)
eval$run <- ifelse(eval$lab == 'device_11',  paste0(gsub('(.*\\.)(\\d)(\\.\\d{2}.*)', '\\2',eval$spectra), eval$run), eval$run)
# For spectra acquired in on device 15 laboratory, there are up to 6 replicates per strain, sometimes having the same label (eg. 1.01.1). The only difference are the last for digits in the spectraname, use these as 'run'
eval$run <- ifelse(eval$lab == 'device_15',  gsub('.*\\.', '', eval$spectra), eval$run)
eval$run <- ifelse(eval$lab == 'device_15',  paste0(gsub('(.*\\.)(\\d)(\\.\\d{2}.*)', '\\2',eval$spectra), eval$run), eval$run)

# for the repetition spectra from CHUV, adding the the 'method' not always makes it unique, add 'replicate' instead. As the 'method' is needed for merge, first add method, then replicate
eval$run <- ifelse(eval$lab == 'device_18',  paste0(gsub('(.*\\.)(\\d)(\\.\\d{2}\\.)(\\d{1})(\\..*)', '\\2',eval$spectra), eval$run), eval$run)
# Add the 'replicate' info only for the ones which where repeated
eval$run <- ifelse(eval$spectra %in% c('device_18.1.01.1.0_A1.1AF70', 'device_18.1.01.2.0_A1.1AF70', 'device_18.2.01.1.0_A2.1AF70', 'device_18.2.01.2.0_A2.1AF70'), paste0(gsub('(.*\\.)(\\d)(\\.\\d{2}\\.)(\\d{1})(\\..*)', '\\4',eval$spectra), '-', eval$run), eval$run)

# The same issue appeaks for spectra acquired in Groningen, 10 spectra have the same 'run' and 'position assignments. These are probably 'Re-measurements'
eval$run <- ifelse(eval$spectra %in% c("device_13.1.42.1.0_H1.c5fc25d2-8481-4892-9c4e-dbb5b60b2bc7","device_13.1.42.2.0_H2.bd991752-c356-4274-a583-7cd8aa1aef8d",
                                       "device_13.1.45.1.0_H7.cba9a00e-2e68-44bd-9144-63b5ea6c8e8e","device_13.1.45.2.0_H8.d81e4cb5-a0e9-409e-9476-87bcf737bf89",
                                       "device_13.1.46.1.0_H9.d48cff06-972f-456b-9e49-651b03661648","device_13.1.46.2.0_H10.e57bf5ee-7b66-4e95-b8b8-e10020e0e82f",
                                       "device_13.2.19.1.0_D3.328e83f7-75fd-43cd-aeb6-c52c8874d7a2","device_13.2.20.2.0_D6.f4a75feb-b23c-47f1-a300-2e7b398ff1ca",
                                       "device_13.2.45.1.0_H7.75f3fc2b-36d0-48ce-9ee9-5b178e04863d","device_13.2.45.2.0_H8.ff86ec78-7ca6-4996-8b0d-dc2c69a3dcb2"), 
                   paste0(eval$run, '-', gsub('(.*\\.)(\\d)(\\.\\d{2}\\.)(\\d{1})(\\..*)', '\\4',eval$spectra)), eval$run)

# For spectra acquired on device_19, remove the 'WASLAB' and 'AcidFormic' tag (not part of the raw 'fid' file and therefor not in the brukerfile)
eval$run <- ifelse(eval$lab == 'device_19', gsub('MALDI 47 WASPLAB ', '',eval$run), eval$run)
eval$run <- ifelse(eval$lab == 'device_19', gsub('\\-AcidFormic', '',eval$run), eval$run)

# check which have a non-unique combination of position, strainnumber and 'run' and add brukercode or run (made unique before) to make unique
eval['to_merge']<-ifelse((!is.na(eval$brukercode) & (!eval$lab %in% c('device_03', 'device_13', 'device_20', 'device_23', 'device_24', 'device_25', 'device_41'))), eval$brukercode, 
                                paste(eval$lab, eval$strainnumber, eval$position, eval$Strain, eval$run))

# check whcih are duplicated
any(duplicated(eval$to_merge)) # now, all have a unique identifier
# dupli_eval<-eval[duplicated(eval$to_merge),] 

# read in brukerreports
dir<-'./Spectra_and_SpeciedID/02_quality_assessment_after_intervention/02_species_ID/microflexBiotyper/02_summarised/'
brukerreport_dirs <- list.files(path = dir,pattern="*.csv")
brukerreport_list <- lapply(paste0(dir, brukerreport_dirs), read.csv, sep=';', stringsAsFactors = FALSE, header=TRUE)

# remove 'spectra' and 'samplename' column from brukerreport and work withe the unique spectraname from the asciis 
for (i in 1:length(brukerreport_dirs)){
  brukerreport_list[[i]]['lab']<-gsub('_brukerreport.*.csv', '', brukerreport_dirs[[i]])
  brukerreport_list[[i]]['samplename']<-NULL
  brukerreport_list[[i]]['spectra']<-NULL
  brukerreport_list[[i]]['run']<-as.character(brukerreport_list[[i]]$run)
}

# convert to character, so that they can be merged
brukerreport_list<-lapply(brukerreport_list, function(df) mutate_at(df, .vars = c("strainnumber",  "col","Strain","species_NGS","genus_NGS"), as.character))

# define which columns are used for merging (these are the ones which are present in all)
cols_merge<-c("lab","strainnumber","col","run","Rank" ,"Score","Species","species_identified","genus_identified", "position",      
              "n_species","n_genera","n_species_over_2","n_genera_over_2","diff","Strain","species_NGS","genus_NGS","correct_genus_identified",
              "correct_species_identified","include")
# summarise all brukerreports to one df
brukerreport_all<-brukerreport_list %>% reduce(full_join, by = cols_merge)
# subset to selected columns and remove empty rows
brukerreport_all<-brukerreport_all[!is.na(brukerreport_all$Rank),c("strainnumber","col","run","position","Rank","Score","Species","lab","species_identified","genus_identified","n_species","n_genera","n_species_over_2","n_genera_over_2","diff","Strain","species_NGS","genus_NGS","correct_genus_identified","correct_species_identified","include", "Method")]

# As in the 'eval' file, for some laboratories the samplename is not unique yet. Make these unique and that the 
# spectra acquired on vitekMS / Shimadzu devices do not need to be made unique, as these can me merged to the VitekMS reports using the sample name
# in 'device_05' the run is not unique, add therefor method
brukerreport_all$run<-ifelse(brukerreport_all$lab == 'device_05', paste0(gsub('(\\d)(\\.\\d{2}.*)', '\\1',brukerreport_all$col), brukerreport_all$run), brukerreport_all$run)

# for the repetition spectra on device_18, adding the the 'method' not always makes it unique, add 'replicate' instead
brukerreport_all$run <- ifelse(brukerreport_all$lab == 'device_18',  paste0(brukerreport_all$Method, 'Souches round 2'), brukerreport_all$run)

# harmoinse 'run' from device_14 and device_45 to parent dir before merge
brukerreport_all$run <- ifelse(brukerreport_all$lab == 'device_14',  paste0('MALDI-TOF-MS-EQA-2_0',gsub('(\\d)(\\.\\d{2}\\.)(\\d{1})', '\\1',brukerreport_all$col)), brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$lab == 'device_45' & grepl('TARGET 1', brukerreport_all$run), 'MALDI-TOF-MS-EQA-2 TARGET1', 
                               ifelse(brukerreport_all$lab == 'device_45' & grepl('TARGET 2', brukerreport_all$run), 'MALDI-TOF-MS-EQA-2 TARGET2', brukerreport_all$run))
#for spectra from device_19, some spectra were put into the directory of another run (in ascii and eval, this will be read), while the fid file which is read in the brukerreports indicates another run. correct for this
brukerreport_all$run <- ifelse(brukerreport_all$lab == 'device_19' & brukerreport_all$col == '1.07.2-0-001', '210414-0926-1011001995', 
                               ifelse(brukerreport_all$lab == 'device_19' & brukerreport_all$col == '2.11.1-0-001', '210414-1111-101100688', 
                                      ifelse(brukerreport_all$lab == 'device_19' & brukerreport_all$col %in% c('2.18.1-0-002','2.18.2-0-003'), '210414-1111-101100688', 
                                             ifelse(brukerreport_all$lab == 'device_19' & brukerreport_all$col %in% c('2.20.1-0-004','2.20.2-0-005'), '210414-1111-101100688',  
                                                    ifelse(brukerreport_all$lab == 'device_19' & brukerreport_all$col %in% c('1.33.2-0-002','1.43.2-0-003'), '210414-0926-1011001995',
                                                           ifelse(brukerreport_all$lab == 'device_19' & brukerreport_all$col %in% c('1.45.1-0-004','1.45.2-0-005'), '210414-0926-1011001995',
                                                                  ifelse(brukerreport_all$lab == 'device_19' & brukerreport_all$col %in% c('2.46.1-0-006','2.46.2-0-007'), '210414-1111-101100688',brukerreport_all$run)))))))
# in device_27 remove double spacing read in brukerreport
brukerreport_all$run <- ifelse(brukerreport_all$lab == 'device_27',  gsub('\\s+', ' ',brukerreport_all$run), brukerreport_all$run)
# harmonise lab naming for device_28 spectra
brukerreport_all$lab <- gsub('device_28_1', 'device_28', brukerreport_all$lab)

# check for which the 'include' column is NA
check<-brukerreport_all[is.na(brukerreport_all$include),] # these are either calibration spectra or empty spectra, its ok to exclude them

# these USB spectra have been acquired unrelated to the EQA, exclude these from further analysis
brukerreport_all$include <- ifelse(grepl('device_17', brukerreport_all$lab) & brukerreport_all$col %in% c("1c80ff26-7888-4463-aaa6-f2b998b76bf2", "13255808-435f-427d-a480-5ea1ce89f6e1", "b9ec6d36-4fe5-4be9-8f0a-69b540f969d6", "f84d83cd-33f3-4828-916f-4feac57f0cdf","faf63369-0e48-46ed-9aa3-00dd01559501", "16e3acb4-af8e-4145-9d15-2c4571e69561", "37ca766a-6e27-49cc-8973-f4eba220cede","22db34b0-2d88-47f8-a14c-01ba78969ea9"), FALSE, brukerreport_all$include) 

# In device_13, not all are unique , as some strains were 'remeasured' (same run and position). check manually which spectrum belongs to which bruker identification and make the 'run' unique.
# For spectra with the same run and position, the spectra with more peaks was assigned to the brukeridentification with the higher score. 
brukerreport_all$run <- ifelse(brukerreport_all$col == '2.19.1' & brukerreport_all$run == '210413-1459-1011028068'& brukerreport_all$position == '0_D3' & brukerreport_all$Score == '1.62', '210413-1459-1011028068-1', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '2.20.2' & brukerreport_all$run == '210413-1459-1011028068' & brukerreport_all$position == '0_D6' & brukerreport_all$Species == 'no peaks found', '210413-1459-1011028068-2', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '1.42.1' & brukerreport_all$run == '210413-1448-1011027967'& brukerreport_all$position == '0_H1'  & brukerreport_all$Score == '1.73', '210413-1448-1011027967-1', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '1.42.2' & brukerreport_all$run == '210413-1448-1011027967'& brukerreport_all$position == '0_H2'  & brukerreport_all$Score == '1.36' & brukerreport_all$lab == 'device_13_2', '210413-1448-1011027967-2', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '1.45.1' & brukerreport_all$run == '210413-1448-1011027967'& brukerreport_all$position == '0_H7'  & brukerreport_all$Score == '1.66', '210413-1448-1011027967-1', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '1.45.2' & brukerreport_all$run == '210413-1448-1011027967'& brukerreport_all$position == '0_H8'  & brukerreport_all$Score == '1.36', '210413-1448-1011027967-2', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '2.45.1' & brukerreport_all$run == '210413-1459-1011028068'& brukerreport_all$position == '0_H7'  & brukerreport_all$Score == '1.58', '210413-1459-1011028068-1', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '2.45.2' & brukerreport_all$run == '210413-1459-1011028068'& brukerreport_all$position == '0_H8'  & brukerreport_all$Score == '1.57', '210413-1459-1011028068-2', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '1.46.1' & brukerreport_all$run == '210413-1448-1011027967'& brukerreport_all$position == '0_H9'  & brukerreport_all$Score == '1.38', '210413-1448-1011027967-1', brukerreport_all$run)
brukerreport_all$run <- ifelse(brukerreport_all$col == '1.46.2' & brukerreport_all$run == '210413-1448-1011027967'& brukerreport_all$position == '0_H10'  & brukerreport_all$Score == '1.78', '210413-1448-1011027967-2', brukerreport_all$run)

# remove the spectra which were identified as wrong genus, these are counted as 'contamination'
brukerreport_incl<-brukerreport_all[brukerreport_all$include == TRUE,]
# remove duplicated spectra
brukerreport_incl<-brukerreport_incl[!duplicated(brukerreport_incl[,]),]

# ensure that strainnumber is two digit
brukerreport_incl['strainnumber']<-ifelse(nchar(brukerreport_incl$strainnumber)==1, as.character(paste0('0', brukerreport_incl$strainnumber)), as.character(brukerreport_incl$strainnumber))
# only assign 'brukercode' if it comes in the right format
brukerreport_incl['brukercode']<-ifelse(grepl('[[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12}', brukerreport_incl$col), gsub('(.*)([[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12})(.*)', '\\2', brukerreport_incl$col), NA)

# set brukerreport NA for Prague and Wales as it is not included in ascii and therefor in the eval file (cannot be used for merging)
brukerreport_incl['brukercode']<-ifelse(brukerreport_incl$lab == 'device_32', NA, as.character(brukerreport_incl$brukercode))
# spectra from lab 'BD' are duplicated, once as M2/MAldiAlt and M3/MaldiNeu --> rename to the same, then remove duplicates


#brukerreport_incl$lab <- gsub('040_University_Hospital_FreiburgiB_M2', '040_University_Hospital_FreiburgiB_MaldiAlt', brukerreport_incl$lab, )
#brukerreport_incl$lab <- gsub('040_University_Hospital_FreiburgiB_M3', '040_University_Hospital_FreiburgiB_MaldiNeu', brukerreport_incl$lab, )

# Remove '_2' tag from Groninge and HUG lab
brukerreport_incl$lab <- gsub('device_13_2', 'device_13', brukerreport_incl$lab)
brukerreport_incl$lab <- gsub('device_19_2', 'device_19', brukerreport_incl$lab)

# remove duplicates
brukerreport_incl <- brukerreport_incl[!duplicated(brukerreport_incl),]

# build column 'to_merge'. this should now be unique
brukerreport_incl['to_merge']<- ifelse(!is.na(brukerreport_incl$brukercode), brukerreport_incl$brukercode, paste(brukerreport_incl$lab, brukerreport_incl$strainnumber, brukerreport_incl$position, brukerreport_incl$Strain, brukerreport_incl$run))

# check merge
check <- brukerreport_incl[brukerreport_incl$to_merge %in% setdiff(brukerreport_incl$to_merge, eval$to_merge),] # these are empty spectra. They were excluded when the asciis where exported, but not from the bruker database species identification. These can be excluded
check <- check[!check$Species %in% c(NA, 'no peaks found'),]

# There can be spectra included, where the ascii is not empty, and still the bruker database assignment says 'no peaks found'. This is possible as different cut-offs for bruker peak picking may have been applied in the flexAnalysis Software connected to the MBT species identification database and which were used to export the asciis. 
# check which duplicated and add 'run' variable to make unique, now none are duplicate
dupli_bruker<-brukerreport_incl[duplicated(brukerreport_incl$to_merge),] # none duplicated anymore
#dupli_bruker<-dupli_bruker[!is.na(dupli_bruker$strainnumber),]
#dupli_bruker_all<-brukerreport_incl[brukerreport_incl$to_merge %in% brukerreport_incl$to_merge[duplicated(brukerreport_incl$to_merge)],]

length(intersect(eval$to_merge, brukerreport_incl$to_merge))

brukerreport_incl<-brukerreport_incl[!is.na(brukerreport_incl$Species),]
setdiff(brukerreport_incl[brukerreport_incl$Species !='no peaks found','to_merge'],eval$to_merge) # all from brukerreports are in eval

# check whether there are still discrepancies before merging
check<-eval[eval$to_merge %in% setdiff(eval$to_merge, brukerreport_incl$to_merge), ]
# do not consider vitek spectra here. there are no discrepancies left
check<-check[!(check$lab %in% c('device_12','device_11','device_35','device_36', 'device_38', 'device_44', 'device_43', 'device_30', 'device_15', 'device_09')),]
# There are no discrepancies left

# remove these columns before merging
brukerreport_incl$brukercode<-NULL
brukerreport_incl$Strain<-NULL
brukerreport_incl$lab<-NULL
brukerreport_incl$position<-NULL
brukerreport_incl$run<-NULL
brukerreport_incl$strainnumber<-NULL

# add 'brukerDB' tag to all brukerID specific colunns efore merging
colnames(brukerreport_incl)[colnames(brukerreport_incl) %in% c("col", "Rank", "Score", "Species", "species_identified", "genus_identified", "n_species", "n_genera", "n_species_over_2", "n_genera_over_2", "diff", "correct_genus_identified", "correct_species_identified")] <- paste0('brukerDB.', c("col", "Rank", "Score", "Species", "species_identified", "genus_identified", "n_species", "n_genera", "n_species_over_2", "n_genera_over_2", "diff", "correct_genus_identified", "correct_species_identified"))

# merge brukerreport and eval
eval2 <- merge(eval, brukerreport_incl, by = 'to_merge', all.x = T)

no_bruker<-eval2[is.na(eval2$brukerDB.Species),]
no_bruker<-no_bruker[!no_bruker$lab %in% c("device_09","device_43","device_07","device_11","device_12","device_15","device_30",	
                                           "device_35","device_44"),]
# now, all have a bruker id
eval2$to_merge<-NULL

# add group information 
groups<-read.csv2('/Users/aline/ESCMID/ESPRIT/04_Strains/Strains_groups_assigned_2021_02.csv', sep=',')
groups['Numbering_Shipment_1']<-ifelse(nchar(groups$Numbering_Shipment_1)==1, as.character(paste0('0', groups$Numbering_Shipment_1)), as.character(groups$Numbering_Shipment_1))

# mere by strainnumber and strain
eval<-merge(eval2, groups, by.x = c('Strain','strainnumber'), by.y = c('Strain','Numbering_Shipment_1'), all.x = TRUE)
eval['spectra']<-gsub('\\#PSD\\=', '', eval$spectra)

# set all below threshold (1.7) to 'no Identification possible)
eval$brukerDB.species_identified<-ifelse(eval$brukerDB.Score < 1.7, 'no ID possible', eval$brukerDB.species_identified)
eval$brukerDB.Species<-ifelse(eval$brukerDB.Score < 1.7, 'no ID possible', eval$brukerDB.Species)

# add vitekMS id
# read VitekMS
VitekMSreport_dirs = list.files(path ='./Spectra_and_SpeciedID/02_quality_assessment_after_intervention/02_species_ID/VitekMS/02_summarised/',pattern="*_VitekMSreport.csv")
VitekMSreport_list = lapply(paste0('./Spectra_and_SpeciedID/02_quality_assessment_after_intervention/02_species_ID/VitekMS/02_summarised/',VitekMSreport_dirs), read.csv, sep=';', stringsAsFactors = FALSE, header=TRUE)

# transform to character before merging
VitekMSreport_list<-lapply(VitekMSreport_list, function(df) mutate_at(df, .vars = c("strainnumber", "Strain","species_NGS","genus_NGS"), as.character))

# define which columns are needed to merge
cols_merge<-c("X","strainnumber","files","position","nb_peaks","max_SN","max_intensities","max_resolution","identifType","species_rank_1","proba_rank_1","score_rank_1","lab","Strain","n_species","diff","diff_proba","species_identified","genus_identified","species_NGS","genus_NGS","correct_genus_identified","correct_species_identified","include")
# summarise all Vitek Reports to one dataframe
VitekMSreport_all<-VitekMSreport_list %>% reduce(full_join, by = cols_merge)

# add column 'samplename_new' to merge to eval
VitekMSreport_all['vitekMSDB.samplename_new'] <- gsub('.txt', '', VitekMSreport_all$files)

# Some spectra where not assigned a species / genus. These are all of 'identifType': NoIdentification / NotEnoughPeaks / TooManyPeaks  / TooNoisy
# Only spectra which have been identified as a wrong genus with high confidence are excluded as regarded as contamination, spectra with one of the above mentioned 'identifType' are not excluded
# Set 'include' therefor to 'TRUE'
VitekMSreport_all['include']<-ifelse(is.na(VitekMSreport_all$include), TRUE, as.character(VitekMSreport_all$include))
# remove spectra where 'included' == 'FALSE, these are regarded as contamination / mixup
VitekMSreport_incl<-VitekMSreport_all[VitekMSreport_all$include!=FALSE,]
# remove first columns (rownames)
VitekMSreport_incl<-VitekMSreport_incl[,2:ncol(VitekMSreport_incl)]
# make the strainnumber two digits 
VitekMSreport_incl$strainnumber<-ifelse(grepl('^[[:digit:]]$', VitekMSreport_incl$strainnumber), paste0('0', VitekMSreport_incl$strainnumber), VitekMSreport_incl$strainnumber)

# add vitekMS tag to columns realated to vitekMS identification vefore merging
colnames(VitekMSreport_incl)[colnames(VitekMSreport_incl) %in% c("files", "position", "nb_peaks", "max_SN", "max_intensities", "max_resolution", "identifType", "species_rank_1", "proba_rank_1", "score_rank_1", "n_species", "diff", "diff_proba", "species_identified", "genus_identified", "correct_genus_identified", "correct_species_identified", "species_NGS","genus_NGS","include")] <- paste0('vitekMSDB.', c("files", "position", "nb_peaks", "max_SN", "max_intensities", "max_resolution", "identifType", "species_rank_1", "proba_rank_1", "score_rank_1", "n_species", "diff", "diff_proba", "species_identified", "genus_identified", "species_NGS", "genus_NGS", "correct_genus_identified", "correct_species_identified", "include"))


# set to 'no ID possible' if identifType where 'NoIdentification', 'NotEnoughPeaks', 'TooNoisy'
VitekMSreport_incl$vitekMSDB.species_identified<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), 'no ID possible', VitekMSreport_incl$vitekMSDB.species_identified)
VitekMSreport_incl$vitekMSDB.species_rank_1<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), 'no ID possible', VitekMSreport_incl$vitekMSDB.species_rank_1)
VitekMSreport_incl$vitekMSDB.proba_rank_1<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), NA, VitekMSreport_incl$vitekMSDB.proba_rank_1)
VitekMSreport_incl$vitekMSDB.score_rank_1<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), NA, VitekMSreport_incl$vitekMSDB.score_rank_1)

# For the LTW spectra, take out the '_02'-tag in the 'eval file before merging
eval$spectra <- ifelse(grepl('device\\_44', eval$spectra), gsub('\\_02$', '', eval$spectra), eval$spectra)

# check if there are spectra for which there are vitek ID, but which are not included in 'eval' file
setdiff(VitekMSreport_incl$vitekMSDB.samplename_new, eval$spectra) # there are none

# merge Vitek reports to eval file
eval<-merge(eval, VitekMSreport_incl, by.x = c('spectra', "Strain","strainnumber","lab"), by.y = c("vitekMSDB.samplename_new","Strain","strainnumber","lab"), all.x = T)
eval$to_merge<-NULL
eval['species_NGS']<-ifelse(!is.na(eval$species_NGS), eval$species_NGS, as.character(eval$vitekMSDB.species_NGS))
eval['genus_NGS']<-ifelse(!is.na(eval$genus_NGS), eval$genus_NGS, as.character(eval$vitekMSDB.genus_NGS))

# set 'include' to false for all which have been identified as another genus with high confidence by the VitekMS Database
eval['include']<-ifelse(!is.na(eval$include), eval$include, as.character(eval$vitekMSDB.include))
eval['include']<-ifelse(eval$spectra %in% setdiff(VitekMSreport_all$vitekMSDB.samplename_new, VitekMSreport_incl$vitekMSDB.samplename_new), FALSE, as.character(eval$include))

# Remove the rows which where 'include' is not 'TRUE'
eval <- eval[eval$include == TRUE,]

# remove these redundant columns
eval$vitekMSDB.species_NGS<-NULL
eval$vitekMSDB.genus_NGS<-NULL
eval$vitekMSDB.include<-NULL

# Add marker based ID
markerID<-read.csv2('./Spectra_and_SpeciedID/02_quality_assessment_after_intervention/02_species_ID/marker_ID/02_summarised/summary_marker_based.csv', sep=',')

# harmonise sample naming
markerID['sample']<-gsub('\\#PSD.*', '', markerID$sample)
markerID['sample']<-gsub(' $', '', markerID$sample)

# make 'strainnumber' two digit
markerID['strainnumber']<-ifelse(grepl('^[[:digit:]]$',markerID$strainnumber), paste0('0', markerID$strainnumber), as.character(markerID$strainnumber))

# remove '02' tag from LTW spectra from the marker based ID file
markerID$sample <- ifelse(grepl('device\\_44', markerID$sample), gsub('\\_02$', '', markerID$sample),markerID$sample)

# harmonise sample naming
eval['spectra']<-gsub('\\#*PSD.*$', '', eval$spectra)
eval['spectra']<-gsub('\\#', '', eval$spectra)
eval['spectra']<-gsub(' $', '', eval$spectra)

# extract which classifier ID. For these specta, it is possible that they are not included in the identification output (although they have been in the input). This is possible as spectra are automatically excluded from the 'classifier' if they are not well calibrated or if they do not match certain quality criteria
classifier_strains <- unique(markerID[!is.na(markerID$profiles_matched), 'strainnumber'])

# check which are missing from markerID. No marker ID, means that these have been excluded for quality reasons (no calibration) 
no_marker<-eval[eval$spectra %in% setdiff(eval$spectra,markerID$sample),]

# all which were identified using the PAPMID DB should be in here, check
no_marker_papmid <- no_marker[!(no_marker$strainnumber %in% classifier_strains),] # none are missing

# check how many were excluded per device
table(no_marker$lab) # all within range 

# check, whether there are some which have a marker ID, but are not in eval 
missing<-markerID[markerID$sample %in% setdiff(markerID$sample,eval$spectra), ] # I checked all of these manually, these where excluded because they where identified as the wrong genus using the bruker /vitek database

# remove these columns, as they are not evaluated
markerID$X.Subunit.Mass.<-NULL
markerID$profiles_matched<-NULL
markerID$Run<-NULL
colnames(markerID)[colnames(markerID) == 'Genus'] <-"genus_identified"

#remove dupliactes if present
markerID<-markerID[!duplicated(markerID),]

# add 'markerID' tag to columns which are associated to marker based ID before merging
colnames(markerID)[colnames(markerID) %in% c("match_count", "species_identified", "genus_identified", "n_species_identified", "n_genera_identified", "correct_genus_identified", "correct_species_identified")] <- paste0('markerDB.', c( "match_count", "species_identified", "genus_identified", "n_species_identified", "n_genera_identified", "correct_genus_identified", "correct_species_identified"))

# merge marker ID to eval
eval<-merge(eval, markerID, by.x = c('spectra',"Strain","strainnumber","species_NGS","genus_NGS"), by.y = c('sample',"Strain","strainnumber","species_NGS","genus_NGS"), all.x = T)
# Add column 'MALDI.device_full' to make with results from first round compatible
eval['MALDI.device_full']<-eval$lab
eval['MALDI.device']<-eval$lab

# in the column 'lab' all spectra acquired in the same laboratory have the same 'lab' value
labs<- read.csv('./Spectra_and_SpeciedID/DeviceID_lab.csv')
eval$lab <- NULL
eval<- merge(eval, labs, by.x = 'MALDI.device_full', by.y = 'MALDI.device', all.x = T)

# remove 'newline' character from the bruker species column before exporting
eval$brukerDB.Species<- sapply(eval$brukerDB.Species,function(x) { gsub("[\r\n]", "", x)})

# add randomly assigned 'device ID' per MALDI-TOF MS dataset
eval['device_id'] <- gsub('device\\_','',eval$MALDI.device_full)

# ensure that this 'device_ID' has two digits
eval$device_id <- ifelse(nchar(eval$device_id) == 1, paste0('0', eval$device_id),eval$device_id)

# remove empty lines
eval$brukerDB.Species<- sapply(eval$brukerDB.Species,function(x) { gsub("[\r\n]", "", x)})
eval <- eval[!is.na(eval$spectra),] # these are all all NA 

# export these
write.table(eval, './Spectra_and_SpeciedID/02_quality_assessment_after_intervention/eval_participating_labs_2.txt', row.names = F, quote = F, dec = '.', sep = ';', eol = "\n")

