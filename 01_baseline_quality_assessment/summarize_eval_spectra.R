library('dplyr')
library('purrr')
library('ggplot2')
library('stringr')
library('tidyr')
library('ggpubr')
library('fastDummies')
library('caret')

# read in spectra eval
dir<-'./Spectra_and_SpeciedID/01_baseline_quality_assessment/01_spectra/03_readout/'
eval_dirs <-list.files(path = dir,pattern="*read_out_sum.csv")
eval_list <-lapply(paste0(dir, eval_dirs), read.csv, sep=';', stringsAsFactors = FALSE, header=TRUE)

# add 'lab' column
for (i in 1:length(eval_dirs)){
  eval_list[[i]]['lab']<-gsub('_read_out_sum.csv','',eval_dirs[[i]])
}

# convert these columns to character
eval_list<-lapply(eval_list, function(df) mutate_at(df, .vars = c("position","frac_peaks_repr"), as.character))
# join all 'eval' files to one dataframe
eval<- eval_list %>% reduce(full_join, by = colnames(eval_list[1]))
# remove NA rows, if there are any
eval<-eval[!is.na(eval$spectra),]
# rename strainnumber column
colnames(eval)[colnames(eval) == 'strain_number'] <- 'strainnumber'

# add 'brukercode' column if it is present in the right format
eval['brukercode']<-ifelse(grepl('[[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12}', eval$spectra), gsub('(.*)([[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12})(.*)', '\\2', eval$spectra), NA)

# make the 'strainnumber' two-digits
eval['strainnumber']<-ifelse(nchar(eval$strainnumber)==1, as.character(paste0('0', eval$strainnumber)), as.character(eval$strainnumber))

# set brukerreport NA for Eurofins as it is not read by brukersoftware and therefor not useful for merge
eval['brukercode']<-ifelse(eval$lab %in% c('device_41'), NA, as.character(eval$brukercode))

# harmonise 'run' for device_45
eval$run <-ifelse(eval$lab == 'device_45', gsub('MALDI-TOF-MS-EQA_1', 'device_45', eval$run), eval$run)

# check which are duplicated and add brukercode / run to make unique
# for spectra acquired on 'device_29' or  'device_45' the brukercode is not read by the brukersoftware and is therefor not useful for merging. Use a combination of run and strain to make this unique
eval['to_merge']<-ifelse(!is.na(eval$brukercode), paste(eval$strainnumber,eval$brukercode, eval$position,eval$Strain,eval$lab), 
                         ifelse(eval$lab == 'device_29', paste(eval$strainnumber,eval$run, eval$position,eval$Strain,eval$lab), 
                                ifelse(eval$lab == 'device_45', paste(eval$strainnumber,eval$run, eval$position,eval$Strain,eval$lab, gsub('(.*\\.)(\\d{1}$)','\\2', eval$spectra)), # for device_45 add additionally the replicate (as 15 spectra are duplicated (same run, same position, same strain))
                                paste(eval$strainnumber,eval$position,eval$Strain,eval$lab))))

# remove empty spectra (n peaks = 0)
eval <- eval[eval$n_peaks > 0,]


# check which are duplicated
dupli_eval<-eval[duplicated(eval$to_merge),]
dupli_eval<-dupli_eval[!(dupli_eval$lab %in% c("device_09","device_11","device_12","device_15","device_30",
                                               "device_35","device_36","device_37","device_43","device_44")), ]
# the ones which are have duplicated 'to_merge' are all from VitekMS. These can be merged by samplename the non-unique 'to_merge' column can be ignored

# read in brukerreports
dir<-'./Spectra_and_SpeciedID/01_baseline_quality_assessment/02_species_ID/microflexBiotyper/02_summarised/'
brukerreport_dirs <- list.files(path = dir,pattern="*.csv")
brukerreport_list <- lapply(paste0(dir, brukerreport_dirs), read.csv, sep=';', stringsAsFactors = FALSE, header=TRUE)

# remove 'spectra' and 'samplename' column from brukerreport and work withe the unique spectraname from the asciis 
for (i in 1:length(brukerreport_dirs)){
  brukerreport_list[[i]]['lab']<-gsub('_brukerreport.*.csv', '', brukerreport_dirs[[i]])
  brukerreport_list[[i]]['samplename']<-NULL
  brukerreport_list[[i]]['spectra']<-NULL
}

brukerreport_list<-lapply(brukerreport_list, function(df) mutate_at(df, .vars = c("strainnumber",  "col","Strain","species_NGS","genus_NGS"), as.character))

# define which columns are used to merge, these are the ones which are present in all files
cols_merge<-c("lab","strainnumber","col","run","Rank" ,"Score","Species","species_identified","genus_identified", "position",      
              "n_species","n_genera","n_species_over_2","n_genera_over_2","diff","Strain","species_NGS","genus_NGS","correct_genus_identified",
              "correct_species_identified","include")
# combine all brukerreports to one dataframe
brukerreport_all<-brukerreport_list %>% reduce(full_join, by = cols_merge)
# remove rows which are NA and select for the columns which are used in the columns needed 
brukerreport_all<-brukerreport_all[!is.na(brukerreport_all$Rank),c("strainnumber","col","run","position","Rank","Score","Species","lab","species_identified","genus_identified","n_species","n_genera","n_species_over_2","n_genera_over_2","diff","Strain","species_NGS","genus_NGS","correct_genus_identified","correct_species_identified","include")]
# remove the newline characters
brukerreport_all$run<-gsub('device_45\n      ','device_45',brukerreport_all$run)
brukerreport_all$run<-gsub('device_45\n      2','device_45 2',brukerreport_all$run)

# check whcih have NA in the 'include' column
check<-brukerreport_all[is.na(brukerreport_all$include),]
# ignore the ones which are empty or calibration spectra
check<-check[!(check$Species == 'no peaks found'| grepl('REFERENZ|Matrix|Negative Control|matrice|matrix|Kontrolle|CTRL|MATRIX|ctl|background control', check$strainnumber)),]
# five spectra are left. 
#   - The 2 spectra acquired on device24/25 are in 'ESCO' run, I assume these are for calibration
#   - The 2 spectra acquired on device 26 are labelled as 'ANA' and no strainnumber can be assigned, as we can not be sure which strain was measured
#   - The 1 spectra acquired on device 31 corresponds to a background control ('BLANK 2')

# remove the ones mentioned above (include = NA) and the ones which were excluded because they have been identified with high confidence as another genus (include == FALSE), as these are regarded as contaminations and will not further be evaluated
brukerreport_incl<-brukerreport_all[brukerreport_all$include == TRUE,]
# remove duplicates
brukerreport_incl<-brukerreport_incl[!duplicated(brukerreport_incl[,]),]

# for Utrecht, strain nr 25 is in folder E3 / E4, but fid file encodes B3 /B4
brukerreport_incl['position']<-ifelse(brukerreport_incl$lab == 'device_14' &
                                        brukerreport_incl$strainnumber == '25' &
                                        brukerreport_incl$position %in% c('0_B3', '0_B4'), gsub('B', 'E', as.character(brukerreport_incl$position)),as.character(brukerreport_incl$position))
# ensure that strainnumber is two digit
brukerreport_incl['strainnumber']<-ifelse(nchar(brukerreport_incl$strainnumber)==1, as.character(paste0('0', brukerreport_incl$strainnumber)), as.character(brukerreport_incl$strainnumber))
# if the brukercode is there in the correct format, use this for merging, its unique
brukerreport_incl['brukercode']<-ifelse(grepl('[[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12}', brukerreport_incl$col), gsub('(.*)([[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12})(.*)', '\\2', brukerreport_incl$col), NA)

# set brukerreport NA for device 06 and 32 as it is not included in ascii and therefor in eval
brukerreport_incl['brukercode']<-ifelse(brukerreport_incl$lab %in% c('device_06','device_32'), NA, as.character(brukerreport_incl$brukercode))

# add number of replicate to make re-measurements performed at sisli eftal hospital unique
brukerreport_incl<- brukerreport_incl %>% group_by(strainnumber,position,Strain,lab, run) %>% mutate(repl_n = row_number())

# If the brukercode is available, use this to merge, if not, use a unique combination of lab, positio and strain
brukerreport_incl['to_merge']<-ifelse(!is.na(brukerreport_incl$brukercode), paste(brukerreport_incl$strainnumber,brukerreport_incl$brukercode, brukerreport_incl$position,brukerreport_incl$Strain,brukerreport_incl$lab), 
                         ifelse(brukerreport_incl$lab == 'device_29', paste(brukerreport_incl$strainnumber,brukerreport_incl$run, brukerreport_incl$position,brukerreport_incl$Strain,brukerreport_incl$lab), 
                                ifelse(brukerreport_incl$lab == 'device_45', paste(brukerreport_incl$strainnumber,brukerreport_incl$run, brukerreport_incl$position,brukerreport_incl$Strain,brukerreport_incl$lab,brukerreport_incl$repl_n), 
                                       paste(brukerreport_incl$strainnumber,brukerreport_incl$position,brukerreport_incl$Strain,brukerreport_incl$lab))))

# for four spectra acquired on device_45, the replicate number is '1', although in the fid file, it states '2'. 
# These are probably swapped, as the fid file with '1' is empty, so is the bruker ID with '2'
# --> first remove empty bruker ID '2', then rename '1' to '2'
brukerreport_incl <- brukerreport_incl[!brukerreport_incl$to_merge == '01 device_45 0_A6 Klebsiella_pneumoniae_602149-19 device_45 2',]
brukerreport_incl <- brukerreport_incl[!brukerreport_incl$to_merge == '21 device_45 0_D10 Bordetella_bronchiseptica_502474-16 device_45 2',]
brukerreport_incl <- brukerreport_incl[!brukerreport_incl$to_merge == '21 device_45 0_D9 Bordetella_bronchiseptica_502474-16 device_45 2',]
brukerreport_incl <- brukerreport_incl[!brukerreport_incl$to_merge == '29 device_452 0_E10 Streptococcus_pseudopneumoniae_610886-17 device_45 2',]

brukerreport_incl$to_merge <- gsub('01 device_45 0_A6 Klebsiella_pneumoniae_602149-19 device_45 1', '01 device_45 0_A6 Klebsiella_pneumoniae_602149-19 device_45 2',brukerreport_incl$to_merge)
brukerreport_incl$to_merge <- gsub('21 device_45 0_D10 Bordetella_bronchiseptica_502474-16 device_45 1', '21 device_45 0_D10 Bordetella_bronchiseptica_502474-16 device_45 2',brukerreport_incl$to_merge)
brukerreport_incl$to_merge <- gsub('21 device_45 0_D9 Bordetella_bronchiseptica_502474-16 device_45 1', '21 device_45 0_D9 Bordetella_bronchiseptica_502474-16 device_45 2',brukerreport_incl$to_merge)
brukerreport_incl$to_merge <- gsub('29 device_452 0_E10 Streptococcus_pseudopneumoniae_610886-17 device_45 1', '29 device_452 0_E10 Streptococcus_pseudopneumoniae_610886-17 device_45 2',brukerreport_incl$to_merge)

brukerreport_incl$repl_n <- NULL

# check which duplicated and add 'run' variable to make unique, where necessary
dupli_bruker<-brukerreport_incl[duplicated(brukerreport_incl$to_merge),]
dupli_bruker<-dupli_bruker[!is.na(dupli_bruker$strainnumber),] # none is duplicated

length(intersect(eval$to_merge, brukerreport_incl$to_merge))

brukerreport_incl<-brukerreport_incl[!is.na(brukerreport_incl$to_merge),] # the rows which are NA in the 'to_merge' column are 'NA' in all other columns too. Remove these
setdiff(brukerreport_incl[brukerreport_incl$Species !='no peaks found',]$to_merge,eval$to_merge) # there are 83 rows in the brukerreports which are 'no peaks found'. Leave this is, as it is possible that there are differences in the peakpicking. The spectra, where no peaks where detected with either peakpicking are excluded during the merge

check<-eval[eval$to_merge %in% setdiff(eval$to_merge, brukerreport_incl[brukerreport_incl$Species !='no peaks found',]$to_merge), ]
check<-check[!(check$lab %in% c('device_12','device_11','device_35', 'device_36', 'device_37', 'device_44', 'device_43', 'device_30', 'device_15', 'device_09')),]
# there is one spectrum acquired on device_01 which cannot be merged. in the brukerreport it is labelled as 'manuell-24' and in the ascii as 'manuell 46'. It is excluded, as it cannot unambiguously me assigned to a strain. 

brukerreport_incl$brukercode<-NULL
brukerreport_incl$Strain<-NULL
brukerreport_incl$lab<-NULL
brukerreport_incl$position<-NULL
brukerreport_incl$run<-NULL
brukerreport_incl$strainnumber<-NULL

# add 'brukerDB' tag to all columns which have to do with the brukerDB before merging
colnames(brukerreport_incl)[colnames(brukerreport_incl) %in% c("col", "Rank", "Score", "Species", "species_identified", "genus_identified", "n_species", "n_genera", "n_species_over_2", "n_genera_over_2", "diff", "correct_genus_identified", "correct_species_identified")] <- paste0('brukerDB.', c("col", "Rank", "Score", "Species", "species_identified", "genus_identified", "n_species", "n_genera", "n_species_over_2", "n_genera_over_2", "diff", "correct_genus_identified", "correct_species_identified"))

# merge brukerreport and eval
eval2 <- merge(eval, brukerreport_incl, by = 'to_merge', all.x = T)

# remove the 'to_merge' column, as it is not required any more
eval2$to_merge<-NULL

# add group information 
groups<-read.csv2('/Users/aline/ESCMID/ESPRIT/04_Strains/Strains_groups_assigned_2021_02.csv', sep=',')
# ensure that the strainnumber is two digits
groups['Numbering_Shipment_1']<-ifelse(nchar(groups$Numbering_Shipment_1)==1, as.character(paste0('0', groups$Numbering_Shipment_1)), as.character(groups$Numbering_Shipment_1))
# merge the 'Group' information to the 'eval' file
eval<-merge(eval2, groups, by.x = c('Strain','strainnumber'), by.y = c('Strain','Numbering_Shipment_1'), all.x = TRUE)
# remove the 'PDD' tag from the spectra column
eval['spectra']<-gsub('\\#PSD\\=', '', eval$spectra)

# set all below threshold (1.7) to 'no Identification possible)
eval$brukerDB.species_identified<-ifelse(eval$brukerDB.Score < 1.7, 'no ID possible', eval$brukerDB.species_identified)
eval$brukerDB.Species<-ifelse(eval$brukerDB.Score < 1.7, 'no ID possible', eval$brukerDB.Species)

# add species identification by the vitekMS DB
# import files
VitekMSreport_dirs = list.files(path ='./Spectra_and_SpeciedID/01_baseline_quality_assessment/02_species_ID/VitekMS/02_summarised/',pattern="*_VitekMSreport.csv")
VitekMSreport_list = lapply(paste0('./Spectra_and_SpeciedID/01_baseline_quality_assessment/02_species_ID/VitekMS/02_summarised/',VitekMSreport_dirs), read.csv, sep=';', stringsAsFactors = FALSE, header=TRUE)

# convert these columns to character before merging
VitekMSreport_list<-lapply(VitekMSreport_list, function(df) mutate_at(df, .vars = c("strainnumber", "samplename_new","Strain","species_NGS","genus_NGS"), as.character))

# define which columns to merge
cols_merge<-c("X","strainnumber","samplename_new","files","nb_peaks","max_SN","max_intensities","max_resolution","identifType","species_rank_1","proba_rank_1","score_rank_1","lab","Strain","n_species","diff","diff_proba","species_identified","genus_identified","species_NGS","genus_NGS","correct_genus_identified","correct_species_identified","include", "position")
# merge all vitekMS files to one dataframe
VitekMSreport_all<-VitekMSreport_list %>% reduce(full_join, by = cols_merge)
# remove hashtag from KLBS2 samplename
#VitekMSreport_all$samplename_new<-gsub('\\#\\.KLBS2\\_', '\\.KLBS2\\_', VitekMSreport_all$samplename_new)

# spectra where 'it could not be 'include' is 'NA' are all of 'identifType' NoIdentification, NotEnoughPeaks or TooNoisy. These are kept in as only spectra are excluded which are identified as an other genus with high confidence, as these are regarded as contaminations. 
VitekMSreport_all['include']<-ifelse(is.na(VitekMSreport_all$include), TRUE, as.character(VitekMSreport_all$include))
# xclude the spectra which were identified as another genus with high confidence
VitekMSreport_incl<-VitekMSreport_all[VitekMSreport_all$include!=FALSE,]
#make sure strainnumber is two digit before merging
VitekMSreport_incl$strainnumber<-ifelse(grepl('^[[:digit:]]$', VitekMSreport_incl$strainnumber), paste0('0', VitekMSreport_incl$strainnumber), VitekMSreport_incl$strainnumber)
# add vitekMS tag to all columns from VitekMS species identification before merging
colnames(VitekMSreport_incl)[colnames(VitekMSreport_incl) %in% c("samplename_new", "files", "position", "nb_peaks", "max_SN", "max_intensities", "max_resolution", "identifType", "species_rank_1", "proba_rank_1", "score_rank_1", "n_species", "diff", "diff_proba", "species_identified", "genus_identified", "correct_genus_identified", "correct_species_identified", "species_NGS","genus_NGS","include")] <- paste0('vitekMSDB.', c("samplename_new", "files", "position", "nb_peaks", "max_SN", "max_intensities", "max_resolution", "identifType", "species_rank_1", "proba_rank_1", "score_rank_1", "n_species", "diff", "diff_proba", "species_identified", "genus_identified", "species_NGS", "genus_NGS", "correct_genus_identified", "correct_species_identified", "include"))

# set all rows which are of identifType 'NoIdentification', 'NotEnoughPeaks', 'TooNoisy' to 'noID possible'
VitekMSreport_incl$vitekMSDB.species_identified<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), 'no ID possible', VitekMSreport_incl$vitekMSDB.species_identified)
VitekMSreport_incl$vitekMSDB.species_rank_1<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), 'no ID possible', VitekMSreport_incl$vitekMSDB.species_rank_1)
VitekMSreport_incl$vitekMSDB.proba_rank_1<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), NA, VitekMSreport_incl$vitekMSDB.proba_rank_1)
VitekMSreport_incl$vitekMSDB.score_rank_1<-ifelse(VitekMSreport_incl$vitekMSDB.identifType %in% c('NoIdentification', 'NotEnoughPeaks', 'TooNoisy'), NA, VitekMSreport_incl$vitekMSDB.score_rank_1)

#remove '_02' tag from LTWAG spectra name
eval['to_merge']<-ifelse(grepl('device_44', eval$spectra), gsub('\\_[[:digit:]]{1,3}$', '', eval$spectra), as.character(eval$spectra))
eval['to_merge']<-ifelse(grepl('device_44', eval$to_merge), gsub('\\_[[:digit:]]{1,3}$', '', eval$to_merge), as.character(eval$to_merge))

# check whether all vitekMS spectra are in 'eval' file
length(intersect(paste(VitekMSreport_incl$vitekMSDB.samplename_new,VitekMSreport_incl$Strain,VitekMSreport_incl$strainnumber,VitekMSreport_incl$lab), paste(eval$to_merge, eval$Strain,eval$strainnumber,eval$lab)))
# check whether there are differences
setdiff(paste(VitekMSreport_incl$vitekMSDB.samplename_new,VitekMSreport_incl$Strain,VitekMSreport_incl$strainnumber,VitekMSreport_incl$lab), paste(eval$to_merge, eval$Strain,eval$strainnumber,eval$lab))
# there are none

# merge vitekMS reports to eval file
eval<-merge(eval, VitekMSreport_incl, by.x = c('to_merge', "Strain","strainnumber","lab"), by.y = c("vitekMSDB.samplename_new","Strain","strainnumber","lab"), all.x = T)
eval$to_merge<-NULL
# for some spectra, the 'species_NGS', 'genus_NGS' and  'include' column are present in only one of the input files. create columns querying both
eval['species_NGS']<-ifelse(!is.na(eval$species_NGS), eval$species_NGS, as.character(eval$vitekMSDB.species_NGS))
eval['genus_NGS']<-ifelse(!is.na(eval$genus_NGS), eval$genus_NGS, as.character(eval$vitekMSDB.genus_NGS))
eval['include']<-ifelse(!is.na(eval$include), eval$include, as.character(eval$vitekMSDB.include))
# then remove the redundant one
eval$vitekMSDB.species_NGS<-NULL
eval$vitekMSDB.genus_NGS<-NULL
eval$vitekMSDB.include<-NULL

# Import the marker based species identification
markerID<-read.csv2('./Spectra_and_SpeciedID/01_baseline_quality_assessment/02_species_ID/marker_ID/02_summarised/markerID_sum_clean_2021_2.csv', sep=',')

# In marker ID file: 
# remove 'PSD' tag from ascii filename
markerID['sample']<-gsub('\\#PSD.*', '', markerID$sample)
# remove everything before lab ID for spectra acquired on device_16
markerID['sample']<-gsub('(.*)(\\.)(device\\_16.*$)', '\\3', markerID$sample)
# remove '$' for ascii spectra
markerID['sample']<-gsub(' $', '', markerID$sample)
# make sure that strainnumber is two digit
markerID['strainnumber']<-ifelse(grepl('^[[:digit:]]$',markerID$strainnumber), paste0('0', markerID$strainnumber), as.character(markerID$strainnumber))

# Include the same harmonisation stept to 'eval' file before merging
#  remove 'PSD' tag from ascii filename
eval['spectra']<-gsub('\\#*PSD.*$', '', eval$spectra)
# remove everything before lab ID for spectra acquired on device_16
eval['spectra']<-gsub('(.*)(\\.)(device\\_16.*$)', '\\3', eval$spectra)
# remove hashtag from samplename
eval['spectra']<-gsub('\\#', '', eval$spectra)
# remove '$' for ascii spectra
eval['spectra']<-gsub(' $', '', eval$spectra)

# check which are missing from markerID. No marker ID, means that these have been excluded for quality reasons, eg. if not enough calibration masses have been detected 
no_marker<-eval[eval$spectra %in% setdiff(eval$spectra,markerID$sample),]
table(no_marker$lab) # all within range 

# only spectra are excluded, if they have gone through a marker 'classifier'. All the ones which went through PAPMID should be here. Check
classifier_strains <- unique(markerID[!is.na(markerID$profiles_matched), 'strainnumber'])
# one spectrum is still missing. #20211001
no_marker_check <- no_marker[!no_marker$strainnumber %in% classifier_strains,]

# check, whether there are some which have a marker ID, but are not in eval (there should not be any)
missing<-markerID[markerID$sample %in% setdiff(markerID$sample,eval$spectra), ]

# create 'id' column to check, whether one onew which are missing here are all missing because they excluded (contamination)
missing['lab']<-gsub('([^\\.]+)(\\..*)', '\\1', missing$sample)
missing['strainnumber']<-gsub('([^\\.]+)(\\.)([[:digit:]]{1,2})(.*)', '\\3', missing$sample)
missing['position']<-ifelse(grepl('.*[[:punct:]]0\\_[[:alpha:]][[:digit:]]{1,2}[[:punct:]]*.*', missing$sample), gsub('(.*[[:punct:]])(0\\_[[:alpha:]][[:digit:]]{1,2})([[:punct:]]*.*)', '\\2', missing$sample), NA)
missing['position']<-ifelse(is.na(missing$position), str_extract(missing$sample,'[[:digit:]][[:alpha:]][[:digit:]]{1,2}$'), missing$position) 
missing['position']<-ifelse(is.na(missing$position), gsub('(.*\\_)([[:alpha:]][[:digit:]])(\\_.*)', '\\2', missing$sample), missing$position)
missing['id']<-paste(missing$strainnumber, missing$position, missing$lab)

# remove empty
missing<-missing[missing$match_count>0,]

# import file listing the ones which were which were excluded (contamination)
check<-read.csv('./Spectra_and_SpeciedID/01_baseline_quality_assessment/02_species_ID/exclude_wrong_genus.csv', sep=';')
wrong_genus<-check
setdiff(missing$id, check$id)# In 'check' the ones which were excluded because they were another genus than what they should be. Identified by Bruker DB and in script 'check_exclude.R'. 

# remove these
missing<-missing[!(missing$id %in% check$id),]
# 15 spectra have no meaning fill samplename, these are ignored 
missing<-missing[nchar(missing$sample)> 2,] 

# the two device_44 spectra which are still in the 'missing' file have also been identified as a wrong genus (Staphylococcus instead of Gardnarella), it is ok that there are excluded. 

# remove the columns indicating which subunits were detected (PAPMID) and which profiles were  matchied (Classifier), as these are not further evaluated
markerID$X.Subunit.Mass.<-NULL
markerID$profiles_matched<-NULL
markerID$Run<-NULL

# Rename 'Genus' to 'genus_identified
colnames(markerID)[colnames(markerID) == 'Genus'] <-"genus_identified"

#remove duplicates (are duplicated for all rown, have probably been run twice through DB)
markerID<-markerID[!duplicated(markerID),]

# add marker ID tag to all columnnames which come from the marker based species identification
colnames(markerID)[colnames(markerID) %in% c("match_count", "species_identified", "genus_identified", "n_species_identified", "n_genera_identified", "correct_genus_identified", "correct_species_identified")] <- paste0('markerDB.', c( "match_count", "species_identified", "genus_identified", "n_species_identified", "n_genera_identified", "correct_genus_identified", "correct_species_identified"))

# merge marker ID to eval
eval<-merge(eval, markerID, by.x = c('spectra',"Strain","strainnumber","species_NGS","genus_NGS"), by.y = c('sample',"Strain","strainnumber","species_NGS","genus_NGS"), all.x = T)
# create column 'AMLDI.device_full'. every measurement dataset receives its on ID
eval['MALDI.device_full']<-eval$lab # includes all specification and also device 35/36/37 and device 38/39/40 for lab BE, where the dataset has been measured by three different people

# create 'MALDI_device' column. IHere, all acquired on the same device have the same MALDI_device
eval['MALDI.device']<-eval$MALDI.device_full
eval$MALDI.device <- gsub('device_36|device_37', 'device_35', eval$MALDI.device)
eval$MALDI.device <- gsub('device_39|device_40', 'device_38', eval$MALDI.device)

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

eval<-eval[!is.na(eval$genus_NGS),] # these are three spectra where the samplename in the brukerreport does not match the samplename in the json file, exclude these

# export this
write.table(eval, './Spectra_and_SpeciedID/01_baseline_quality_assessment/eval_participating_labs_1.txt', row.names = F, quote = F, dec = '.', sep = ';', eol = "\n")

