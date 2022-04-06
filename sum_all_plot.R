library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(stringr)


#import sumary file from the baseline quality assessment
sum_round_1 <- read.csv2('./Spectra_and_SpeciedID/01_baseline_quality_assessment/eval_participating_labs_1.txt')
# remove empty spectra
sum_round_1 <- sum_round_1[!sum_round_1$n_peaks == '0',]

# extract which are bruker and which are vitekMS / shimadzu devices, by checking whach DB was used for species Identification
bruker.devices<-as.character(unique(sum_round_1[!is.na(sum_round_1$brukerDB.correct_species_identified), 'MALDI.device_full']))
VitekMS.devices<-as.character(unique(sum_round_1[!is.na(sum_round_1$vitekMSDB.correct_species_identified), 'MALDI.device_full']))
as.character(unique(sum_round_1[!is.na(sum_round_1$vitekMSDB.correct_species_identified), 'device_id']))

# define shimadzu divices
Shimadzu.devices<-c("device_35","device_36","device_37","device_43","device_44")

#import summary file from the spectra acquired after the intervention
sum_round_2 <- read.csv2('./Spectra_and_SpeciedID/02_quality_assessment_after_intervention/eval_participating_labs_2.txt')
# remove 'to_merge' column, not needed any more
sum_round_2$to_merge <- NULL

#check whether there are differences
setdiff(sum_round_1$MALDI.device, sum_round_2$MALDI.device) # these 4 laboratories have not participated in the quality assessment after the intervention
setdiff(sum_round_2$MALDI.device, sum_round_1$MALDI.device) # all which participated in the intervention, have already participated in baseline intervention

setdiff(colnames(sum_round_1), colnames(sum_round_2)) # check difference in columns between the two files
setdiff(colnames(sum_round_2), colnames(sum_round_1))

# Add method column. call method from round 1, method 0 (routinely acquired)
sum_round_1['method'] <- '0'
sum_round_1$MALDI.device.1 <- NULL
# add replicate (call replicate_eff, as it describes the number of *effective* replicates rather that what the spectra were called)
sum_round_1 <- sum_round_1 %>% group_by(strainnumber, MALDI.device) %>% mutate(replicate_eff = row_number())

# Same for second round add replicate (call repliacate_eff, as it describes the number of *effective* replicates rather that what the spectra were called)
sum_round_2 <- sum_round_2 %>% group_by(strainnumber, MALDI.device, method) %>% mutate(replicate_eff = row_number())

# remove columns which are not there for both
sum_round_1 <- sum_round_1[, intersect(colnames(sum_round_1), colnames(sum_round_2))]
sum_round_2 <- sum_round_2[, intersect(colnames(sum_round_1), colnames(sum_round_2))]

# merge the two datasets
sum <- merge(sum_round_1, sum_round_2, by = intersect(colnames(sum_round_1), colnames(sum_round_2)), all = T)

# summarise baseline quality assessment per device. Summarise the following endpoints:
# how many ribosomal subunits were found, what the reproducibility between 
# what the reproducibility between two technical replicates and
# mean measurement error

# safe these as numeric
sum$frac_peaks_repr_2spectra <- as.numeric(as.character(sum$frac_peaks_repr_2spectra))
sum$mean_dist_ppm <- as.numeric(as.character(sum$mean_dist_ppm))

#calculate interquartile ranges and median per device
sum_device <- sum %>% group_by(method, device_id) %>%
  summarise(lower_iqr_n_ribos = quantile(n_ribos_detected, probs = 0.25), 
            median_n_ribos = quantile(n_ribos_detected, probs = 0.5), 
            higher_iqr_n_ribos = quantile(n_ribos_detected, probs = 0.75), 
            lower_iqr_frac = quantile(frac_peaks_repr_2spectra, probs = 0.25), 
            median_frac = quantile(frac_peaks_repr_2spectra, probs = 0.5), 
            higher_iqr_frac = quantile(frac_peaks_repr_2spectra, probs = 0.75), 
            lower_iqr_m_error = quantile(mean_dist_ppm, probs = 0.25, na.rm = T), 
            median_m_error = quantile(mean_dist_ppm, probs = 0.5, na.rm = T), 
            higher_iqr_m_error = quantile(mean_dist_ppm, probs = 0.75, na.rm = T))
# check values for baseline quality assessment
# View(sum_device[sum_device$method == '0',])

# check the median over all device
quantile(sum$n_ribos_detected)
quantile(sum$mean_dist_ppm, na.rm = T)

# compare gram positives and gram negatives in terms of number of ribosomal subunits detected and reproducibility between technical replicates 
sum$Gram <- gsub('variable', 'positive', sum$Gram)
sum$Gram <- gsub('Positive', 'positive', sum$Gram)
sum$Gram <- gsub('Negative', 'negative', sum$Gram)
sum_gram <- sum %>% group_by(method, Gram) %>%
  summarise(lower_iqr_n_ribos = quantile(n_ribos_detected, probs = 0.25), 
            median_n_ribos = quantile(n_ribos_detected, probs = 0.5), 
            higher_iqr_n_ribos = quantile(n_ribos_detected, probs = 0.75), 
            lower_iqr_frac = quantile(frac_peaks_repr_2spectra, probs = 0.25), 
            median_frac = quantile(frac_peaks_repr_2spectra, probs = 0.5), 
            higher_iqr_frac = quantile(frac_peaks_repr_2spectra, probs = 0.75) )


# use unpaired wilcoxon rank test, as the spectra are from different strains
wilcox.test(sum[sum$method == '0' & sum$Gram == 'positive', 'n_ribos_detected'], 
            sum[sum$method == '0' & sum$Gram == 'negative', 'n_ribos_detected'], 
            paired = F)
wilcox.test(sum[sum$method == '0' & sum$Gram == 'positive', 'frac_peaks_repr_2spectra'], 
            sum[sum$method == '0' & sum$Gram == 'negative', 'frac_peaks_repr_2spectra'], 
            paired = F)


# split up streptococcus in viridans and non-viridans
sum['Group']<-as.factor(ifelse(grepl('infantis|pseudopneumoniae|pneumoniae|gordonii', as.character(sum$Strain)) & sum$genus_NGS == "Streptococcus", 'Viridans Streptococci', 
                                     ifelse(grepl('equinus|gallolyticus|lutetiensis|dysgalactiae', as.character(sum$Strain)) & sum$genus_NGS == "Streptococcus", 'Other Streptococci', as.character(sum$Group))))


# set all bruker ID which have a score lower than 1.7 to 'no ID'
sum["brukerDB.correct_species_identified"]<-ifelse(sum$brukerDB.Score<1.7, 'no ID', sum$brukerDB.correct_species_identified)
sum['brukerDB.correct_species_identified']<-ifelse(sum$brukerDB.correct_species_identified %in% c('1', 'TRUE'), TRUE, 
                                                         ifelse(sum$brukerDB.correct_species_identified %in% c('1', 'FALSE'), FALSE, 
                                                                ifelse(sum$brukerDB.correct_species_identified == 'no ID', 'no ID', NA)))

# marker based db: set all hits with a match count lower than a certain threshold to 'no ID'
# Use the following thresholds: 
# E. coli / Shigella: 15
# Enterobacter: 20
# Klebsiella: 15
# S.aureus 7: 
# all other: 10
sum['species_NGS']<-gsub('([[:alpha:]]+)(\\_)([[:alpha:]]+)(\\_)(.*)', '\\3', sum$Strain)
sum['genus_NGS']<-gsub('([[:alpha:]]+)(\\_)([[:alpha:]]+)(\\_)(.*)', '\\1', sum$Strain)
#set threshold
sum['marker.theshold']<-ifelse(sum$genus_NGS %in% c('Escherichia', 'Shigella', 'Klebsiella'), 15,
                                     ifelse(sum$genus_NGS == 'Enterobacter', 20,
                                            ifelse(sum$genus_NGS == 'Staphylococcus', 7, 10)))
# set all values to 'not-odentified' if the match count is below 0
sum['markerDB.match_count']<-ifelse(is.na(sum$markerDB.match_count), 0, as.numeric(as.character(sum$markerDB.match_count)))                                  
sum["markerDB.species_identified"]<-ifelse(as.numeric(sum$markerDB.match_count) < sum$marker.theshold, 'no identification possible',sum$markerDB.species_identified)
sum["markerDB.Genus"]<-ifelse(as.numeric(sum$markerDB.match_count) < sum$marker.theshold, 'no identification possible',sum$markerDB.genus_identified)
#sum["markerDB.match_count"]<-ifelse(as.numeric(sum$markerDB.match_count) < sum$marker.theshold, 0,sum$markerDB.match_count)
sum["markerDB.n_species_identified"]<-ifelse( as.numeric(sum$markerDB.match_count) < sum$marker.theshold, 0,sum$markerDB.n_species_identified)
sum["markerDB.n_genera_identified"]<-ifelse( as.numeric(sum$markerDB.match_count) < sum$marker.theshold, 0,sum$markerDB.n_genera_identified)
sum["markerDB.correct_species_identified"]<-ifelse( as.numeric(sum$markerDB.match_count) < sum$marker.theshold, 'no ID',sum$markerDB.correct_species_identified)
sum$marker.theshold<-NULL
sum['markerDB.correct_single_species']<-ifelse(sum$markerDB.correct_species_identified == TRUE & sum$markerDB.n_species_identified == '1', TRUE, FALSE)

# Vitek MS
# rename vitekMSDB.identifType column
sum['VitekMS.DBidentifType']<-sum$vitekMSDB.identifType
sum$vitekMSDB.identifType<-NULL

# Set all vitekMS ID who have a indentification type of 'no ID, 'not enough peaks', 'too noisy' all to 'no ID'
sum['VitekMS.DBcorrect_species_identified']<-ifelse(sum$VitekMS.DBidentifType %in% c('NoIdentification','NotEnoughPeaks','TooNoisy'), 'noID',
                                                          ifelse(sum$vitekMSDB.correct_species_identified == 0, FALSE,
                                                                 ifelse(sum$vitekMSDB.correct_species_identified == 1, TRUE, NA)))

# add the 'mean intensity peaks' ny deviding total int / n_peaks
sum['mean_intensity_peaks']<-as.numeric(as.character(sum$total_intensity_peaks)) / as.numeric(as.character(sum$n_peaks))
# correct for typo
colnames(sum)[colnames(sum) == 'median_intensity'] <- 'median_intensity'

# add which maldi has been used
sum['MALDI.device.type']<-ifelse(sum$MALDI.device_full %in% bruker.devices, 'MBT Biotyper', 'VitekMS / Axima Confidence')
sum['MALDI.device.type']<-factor(sum$MALDI.device.type, levels = c('MBT Biotyper', 'VitekMS / Axima Confidence'))

## Summarise per strain, protocol and database how many correct identification
sum_id <- sum[sum$method %in% c(0,1,2),] %>% # leave out 'method 3'
  group_by(strainnumber, Strain, method) %>%
  summarise(correct_bruker =  sum(grepl('TRUE', brukerDB.correct_species_identified)),
            all_bruker = sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            perc_correct_bruker = sum(grepl('TRUE', brukerDB.correct_species_identified))/ sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
            correct_vitek = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified)), 
            all_vitek = sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)), 
            perc_correct_vitek = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified)) / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)), 
            correct_marker = sum(grepl('TRUE', markerDB.correct_species_identified)), 
            all_marker = sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)), 
            perc_correct_marker = sum(grepl('TRUE', markerDB.correct_species_identified)) / sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)))

write.csv2(sum_id, './01_Spectra/04_outputs/sum_id.csv', quote = F, row.names = F)

### plot baseline quality assessment results
# plot the mean measurement error [ppm] per device
dist_ppm<-ggplot(sum[sum$method == '0',], aes(x=reorder(device_id, mean_dist_ppm,na.rm = TRUE, FUN = median), y=mean_dist_ppm, col = MALDI.device.type)) +
  geom_boxplot(lwd=0.2, outlier.size = 0.5) +  scale_color_manual(values=c("goldenrod3","turquoise")) +
  ylab('Mean measurement error [ppm]') +
  xlab('Device ID') + theme(axis.text.x = element_text(angle = 60,  size = 7), axis.text.y = element_text(size = 7)) +
  geom_hline(aes(yintercept = quantile(sum[sum$method == '0','mean_dist_ppm'], na.rm = T, probs = 0.5)), linetype="solid", color = "red", size=0.4) +
  geom_hline(aes(yintercept = quantile(sum[sum$method == '0','mean_dist_ppm'], na.rm = T, probs = 0.25)), linetype="dashed", color = "red", size=0.4) +
  geom_hline(aes(yintercept = quantile(sum[sum$method == '0','mean_dist_ppm'], na.rm = T, probs = 0.75)), linetype="dashed", color = "red", size=0.4)

dist_ppm

pdf('./mean_dist_labs.pdf', height = 3, width = 5.5)
dist_ppm + theme(axis.text.x = element_text(angle = 60,  size = 7, hjust = 1), axis.text.y = element_text(size = 7), strip.text = element_text(size=3),legend.position = 'bottom')
dev.off()

# In order to plot the mass spectral quality acquired of the baseline quality evaluation, transform the sum[sum$method == '0',] dataframe to 'long' format and add proportions of species identification 
# calculate proportions of species identification evaluations per device (ie what fraction of spectra was correctly identified)
eval2 <- sum[sum$method == '0',] %>%
  group_by(MALDI.device_full) %>%
                  #bruker ID summary
  dplyr::summarize(prop.bruker.correct.single.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
                   prop.bruker.correct.multi.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 >= 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
                   prop.bruker.correct.under.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
                   prop.noID.bruker = sum(grepl('no ID', brukerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
                   prop.wong.bruker.over.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score >= 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
                   prop.wong.bruker.under.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score < 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
                   # marker ID summary
                   prop.correct.marker.single = sum(grepl('TRUE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),             
                   prop.wrong.marker.single = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('FALSE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   prop.correct.marker.multi = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   prop.noID.marker = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('no ID', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   # VitekMS ID summary
                   prop.correct.VitekMS.single = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
                   prop.correct.VitekMS.lowdiscriminatory = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
                   prop.wrong.VitekMS.single = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
                   prop.wrong.VitekMS.lowdiscriminatory = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
                   prop.noID.VitekMS = sum(grepl('noID', VitekMS.DBcorrect_species_identified)) /sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)), 
                   n_correct_single = sum(grepl('TRUE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)),
                   n_all_spectra_ID = sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified))) 

eval2_sum_correct<-eval2
# add proportion of 'marker correct' referring to the ones which have been identified on single species level AND the ones which which include the correct species in their multispecies ID
eval2_sum_correct['marker_correct']<-eval2$prop.correct.marker.single + eval2$prop.correct.marker.multi
# remove these columns, not needed for further analyses
eval2$n_correct_single<-NULL
eval2$n_all_spectra_ID<-NULL

# check prop of correct species identification per strain for the baseline quality assessment
eval2_strain <- sum[sum$method == '0',] %>%
  group_by(Strain, strainnumber) %>%
  dplyr::summarize(prop.correct.marker.single = sum(grepl('TRUE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),             
                   prop.wrong.marker.single = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('FALSE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   prop.correct.marker.multi = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   prop.noID.marker = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('no ID', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   n_correct_single = sum(grepl('TRUE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)),
                   n_all_spectra_ID = sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified))) 
# add proportion of 'marker correct' referring to the ones which have been identified on single species level AND the ones which which include the correct species in their multispecies ID
eval2_strain['marker_correct']<-eval2_strain$prop.correct.marker.single + eval2_strain$prop.correct.marker.multi

# the eval2_strain will not be used for further analysis, continue with eval2 which includes proportions per device

# merge these proportions back to the dataframe
eval<-merge(sum[sum$method == '0',] , eval2, by = intersect(colnames(sum[sum$method == '0',]), colnames(eval2)), all.x  = T)
#correct for typo
names(eval)[names(eval) == 'meadian_intensity'] <- 'median_intensity'
#convert to 'long format
eval_prop_long<-gather(eval2, "classification_eval_prop", "value", -c(MALDI.device_full))
#Gather all endpoints to one column in order to plot all at once using facet_wrap
eval_long<-gather(eval[,c('max_mass_ribo', 'spectra',"n_peaks", "mean_intensity_peaks","highest_peaks", "mass_90", "n_high_peaks", "frac_peaks_repr_2spectra", "n_ribos_detected", "n_predicted_su", "rel_amount_su_detected", "mean_dist", "mean_dist_ppm", "mean_intensity","median_intensity", "total_intensity_peaks","brukerDB.Score", "brukerDB.n_species", "brukerDB.n_genera", "brukerDB.n_species_over_2", "brukerDB.n_genera_over_2", "brukerDB.diff", "brukerDB.correct_genus_identified", "brukerDB.correct_species_identified",  "markerDB.match_count",  "markerDB.n_genera_identified",  "markerDB.correct_genus_identified", "markerDB.correct_species_identified", "markerDB.correct_single_species")], "Endpoint", "value", -spectra)
# merge back to the dataframe
eval_long<-merge(eval[,c("spectra","Strain", "strainnumber", "Shape", "Gram", "aerobic", "Family", "Group", "spectra", "species_NGS", "genus_NGS", 'MALDI.device_full', 'device_id')], eval_long, by='spectra')
# convert TRUE/FALSE to 1/0 
eval_long['value']<-ifelse(eval_long$value == 'TRUE', 1,
                           ifelse(eval_long$value == 'FALSE', 0, eval_long$value))
# convert ',' to '.' to be able to interpret as numerical
eval_long["value"]<-gsub(',', '.', eval_long$value)
eval_long["value"]<-as.numeric(as.character(eval_long$value)) 
# add a tag which DB was used
eval_prop_long['DB']<-ifelse(grepl('bruker', eval_prop_long$classification_eval_prop), 'bruker', 
                             ifelse(grepl('marker', eval_prop_long$classification_eval_prop), 'marker', 
                                    ifelse(grepl('VitekMS', eval_prop_long$classification_eval_prop), 'VitekMS', NA)))
# create 'eval' column
eval_prop_long['eval']<-eval_prop_long$classification_eval_prop

# check whether there is a correlation between (i) the fraction of correctly identified spectra (single and multi-species identification together) and (ii) the median number of detected marker masses per device
sum_round_1_median_marker <- sum_round_1 %>% 
  group_by(MALDI.device_full) %>%
  summarise(median_nr_marker = median(n_ribos_detected))

corr_check <- merge(eval2[, c('MALDI.device_full','prop.correct.marker.single','prop.correct.marker.multi')], sum_round_1_median_marker, by = 'MALDI.device_full')
corr_check['correct_marker_prop'] <- corr_check$prop.correct.marker.single + corr_check$prop.correct.marker.multi

cor(corr_check$median_nr_marker, corr_check$correct_marker_prop, method = "pearson")

# define colours which species evaluation willbe drawn in
#"prop.bruker.correct.single.over.two"      #009E73
#"prop.correct.marker.single"               #009E73
#"prop.correct.VitekMS.single"              #009E73
#"prop.bruker.correct.multi.over.two"       "#0072B2"
#"prop.correct.marker.multi"                "#0072B2"
#"prop.correct.VitekMS.lowdiscriminatory"  "#56B4E9"
#"prop.bruker.correct.under.two"           "#56B4E9"
#"prop.noID.bruker"                         "#999999"     
#"prop.noID.marker.single"                  "#999999"
#"prop.noID.VitekMS"                        "#999999"
#"prop.wong.bruker.under.two"               "#E69F00"
#"prop.wrong.VitekMS.lowdiscriminatory"     "#E69F00"
#"prop.wong.bruker.over.two"                "#D55E00"
#"prop.wrong.marker.single"                 "#D55E00"
#"prop.wrong.VitekMS.single"                "#D55E00"  

#for overall figure, exclusively include marker ID. Thus remove 'Vitek' and 'Bruker' DB evaluations
eval_prop_long<-eval_prop_long[eval_prop_long$DB == 'marker',]
# define the ordering of the colour stack by defining the factor levels
eval_prop_long$eval<- factor(eval_prop_long$eval, levels = c("prop.wrong.marker.single", "prop.noID.marker", 
                                                             "prop.correct.marker.multi", 
                                                             "prop.correct.marker.single"))
# Plot the number of ribosomal subunits and the evaluation of the marker based species identification per device
plot_data<-eval_long[eval_long$Endpoint %in% c("n_ribos_detected","median_intensity","total_intensity_peaks", "frac_peaks_repr_2spectra"),]
# add column Endpoint_relative including the relative amount of subunits detected (will be used when comparing between different bacterial groups)
plot_data_rel<-eval_long[eval_long$Endpoint  == "rel_amount_su_detected",c("spectra","value", "Endpoint")]
# tag these 'relative' columns with 'rel' tag
colnames(plot_data_rel)<-c("spectra","value.rel","Endpoint.rel")
# merge back to the original data 
plot_data<-merge(plot_data, plot_data_rel, by.x="spectra", by.y = "spectra", all =T)
# only use relative value for the number of ribosomal subunits detected when comparing between different bacterial groups
plot_data$Endpoint.rel <- ifelse(plot_data$Endpoint == 'n_ribos_detected', plot_data$Endpoint.rel, plot_data$Endpoint)
plot_data$value.rel <- ifelse(plot_data$Endpoint == 'n_ribos_detected', plot_data$value.rel, plot_data$value)
# rename 'DB' column to endpoint
names(eval_prop_long)[names(eval_prop_long) == 'DB'] <- 'Endpoint'
eval_prop_long$classification_eval_prop<-NULL

# add device id first
eval_prop_long<-merge(plot_data[!duplicated(plot_data$MALDI.device_full),c("MALDI.device_full", "device_id")], eval_prop_long, by = "MALDI.device_full", all.y  = TRUE)
plot_data<-merge(plot_data, eval_prop_long, by = intersect(colnames(eval_prop_long), colnames(plot_data)), all = TRUE)

#take out replicate datasets acquired at laboratory BE as these were acquired on the same device
plot_data<-plot_data[!plot_data$MALDI.device_full %in% c("device_36","device_37","device_39","device_40"),]
#add which MALDI type have been used
plot_data$MALDI.type<-ifelse(plot_data$MALDI.device_full %in% VitekMS.devices, 'VitekMS', 
                             ifelse(plot_data$MALDI.device_full %in% bruker.devices,'MicroflexBruker', NA))

# for calculating the stats on fractions of peaks which could be reproduces between two technical replica, set value to NA for device_30 & device_23 (only one technical replicate acquired)
plot_data$value <- ifelse(plot_data$MALDI.device_full %in% c('device_23','device_30') & plot_data$Endpoint == "frac_peaks_repr_2spectra", NA, plot_data$value)
# add overall quantiles to plot
plot_data <- plot_data %>% group_by(Endpoint) %>%
  mutate(median = median(value,na.rm = TRUE), lowerIQR = quantile(value, probs = 0.25,na.rm = TRUE)[[1]], upperIQR = quantile(value, probs = 0.75,na.rm = TRUE)[[1]]) %>% ungroup()
# convert to dataframe
plot_data<-as.data.frame(plot_data)
# set qantiles to NA for species identification parts of the plot
plot_data$median <- ifelse(plot_data$Endpoint %in% c('marker', 'bruker', 'VitekMS'), NA, plot_data$median)
plot_data$lowerIQR <- ifelse(plot_data$Endpoint %in% c('marker', 'bruker', 'VitekMS'), NA, plot_data$lowerIQR)
plot_data$upperIQR <- ifelse(plot_data$Endpoint %in% c('marker', 'bruker', 'VitekMS'), NA, plot_data$upperIQR)
# calculate median per device and order according to the median of ribosomal subunuts detected per lab
plot_data<- plot_data %>% group_by(Endpoint, device_id) %>% mutate(median_per_lab = median(value, na.rm = T), lowerIQR_per_lab = quantile(value, probs = 0.25,na.rm = TRUE)[[1]], upperIQR_per_lab = quantile(value, probs = 0.75,na.rm = TRUE)[[1]]) %>% ungroup() 
plot_data<- plot_data %>% group_by(Endpoint, device_id) %>% mutate(order = median(value, na.rm = T)) 
plot_data$order <- ifelse(grepl('n_ribos_detected',plot_data$Endpoint), plot_data$order, NA)

# Plot the relative number of ribosomal subunits detected per group for spectra acquired in all laboratories
plot_data.groups<-plot_data[!plot_data$Endpoint.rel %in% c('marker','median_intensity','total_intensity_peaks'),]
# remove lines which are missing a Group assignment (ID)
plot_data.groups<-plot_data.groups[!is.na(plot_data.groups$Group),]
plot_data.groups<-plot_data.groups[!is.na(plot_data.groups$Endpoint.rel),]

# Define labels for endpoints
plot_data.groups$Endpoint.rel<-factor(plot_data.groups$Endpoint.rel, levels =c("rel_amount_su_detected","frac_peaks_repr_2spectra"),
                                      labels = str_wrap(c("% marker masses", "Reproducibility"), width = 15))
# define how the bacterial groups should be labelled
group_labels_strep_together <- c("Enterobacteriaceae","Burkholderia", "Bordetella","Gram negative Anaerobes", 
                                 "Listeria", "Streptococci", "Staphylococcus", "Actinobacteria")
group_labels_strep_together <- str_wrap(group_labels_strep_together, width = 15)
# add 'Streptococci' groups to one
plot_data.groups$Group<-ifelse(grepl('Streptococci', plot_data.groups$Group), 'Streptococci', as.character(plot_data.groups$Group))
# define order in the plot
plot_data.groups$Group<-factor(plot_data.groups$Group, levels = c("Enterobacteriaceae",  "Burkholderia", "Bordetella", "Anaerob_Gram_negative", 
                                                                  "Listeria","Streptococci", "Staphylococcus", "Actinobacteria"),
                               labels = group_labels_strep_together)

# plot
group <- ggplot(plot_data.groups, aes(x=Group, y=as.numeric(as.character(value.rel)), ymin=0,ymax=value))+facet_grid(Endpoint.rel ~ .)  + ylim(0,1) +
  ylab('') +xlab('') + theme(axis.text.x = element_text(angle = 60,  size = 7, hjust = 1), axis.text.y = element_text(size = 7), strip.text = element_text(size=9))  + geom_boxplot(lwd=0.2, outlier.size = 0.5) 
group

# export
pdf('./participating_labs_overall_groups.pdf', height = 3.5, width = 2)
group
dev.off()

#plot the number of ribosomal subunits per lab
# for this exclude the identification by Vitek and microflex, only consider marker ID as this includes all species and is independant of the device
plot_data_n_ribos<-plot_data[!grepl('bruker|VitekMS', plot_data$eval) & !grepl('intensity', plot_data$Endpoint),]
plot_data_n_ribos<-plot_data_n_ribos[!is.na(plot_data_n_ribos$MALDI.device_full),]

# define how the endpoints should be labelled
labels_endpoints<-c("# marker masses", "Reproducibility", "Marker based ID")
# introduce line breaks if necessary
labels_endpoints<-str_wrap(labels_endpoints, width = 13)
# define the order of the endpoints
plot_data_n_ribos$Endpoint<- factor(plot_data_n_ribos$Endpoint, levels = c("n_ribos_detected","frac_peaks_repr_2spectra", "marker"), 
                                   labels = labels_endpoints)

# Exclude reproducibility
plot_data_n_ribos <- plot_data_n_ribos[plot_data_n_ribos$Endpoint %in% c("# marker\nmasses","Marker based\nID"),]
plot_data_n_ribos$Endpoint<- factor(plot_data_n_ribos$Endpoint, levels = c("# marker\nmasses","Marker based\nID"))

# Add which MALDI-types have been used and colozr accordingly
plot_data_n_ribos['MALDI.device.type']<- ifelse(plot_data_n_ribos$MALDI.device_full %in% bruker.devices, 'MBT Biotyper', 'VitekMS / Axima Confidence')
plot_data_n_ribos['MALDI.device.type'] <- factor(plot_data_n_ribos$MALDI.device.type, levels = c('MBT Biotyper', 'VitekMS / Axima Confidence'))

#plot
f1 <- ggplot(plot_data_n_ribos, aes(x=reorder(device_id, order, na.rm = TRUE), y=as.numeric(as.character(value)), ymin=0,ymax=value))+facet_grid(Endpoint ~ ., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('Device ID')  
f2 <- f1+geom_boxplot(data = plot_data_n_ribos[grepl('.*marker.*masses*',plot_data_n_ribos$Endpoint),], lwd=0.2, outlier.size = 0.5, aes(col = MALDI.device.type)) +  scale_color_manual(values=c("goldenrod3","turquoise")) #+ geom_point(data=ylim_points,x=NA)
f5 <- f2+geom_col(data = plot_data_n_ribos[grepl('.*Marker.*ID.*',plot_data_n_ribos$Endpoint),], aes(fill=eval), width = 0.7)
f7 <- f5+scale_fill_manual('Identification',
                                  breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker","prop.wrong.marker.single" ),
                                  values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
  
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('') +   geom_hline(aes(yintercept = median, group = Endpoint), linetype="solid", color = "red", size=0.4) +
  geom_hline(aes(yintercept = lowerIQR, group = Endpoint), linetype="dashed", color = "red", size=0.4) +
  geom_hline(aes(yintercept = upperIQR, group = Endpoint), linetype="dashed", color = "red", size=0.4) 

f7 <- f7  + theme(axis.text.x = element_text(angle = 60,  size = 7, hjust = 1), axis.text.y = element_text(size = 7), strip.text = element_text(size=3)) 
f7 

#export
pdf('./overall_labs_quality.pdf', width=7, height = 4)
f7 + theme(legend.position = 'none')
dev.off()

# assess how the measurement error has changed between the different sample preperation protocols
# convert to numeric
sum$mean_dist_ppm <- as.numeric(as.character(sum$mean_dist_ppm))

# exclude devices which have not acquired spectra using 'method 1'
sum_err <- sum %>%
  group_by(device_id) %>%
  mutate(check=any(grepl('1', method))) %>% 
  filter(check == TRUE) %>% 
  select(-check)

# order devices by error in acquired in the baseline quality assessment 
device_id_order_ppm <- sum_err[sum_err$method == '0', c('mean_dist_ppm', 'device_id')]
device_id_order_ppm <- reorder(device_id_order_ppm$device_id, device_id_order_ppm$mean_dist_ppm,na.rm = TRUE, FUN = median)
sum_err$device_id <- factor(sum_err$device_id, levels = levels(device_id_order_ppm))
# plot. INclude 500ppm as line as this is often used as cut-off for marker based identification tools, eg in the bruker subtyping module
error <- ggplot(sum_err[sum_err$method %in% c(0,1) ,], aes(x=device_id, y=mean_dist_ppm, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  geom_hline(yintercept=500) +
  ylab('Mean measurement error [ppm]') +
  xlab('Device') +
  scale_fill_manual(values=c('white', 'blue'))
error

# export
pdf('./Round2_measurement_error.pdf', height= 2.8, width = 12)
error
dev.off()

# define endpoints
endpoints <- c("n_peaks","highest_peaks","mass_90" ,"n_high_peaks","total_intensity_peaks",               
               "frac_peaks_repr","frac_peaks_repr_2spectra",     
               "n_ribos_detected","n_predicted_su","rel_amount_su_detected"  ,
               "mean_dist", "mean_dist_ppm","mean_intensity" ,         
               "median_intensity","max_mass_ribo")

# check in how many labs the number of ribos is significantly higher / lower with methods 1 and 2 compared to method 0
# in order to do so, check that there is the same amount of samples per lab / strain / method
check_n <- as.data.frame(table(sum$strainnumber, sum$MALDI.device, sum$method))
colnames(check_n) <- c('strainnumber', 'device', 'method', 'n_spectra')
#correct for typo
names(sum)[names(sum) == 'meadian_intensity'] <- 'median_intensity'
# convert all endpoints to numeric
sum <- sum %>% mutate_at(endpoints, ~as.numeric(as.character(.))) %>%
  mutate_at(endpoints[!grepl('dist', endpoints)], ~replace(., is.na(.), 0)) # if no ribos are found, some endpoints are NA. replace these with 0, except for measurement error, which is truly NA

# convert to long format
sum_long <- pivot_longer(
  sum,
  all_of(endpoints),
  names_to = "Endpoint")

# remove the rows where no ribos are detected and the measurement error can therefor not be assessed (only the rows referrring to the measrurement error are removed, not the ones referring to the number of ribosomal subunits detected)
sum_long <- sum_long[!is.na(sum_long$value), ]

# choose methods 0 an 1 to compare
sum_long_0_1 <- sum_long[sum_long$method %in% c(0,1),]

# chose replicates which are present in both datasets
sum_long_0_1['spectra_in_01'] <- paste(sum_long_0_1$MALDI.device, sum_long_0_1$strainnumber, sum_long_0_1$replicate_eff, sum_long_0_1$Endpoint)
spectra_0_1 <- intersect(sum_long_0_1[sum_long_0_1$method == '0', ]$spectra_in_01, sum_long_0_1[sum_long_0_1$method == '1', ]$spectra_in_01)
sum_long_0_1 <- sum_long_0_1[sum_long_0_1$spectra_in_01 %in% spectra_0_1,]

# convert method to factors for it not to be interpreted as number
sum_long_0_1$method <- as.factor(sum_long_0_1$method)

# calculate median per divice, method and endpoint
sum_long_0_1_median <- sum_long_0_1 %>%   
  group_by(method,device_id, Endpoint) %>%
  summarise(median = median(value)) %>%
  pivot_wider(names_from = method, values_from = median, names_prefix = "median_method_")

# use a paired wicoxon rank test to compare between the methids (same strains measured)
stats_0_1 <- sum_long_0_1 %>% 
  group_by(device_id, Endpoint) %>%
  rstatix::wilcox_test(value ~ method, paired = T) %>%
  add_significance() %>% 
  add_xy_position(x = "device_id", dodge = 0.8, scales = "free")

# add median values to check which are sign increased or decreased
stats_0_1 <- merge(stats_0_1, sum_long_0_1_median, by = c('device_id', 'Endpoint'), all.x = T)
stats_0_1 <- stats_0_1 %>% 
  mutate(median_diff = (median_method_1 - median_method_0))

# check in how many labs the number of ribos increased 
stats_0_1_n_ribos <- stats_0_1[stats_0_1$Endpoint == "n_ribos_detected",]
stats_0_1_n_ribos['increased_1']<-ifelse(stats_0_1_n_ribos$p< 0.05 & stats_0_1_n_ribos$median_diff > 0, 'increase', 
                                         ifelse(stats_0_1_n_ribos$p< 0.05 & stats_0_1_n_ribos$median_diff < 0, 'decrease', 'no significant change'))
table(stats_0_1_n_ribos['increased_1'])
table(sort(stats_0_1_n_ribos[stats_0_1_n_ribos$median_method_0 <15,]$increased_1))

# check reproducibility
stats_0_1_frac2sp <- stats_0_1[stats_0_1$Endpoint == "frac_peaks_repr_2spectra",]
stats_0_1_frac2sp['increased_1']<-ifelse(stats_0_1_frac2sp$p< 0.05 & stats_0_1_frac2sp$median_diff > 0, 'increase', 
                                         ifelse(stats_0_1_frac2sp$p< 0.05 & stats_0_1_frac2sp$median_diff < 0, 'decrease', 'no significant change'))
table(stats_0_1_frac2sp['increased_1'])

# check measurement error
stats_0_1_n_mean_dist_ppm <- stats_0_1[stats_0_1$Endpoint == "mean_dist_ppm",]
stats_0_1_n_mean_dist_ppm['increased_1']<-ifelse(stats_0_1_n_mean_dist_ppm$p< 0.05 & stats_0_1_n_mean_dist_ppm$median_diff > 0, 'increase', 
                                                ifelse(stats_0_1_n_mean_dist_ppm$p< 0.05 & stats_0_1_n_mean_dist_ppm$median_diff < 0, 'decrease', 'no significant change'))
table(stats_0_1_n_mean_dist_ppm['increased_1'])
stats_0_1_n_mean_dist_ppm['y.position'] <- 800

# check total intensity
stats_0_1_n_tot_int <- stats_0_1[stats_0_1$Endpoint == "total_intensity_peaks",]
stats_0_1_n_tot_int['increased_1']<-ifelse(stats_0_1_n_tot_int$p< 0.05 & stats_0_1_n_tot_int$median_diff > 0, 'increase', 
                                                 ifelse(stats_0_1_n_tot_int$p< 0.05 & stats_0_1_n_tot_int$median_diff < 0, 'decrease', 'no significant change'))
table(stats_0_1_n_tot_int['increased_1'])

# check median int ribos
stats_0_1_median_int_ribos <- stats_0_1[stats_0_1$Endpoint == "median_intensity",]
stats_0_1_median_int_ribos['increased_1']<-ifelse(stats_0_1_median_int_ribos$p< 0.05 & stats_0_1_median_int_ribos$median_diff > 0, 'increase', 
                                           ifelse(stats_0_1_median_int_ribos$p< 0.05 & stats_0_1_median_int_ribos$median_diff < 0, 'decrease', 'no significant change'))
table(stats_0_1_median_int_ribos['increased_1'])

## plot measurement error incl.  p-values
# order by the median measurement error recorded at baseline quality assessment
device_id_order_ppm <- sum_err[sum_err$method == '0', c('mean_dist_ppm', 'device_id')]
device_id_order_ppm <- reorder(device_id_order_ppm$device_id, device_id_order_ppm$mean_dist_ppm,na.rm = TRUE, FUN = median)
sum_long_0_1_err <- sum_long_0_1[sum_long_0_1$Endpoint == 'mean_dist_ppm',]
sum_long_0_1_err$device_id <- factor(sum_long_0_1_err$device_id, levels = levels(device_id_order_ppm))

# plot, colour by method
error <- ggplot(sum_long_0_1_err, aes(x=device_id, y=value)) + 
  geom_boxplot(aes(fill = method)) + 
  geom_hline(yintercept=500) +
  ylab('Mean measurement error [ppm]') +
  geom_hline(yintercept=500) +
  ylab('Mean measurement error [ppm]') +
  xlab('Device') +
  scale_fill_manual(values=c('white', 'blue')) + theme(legend.position = 'bottom')

# format df to plot p-values
stats_0_1_error <- stats_0_1[stats_0_1$Endpoint == 'mean_dist_ppm',]
order_df <- data.frame(device_id = levels(device_id_order_ppm), order = 1:length(levels(device_id_order_ppm)))
order_df_ribos <- order_df
stats_0_1_error <- merge(stats_0_1_error, order_df, by= 'device_id')
stats_0_1_error$x <- stats_0_1_error$order
stats_0_1_error['xmin'] <- stats_0_1_error$x - 0.2
stats_0_1_error['xmax'] <- stats_0_1_error$x + 0.2

stats_0_1_error$y.position <- 700
stats_0_1_error['ymin'] <- 680
stats_0_1_error['ymax'] <- 720
stats_0_1_error['p.formated']<- p_format(p_round(stats_0_1_error$p), digits = 2)

# add p-values to plot
error <- error + 
  stat_pvalue_manual(stats_0_1_error, label = "{p.formated}\n{p.signif}",  bracket.nudge.y = 1, label.size = 2) 

# export
pdf('./Round2_measurement_error.pdf', height= 3.5, width = 12)
error
dev.off()

# Plot the number of detected subunits between methid 0 and method 1
# define order of devices, order by the median number of ribos dected in baseline quality assessment, excluding devices for which no spectra from mathod 1 are available
device_id_order_n_ribos <- sum[sum$method == '0' & !(sum$device_id %in% c("8", "7", "16","29","31")),c('n_ribos_detected', 'device_id')]
device_id_order_n_ribos <- reorder(device_id_order_n_ribos$device_id, device_id_order_n_ribos$n_ribos_detected,na.rm = TRUE, FUN = median)

# choose the endpoint
sum_long_0_1_n_ribos <- sum_long_0_1[sum_long_0_1$Endpoint == 'n_ribos_detected',]
# convert factors to characters
sum_long_0_1_n_ribos <- sum_long_0_1_n_ribos %>% mutate_if(is.factor, as.character)

# order as previously defined
sum_long_0_1_n_ribos$device_id <- factor(sum_long_0_1_n_ribos$device_id, levels = levels(device_id_order_n_ribos))
sum_long_0_1_n_ribos<-sum_long_0_1_n_ribos[!is.na(sum_long_0_1_n_ribos$device_id),]

# plot 
ribos_0_1 <- ggplot(sum_long_0_1_n_ribos, aes(x=device_id, y=value, drop=FALSE)) + 
  geom_boxplot(aes(fill = method)) + 
  ylab('Number of ribosomal marker masses detected') +
  xlab('Device') +
  scale_fill_manual(values = c('white', 'blue'), drop=FALSE) + 
  ylim(0,38) +
  theme(legend.position = 'bottom')

# format df to add p-value
stats_0_1_nribos <- stats_0_1[stats_0_1$Endpoint == 'n_ribos_detected',]
order_df <- data.frame(device_id = levels(device_id_order_n_ribos), order = 1:length(levels(device_id_order_n_ribos)))
stats_0_1_nribos <- merge(stats_0_1_nribos, order_df, by= 'device_id')
stats_0_1_nribos$x <- stats_0_1_nribos$order
stats_0_1_nribos['xmin'] <- stats_0_1_nribos$x - 0.2
stats_0_1_nribos['xmax'] <- stats_0_1_nribos$x + 0.2
stats_0_1_nribos['p.formated']<- p_format(p_round(stats_0_1_nribos$p), digits = 2)
stats_0_1_nribos$y.position <- 32
# add p-value to plot
ribos_0_1 <- ribos_0_1 + 
  stat_pvalue_manual(stats_0_1_nribos, label = "{p.formated}\n{p.signif}", bracket.nudge.y = 1, label.size = 2) 

# export
pdf('./Round2_n_ribos_01.pdf', height= 3.5, width = 12)
ribos_0_1
dev.off()

# check reproducibility
sum_long_0_1[sum_long_0_1$Endpoint == 'frac_peaks_repr_2spectra', ] %>% 
  group_by(method) %>%
  summarise_at(vars(value),
               list(min=min, Q1=~quantile(., probs = 0.25, na.rm = T),
                    median=median, Q3=~quantile(., probs = 0.75, na.rm = T),
                    max=max)) # overall no change can be observed

wilcox.test(sum_long_0_1[sum_long_0_1$Endpoint == 'frac_peaks_repr_2spectra' & sum_long_0_1$method == '0', ]$value, 
            sum_long_0_1[sum_long_0_1$Endpoint == 'frac_peaks_repr_2spectra' & sum_long_0_1$method == '1', ]$value, 
            paired = F) # overall no change can be observed


##### compare method 0 to 2
# choose methods 0 and 2
sum_long_0_2 <- sum_long[sum_long$method %in% c(0,2),]

# chose replicates which are present in both datasets
sum_long_0_2['spectra_in_02'] <- paste(sum_long_0_2$MALDI.device, sum_long_0_2$strainnumber, sum_long_0_2$replicate_eff, sum_long_0_2$Endpoint)
spectra_0_2 <- intersect(sum_long_0_2[sum_long_0_2$method == '0', ]$spectra_in_02, sum_long_0_2[sum_long_0_2$method == '2', ]$spectra_in_02)
sum_long_0_2 <- sum_long_0_2[sum_long_0_2$spectra_in_02 %in% spectra_0_2,]

# remove devices for which there are no spectra available acquired with method 2
sum_long_0_2 <- sum_long_0_2[!sum_long_0_2$device_id %in% c("44", "41", "31", "29", "16", "07", "08"),]

# convert method to factors for it not to be interpreted as number
sum_long_0_2$method <- as.factor(sum_long_0_2$method)

# compare the number of ribosomal subunits detected between method 0 and 2
sum_long_0_2_n_ribos <- sum_long_0_2[sum_long_0_2$Endpoint == 'n_ribos_detected',]
# convert all factors to character
sum_long_0_2_n_ribos <- sum_long_0_2_n_ribos %>% mutate_if(is.factor, as.character)

# order by the number of ribosomal subunits detected at baseline quality assessment
sum_long_0_2_n_ribos$device_id <- factor(sum_long_0_2_n_ribos$device_id, levels = levels(device_id_order_n_ribos))

# plot
ribos_0_2 <- ggplot(sum_long_0_2_n_ribos, aes(x=device_id, y=value, drop=FALSE)) + 
  geom_boxplot(aes(fill = method)) + 
  ylab('Number of ribosomal marker masses detected') +
  xlab('Device') +
  scale_fill_manual(values = c('white', 'orange')) + 
  ylim(0,38) +
  theme(legend.position = 'bottom')

# calculate stats (use wilcoxon rank test paired as here the same samples are compared)
stats_0_2 <- sum_long_0_2 %>% 
  group_by(device_id, Endpoint) %>%
  rstatix::wilcox_test(value ~ method, paired = T) %>%
  add_significance() %>% 
  add_xy_position(x = "device_id", dodge = 0.8, scales = "free")

# compare medians in order to assess whether there was an increase or a decrease
sum_long_0_2_median <- sum_long_0_2 %>%   
  group_by(method,device_id, Endpoint) %>%
  summarise(median = median(value)) %>%
  pivot_wider(names_from = method, values_from = median, names_prefix = "median_method_")

# add median values to check which are sign increased or decreased
stats_0_2 <- merge(stats_0_2, sum_long_0_2_median, by = c('device_id', 'Endpoint'), all.x = T)
stats_0_2 <- stats_0_2 %>% 
  mutate(median_diff = (median_method_2 - median_method_0))

# assess the number of ribosomal subunits detected
stats_0_2_n_ribos <- stats_0_2[stats_0_2$Endpoint == 'n_ribos_detected',]
# check how many have increased
stats_0_2_n_ribos['increased_2']<-ifelse(stats_0_2_n_ribos$p< 0.05 & stats_0_2_n_ribos$median_diff > 0, 'increase', 
                                         ifelse(stats_0_2_n_ribos$p< 0.05 & stats_0_2_n_ribos$median_diff < 0, 'decrease', 'no significant change'))
table(stats_0_2_n_ribos['increased_2'])

# check reproducibility
stats_0_2_frac2sp <- stats_0_2[stats_0_2$Endpoint == "frac_peaks_repr_2spectra",]
stats_0_2_frac2sp['increased_2']<-ifelse(stats_0_2_frac2sp$p< 0.05 & stats_0_2_frac2sp$median_diff > 0, 'increase', 
                                         ifelse(stats_0_2_frac2sp$p< 0.05 & stats_0_2_frac2sp$median_diff < 0, 'decrease', 'no significant change'))
table(stats_0_2_frac2sp['increased_2'])

# define 'order_df' and format df to add p-value on plot
order_df <- data.frame(device_id = levels(device_id_order_n_ribos), order = 1:length(levels(device_id_order_n_ribos)))
stats_0_2_n_ribos <- merge(stats_0_2_n_ribos, order_df, by= 'device_id')
stats_0_2_n_ribos$x <- stats_0_2_n_ribos$order
stats_0_2_n_ribos['xmin'] <- stats_0_2_n_ribos$x - 0.2
stats_0_2_n_ribos['xmax'] <- stats_0_2_n_ribos$x + 0.2
stats_0_2_n_ribos['p.formated']<- p_format(p_round(stats_0_2_n_ribos$p), digits = 2)
stats_0_2_n_ribos$y.position <- 32
# add p-value
ribos_0_2 <- ribos_0_2 + 
  stat_pvalue_manual(stats_0_2_n_ribos, label = "{p.formated}\n{p.signif}", vjust = -1, bracket.nudge.y = 1, label.size = 2) 

# export
pdf('./Round2_n_ribos_02.pdf', height= 3.5, width = 12)
ribos_0_2
dev.off()

## compare methods 1 to 2
sum_long <- sum_long[!is.na(sum_long$value), ]
sum_long_1_2 <- sum_long[sum_long$method %in% c(1,2),]

# chose replicates which are present in both datasets
sum_long_1_2['spectra_in_12'] <- paste(sum_long_1_2$MALDI.device, sum_long_1_2$strainnumber, sum_long_1_2$replicate_eff, sum_long_1_2$Endpoint)
spectra_1_2 <- intersect(sum_long_1_2[sum_long_1_2$method == '1', ]$spectra_in_12, sum_long_1_2[sum_long_1_2$method == '2', ]$spectra_in_12)

sum_long_1_2 <- sum_long_1_2[sum_long_1_2$spectra_in_12 %in% spectra_1_2,]
# convert method to factors for it not to be interpreted as number
sum_long_1_2$method <- as.factor(sum_long_1_2$method)
# remove the devices which are not present in both datasets 
sum_long_1_2 <- sum_long_1_2[!(sum_long_1_2$device_id %in%c('44', '41','31', '29', '16', '7', '8')),]

# compare the number of detected marker masses between the methods 0 and 1
sum_long_1_2_n_ribos <- sum_long_1_2[sum_long_1_2$Endpoint == 'n_ribos_detected',]
sum_long_1_2_n_ribos$device_id <- factor(sum_long_1_2_n_ribos$device_id, levels = levels(device_id_order_n_ribos))

# plot
ribos_1_2 <- ggplot(sum_long_1_2_n_ribos, aes(x=device_id, y=value, drop=FALSE)) + 
  geom_boxplot(aes(fill = method)) + 
  ylab('Number of ribosomal marker masses detected') +
  xlab('Device') +
  scale_fill_manual(values = c('white', 'orange')) + 
  ylim(0,38) +
  theme(legend.position = 'bottom')

# stats
stats_1_2 <- sum_long_1_2 %>% 
  group_by(device_id, Endpoint) %>%
  rstatix::wilcox_test(value ~ method, paired = T) %>%
  add_significance("p") %>% 
  add_xy_position(x = "device_id", dodge = 0.8, scales = "free")

# claclulate media
sum_long_1_2_median <- sum_long_1_2 %>%   
  group_by(method,device_id, Endpoint) %>%
  summarise(median = median(value)) %>%
  pivot_wider(names_from = method, values_from = median, names_prefix = "median_method_")

# add median values to check which are sign increased or decreased
stats_1_2 <- merge(stats_1_2, sum_long_1_2_median, by = c('device_id', 'Endpoint'), all.x = T)
stats_1_2 <- stats_1_2 %>% 
  mutate(median_diff = (median_method_2 - median_method_1))

# choose endpoint
stats_1_2_n_ribos <- stats_1_2[stats_1_2$Endpoint == 'n_ribos_detected',]
stats_1_2_n_ribos['increased_2']<-ifelse(stats_1_2_n_ribos$p< 0.05 & stats_1_2_n_ribos$median_diff > 0, 'increase', 
                                         ifelse(stats_1_2_n_ribos$p< 0.05 & stats_1_2_n_ribos$median_diff < 0, 'decrease', 'no significant change'))
table(stats_1_2_n_ribos['increased_2'])

# check reproducibility
stats_1_2_frac2sp <- stats_1_2[stats_1_2$Endpoint == "frac_peaks_repr_2spectra",]
stats_1_2_frac2sp['increased_2']<-ifelse(stats_1_2_frac2sp$p< 0.05 & stats_1_2_frac2sp$median_diff > 0, 'increase', 
                                         ifelse(stats_1_2_frac2sp$p< 0.05 & stats_1_2_frac2sp$median_diff < 0, 'decrease', 'no significant change'))
table(stats_1_2_frac2sp['increased_2'])

## compare the effects of each sample preperation per group
# calculate the stats between method 0 and 1
stats_0_1_group <- sum_long_0_1 %>% 
  group_by(Group, Endpoint) %>%
  rstatix::wilcox_test(value ~ method, paired = T) %>%
  add_significance() %>% 
  add_xy_position(x = "Group", dodge = 0.8, scales = "free")

# calculate the stats between method 0 and 2
stats_0_2_group <- sum_long_0_2 %>% 
  group_by(Group, Endpoint) %>%
  rstatix::wilcox_test(value ~ method, paired = T) %>%
  add_significance() %>% 
  add_xy_position(x = "Group", dodge = 0.8, scales = "free")

# and between method 1 and 2
stats_1_2_group <- sum_long_1_2 %>% 
  group_by(Group, Endpoint) %>%
  rstatix::wilcox_test(value ~ method, paired = T) %>%
  add_significance() %>% 
  add_xy_position(x = "Group", dodge = 0.8, scales = "free")

# format the df and stats to plot comparison 0 and 1
stats_0_1_group_nribos <- stats_0_1_group[stats_0_1_group$Endpoint == 'n_ribos_detected',]
stats_0_1_group_nribos['p.formated']<- p_format(p_round(stats_0_1_group_nribos$p), digits = 2)
stats_0_1_group_nribos$y.position <- 32
sum_long_0_1_n_ribos <- sum_long_0_1_n_ribos[!is.na(sum_long_0_1_n_ribos$Group),]
# plot
ribos_0_1_group <- ggplot(sum_long_0_1_n_ribos, aes(x=Group, y=value)) + 
  geom_boxplot(aes(fill = method)) + 
  ylab('Number of ribosomal marker masses detected') +
  xlab('') +
  scale_fill_manual(values = c('white', 'blue')) + 
  ylim(0,38) +
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# add p-values
ribos_0_1_group <- ribos_0_1_group + 
  stat_pvalue_manual(stats_0_1_group_nribos, label = "{p.formated}\n{p.signif}", vjust = -1, bracket.nudge.y = 1, label.size = 2) 
ribos_0_1_group 


# format the df and stats to plot comparison 0 and 2
stats_0_2_group_nribos <- stats_0_2_group[stats_0_2_group$Endpoint == 'n_ribos_detected',]
stats_0_2_group_nribos['p.formated']<- p_format(p_round(stats_0_2_group_nribos$p), digits = 2)
stats_0_2_group_nribos$y.position <- 32
sum_long_0_2_n_ribos <- sum_long_0_2_n_ribos[!is.na(sum_long_0_2_n_ribos$Group),]
# plot
ribos_0_2_group <- ggplot(sum_long_0_2_n_ribos, aes(x=Group, y=value)) + 
  geom_boxplot(aes(fill = method)) + 
  ylab('Number of ribosomal marker masses detected') +
  xlab('') +
  scale_fill_manual(values = c('white', 'orange')) + 
  ylim(0,38) +
  theme(axis.text.x = element_text(angle = 60, hjust=1))
# add p-value
ribos_0_2_group <- ribos_0_2_group + 
  stat_pvalue_manual(stats_0_2_group_nribos, label = "{p.formated}\n{p.signif}", vjust = -1, bracket.nudge.y = 1, label.size = 2) 
 
####### compare 0 to 1 to 2 also including species identification
# set to 'not ID' when under certain threshold
sum['mean_intensity_peaks']<-sum$total_intensity_peaks / sum$n_peaks

# split up streptococcus in viridans and non-viridans
sum['Group']<-as.factor(ifelse(grepl('infantis|pseudopneumoniae|pneumoniae|gordonii', as.character(sum$Strain)) & sum$genus_NGS == "Streptococcus", 'Viridans Streptococci', 
                                     ifelse(grepl('equinus|gallolyticus|lutetiensis|dysgalactiae', as.character(sum$Strain)) & sum$genus_NGS == "Streptococcus", 'Other Streptococci', as.character(sum$Group))))

# gather to long format
sum_long <- pivot_longer(
  sum,
  endpoints,
  names_to = "Endpoint")

# remove the ones were the value is 'NA' as no ribosomal markers could be detected (then, eg., the mean_dist cannot be assessed)
sum_long <- sum_long[!is.na(sum_long$value), ]
sum_long_0_1_2 <- sum_long

# chose replicates which are present in all three datasets
sum_long_0_1_2['spectra_in_012'] <- paste(sum_long_0_1_2$MALDI.device, sum_long_0_1_2$strainnumber, sum_long_0_1_2$replicate_eff, sum_long_0_1_2$Endpoint)
sum_long_0_1_2 <-  sum_long_0_1_2[sum_long_0_1_2$spectra_in_012 %in% intersect(spectra_0_1, spectra_0_2) & sum_long_0_1_2$method %in% c(0,1,2),]

# calculate quantiles per group
sum_per_group <- sum_long_0_1_2[sum_long_0_1_2$Endpoint %in% c('n_ribos_detected', 'frac_peaks_repr_2spectra'), ] %>% group_by(method, Group, Endpoint) %>%
  summarise(lower_iqr = quantile(value, probs = 0.25), 
            median = quantile(value, probs = 0.5), 
            higher_iqr = quantile(value, probs = 0.75))

# convert method to factors for it not to be interpreted as number
sum_long_0_1_2$method <- factor(sum_long_0_1_2$method, levels = c(0,1,2))

# choose for the number of marker mass endpoint
sum_long_0_1_2_n_ribos <- sum_long_0_1_2[sum_long_0_1_2$Endpoint == 'n_ribos_detected',]

### add proportion of correctly identified by device, group and method

# add proportions
eval2 <- sum_long_0_1_2 %>%
  group_by(MALDI.device_full, method, Group) %>%
  dplyr::summarize(prop.correct.marker.single = sum(grepl('TRUE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),             
                   prop.wrong.marker.single = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('FALSE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   prop.correct.marker.multi = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
                   prop.noID.marker = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('no ID', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified))) 

eval<-merge(sum, eval2, by = intersect(colnames(sum), colnames(eval2)), all.x  = T)
# correct for typo
colnames(eval)[colnames(eval) == 'median_intensity'] <- 'median_intensity'
# convert to long format
eval_prop_long<-gather(eval2, "classification_eval_prop", "value", -c(MALDI.device_full, method, Group))
#Gather all endpoints to one column in order to plot all at once using facet_wrap
eval_long<-gather(eval[,c('max_mass_ribo', 'spectra',"n_peaks", "mean_intensity_peaks","highest_peaks", "mass_90", "n_high_peaks", "frac_peaks_repr_2spectra", "n_ribos_detected", "n_predicted_su", "rel_amount_su_detected", "mean_dist", "mean_dist_ppm", "mean_intensity","median_intensity", "total_intensity_peaks","brukerDB.Score", "brukerDB.n_species", "brukerDB.n_genera", "brukerDB.n_species_over_2", "brukerDB.n_genera_over_2", "brukerDB.diff", "brukerDB.correct_genus_identified", "brukerDB.correct_species_identified",  "markerDB.match_count",  "markerDB.n_genera_identified",  "markerDB.correct_genus_identified", "markerDB.correct_species_identified", "markerDB.correct_single_species")], "Endpoint", "value", -spectra)
eval_long<-merge(eval[,c("spectra","Strain", "strainnumber", "Group", "spectra", "species_NGS", "genus_NGS", 'MALDI.device_full', 'device_id', 'method')], eval_long, by='spectra')
eval_long['value']<-ifelse(eval_long$value == 'TRUE', 1,
                           ifelse(eval_long$value == 'FALSE', 0, eval_long$value))
eval_long["value"]<-gsub(',', '.', eval_long$value)
eval_long["value"]<-as.numeric(as.character(eval_long$value)) 

# all identifications evaluated here were identified using the marker based database
eval_prop_long['DB']<-'marker'
eval_prop_long['eval']<-eval_prop_long$classification_eval_prop

#for overall figure, exclusively include marker ID
eval_prop_long<-eval_prop_long[eval_prop_long$DB == 'marker',]
eval_prop_long$eval<- factor(eval_prop_long$eval, levels = c("prop.wrong.marker.single", "prop.noID.marker", 
                                                             "prop.correct.marker.multi", 
                                                             "prop.correct.marker.single"))

# select for the endpoints to plot
plot_data<-eval_long[eval_long$Endpoint %in% c("n_ribos_detected","median_intensity","total_intensity_peaks", "frac_peaks_repr_2spectra"),]

# rename DB to endpoint
names(eval_prop_long)[names(eval_prop_long) == 'DB'] <- 'Endpoint'
eval_prop_long$classification_eval_prop<-NULL

# add device id first
eval_prop_long<-merge(plot_data[!duplicated(plot_data$MALDI.device_full),c("MALDI.device_full", "device_id")], eval_prop_long, by = "MALDI.device_full", all.y  = TRUE)
plot_data<-merge(plot_data, eval_prop_long, by = intersect(colnames(eval_prop_long), colnames(plot_data)), all = TRUE)

# define how endpoints should be labelled
labels_endpoints<-c("# marker\nmasses", "Relative int. marker masses (log10)", "Total int. peaks (log10)", "Reproducibility", "Marker based ID")
labels_endpoints<-str_wrap(labels_endpoints, width = 13)
# relabel
plot_data$Endpoint<- factor(plot_data$Endpoint, levels = c("n_ribos_detected", "median_intensity", "total_intensity_peaks", "frac_peaks_repr_2spectra", "marker"), 
                                   labels = labels_endpoints)
# convert device id to factors
plot_data$device_id <- as.factor(plot_data$device_id)
# only plot the number of marker masses, the reproducibility and the the marker based species ID
plot_data_short <- plot_data[grepl('\\#.*marker.*|.*eproducibility|.*Marker.*ID.*',plot_data$Endpoint),]

# do only display the laboratories who have acquired spectra with methods 0, 1 and 2 
plot_data_short <- plot_data_short[plot_data_short$device_id %in% sum_long_0_1_2$device_id & plot_data_short$method %in% c(0,1,2),]
plot_data_short <- merge(plot_data_short, order_df_ribos, by = 'device_id')

# convert the ID column onto one value per device which sums up to 1. (multiply each group value by the number of spectra in each group and devide the overall results through the total number of spectra per device in order to account for the fact that for each group a different amount of spectra was acquired)
plot_data_short <- plot_data_short %>% 
  group_by(device_id, method) %>%
  mutate(n_distinct_groups = n_distinct(Group, na.rm = TRUE),
         n_distinct_spectra_per_device = n_distinct(spectra, na.rm = TRUE)) %>%
  group_by(Group, device_id, method) %>%
  mutate(n_distinct_spectra_per_group = n_distinct(spectra, na.rm = TRUE)) %>%
  group_by(Group, device_id, Endpoint, method) %>%
  mutate(value_g_norm = (as.numeric(as.character(value))*n_distinct_spectra_per_group)/n_distinct_spectra_per_device)

# use this normalised value only for the marker ID, the others are not affected by differet weightig of groups
plot_data_short$value_g <- ifelse(grepl('.*Marker.*ID.*',plot_data_short$Endpoint), plot_data_short$value_g_norm, plot_data_short$value)

## also, only display the devices for which spectra from all groups are available acquired with all three methods
incomplete <- as.data.frame(unique(plot_data_short[plot_data_short$n_distinct_groups != max(plot_data_short$n_distinct_groups),'device_id']))
plot_data_short <- plot_data_short[!plot_data_short$device_id %in% union(incomplete$device_id, c(28,44,41)),]
plot_data_short$device_id <- factor(plot_data_short$device_id, levels = levels(device_id_order_n_ribos)[levels(device_id_order_n_ribos) %in% plot_data_short$device_id])

#plot
f1 <- ggplot(plot_data_short, aes(x=method, y=as.numeric(as.character(value_g)), ymin=0,ymax=value))+facet_grid(Endpoint ~ ., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')  
f2 <- f1+geom_boxplot(data = plot_data_short[grepl('\\#.*marker.*',plot_data_short$Endpoint),], lwd=0.2, outlier.size = 0.5, outlier.stroke = 0.2) 
f4 <- f2+geom_boxplot(data = plot_data_short[grepl('.*eproducibility.*',plot_data_short$Endpoint),], lwd=0.2, outlier.size = 0.5, outlier.stroke = 0.2)
f5 <- f4+geom_col(data = plot_data_short[grepl('.*Marker.*ID.*',plot_data_short$Endpoint),], aes(fill=eval), width = 0.7)
f7<- f5+scale_fill_manual('Identification',
                          breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker","prop.wrong.marker.single" ),
                          values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
  
  
  theme(axis.text.x = element_text(angle = 0,  size = 7), axis.text.y = element_text(size = 7), legend.position = 'bottom') +
  ylab('') +xlab('') + facet_grid(Endpoint~device_id,  scales="free_y") + theme(panel.spacing = unit(0.2, "lines"))

f7 
# export 
pdf('./Round2_sum_fig_per_lab.pdf', width = 7, height = 4.5)
f7
dev.off()

# plot the same but split by group
# convert the ID column onto one value per group which sums up to 1. (multiply each device value by the number of spectra in each device and divide the overall results through the total number of spectra per group in order to account for the fact that for each device a different amount of spectra was acquired)
# remove previous calculations
plot_data_short$value_g_norm <- NULL
plot_data_short$value_g <- NULL
plot_data_short$n_distinct_devices <- NULL
plot_data_short$n_distinct_spectra_per_device_and_group <- NULL
plot_data_short$n_distinct_spectra_per_device <- NULL
plot_data_short$n_distinct_spectra_per_group <- NULL

# calculate weighted sum of proportions
plot_data_short <- plot_data_short %>% 
  group_by(Group, method) %>%
  mutate(n_distinct_devices = n_distinct(device_id, na.rm = TRUE),
         n_distinct_spectra_per_group = n_distinct(spectra, na.rm = TRUE)) %>%
  group_by(device_id, Group, method) %>%
  mutate(n_distinct_spectra_per_device_and_group = n_distinct(spectra, na.rm = TRUE)) %>%
  group_by(Group, device_id, Endpoint, method) %>%
  mutate(value_g_norm = (as.numeric(as.character(value))*n_distinct_spectra_per_device_and_group)/ n_distinct_spectra_per_group)

# use this normalised value only for the marker ID, the others are not affected by differet weightig of groups
plot_data_short['value_g'] <- ifelse(grepl('.*Marker.*ID.*',plot_data_short$Endpoint), plot_data_short$value_g_norm, plot_data_short$value)

# plot
f1 <- ggplot(plot_data_short, aes(x=method, y=as.numeric(as.character(value_g)), ymin=0,ymax=value))+facet_grid(Endpoint ~ ., scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')  
f2 <- f1+geom_boxplot(data = plot_data_short[grepl('\\#.*marker.*',plot_data_short$Endpoint),], lwd=0.2, outlier.size = 0.5, outlier.stroke = 0.2) 
f4 <- f2+geom_boxplot(data = plot_data_short[grepl('.*eproducibility.*',plot_data_short$Endpoint),], lwd=0.2, outlier.size = 0.5, outlier.stroke = 0.2)
f5 <- f4+geom_col(data = plot_data_short[grepl('.*Marker.*ID.*',plot_data_short$Endpoint),], aes(fill=eval), width = 0.7)
f7<- f5+scale_fill_manual('Identification',
                          breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker","prop.wrong.marker.single" ),
                          values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
  
  
  theme(axis.text.x = element_text(angle = 0,  size = 7), axis.text.y = element_text(size = 7), legend.position = 'bottom') +
  ylab('') +xlab('') + facet_grid(Endpoint~Group,  scales="free_y") + theme(panel.spacing = unit(0.2, "lines"))

# export
pdf('./Round2_sum_fig_per_group_no_lines.pdf', width = 4, height = 4.2)
f7
dev.off()

