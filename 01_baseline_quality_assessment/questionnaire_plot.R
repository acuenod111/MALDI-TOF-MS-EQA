# load packages
library('dplyr')
library('tidyr')
library('tidyverse')
library('cowplot')
library('purrr')
library('sjPlot')
library('ggpubr')
library('rstatix')


#### Define functions
## This function 'ID_eval' counts the proportion of spectra which have been identified correctly / multispecies correctly / not identified / wrongly identified with three different spectra identification databases (microflex Biotyper database, VitekMS database or a marker based spectral identification)
## The input to this function are 
#  (i) the eval file and 
#  (ii) one or two grouping variables (columns of the 'eval' file), by which the spectra identification should be summerised (eg. by bacterial group, by center, or by answer to a specific questionnaire question)

ID_eval <- function (eval_file, grouping_variable1, grouping_variable2) {
  myenc1 <- enquo(grouping_variable1)
  myenc2 <- enquo(grouping_variable2)
    ID_eval_summary <- eval_file %>%
    group_by(!!myenc1, !!myenc2) %>%
    dplyr::summarize(prop.bruker.correct.single.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
              prop.bruker.correct.multi.over.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score >= 2 & brukerDB.n_species_over_2 >= 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
              prop.bruker.correct.under.two = sum(grepl('TRUE', brukerDB.correct_species_identified) & brukerDB.Score < 2) / sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)), 
              prop.noID.bruker = sum(grepl('no ID', brukerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
              prop.wong.bruker.over.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score >= 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
              prop.wong.bruker.under.two = sum(grepl('FALSE', brukerDB.correct_species_identified) & brukerDB.Score < 2) /sum(grepl('TRUE|FALSE|no ID', brukerDB.correct_species_identified)),  
            
              prop.correct.marker.single = sum(grepl('TRUE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),             
              prop.wrong.marker.single = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('FALSE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
              prop.correct.marker.multi = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('TRUE', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
              prop.noID.marker = sum(grepl('FALSE', markerDB.correct_single_species) & grepl('no ID', markerDB.correct_species_identified)) /sum(grepl('TRUE|FALSE|no ID', markerDB.correct_species_identified)),
            
              prop.correct.VitekMS.single = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
              prop.correct.VitekMS.lowdiscriminatory = sum(grepl('TRUE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
              prop.wrong.VitekMS.single = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'SingleChoice') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
              prop.wrong.VitekMS.lowdiscriminatory = sum(grepl('FALSE', VitekMS.DBcorrect_species_identified) & VitekMS.DBidentifType == 'LowDiscrimination') / sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)),
              prop.noID.VitekMS = sum(grepl('noID', VitekMS.DBcorrect_species_identified)) /sum(grepl('TRUE|FALSE|noID', VitekMS.DBcorrect_species_identified)))
  return(ID_eval_summary)}

## This function brings together the summarised proposrtions (see above) with the 'eval' file which contains the compiled quality feature data per spectrum
# The inputs for this function are
# (i) the 'eval_file', one row per spectrum, each column representing a spectral feature, used for quality assessment (eg. the number of ribosomal subunits detected, total intensity, measurement error, etc.)
# (ii) the 'ID_eval_out' file (output of the above function 'ID_eval')
# (iii) 1-2 grouping variables by which the spectra should be grouped. (eg. group all spectra together of laboratories which gave the same answer to a particular question in the questionnaire, group all spectra together which have been acquired on the same device, etc.)
# the outputted dataframe ('plot_data') can serve as input files for the plotting functions (defined below)

summarise_prop_and_eval<-function(eval_file, ID_eval_out, grouping_variable1, grouping_variable2='') {
  eval<-merge(eval_file, ID_eval_out, by = intersect(colnames(eval_file), colnames(ID_eval_out)), all.x  = T)
  eval_prop_long<-gather(ID_eval_out, "classification_eval_prop", "value", -c(!!grouping_variable1, !!grouping_variable2))
  
  #Gather all endpoints to one column in order to plot all at once using facet_wrap
  eval_long<-gather(eval[,c('max_mass_ribo', 'spectra',"n_peaks", "mean_intensity_peaks","highest_peaks", "mass_90", "n_high_peaks", "frac_peaks_repr_2spectra", "n_ribos_detected", "n_predicted_su", "rel_amount_su_detected", "mean_dist", "mean_dist_ppm", "mean_intensity","median_intensity", "total_intensity_peaks","brukerDB.Score", "brukerDB.n_species", "brukerDB.n_genera", "brukerDB.n_species_over_2", "brukerDB.n_genera_over_2", "brukerDB.diff", "brukerDB.correct_genus_identified", "brukerDB.correct_species_identified",  "markerDB.match_count",  "markerDB.n_genera_identified",  "markerDB.correct_genus_identified", "markerDB.correct_species_identified", "markerDB.correct_single_species")], "Endpoint", "value", -spectra)
  eval_long<-merge(eval[,c("spectra","Strain", "strainnumber", "Shape", "Gram", "aerobic", "Family", "Group", "spectra", "species_NGS", "genus_NGS", 'MALDI.device_full', grouping_variable1, grouping_variable2)], eval_long, by='spectra')
  eval_long['value']<-ifelse(eval_long$value == 'TRUE', 1,
                           ifelse(eval_long$value == 'FALSE', 0, eval_long$value))
  eval_long["value"]<-gsub(',', '.', eval_long$value)
  eval_long["value"]<-as.numeric(as.character(eval_long$value)) 
  eval_prop_long['DB']<-ifelse(grepl('bruker', eval_prop_long$classification_eval_prop), 'bruker', 
                             ifelse(grepl('marker', eval_prop_long$classification_eval_prop), 'marker', 
                                    ifelse(grepl('VitekMS', eval_prop_long$classification_eval_prop), 'VitekMS', NA)))

  eval_prop_long['eval']<-eval_prop_long$classification_eval_prop
  eval_prop_long<-eval_prop_long[eval_prop_long$DB == 'marker',]
  eval_prop_long$eval<- factor(eval_prop_long$eval, levels = c("prop.wrong.marker.single", "prop.noID.marker", "prop.correct.marker.multi", "prop.correct.marker.single"))
  # only include n_ribos detected and reproducibility for plotting
  plot_data<-eval_long[eval_long$Endpoint %in% c("n_ribos_detected","frac_peaks_repr_2spectra"),]

  # Add column 'rel_amount_su_detected' displaying the relative amount of subunits detected. 
  # This is important when we are comparing between different phylogenetic groups as not all strains encode the same number of ribosomal subunits within the mass range of MALDI-TOF MS. 
  # When comparing between different measurements of the same strains, the absolut number of detected ribosomal subunits can be used
  plot_data_rel<-eval_long[eval_long$Endpoint  == "rel_amount_su_detected",c("spectra","value", "Endpoint")]
  colnames(plot_data_rel)<-c("spectra","value.rel","Endpoint.rel")
  plot_data<-merge(plot_data, plot_data_rel, by.x="spectra", by.y = "spectra", all =T)
  plot_data$Endpoint.rel <- ifelse(plot_data$Endpoint == 'n_ribos_detected', plot_data$Endpoint.rel, plot_data$Endpoint)
  plot_data$value.rel <- ifelse(plot_data$Endpoint == 'n_ribos_detected', plot_data$value.rel, plot_data$value)
  
  # add 'lab' ID; the first three digits in each spectrum name correspond to a laboratory
  plot_data['lab']<-substring(plot_data$spectra, 1,3)
  # rename DB to endpoint
  names(eval_prop_long)[names(eval_prop_long) == 'DB'] <- 'Endpoint'
  eval_prop_long$classification_eval_prop<-NULL
  # merge the two dataframes together
  plot_data<-merge(plot_data, eval_prop_long, by = intersect(colnames(eval_prop_long), colnames(plot_data)), all = TRUE)
  plot_data$Endpoint<- factor(plot_data$Endpoint, levels = c("n_ribos_detected", "frac_peaks_repr_2spectra", "marker"), 
                              labels = labels_endpoints)
  plot_data<-plot_data[!is.na(plot_data[grouping_variable1]),]
  plot_data<-plot_data[!is.na(plot_data[grouping_variable2]),]
  return(plot_data)}

# This function plots the 'plot_data' (summarized by the above defined function): 
# The numerical spectral quality factors (number of detected marker masses and reproducibility) are displayed as boxplots and the evaluation of the species identification through a marker based approach as stacked colored barplots
# Inputs to this function are:
# (i) Spectral evaluation data, summarised as 'plot_data' dataframe by the above defined function
# (ii) A grouping variable according to which the spectra are split (eg. group all spectra of laboratories which have given the same answer to a particular question in the questionnaire). Grouped spectral data are displayed as sperate columns in the plot
plot_plot_data<-function(plot_data, grouping_variable){
  grouping_variable <- enquo(grouping_variable)
  f1 <- ggplot(plot_data, aes(x=!!grouping_variable, y=as.numeric(as.character(value)), ymin=0,ymax=value))+facet_grid(Endpoint ~ MALDI.system, scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('') +xlab('')  
  f2 <- f1+geom_boxplot(data = plot_data[grepl('N.*marker.*',plot_data$Endpoint),])
  f4 <- f2+geom_boxplot(data = plot_data[grepl('.*eproducibility.*',plot_data$Endpoint),])
  f5 <- f4+geom_col(data = plot_data[grepl('.*Marker.*ID.*',plot_data$Endpoint),], aes(fill=eval), width = 0.7)
  f7 <- f5+scale_fill_manual('Identification',
                           breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker","prop.wrong.marker.single" ),
                           values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('') 
  return(f7)}

# Same as the above function with the addition that one grouping variable (eg. the answer to a question) can behighlighted (eg. gree background for factors, which have a positive effect on mass spectral quality)
# Inputs to this function are:
# (i) Spectral evaluation data, summarised as 'plot_data' dataframe by the above defined function
# (ii) A grouping variable according to which the spectra are split (eg. group all spectra of laboratories which have given the same answer to a particular question in the questionnaire). Grouped spectral data are displayed as sperate columns in the plot
# (iii) Hex code of the highlighting colour
plot_plot_data_highlight<-function(plot_data, grouping_variable, highlight_color){
  gp_enc <- deparse(substitute(grouping_variable))
  avals1 = c(1, rep(1, length(unique(plot_data[[gp_enc]]))))
  avalsHex = paste0("#000000", toupper(as.hexmode(round(avals1*255))))
  grouping_variable <- enquo(grouping_variable)
  f1 <- ggplot(plot_data, aes(x=as.factor(!!grouping_variable), y=as.numeric(as.character(value)), ymin=0,ymax=value)) + annotate("rect", xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = highlight_color, alpha = .3) +facet_grid(Endpoint ~ MALDI.system, scales = 'free_y') + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylab('') +xlab('')   
  f2 <- f1+geom_boxplot(data = plot_data[grepl('N.*marker.*',plot_data$Endpoint),], aes(color=!!grouping_variable, alpha = !!grouping_variable), lwd=0.3) +  scale_alpha_manual(values = avals1) + scale_colour_manual(values = avalsHex)
  f4 <- f2+geom_boxplot(data = plot_data[grepl('.*eproducibility.*',plot_data$Endpoint),], aes(color=!!grouping_variable, alpha =!!grouping_variable), lwd=0.3) 
  f5 <- f4+geom_col(data = plot_data[grepl('.*Marker.*ID.*',plot_data$Endpoint),], aes(fill=eval, alpha =!!grouping_variable), width = 0.7)
  f7 <- f5+scale_fill_manual('Identification',
                             breaks = c("prop.correct.marker.single","prop.correct.marker.multi", "prop.noID.marker","prop.wrong.marker.single" ),
                             values = c("#009E73", "#0072B2", "#999999","#D55E00")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('') +xlab('') 
  return(f7)}


#### Define how to label endpoints
# labels_endpoints<-c("Number of phylogenetic markers detected", "Median relative intensity of the phylogenetic markers", "Total intensity of all peaks (log10)", "Fraction of reproducibly detected peaks", "Evaluation of a marker based species identification")

# Exclusively plot the Number of marker masses, the Reproducibility (between 2 spectra) and the evaluation of the marker based ID
labels_endpoints<-c("No. marker masses", "Reproducibility", "Marker based ID")

# Introduce line breaks for very long labels
labels_endpoints<-str_wrap(labels_endpoints, width = 15)

# Read in the 'eval' file, has been exported from "summarize_eval_participating_labs.R"
eval<-read.csv("./Spectra_and_SpeciedID/01_baseline_quality_assessment/eval_participating_labs_1.txt", header = T, sep=';')
# Exclude empty spectra from further analysis
eval <- eval[eval$n_peaks > 0, ]

# correct for typo in column name
colnames(eval)[colnames(eval) == 'meadian_intensity'] <- 'median_intensity'

# Spectra which have a microflex Biotyper species ID have been acquired on Bruker devices
bruker.devices<-as.character(unique(eval[!is.na(eval$brukerDB.correct_species_identified), 'MALDI.device_full']))
# Spectra which have a VitekMS species ID have been acquired on a VitekMS / Shimadzu device. 
# 'VitekMS.devices' include both of these, 'Shimadzu.devices' only the spectra which have been acquired on shimadzu devices and 'VitekMS.devices.exclusive' exclusively the ones which have been acquired on VitekMS devices
VitekMS.devices<-as.character(unique(eval[!is.na(eval$vitekMSDB.correct_species_identified), 'MALDI.device_full']))
Shimadzu.devices<-c("device_35","device_36","device_37","device_44","device_45")

# Spectra acquired on VitekMS devices are output as peaklist, whereas for spectra acquired on Shimadzu devices, the peak picking has been performed via laungepad software
VitekMS.devices.exclusive<-setdiff(VitekMS.devices,Shimadzu.devices)

# Add a column 'mean inensity per peak'
eval['mean_intensity_peaks']<-eval$total_intensity_peaks / eval$n_peaks

# split up the group of strains 'streptococcus' in 'non-viridans streptococci' and 'viridans streptococci' as the latter ones are of special interest in clinical diagnostics
eval['Group']<-as.factor(ifelse(grepl('infantis|pseudopneumoniae|pneumoniae|gordonii', as.character(eval$Strain)), 'Viridans Streptococci', 
                                     ifelse(grepl('equinus|gallolyticus|lutetiensis|dysgalactiae', as.character(eval$Strain)), 'Other Streptococci', as.character(eval$Group))))


# For the evaluation of identification with the MBT Bruker database, set everything to 'no identification possible' which has a score lower than 1.7
eval["brukerDB.correct_species_identified"]<-ifelse(eval$brukerDB.Score<1.7, 'no ID', eval$brukerDB.correct_species_identified)
eval['brukerDB.correct_species_identified']<-ifelse(eval$brukerDB.correct_species_identified == '1', TRUE, 
                                                         ifelse(eval$brukerDB.correct_species_identified == '0', FALSE, 
                                                                ifelse(eval$brukerDB.correct_species_identified == 'no ID', 'no ID', eval$brukerDB.correct_species_identified)))

# For the evaluation of identification with the marker based, set everything to 'no identification possible' which has a score lower than a group specific threshold
# Use the following thresholds: 
# E. coli / Shigella: 15
# Enterobacter: 20
# Klebsiella: 15
# S.aureus 7: 
# all other: 10

# Extract the species as and genus as assigned by WGS
eval['species_NGS']<-gsub('([[:alpha:]]+)(\\_)([[:alpha:]]+)(\\_)(.*)', '\\3', eval$Strain)
eval['genus_NGS']<-gsub('([[:alpha:]]+)(\\_)([[:alpha:]]+)(\\_)(.*)', '\\1', eval$Strain)

# set the marker per bacterial group as defined above
eval['marker.theshold']<-ifelse(eval$genus_NGS %in% c('Escherichia', 'Shigella', 'Klebsiella'), 15,
                                     ifelse(eval$genus_NGS == 'Enterobacter', 20,
                                            ifelse(eval$genus_NGS == 'Staphylococcus', 7, 10)))

# set all identifications with a score lower tham the threshold to 'no identification possible'
eval["markerDB.species_identified"]<-ifelse(as.numeric(eval$markerDB.match_count) < eval$marker.theshold, 'no identification possible',eval$markerDB.species_identified)
eval["markerDB.Genus"]<-ifelse( as.numeric(eval$markerDB.match_count) < eval$marker.theshold, 'no identification possible',eval$markerDB.species_identified)

# If no identification was possible 0 genera and 0 species were identified
eval["markerDB.n_species_identified"]<-ifelse( as.numeric(eval$markerDB.match_count) < eval$marker.theshold, 0,eval$markerDB.n_species_identified)
eval["markerDB.n_genera_identified"]<-ifelse( as.numeric(eval$markerDB.match_count) < eval$marker.theshold, 0,eval$markerDB.n_genera_identified)
eval["markerDB.correct_species_identified"]<-ifelse( as.numeric(eval$markerDB.match_count) < eval$marker.theshold, 'no ID',eval$markerDB.correct_species_identified)
eval$marker.theshold<-NULL

# introduce column 'correct_single_species' which states whether the correct species was identified as only species
eval['markerDB.correct_single_species']<-ifelse(eval$markerDB.correct_species_identified == TRUE & eval$markerDB.n_species_identified == '1', TRUE, FALSE)

# Vitek MS
# Rename identification type
eval['VitekMS.DBidentifType']<-eval$vitekMSDB.identifType
eval$vitekMSDB.identifType<-NULL

# set all VitekMS identifications with an 'Identification type' = 'NoIdentification','NotEnoughPeaks','TooNoisy' to 'no ID'
eval['VitekMS.DBcorrect_species_identified']<-ifelse(eval$VitekMS.DBidentifType %in% c('NoIdentification','NotEnoughPeaks','TooNoisy'), 'noID',
                                                          ifelse(eval$vitekMSDB.correct_species_identified == 0, FALSE,
                                                                 ifelse(eval$vitekMSDB.correct_species_identified == 1, TRUE, NA)))

# Import questionnaire
questionnaire_raw<-read.csv2('./Table\ S3_questionnaire.csv', sep=',',na.strings=c(""," ","NA"))
questionnaire_raw<-questionnaire_raw[!is.na(questionnaire_raw$lab),]
# save all columns as character
questionnaire_raw<-questionnaire_raw %>%
  mutate_all(as.character)

# remove lab name and when the questionnaire has been filled out
questionnaire_raw$What.is.the.name.of.your.laboratory.<-NULL
questionnaire_raw$Zeitstempel<-NULL
questionnaire<-questionnaire_raw

# Add MALDI.device to merge, 'MALDI.system' (bruker or shimadzu/vitekMS column) and MALDI.type (exact version)
# The lab 'AL' the in the questionnaire using 2 devices (smart and sirius), but only one set of spectra has been received, measured on a MBT sirius
questionnaire['Which.MALDI.TOF.MS.system.do.you.use.in.for.this.ringtrial.']<-ifelse(grepl('AL', questionnaire$lab), 'Microflex Biotyper by Bruker Daltonics', questionnaire$Which.MALDI.TOF.MS.system.do.you.use.in.for.this.ringtrial.)

# If a lab uses more than one MALDI device in this EQA, reformat so that there is one row per device 
two.devices<-questionnaire[questionnaire$lab %in% c('AJ', 'AT', 'AV', 'BD', 'BE'),]
device_1<-two.devices
device_1['MALDI.device']<-ifelse(grepl('AJ', device_1$lab), 'device_10',
                                 ifelse(grepl('AT', device_1$lab), 'device_21',
                                        ifelse(grepl('AV', device_1$lab), 'device_24',
                                               ifelse(grepl('BD', device_1$lab), 'device_33',
                                                      ifelse(grepl('BE', device_1$lab), 'device_35', NA)))))

device_2<-two.devices
device_2['MALDI.device']<-ifelse(grepl('AJ', device_2$lab), 'device_11',
                                 ifelse(grepl('AT', device_2$lab), 'device_22',
                                        ifelse(grepl('AV', device_2$lab), 'device_25',
                                               ifelse(grepl('BD', device_2$lab), 'device_34',
                                                      ifelse(grepl('BE', device_2$lab), 'device_38', NA)))))

two.devices<-rbind(device_1, device_2)

# remove these (have two devices and have just been added as such in previous step)
questionnaire<-questionnaire[!grepl('AJ|AT|AV|BD|BE', questionnaire$lab),]
# add two devices
questionnaire<-merge(questionnaire, two.devices, by=intersect(colnames(questionnaire), colnames(two.devices)), all=TRUE)
questionnaire['MALDI.device']<-ifelse(is.na(questionnaire$MALDI.device), paste0('device_', questionnaire$Device.ID), questionnaire$MALDI.device)

setdiff(questionnaire$MALDI.device, eval$MALDI.device)

# for some laboratories, add manually which device have been used
questionnaire['MALDI.system']<-ifelse(questionnaire$Which.MALDI.TOF.MS.system.do.you.use.in.for.this.ringtrial. !='Other', questionnaire$Which.MALDI.TOF.MS.system.do.you.use.in.for.this.ringtrial,
                                      ifelse(questionnaire$MALDI.device %in% c('device_03', 'device_07', 'device_21', 'device_22', 'device_38'), 'Microflex Biotyper by Bruker Daltonics', 
                                             ifelse(questionnaire$MALDI.device %in% c('device_35', 'device_43', 'device_44'), 'Shimadzu/VitekMS', NA)))
# summarize more broad into 'MALDI systems'
questionnaire['MALDI.system']<-gsub('VitekMS by Biomérieux', 'Shimadzu/VitekMS', questionnaire$MALDI.system)
questionnaire['MALDI.system']<-gsub('Microflex Biotyper by Bruker Daltonics', 'microflexBiotyper', questionnaire$MALDI.system)
questionnaire['MALDI.system']<-ifelse(questionnaire$lab == 'device_26', 'microflexBiotyper', questionnaire$MALDI.system)
questionnaire['MALDI.system']<-ifelse(questionnaire$MALDI.device == 'device_11', 'Shimadzu/VitekMS', questionnaire$MALDI.system)

# Add which 'type' of MALDI has been used, also distinguishing between the different MBT methods and between VitekMS and Shimadzu devices
questionnaire['MALDI.type']<-ifelse(!is.na(questionnaire$Which.Microflex.Biotyper.System.do.you.use.for.this.ringtrial.), questionnaire$Which.Microflex.Biotyper.System.do.you.use.for.this.ringtrial.,
                                    ifelse(questionnaire$MALDI.device %in% c("device_03","device_07"), NA, 
                                           ifelse(questionnaire$MALDI.device %in% c("device_12","device_09", "device_15", "device_30", 'device_11') , 'VitekMS by Biomérieux',
                                                  ifelse(questionnaire$MALDI.device %in% c("device_14", 'device_07', "device_26") , 'MBT Sirius',
                                                         ifelse(questionnaire$MALDI.device %in% c("device_21") , "microflex LT/SH",
                                                                ifelse(questionnaire$MALDI.device %in% c("device_35","device_44","device_433"), "Shimadzu (Axima Confidence)", NA))))))

questionnaire['MALDI.type']<- ifelse(questionnaire$MALDI.device == 'device_11', 'VitekMS by Biomérieux', questionnaire$MALDI.type)

# summarize which MALDI to single variable
questionnaire$MALDI.type <- gsub('MBT Sirius','microflex LT/SH “smart”',questionnaire$MALDI.type)
questionnaire$Which.MALDI.TOF.MS.system.do.you.use.in.for.this.ringtrial.<-NULL                                             
questionnaire$Which.MALDI.TOF.MS.System.do.you.use.<-NULL

## remove freetext question for statistical analysis
questionnaire$For.which.species.have.you.extended.the.database.of.your.MALDI.TOF.MS.system.<-NULL
questionnaire$Which.strains.do.you.use.to.control.the.quality.of.your.MALDI.TOF.MS.measurements.<-NULL
questionnaire$How.do.you.prepare.acid.fast.bacteria.<-NULL
questionnaire$For.which.species.do.you.confirm..perform.further.tests.<-NULL
questionnaire$For.which.species.do.you.routinely.perform.more.elaborate.sample.preparation.protocols..such.as.protein.exctraction.using.EtOH....multiple.answers.possible.<-NULL
questionnaire$Strains.are.generally.reported.on.genus.level..if.the.score.is.higher.than.<-NULL
questionnaire$Strains.are.generally.reported.on.species.level..if.the.score.is.higher.than.<-NULL
questionnaire$Where.do.you.store.the.target.plates.<-NULL # question not clear enough
questionnaire$How.else.did.you.process.the.received.strains.<-NULL
questionnaire$If.two.species.have.a.score.over.the.species.threshold..you.report<-NULL #most frequent answer (>60%) is 'none of the above'

# Convert the answers to factors (binary or few levels only)
questionnaire['FA.applied.routinely.bacteria']<-ifelse(grepl('For all microorganisms|unknown bacteria|for all bacteria',questionnaire$For.which.samples.do.you.use.formic.acid..multiple.answers.possible..),'yes','no')
questionnaire['Regular.quality.assessment']<-ifelse(questionnaire$Do.you.regularly.assess.the.quality.of.your.MALDI.TOF.MS.measurements.by.e.g..measuring.one.or.more.reference.strain.and.assessing.the.scores.or.the.number.of.detected.marker.masses. == 'No', 'no', 
                                                    ifelse(is.na(questionnaire$Do.you.regularly.assess.the.quality.of.your.MALDI.TOF.MS.measurements.by.e.g..measuring.one.or.more.reference.strain.and.assessing.the.scores.or.the.number.of.detected.marker.masses.),NA,'yes'))
questionnaire['MALDI.station']<-ifelse(grepl('Yes|we have 6-8 people trained and competent on a weekly rota system',questionnaire$Is.there.a.MALDI.TOF.MS.workstation.in.your.laboratory...For.a.certain.period..one.or.more.member.of.staff.is.responsible.for.all.MALDI.TOF.MS.measurements.), 'yes', 
                                       ifelse(is.na(questionnaire$Is.there.a.MALDI.TOF.MS.workstation.in.your.laboratory...For.a.certain.period..one.or.more.member.of.staff.is.responsible.for.all.MALDI.TOF.MS.measurements.), NA, 'no'))
questionnaire['samplepre_training']<-ifelse(questionnaire$Do.you.perform.training.for.MALDI.TOF.MS.sample.processing. %in% c('No', 'will be the case in future'), 'no', 
                                            ifelse(is.na(questionnaire$Do.you.perform.training.for.MALDI.TOF.MS.sample.processing.), NA, 'yes'))
questionnaire['regular.hardware.service']<-ifelse(questionnaire$Is.there.a.regular.hardware.service.on.your.MALDI.TOF.MS.device.by.the.provider. %in% c('No', 'will be the case in future'), 'no', 
                                                  ifelse(is.na(questionnaire$Is.there.a.regular.hardware.service.on.your.MALDI.TOF.MS.device.by.the.provider.), NA, 'yes'))
questionnaire['regular.hardware.service']<-ifelse(questionnaire$regular.hardware.service == '', NA, questionnaire$regular.hardware.service)

questionnaire['How.clean.target.plate']<-ifelse(grepl('(B|b)ruker|see SOP|cleaning protocol according to manufacturer', questionnaire$How.do.you.clean.the.target.plate.), 'TFA', 
                                                ifelse(grepl('((E|e)thanol)*.*(TFA|trifluor)', questionnaire$How.do.you.clean.the.target.plate.), 'TFA',
                                                       ifelse(grepl('(TFA|(t|T)rifluor).*((E|e)thanol)*', questionnaire$How.do.you.clean.the.target.plate.), 'TFA',
                                                              ifelse(grepl('methan(o)*l.*aceton', questionnaire$How.do.you.clean.the.target.plate.), 'Methanole-Acetone',
                                                                     ifelse(grepl('Disposable', questionnaire$Which.target.plates.do.you.use.),'Disposable', 'Other')))))
# remove these columns including freetext, they have just been summarised
questionnaire$How.do.you.clean.the.target.plate.<-NULL
questionnaire$Is.there.a.regular.hardware.service.on.your.MALDI.TOF.MS.device.by.the.provider.<-NULL
questionnaire$Do.you.perform.training.for.MALDI.TOF.MS.sample.processing.<-NULL
questionnaire$For.which.samples.do.you.use.formic.acid..multiple.answers.possible..<-NULL
questionnaire$Do.you.regularly.assess.the.quality.of.your.MALDI.TOF.MS.measurements.by.e.g..measuring.one.or.more.reference.strain.and.assessing.the.scores.or.the.number.of.detected.marker.masses.<-NULL
questionnaire$Is.there.a.MALDI.TOF.MS.workstation.in.your.laboratory...For.a.certain.period..one.or.more.member.of.staff.is.responsible.for.all.MALDI.TOF.MS.measurements.<-NULL

# Harmonise answers to a few levels and shorten answers to display on figures
questionnaire$What.was.the.purpose.of.extending.the.database.of.your.MALDI.TOF.MS.system.<-as.factor(ifelse(grepl('species (identification|ID)', questionnaire$What.was.the.purpose.of.extending.the.database.of.your.MALDI.TOF.MS.system.), 'Routine ID',
                                                                                                  ifelse(grepl('Research', questionnaire$What.was.the.purpose.of.extending.the.database.of.your.MALDI.TOF.MS.system.),  'Research', NA)))
questionnaire['room.sun']<-factor(ifelse(!is.na(questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), 
                                  ifelse(grepl('the sun shines on the MALDI',questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..),  'sun', 'no.sun'), NA), levels = c('sun','no.sun'))
questionnaire['room.aircondition']<-factor(ifelse(!is.na(questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), 
                                           ifelse(grepl('the room is air conditioned',questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), 'aircon', 'no.aircon'), NA), levels = c('aircon', 'no.aircon'))

# for questions where multiple nominal categories can be selected, create contrast for every selected (1) vs. not selected (-1)
# which DB are used
questionnaire['which.VitekMSDB.v2']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.VitekMS.system), 
                                            ifelse(questionnaire$Which.databases.are.installed.on.your.VitekMS.system. == 'VitekMS Version 2','yes', 'no'), NA),  levels = c('yes', 'no'))
questionnaire['which.VitekMSDB.v4']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.VitekMS.system), 
                                            ifelse(questionnaire$Which.databases.are.installed.on.your.VitekMS.system. == 'VitekMS Version 4','yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.VitekMSDB.Saramis']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.VitekMS.system), 
                                                 ifelse(grepl('Saramis', questionnaire$Which.databases.are.installed.on.your.VitekMS.system.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.RUO']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                               ifelse(grepl('RUO',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.IVD']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                               ifelse(grepl('IVD',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.Security.Related']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                                            ifelse(grepl('Security\\-Related',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.Filamentous.Fungi']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                                             ifelse(grepl('Filamentous Fungi',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.Mycobacteria']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                                        ifelse(grepl('Mycobacteria',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.Sepsityper']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                                      ifelse(grepl('Sepsityper',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.in.house']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                                    ifelse(grepl('((O|o)wn)|(in\\-house)|(homemade)|(Uniklinik)',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.microflexDB.Subtyping']<-factor(ifelse(!is.na(questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..), 
                                                     ifelse(grepl('Subtyping',questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire$Which.databases.are.installed.on.your.Microflex.Biotyper.system..more.than.one.answer.possible..<-NULL

# which quality parameters are used
questionnaire['which.quality.parameter.score']<-factor(ifelse(!is.na(questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.), 
                                                       ifelse(grepl('The score assigned by the MALDI-TOF MS Software',questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.quality.parameter.marker']<-factor(ifelse(!is.na(questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.), 
                                                        ifelse(grepl('The number of predefined marker masses detected',questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.quality.parameter.total_n_peaks']<-factor(ifelse(!is.na(questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.), 
                                                               ifelse(grepl('The total amount of peaks detected',questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.quality.parameter.highest_peak']<-factor(ifelse(!is.na(questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.), 
                                                              ifelse(grepl('The peak with the highest m.z value detected',questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.quality.parameter.Time.per.spot']<-factor(ifelse(!is.na(questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.), 
                                                               ifelse(grepl('Time',questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['which.quality.parameter.signal.to.noise']<-factor(ifelse(!is.na(questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.), 
                                                                 ifelse(grepl('(s|S)ignal.*(n|N)oise',questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire$Which.parameters.do.you.examine.when.assessing.the.quality.of.your.MALDI.TOF.MS.measurements...multiple.options.possible.<-NULL

# what reason to call helpdesk
questionnaire['reason.call.helpdesk.Software']<-factor(ifelse(!is.na(questionnaire$What.are.reasons.for.which.you.call.the.helpdesk.), 
                                                       ifelse(grepl('Software',questionnaire$What.are.reasons.for.which.you.call.the.helpdesk.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['reason.call.helpdesk.Hardware']<-factor(ifelse(!is.na(questionnaire$What.are.reasons.for.which.you.call.the.helpdesk.), 
                                                       ifelse(grepl('Hardware',questionnaire$What.are.reasons.for.which.you.call.the.helpdesk.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['reason.call.helpdesk.Calibration']<-factor(ifelse(!is.na(questionnaire$What.are.reasons.for.which.you.call.the.helpdesk.), 
                                                          ifelse(grepl('Calibration',questionnaire$What.are.reasons.for.which.you.call.the.helpdesk.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire$What.are.reasons.for.which.you.call.the.helpdesk.<-NULL

# How formic acid is applied
questionnaire['how.apply.FA.manually.single.step']<-factor(ifelse(!is.na(questionnaire$How.do.you.apply.formic.acid.and.matrix..multiple.answers.possible..), 
                                                           ifelse(grepl('Manually using a single step pipette',questionnaire$How.do.you.apply.formic.acid.and.matrix..multiple.answers.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['how.apply.FA.manually.multi.step']<-factor(ifelse(!is.na(questionnaire$How.do.you.apply.formic.acid.and.matrix..multiple.answers.possible..), 
                                                          ifelse(grepl('Manually using a multistep pipette',questionnaire$How.do.you.apply.formic.acid.and.matrix..multiple.answers.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['how.apply.FA.Galaxy']<-factor(ifelse(!is.na(questionnaire$How.do.you.apply.formic.acid.and.matrix..multiple.answers.possible..), 
                                             ifelse(grepl('Galaxy',questionnaire$How.do.you.apply.formic.acid.and.matrix..multiple.answers.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire$How.do.you.apply.formic.acid.and.matrix..multiple.answers.possible..<-NULL

# which media are routinely used
questionnaire['agar.sheep.blood']<-factor(ifelse(!is.na(questionnaire$What.agar.plates.are.the.microbial.colonies.which.are.measured.on.the.MALDI.TOF.MS.grown.on.), 
                                          ifelse(grepl('5\\% Sheep agar plates',questionnaire$What.agar.plates.are.the.microbial.colonies.which.are.measured.on.the.MALDI.TOF.MS.grown.on.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['agar.chrom']<-factor(ifelse(!is.na(questionnaire$What.agar.plates.are.the.microbial.colonies.which.are.measured.on.the.MALDI.TOF.MS.grown.on.), 
                                    ifelse(grepl('Chrom agar plates',questionnaire$What.agar.plates.are.the.microbial.colonies.which.are.measured.on.the.MALDI.TOF.MS.grown.on.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['agar.other']<-factor(ifelse(!is.na(questionnaire$What.agar.plates.are.the.microbial.colonies.which.are.measured.on.the.MALDI.TOF.MS.grown.on.), 
                                    ifelse(grepl('Other',questionnaire$What.agar.plates.are.the.microbial.colonies.which.are.measured.on.the.MALDI.TOF.MS.grown.on.),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire$What.agar.plates.are.the.microbial.colonies.which.are.measured.on.the.MALDI.TOF.MS.grown.on.<-NULL

# How to prevent dust to enter the machine
questionnaire['dust.gloves']<-factor(ifelse(!is.na(questionnaire$How.do.you.prevent.dust.to.come.in.to.the.flight.tube...more.than.one.answer.possible..), 
                                     ifelse(grepl('Wear gloves when inserting the target plate',questionnaire$How.do.you.prevent.dust.to.come.in.to.the.flight.tube...more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['dust.other']<-factor(ifelse(!is.na(questionnaire$How.do.you.prevent.dust.to.come.in.to.the.flight.tube...more.than.one.answer.possible..), 
                                    ifelse(grepl('Other',questionnaire$How.do.you.prevent.dust.to.come.in.to.the.flight.tube...more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire['dust.tissue']<-factor(ifelse(!is.na(questionnaire$How.do.you.prevent.dust.to.come.in.to.the.flight.tube...more.than.one.answer.possible..), 
                                     ifelse(grepl('Remove dust with dust tissue before inserting target plate',questionnaire$How.do.you.prevent.dust.to.come.in.to.the.flight.tube...more.than.one.answer.possible..),'yes', 'no'), NA), levels = c('yes', 'no'))
questionnaire$How.do.you.prevent.dust.to.come.in.to.the.flight.tube...more.than.one.answer.possible..<-NULL

# Shorten answers for display reasons
questionnaire$How.did.you.process.the.received.strains.<-gsub('Streak eSwabs on agar plate proceed with MALDI-TOF MS measurements', 'direct.streak.out.from.eswab', questionnaire$How.did.you.process.the.received.strains.)
questionnaire$How.did.you.process.the.received.strains.<-gsub('Streak eSwabs on agar plate, subculture the strains once on another agar plate and proceed with MALDI-TOF MS measurements', 'streak.out.from.eswab.subcultured.once', questionnaire$How.did.you.process.the.received.strains.)
questionnaire$How.did.you.process.the.received.strains.<-gsub('Freeze the strains, streak out on agar plate and proceed with MALDI-TOF MS measurements', 'freeze.streak.out', questionnaire$How.did.you.process.the.received.strains.)
questionnaire$How.did.you.process.the.received.strains.<-gsub('Freeze the strains, streak out on agar plate, subculture the strains once on another agar and proceed with MALDI-TOF MS measurements', 'freeze.streak.out.subculture.once', questionnaire$How.did.you.process.the.received.strains.)

# Summarise 'BL21' and 'ATCC' calibration strains to 'Other'
questionnaire$Which.strain.do.you.use.to.calibrate.your.MALDI.TOF.MS.device.<-as.factor(ifelse(grepl('BL21|ATCC', questionnaire$Which.strain.do.you.use.to.calibrate.your.MALDI.TOF.MS.device.), 'Other', as.character(questionnaire$Which.strain.do.you.use.to.calibrate.your.MALDI.TOF.MS.device.)))

# For ordinal variables (variables, which have a natural order), encode these in this order
# How often the room is whipped
questionnaire['room.whiped']<-as.factor(ifelse(grepl('The floor is wiped moist once a week',questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), 'weekly',
                                               ifelse(grepl('The floor is wiped moist once a day or more often',questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), 'daily', 
                                                      ifelse(grepl('The floor is wiped moist less than once a week',questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), 'less.than.weekly','missing'))))
# How many people work in the same room
questionnaire['room.people']<-as.factor(ifelse(grepl('on average work 16 \\- 30 people',questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), '16-30',
                                               ifelse(grepl('on average work 0 \\- 15 people',questionnaire$The.MALDI.is.located.in.a.room.where..multiple.answers.possible..), '0-15', 'missing')))
# When was the last time the software was updated
questionnaire$When.was.the.last.software.update..e.g..database..on.your.MALDI.TOF.MS.device.<-as.factor(gsub('Between six and 12 months ago', 'Less than 12 months ago', questionnaire$When.was.the.last.software.update..e.g..database..on.your.MALDI.TOF.MS.device.))

# Merge the eval dataframe (including mass spectral quality features for each spectrum) to the questionnaire (answers summarised and cleaned up)
eval_q<-merge(eval, questionnaire, by = c(intersect(colnames(eval), colnames(questionnaire))))
eval_stat_q<-eval_q

# set all missing answers to NA
eval_stat_q[70:138] <- lapply(eval_stat_q[70:138], function(x) gsub("missing", NA, x))

# only include dataset per device (at Mabritec for the first round of the EQA, three different people measure the same strains. Consider only 1 dataset per device)
eval_stat_q<-eval_stat_q[!grepl('device_36|device_37|device_39|device_40',eval_stat_q$MALDI.device_full),]

# rename calibration column, shorten
colnames(eval_stat_q)[colnames(eval_stat_q) == 'How.often.do.you.calibrate.the.your.MALDI.TOF.MS.system.by.measuring.a.reference.strain.for.which.the.MALDI.TOF.MS.Software.includes.marker.masses..E.coli.BST...DH5alpha..'] <- 'calibration.freq'

# plot how often the device is calibrated vs. the mean measurement ppm detected [error]
eval_stat_q$calibration.freq <- factor(eval_stat_q$calibration.freq, 
         levels = c('With every run', 'Once per day', 'Once per week', 'Once per month'))
eval_stat_q['cal.run']<-ifelse(grepl('run', eval_stat_q$calibration.freq), 'With every run', 'less often')
eval_stat_q$cal.run<-factor(eval_stat_q$cal.run, levels = c('With every run', 'less often'))
# At two participating laboratory, only one spectrum per strain was acquired. The reproducibily between technical 
eval_stat_q$frac_peaks_repr_2spectra<-ifelse(grepl('.*device\\_30.*|.*device\\_42.*', eval_stat_q$MALDI.device), NA, eval_stat_q$frac_peaks_repr_2spectra)

# Calibration 
# calculate p-values between 'with every run' and less often
cal.run.pval <- eval_stat_q %>%
  wilcox_test(mean_dist_ppm ~ cal.run, paired = F)
cal.run.pval <- cal.run.pval %>% add_xy_position(x = "cal.run")

avals1 = c(1, rep(0.4, length(eval_stat_q$calibration.freq)))
avalsHex = paste0("#000000", toupper(as.hexmode(round(avals1*255))))

# Plot the measurement errors vs. calibrating frequencies, highlighting in green 'with every run', as we find a significantly lower measurement error for devices which are calibrated that often
cal.participating <- ggplot(eval_stat_q, aes(x=calibration.freq, y= mean_dist_ppm)) + annotate("rect", xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = '#2A9D8F', alpha = .3) +
  geom_boxplot(aes(color=calibration.freq, alpha = calibration.freq)) +  scale_alpha_manual(values = avals1) + scale_colour_manual(values = avalsHex) + 
  stat_pvalue_manual(cal.run.pval, tip.length = 0, step.increase = 0.2) + ylab('Mean measurement error [ppm]') + xlab('') + ggtitle('Calibration frequency') + theme(legend.position = 'none')

pdf('./questionnaire_cal.pdf', width = 4.5, height = 3.75)
cal.participating
dev.off()

# Summarise spectra evaluation per device
eval_per_device<-eval %>% group_by(MALDI.device) %>%
  summarize(mean_dist_ppm = mean(mean_dist_ppm, na.rm = T), n_ribos_detected = mean(n_ribos_detected), total_intensity_peaks = mean(total_intensity_peaks), median_intensity = mean(median_intensity, na.rm = T), frac_peaks_repr_2spectra = mean(frac_peaks_repr_2spectra))
eval_q_sum<-merge(eval_per_device, questionnaire, by = c(intersect(colnames(eval_per_device), colnames(questionnaire))))
# At two participating laboratory, only one spectrum per strain was acquired. The reproducibily between technical 
eval_q_sum$frac_peaks_repr_2spectra<-ifelse(grepl('.*device\\_30.*|.*device\\_42.*', eval_q_sum$MALDI.device), NA, eval_q_sum$frac_peaks_repr_2spectra)

# remove columns which are the same for all devices (same answer from all laboratories)
eval_stat_q<-Filter(function(x) length(unique(na.omit(x)))!=1, eval_stat_q)

# If no ribosomal subunit was detected, 'max mass ribo' the highest mass of a ribosomal subunit is set to 3'000, the lowest possible mass in a recorded mass spectrum. 
eval_stat_q$max_mass_ribo<-ifelse(is.na(eval_stat_q$max_mass_ribo), 3000, as.numeric(gsub(',', '.', eval_stat_q$max_mass_ribo)))
# If no ribosomal subunit was detected, 'mean intensity', referring to the mean intensity of all ribosomal subunits is set to 0
eval_stat_q$mean_intensity<-ifelse(is.na(eval_stat_q$mean_intensity), 0, as.numeric(gsub(',', '.', eval_stat_q$mean_intensity)))
# reorder MALDI systems
eval_stat_q$MALDI.system<-factor(eval_stat_q$MALDI.system, levels = c('microflexBiotyper', 'Shimadzu/VitekMS'), labels = c('MBT', 'VitekMS'))

# check median, IQR and p-values for the most important changes. Use unpaired wilcoxon rank tests, as the number of samples might differ between the comparing groups although the same strains are compared (as the number of laboratories answering the same thing might differ per question)
# how the received strains were processed
eval_stat_q['how_processed_stats'] <- ifelse(eval_stat_q$How.did.you.process.the.received.strains. == "freeze.streak.out", "freeze.streak.out", 
                                             ifelse(is.na(eval_stat_q$How.did.you.process.the.received.strains.), NA, 'Other'))
eval_stat_q %>% 
  group_by(how_processed_stats) %>%
  summarise_at(vars(n_ribos_detected),
               list(min=min, Q1=~quantile(., probs = 0.25, na.rm = T),
                    median=median, Q3=~quantile(., probs = 0.75, na.rm = T),
                    max=max))
  
wilcox.test(eval_stat_q[eval_stat_q$how_processed_stats == 'freeze.streak.out', 'n_ribos_detected'], 
            eval_stat_q[eval_stat_q$how_processed_stats == 'Other', 'n_ribos_detected'], 
            paired = F)
eval_stat_q$how_processed_stats <- NULL

# Which target plates are used
eval_stat_q %>% 
  group_by(Which.target.plates.do.you.use.) %>%
  summarise_at(vars(n_ribos_detected),
               list(min=min, Q1=~quantile(., probs = 0.25, na.rm = T),
                    median=median, Q3=~quantile(., probs = 0.75, na.rm = T),
                    max=max))

wilcox.test(eval_stat_q[eval_stat_q$Which.target.plates.do.you.use. == 'Steel target', 'n_ribos_detected'], 
            eval_stat_q[eval_stat_q$Which.target.plates.do.you.use. == 'Disposable', 'n_ribos_detected'], 
            paired = F)

# How the target plates are cleaned
eval_stat_q['how_cleaning_stats'] <- ifelse(eval_stat_q$How.clean.target.plate  == "Methanole-Acetone", "Methanole-Acetone", 
                                            ifelse(is.na(eval_stat_q$How.clean.target.plate), NA, 'Other'))
eval_stat_q %>% 
  group_by(how_cleaning_stats) %>%
  summarise_at(vars(n_ribos_detected),
               list(min=min, Q1=~quantile(., probs = 0.25, na.rm = T),
                    median=median, Q3=~quantile(., probs = 0.75, na.rm = T),
                    max=max))

wilcox.test(eval_stat_q[eval_stat_q$how_cleaning_stats == 'Methanole-Acetone', 'n_ribos_detected'], 
            eval_stat_q[eval_stat_q$how_cleaning_stats == 'Other', 'n_ribos_detected'], 
            paired = F)
eval_stat_q$how_cleaning_stats <- NULL

# How often a hardware service is performed
eval_stat_q %>% 
  group_by(regular.hardware.service) %>%
  summarise_at(vars(n_ribos_detected),
               list(min=min, Q1=~quantile(., probs = 0.25, na.rm = T),
                    median=median, Q3=~quantile(., probs = 0.75, na.rm = T),
                    max=max))

wilcox.test(eval_stat_q[eval_stat_q$regular.hardware.service == 'yes', 'n_ribos_detected'], 
            eval_stat_q[eval_stat_q$regular.hardware.service == 'no', 'n_ribos_detected'], 
            paired = F)

# Whether the participating laboratory works with a MALDI workstation (for a certain period of time, one person perocessed all sample for MALDI measurements)
eval_stat_q %>% 
  group_by(MALDI.station) %>%
  summarise_at(vars(n_ribos_detected),
               list(min=min, Q1=~quantile(., probs = 0.25, na.rm = T),
                    median=median, Q3=~quantile(., probs = 0.75, na.rm = T),
                    max=max))

wilcox.test(eval_stat_q[eval_stat_q$MALDI.station == 'yes', 'n_ribos_detected'], 
            eval_stat_q[eval_stat_q$MALDI.station == 'no', 'n_ribos_detected'], 
            paired = F)

# How long the matrix is kept in the routine workflow for 
eval_stat_q['how_long_matrix_stats'] <- ifelse(eval_stat_q$How.long.do.you.keep.the.same.matrix.in.the.workflow %in% c("maximal 1 day", "maximal a week"), "no longer than a week", 
                                             ifelse(is.na(eval_stat_q$How.long.do.you.keep.the.same.matrix.in.the.workflow.), NA, 'Other'))

eval_stat_q %>% 
  group_by(how_long_matrix_stats) %>%
  summarise_at(vars(n_ribos_detected),
               list(min=min, Q1=~quantile(., probs = 0.25, na.rm = T),
                    median=median, Q3=~quantile(., probs = 0.75, na.rm = T),
                    max=max))

wilcox.test(eval_stat_q[eval_stat_q$how_long_matrix_stats == 'no longer than a week', 'n_ribos_detected'], 
            eval_stat_q[eval_stat_q$how_long_matrix_stats == 'Other', 'n_ribos_detected'], 
            paired = F)

eval_stat_q$how_long_matrix_stats <- NULL

# Using the functions defined at the beginning of this script, first summarise the mass spectral features, grouped by answers to questions, associated with the biggest changes in mass spectral quality and MALD system and plot them, highlighting in green the 'better' practice or in red the 'worde' practice
# summarise by how the strains were processed for these measurements
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = How.did.you.process.the.received.strains., grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "How.did.you.process.the.received.strains.", grouping_variable2 = "MALDI.system")
# a warning appears: NA introduced by coersion. There are 91 rows who have a NA 'value' entry. These are the 'frac_peaks_repr_2spectra' from device_30 and device_42, were only one technical replicate was measured in baseline quality assessment
# Define the order in which the different procedures are plotted
plot_data$How.did.you.process.the.received.strains.<-factor(plot_data$How.did.you.process.the.received.strains., c("freeze.streak.out", "freeze.streak.out.subculture.once", "direct.streak.out.from.eswab", "streak.out.from.eswab.subcultured.once"), 
                                                              labels = c('procedure A', 'procedure B', 'procedure C', 'procedure D'))
# Define the margins where to plot the red highlight
plot_data['highlight']<-ifelse(plot_data$How.did.you.process.the.received.strains. == "procedure A", 1, 0.4)

# plot
f7_how_processed <- plot_plot_data_highlight(plot_data = plot_data, grouping_variable = How.did.you.process.the.received.strains., highlight_color = "#E76F51")
# Remove legend, add title
f7_how_processed <- f7_how_processed + theme(legend.position = 'none') + ggtitle('Processing of the strains')
# add p-values
f7_how_processed <- f7_how_processed + stat_compare_means(method = 'wilcox', paired = F, label = 'p.signif', label.y.npc = 0.9, ref.group = 'procedure A', alpha = c(rep(1, 10),rep(0, 5)))

# Repeat this procedure for multiple answers and inspect plots
# How many samples are measured in the laboratory
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = How.many.MALDI.TOF.MS.measurements.per.day.have.been.performed.on.average.within.the.last.month., grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "How.many.MALDI.TOF.MS.measurements.per.day.have.been.performed.on.average.within.the.last.month.", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = How.many.MALDI.TOF.MS.measurements.per.day.have.been.performed.on.average.within.the.last.month.)
f7_how 

# When are the samples measured after they have been applied onto the target plate
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = When.are.the.MALDI.TOF.MS.targets.measured., grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "When.are.the.MALDI.TOF.MS.targets.measured.", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = When.are.the.MALDI.TOF.MS.targets.measured.)
f7_how

# Which MALDI devices was used
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = MALDI.type, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "MALDI.type", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = MALDI.type)
f7_how

# Which target plates where used for these measurements
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = Which.target.plates.do.you.use., grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "Which.target.plates.do.you.use.", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
plot_data$Which.target.plates.do.you.use.<-factor(plot_data$Which.target.plates.do.you.use., levels = c("Steel target", "Disposable"))
f7_which_target <- plot_plot_data_highlight(plot_data = plot_data, grouping_variable = Which.target.plates.do.you.use., highlight_color = "#2A9D8F")
f7_which_target <- f7_which_target + ggtitle('Target plates used') + stat_compare_means(method = 'wilcox', paired = F, label = 'p.signif', label.y.npc = 0.9, ref.group = 'Steel target', alpha = c(rep(1, 4), rep(0, 2)))

# How these target plates are routinely cleaned
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = How.clean.target.plate, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "How.clean.target.plate", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
plot_data <- plot_data[!plot_data$How.clean.target.plate == 'Disposable', ]
plot_data$How.clean.target.plate<-factor(plot_data$How.clean.target.plate, levels = c("Methanole-Acetone", "TFA","Other"))
f7_how_clean <- plot_plot_data_highlight(plot_data = plot_data, grouping_variable = How.clean.target.plate, highlight_color = "#2A9D8F")
f7_how_clean <- f7_how_clean + ggtitle('Cleaning protocol for steel targets') + stat_compare_means(method = 'wilcox', paired = F, label = 'p.signif', label.y.npc = 0.9, ref.group = 'Methanole-Acetone', alpha = c(rep(1, 6), rep(0, 3)))

targetplot<-cowplot::plot_grid(f7_how_processed + theme(axis.text=element_text(size=8), strip.text.y = element_text(size = 5.5), axis.title=element_text(size=8), legend.position = 'none'),
                              f7_which_target+ theme(axis.text=element_text(size=8), strip.text.y = element_text(size = 5.5), axis.title=element_text(size=8), legend.position = 'none'),
                               f7_how_clean+ theme(axis.text=element_text(size=8), strip.text.y = element_text(size = 5.5), axis.title=element_text(size=8), legend.position = 'none'),
                               rel_widths = c(0.2, 0.14, 0.18), ncol = 3, align = 'h')

#pdf('./questionnaire_targets_2.pdf', width = 8.28, height = 4)
pdf('./questionnaire_targets_2b.pdf', width = 8.28, height = 6)
targetplot
dev.off()

# Whether wooden toothpicks or plastic inoculation needles are used to transfer bacterial material from the media onto the target plate
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = What.do.you.use.for.applying.the.bacterial.colony.onto.the.MALDI.TOF.MS.target.plate., grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "What.do.you.use.for.applying.the.bacterial.colony.onto.the.MALDI.TOF.MS.target.plate.", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = What.do.you.use.for.applying.the.bacterial.colony.onto.the.MALDI.TOF.MS.target.plate.)
f7_how

# Whether the mass spectral quality is routinely assessed
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = Regular.quality.assessment, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "Regular.quality.assessment", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = Regular.quality.assessment)
f7_how

# Whether there are regular hardware services performed to maintain the device
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = regular.hardware.service, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "regular.hardware.service", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
plot_data$regular.hardware.service<-factor(plot_data$regular.hardware.service, levels = c("yes", "no"))
f7_hardware <- plot_plot_data_highlight(plot_data = plot_data, grouping_variable = regular.hardware.service, highlight_color = "#2A9D8F")
f7_hardware <- f7_hardware + ggtitle('Regular hardware services') + stat_compare_means(method = 'wilcox', paired = F, label = 'p.signif', label.y.npc = 0.9, ref.group = 'yes', alpha = c(rep(1, 2), rep(0, 1)))
# a warning will be displayed stating tha the computation failed. this comes from the fact that there is no two groups to compare in the vitekMS group, all laboratories stated regularly performing a hardware service. 

# Whether the technicians operating the MALDI undergo regulat sample preparation training
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = samplepre_training, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "samplepre_training", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = samplepre_training)
f7_how

# Whether strains measured in routine diagnostics are grown on CHROM agar or not
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = agar.chrom, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "agar.chrom", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = agar.chrom)
f7_how

# Whether strains measured in routine diagnostics are grown on Sheep Blood agar or not
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = agar.sheep.blood, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "agar.sheep.blood", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = agar.sheep.blood)
f7_how

# Whether strains measured in routine diagnostics are grown on 'Other' (other than sheep blood and CHROM agar) agar or not
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = agar.other, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "agar.other", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = agar.other)
f7_how

# How long the matrix is kept in the routine workflow for
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = How.long.do.you.keep.the.same.matrix.in.the.workflow., grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "How.long.do.you.keep.the.same.matrix.in.the.workflow.", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
plot_data$How.long.do.you.keep.the.same.matrix.in.the.workflow.<-factor(plot_data$How.long.do.you.keep.the.same.matrix.in.the.workflow., levels = c('maximal 1 day', 'maximal a week', 'maximal a month'), labels = c('max. 1 day','max. 1 week','max. 1 month'))
f7_how_long_matrix <- plot_plot_data_highlight(plot_data = plot_data, grouping_variable = How.long.do.you.keep.the.same.matrix.in.the.workflow., highlight_color = "#2A9D8F")
f7_how_long_matrix <- f7_how_long_matrix + ggtitle('Duration of matrix usage') + stat_compare_means(method = 'wilcox', paired = F, label = 'p.signif', label.y.npc = 0.9, ref.group = 'max. 1 day', alpha = c(rep(1, 8), rep(0, 4)))

# Whether routinely, unknown samples are overlaid with formic acid or not
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = FA.applied.routinely.bacteria, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "FA.applied.routinely.bacteria", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
f7_how<-plot_plot_data(plot_data = plot_data, grouping_variable = FA.applied.routinely.bacteria)
f7_how

# Whether the participating laboratory works with a MALDI workstation (for a certain period of time, one person perocessed all sample for MALDI measurements)
eval2<-ID_eval(eval_file = eval_stat_q, grouping_variable1 = MALDI.station, grouping_variable2 = MALDI.system)
plot_data<-summarise_prop_and_eval(eval_file = eval_stat_q, ID_eval_out = eval2, grouping_variable1 = "MALDI.station", grouping_variable2 = "MALDI.system")
# remove the 91 rows with NA values which correspond to the frac_peaks_repr_2spectra fro device 30 and device 42 were only one technical replicte was measured for baseline quality assessment
plot_data <- plot_data[!is.na(plot_data$value),]
plot_data$MALDI.station<-factor(plot_data$MALDI.station, levels = c("yes", "no"))
f7_MALDI_station <- plot_plot_data_highlight(plot_data = plot_data, grouping_variable = MALDI.station, highlight_color = "#2A9D8F")
f7_MALDI_station <- f7_MALDI_station + ggtitle('Rotating MALDI workstation')+ stat_compare_means(method = 'wilcox', paired = F, label = 'p.signif', label.y.npc = 0.9, ref.group = 'yes', alpha = c(rep(1, 4), rep(0, 2)))

# summarise th answers associated with the biggest changes in mass spectral quality into a plot
metaplot<-cowplot::plot_grid(f7_hardware + theme(axis.text=element_text(size=8), strip.text.y = element_text(size = 6), axis.title=element_text(size=8), legend.position = 'none'),
                             f7_MALDI_station + theme(axis.text=element_text(size=8), strip.text.y = element_text(size = 6), axis.title=element_text(size=8), legend.position = 'none'),
                             f7_how_long_matrix + theme(axis.text=element_text(size=8), strip.text.y = element_text(size = 6), axis.title=element_text(size=8), legend.position = 'none'),
                             rel_widths = c(0.16, 0.16, 0.2), ncol = 3, align = 'h')
# a warning will be displayed stating the computation failed. this comes from the fact that there is no two groups to compare in the vitekMS group, all laboratories stated regularly performing a hardware service. 

#pdf('./questionnaire_metaplot_2.pdf', width = 8.28, height = 4)
pdf('./questionnaire_metaplot_2b.pdf', width = 8.28, height = 6)
metaplot
dev.off()

