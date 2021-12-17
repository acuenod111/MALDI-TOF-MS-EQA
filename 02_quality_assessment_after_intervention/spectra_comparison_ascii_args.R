
######## This script summarises the following information for the test spectra acquired: total number of peaks, number of predicted masses detected (absolute and relative) (every spectra is exclusively queried for the subunits of the relative genome), average distance of predicted to detected, highest mass
# The following input arguments are required:
# (i) directory to ascii 
# (ii) prefix of the lab, eg. 'device_01'
## same for all (iii) path to the id file (Species ID either from from the microflex biotyper DB or from the VitekMS DB)
## same for all (iv) path to the file 'predicted_masses.csv' including all predicted masses and produced by the script 'import_GAP.R'
## same for all (v) path to the file 'Strains_numbering_first_Shihpment.csv'
# (vi) name and path of the summary output file (1 row per spectra)
# (vii) name and file to the to the per subunits output file (one row per detected subunit)
# When running the script, a warning message will be evoked through MALDIQuant which is triggered by empty spectra. This can be ignored. 


library('MALDIquant')
library('MALDIquantForeign')
library('stringr')
library('dplyr')
library('purrr')
library('tidyr')

# define arguments
args = commandArgs(trailingOnly=TRUE)

# import peak lists in ascii format
asciis_dir<-list.files(
  path = args[1],        # directory to search within
  pattern = "*.txt$", 
  recursive = TRUE,         
  full.names = TRUE         
)

asciis = lapply(asciis_dir, read.delim)

# define empty lists for the variables which should be extracted
spectra_name<-list()
peaks <- list()
peaks2<-list()
spectra_name_run<-list()
run_raw<-list()
run<-list()
# loop though each file, add normalised intensity and reformat to dataframe
for (i in 1:length(asciis)){
  spectra_name<-as.character(asciis[[i]][1,])
  spectra_name<-gsub('.*Sample\\=','',spectra_name)
  spectra_name<-gsub(' PSD.*\\#$', '', spectra_name)
  spectra_name<-gsub('\\#PSD\\=\\#$', '', spectra_name)
  run_raw<- gsub(' $', '', gsub('COM\\=','', asciis[[i]][4,]))
  run<-ifelse(grepl('\\\\', run_raw), str_split(run_raw, '\\\\')[[1]][length(str_split(run_raw, '\\\\')[[1]])-3], gsub('\\_[[:digit:]][[:alpha:]][[:digit:]]{1,2}$','',run_raw))
  spectra_name_run<-append(spectra_name_run, paste(spectra_name, run, sep=';'))
  peaks[[spectra_name]]<-separate(as.data.frame(asciis[[i]][-c(1:6),]), col="asciis[[i]][-c(1:6), ]", into = c("mass", "intensity_raw"), sep=" ") %>% unique() %>%
    mutate(mass =  as.numeric(as.character(mass))) %>%
    mutate(intensity_norm =  as.numeric(intensity_raw)/median(as.numeric(intensity_raw)))
  # set mass range to 3'000 - 20'000
  peaks2[[spectra_name]]<-createMassPeaks(peaks[[spectra_name]]$mass[peaks[[spectra_name]]$mass > 3000 & peaks[[spectra_name]]$mass < 20000], peaks[[spectra_name]]$intensity_norm[peaks[[spectra_name]]$mass > 3000 & peaks[[spectra_name]]$mass < 20000])
}
# merge all to the same dataframe
spectra_name_run_df<-do.call(rbind.data.frame, spectra_name_run)
spectra_name_run_df<-str_split_fixed(spectra_name_run_df[,1], ";", 2)
# split run and spectra into two columns
colnames(spectra_name_run_df)<-c('spectra', 'run')
spectra_name_run_df<-as.data.frame(spectra_name_run_df)

#extract the number of peaks, the samplename, the highest peak recorded, the peak recorded at the 90th percentile, the number of peaks above 10'000 and the total intensity per spectrum
n_peaks<-list()
sample_name<-list()
highest_peaks<-list()
sample_with_position<-list()
mass_90<-list()
n_high_peaks<-list()
total_intensity_peaks<-list()
for (i in 1:length(peaks)){
   n_peaks[i]<-length(peaks2[[i]]@mass)
   highest_peaks[i]<-ifelse(n_peaks[i] != 0, max(peaks2[[i]]@mass), NA)
   mass_90[i]<-ifelse(n_peaks[i] != 0, quantile(as.numeric(peaks2[[i]]@mass),probs = .90), NA)
   n_high_peaks[i]<-sum(peaks2[[i]]@mass>10000)
   total_intensity_peaks[i]<-sum(as.numeric(peaks[[i]]$intensity_raw))
}
# merge to one df
n_peaks<-data.frame(n_peaks = as.character(n_peaks), highest_peaks = as.character(highest_peaks), mass_90 = as.character(mass_90), spectra = names(peaks), n_high_peaks = as.character(n_high_peaks), total_intensity_peaks = as.character(total_intensity_peaks))
# add strainnumber and position
n_peaks['strain_number']<-ifelse(grepl(paste0(args[2],'\\.(1|2|3)\\.\\d{1,2}.*'),n_peaks$spectra), gsub(paste0('(',args[2],'\\.(1|2|3)\\.)','(\\d{1,2})(.*)'),'\\3', n_peaks$spectra), NA)
n_peaks['method']<-ifelse(grepl(paste0(args[2],'\\.(1|2|3)\\.\\d{1,2}.*'),n_peaks$spectra), gsub(paste0('(',args[2],'\\.)((1|2|3))','(\\.\\d{1,2})(.*)'),'\\2', n_peaks$spectra), NA)
n_peaks['position']<-ifelse(grepl('.*[[:punct:]]0\\_[[:alpha:]][[:digit:]]{1,2}[[:punct:]]*.*', n_peaks$spectra), gsub('(.*[[:punct:]])(0\\_[[:alpha:]][[:digit:]]{1,2})([[:punct:]]*.*)', '\\2', n_peaks$spectra), NA)
n_peaks['position']<-ifelse(is.na(n_peaks$position) & grepl('(.*\\.)([[:digit:]][[:alpha:]][[:digit:]]{1,2})(\\..*)', n_peaks$spectra), gsub('(.*\\.)([[:digit:]][[:alpha:]][[:digit:]]{1,2})(\\..*)', '\\2',n_peaks$spectra), 
                            ifelse(is.na(n_peaks$position) & !grepl('(.*\\.)([[:digit:]][[:alpha:]][[:digit:]]{1,2})(\\..*)', n_peaks$spectra), NA, n_peaks$position))
n_peaks['position']<-ifelse(is.na(n_peaks$position), gsub('(.*\\.)([[:digit:]]*[[:alpha:]][[:digit:]])(\\..*)', '\\2', n_peaks$spectra), n_peaks$position)
# add run
n_peaks<-merge(spectra_name_run_df, n_peaks, by='spectra', all.x = T)
# add brukercode, if there
n_peaks['brukercode']<-ifelse(grepl('[[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12}', n_peaks$spectra), gsub('(.*)([[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12})(.*)', '\\2', n_peaks$spectra), NA)
n_peaks['brukercode']<-if(grepl('043_Eurofins', args[2])){NA}else{as.character(n_peaks$brukercode)}
n_peaks['to_exclude']<-ifelse(!is.na(n_peaks$brukercode), paste(n_peaks$strain_number,n_peaks$brukercode, n_peaks$position), 
                                  paste(n_peaks$method, n_peaks$strain_number,n_peaks$position))

#exclude background control
n_peaks<-n_peaks[!is.na(n_peaks$strain_number),]
n_peaks['strain_number']<-ifelse(nchar(n_peaks$strain_number)==1, as.character(paste0('0', n_peaks$strain_number)), as.character(n_peaks$strain_number))
n_peaks['spectra']<-gsub('\\#PSD\\=\\#$', '', as.character(n_peaks$spectra))

#exclude mixup / contaminations. for this import the species identification of the microflex Bruker / the VitekMS database
path_to_id<-args[3]
id_path <- list.files(path = path_to_id, pattern = paste0(args[2], "\\_.*","\\.csv"))

# if no bruker ID / VitekMS ID yet, skip this step
if (length(id_path) > 0){
  if(length(id_path) > 1){
    id_list<-lapply(paste0(path_to_id,id_path), read.csv, sep=';', stringsAsFactors = FALSE, header=TRUE)
    vars <- if(grepl('bruker', args[3])){c("strainnumber", "col","Strain","species_NGS","genus_NGS")} else{ c("strainnumber", "Strain","species_NGS","genus_NGS")}
    id_list<-lapply(id_list, function(df) mutate_at(df, .vars = vars, as.character))
    id<-id_list %>% reduce(full_join, by = intersect(colnames(id_list[1]), colnames(id_list[2])))
    id<-id[!is.na(id$strainnumber),]
  } else{
    id<-read.csv(paste0(path_to_id,id_path), sep =';')
  }
  if(any(id[!is.na(id$include), 'include']==FALSE)){
    id_exclude<-id[id$include == FALSE,]
    id_exclude<-id_exclude[!is.na(id_exclude$strainnumber),]
    if(grepl('brukerreport', args[3])){
      id_exclude['strainnumber']<-ifelse(nchar(as.character(id_exclude$strainnumber))==1, as.character(paste0('0', id_exclude$strainnumber)), as.character(id_exclude$strainnumber))
      id_exclude['brukercode']<-ifelse(grepl('[[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12}', id_exclude$col), gsub('(.*)([[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[:alnum:]]{12})(.*)', '\\2', id_exclude$col), NA)
      id_exclude['to_exclude']<-ifelse(!is.na(id_exclude$brukercode), paste(id_exclude$strainnumber,id_exclude$brukercode, id_exclude$position), 
                                           paste(id_exclude$method, id_exclude$strainnumber,id_exclude$position))
      n_peaks<-n_peaks[!(n_peaks$to_exclude %in% id_exclude$to_exclude),]
      peaks2<-peaks2[!(n_peaks$to_exclude %in% id_exclude$to_exclude)]
    } else { 
      id_exclude['samplename_new'] <- gsub('.txt','', id_exclude$files)
      n_peaks<-n_peaks[!(n_peaks$spectra %in% id_exclude$samplename_new),]
      peaks2<-peaks2[!(n_peaks$spectra %in% id_exclude$samplename_new)]}}
  } else{
    print('no brukerreport / Vitekreport yet, this set of spectra might still include contaminations / mixup')
  }

n_peaks$brukercode<-NULL
n_peaks$to_exclude<-NULL

# extract which ribosomal subunits have been detected, which have been detected and what was the average distance to detected ribos
# read GAP
ribos<-read.csv2(args[4])
numbering<-read.csv2(args[5], sep=',')

ribos<-merge(ribos, numbering, by = 'Strain', all.x = T)
ribos<-ribos[ribos$within_mass_range=='TRUE',]
ribos['Numbering_Shipment_1']<-ifelse(nchar(ribos$Numbering_Shipment_1)==1, paste0('0', ribos$Numbering_Shipment_1), ribos$Numbering_Shipment_1)

# set error
error <- 800
ppm <- error / 1000000

# define variables which are substracted as empty lists
subunit_detected_all_in_one_list<-list()
subunit_detected_all_in_one_each_su<-data.frame()
mass_detected_all_in_one_each_su<-data.frame()
intensity_detected_all_in_one_each_su<-data.frame()
n_predicted_su_all <- list()
frac_peaks_repr<-list()
frac_peaks_repr_2<-list()

# for each strain, loop through all technical replicates and look for peaks within error range of the predicted mass
for (strain in unique(ribos$Strain)){ # Loop through all strains 
  strain_number <- unique(ribos[ribos$Strain == strain, 'Numbering_Shipment_1'])
  peaks_of_interest<-if (any(n_peaks$strain_number == strain_number)){ #only compare the spectra to the predicted ribos of that strain # add if statement, if strain is missing, no error
    peaks2[names(peaks2) %in% n_peaks[n_peaks$strain_number == strain_number,'spectra']] 
  } else {
    next
  }
  n_peak_of_interest<-if (any(n_peaks$strain_number == strain_number)){ #only compare the spectra to the predicted ribos of that strain # add if statement, if strain is missing, no error
    n_peaks[n_peaks$strain_number == strain_number,] 
  } else {
    next
  }
  subunit_detected_per_strain<-list()
  
  for (i in 1:length(peaks_of_interest)){
    subunit_detected<- list()
    mass_detected<- list()
    intensity_detected<- list()
    sample_with_position<-names(peaks_of_interest)[[i]]
    # add % subunits detected (not all strains have the same number of predicted masses)
    n_predicted_su <- length(ribos[ribos$Strain == strain, 'Subunit'])
    
    for (subunit in ribos[ribos$Strain == strain, 'Subunit']){ # check for the presence of a peak for each subunit (mass +/- error range)
      mass <- as.numeric(as.character(ribos[ribos$Strain == strain & ribos$Subunit == subunit, 'Mass']))
      subunit_temp<-if (any(peaks_of_interest[[i]]@mass<(mass+(mass*ppm)) &  peaks_of_interest[[i]]@mass>(mass-(mass*ppm)))){
        paste(subunit, mass) # If nothing  is detected include NA, so that if non is detected, still someting can be appended
      } else {
        NA
      }
      subunit_detected <- append(subunit_detected, subunit_temp)
      subunit_detected<- unlist(unique(subunit_detected[!is.na(subunit_detected)])) # remove where not detected
      
      mass_detected_temp <- if (!is.na(subunit_temp)){
        as.character(peaks_of_interest[[i]]@mass[peaks_of_interest[[i]]@mass<(mass+(mass*ppm)) &  peaks_of_interest[[i]]@mass>(mass-(mass*ppm))])
      } else {
        NA
      }
      mass_detected_temp <- paste(subunit, mass_detected_temp)
      mass_detected <- append(mass_detected, mass_detected_temp)
      mass_detected <- mass_detected[!grepl('NA$',mass_detected)] # remove where not detected
      
      intensity_detected_temp <- if (!is.na(subunit_temp)){
        as.character(peaks_of_interest[[i]]@intensity[peaks_of_interest[[i]]@mass<(mass+(mass*ppm)) &  peaks_of_interest[[i]]@mass>(mass-(mass*ppm))])
      } else {
        NA
      }
      intensity_detected_temp <- paste(subunit, intensity_detected_temp)
      intensity_detected <- append(intensity_detected, intensity_detected_temp)
      intensity_detected <- intensity_detected[!grepl('NA$',intensity_detected)] # remove where not detected
      
    }
    
    subunit_detected_all_in_one_list[[sample_with_position]] <- paste(subunit_detected, collapse = ',') # summarise
    if(length(peaks_of_interest[[i]]@mass)> 0){
                frac_peaks_repr[[sample_with_position]]<-length(filterPeaks(binPeaks(peaks_of_interest, method='relaxed', tolerance = ppm), minFrequency=0.51)[[i]]@mass) / length(peaks_of_interest[[i]]@mass)
              } else {
                frac_peaks_repr[[sample_with_position]]<-0
              }
              if(length(peaks_of_interest[[i]]@mass) > 0 & length(peaks_of_interest[[length(peaks_of_interest)+1 -i]])){ # read out the reproducibility of two random spectra of the same technical replicate
                frac_peaks_repr_2[[sample_with_position]]<-length(filterPeaks(binPeaks(peaks_of_interest[c(i, length(peaks_of_interest)+1-i)], method='relaxed', tolerance = ppm), minFrequency=0.51)[[1]]@mass) / length(peaks_of_interest[[i]]@mass)
              }
              else{
                frac_peaks_repr_2[[sample_with_position]]<-0
              }
    n_predicted_su_all[strain] <- n_predicted_su
    subunit_detected_all_in_one_each_su <- rbind(subunit_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(subunit_detected), " "))) %>% mutate(spectra = sample_with_position))
    mass_detected_all_in_one_each_su <- rbind(mass_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(mass_detected), " "))) %>% mutate(spectra = sample_with_position))
    intensity_detected_all_in_one_each_su <- rbind(intensity_detected_all_in_one_each_su,as.data.frame(do.call(rbind, strsplit(as.character(intensity_detected), " "))) %>% mutate(spectra = sample_with_position))
  }
}

# summarise
colnames(subunit_detected_all_in_one_each_su)<-c('Subunit', 'predicted_mass', 'spectra')
colnames(mass_detected_all_in_one_each_su)<-c('Subunit', 'detected_mass', 'spectra')
colnames(intensity_detected_all_in_one_each_su)<-c('Subunit', 'intensity', 'spectra')

# merge all subunits
read_out_each_su<-merge(subunit_detected_all_in_one_each_su, mass_detected_all_in_one_each_su, by=c(intersect(colnames(subunit_detected_all_in_one_each_su), colnames(mass_detected_all_in_one_each_su))))
read_out_each_su<-merge(read_out_each_su, mass_detected_all_in_one_each_su, by=c(intersect(colnames(read_out_each_su), colnames(mass_detected_all_in_one_each_su))))
read_out_each_su<-merge(read_out_each_su, intensity_detected_all_in_one_each_su, by=c(intersect(colnames(read_out_each_su), colnames(intensity_detected_all_in_one_each_su))))

# add 'dist' columns, specifiying what the distance between the predicted mass of the subunit and the detected peak was
read_out_each_su['dist']<-abs(as.numeric(as.character(read_out_each_su$predicted_mass)) - as.numeric(as.character(read_out_each_su$detected_mass)))
read_out_each_su['dist_ppm']<-(read_out_each_su$dist / as.numeric(as.character(read_out_each_su$predicted_mass))) *1000000

# If multiple peaks per subunit detected, only keep the one with higher intensity
read_out_each_su<-read_out_each_su %>% group_by(spectra,Subunit) %>%
  filter(intensity == max(intensity)) %>%
  filter(dist_ppm == min(dist_ppm))

# summarise values for each spectrum
read_out_each_su_sum<-read_out_each_su %>% group_by(spectra) %>%
  summarize(mean_dist = mean(as.numeric(as.character(dist))), mean_dist_ppm = mean(as.numeric(as.character(dist_ppm))), mean_intensity = mean(as.numeric(as.character(intensity))), meadian_intensity =  median(as.numeric(as.character(intensity))), max_mass_ribo = max(as.numeric(as.character(predicted_mass))))

# Add fraction of reproducibly detected peaks to the 'n_peaks' dataframe
frac_peaks_repr_df<-do.call(rbind, frac_peaks_repr)
colnames(frac_peaks_repr_df)<-'frac_peaks_repr'
frac_peaks_repr_df<- cbind(spectra = rownames(frac_peaks_repr_df), data.frame(frac_peaks_repr_df, row.names=NULL))
frac_peaks_repr_df_2spectra<-do.call(rbind, frac_peaks_repr_2)
colnames(frac_peaks_repr_df_2spectra)<-'frac_peaks_repr_2spectra'
frac_peaks_repr_df_2spectra<- cbind(spectra = rownames(frac_peaks_repr_df_2spectra), data.frame(frac_peaks_repr_df_2spectra, row.names=NULL))
frac_peaks_repr_df<- merge(frac_peaks_repr_df, frac_peaks_repr_df_2spectra, by = 'spectra')
# merge this to the n_peaks df
n_peaks<-merge(n_peaks, frac_peaks_repr_df, by='spectra')
# add the information on the ribosomal subunits which have been detected
ribos_detected<-do.call(rbind, subunit_detected_all_in_one_list)
colnames(ribos_detected)<-'subunits_detected'
ribos_detected<- cbind(spectra = rownames(ribos_detected), data.frame(ribos_detected, row.names=NULL))
# extract strainnumber from filename
ribos_detected['strain_number']<-ifelse(grepl(paste0(args[2],'\\.(1|2|3)\\.\\d{1,2}.*'),ribos_detected$spectra), gsub(paste0('(',args[2],'\\.(1|2|3)\\.)','(\\d{1,2})(.*)'),'\\3', ribos_detected$spectra), NA)
ribos_detected['strain_number']<-gsub('(^\\d{1}$)','0\\1', ribos_detected$strain_number)
numbering['Numbering_Shipment_1']<-gsub('(^\\d{1}$)','0\\1', numbering$Numbering_Shipment_1)
ribos_detected<-merge(ribos_detected, numbering, by.x = 'strain_number', by.y = 'Numbering_Shipment_1', all.x = T)

# Count subunits detected. There is one more subunit detected that there are commas, so add 1 to count the subunits detected
ribos_detected['n_ribos_detected']<-ifelse(ribos_detected$subunits_detected == "",0,str_count(ribos_detected$subunits_detected, ',') + 1)
n_predicted_su_all<-t(as.data.frame(n_predicted_su_all))
colnames(n_predicted_su_all)<-'n_predicted_su'
rownames(n_predicted_su_all)<-gsub('\\.', '-', rownames(n_predicted_su_all))
rownames(n_predicted_su_all)<-gsub('\\-TS\\-', '\\(TS\\)', rownames(n_predicted_su_all))
ribos_detected<-merge(ribos_detected, n_predicted_su_all, by.x = 'Strain', by.y = "row.names", all=T)
ribos_detected['rel_amount_su_detected']<-ribos_detected$n_ribos_detected / ribos_detected$n_predicted_su

# merge
read_out<-merge(n_peaks, ribos_detected, by=intersect(colnames(n_peaks), colnames(ribos_detected)), all.x = T)
read_out<-merge(read_out, read_out_each_su_sum, by = intersect(colnames(read_out), colnames(read_out_each_su_sum)), all.x = T)

#export
write.table(read_out, args[6], row.names=FALSE, quote = FALSE, dec = '.', sep=';')
write.table(read_out_each_su, args[7], row.names=FALSE, quote = FALSE, dec = '.', sep=';')
