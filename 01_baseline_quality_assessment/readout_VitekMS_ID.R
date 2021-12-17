# Read out Vitek MS
# add 'include' column, add, how many species identified, spectraname and strainnumber and merge with group
## This summarises the the VitekMS species identification (one line per spectrum). It requires the following four arguments 
#(i) input path to the VitekMS .csv file file
#(ii) input to the .csv file containing the translation from sampleID to strainnumber
#(iii) Input file which which strainnumber is which species per NGS.
#(iv) output file to the csv file

args = commandArgs(trailingOnly=TRUE)

# read in VitekMS report
vitekreport<-read.csv2(args[1], sep=';')
#exclude empty lines
vitekreport<-vitekreport[!(vitekreport$files=='' |grepl('MICROFLEX', vitekreport$files)),]
#add samplename
vitekreport['samplename_old']<-gsub('\\.mz.+$|\\.txt$', '', vitekreport$files)
vitekreport['samplename_old']<-if(grepl('device_12', args[2])){gsub('([[:digit:]]{2}\\-[[:digit:]]{2}\\-[[:digit:]]{2}\\_)(micro\\_.*)', '\\2', as.character(vitekreport$samplename_old))
} else {
  if(grepl('device_11', args[2])){gsub('([[:digit:]]{2}\\-[[:digit:]]{2}\\-[[:digit:]]{2}\\_Labo1(\\_R|\\_E))(.*)', '\\3', as.character(vitekreport$samplename_old))}
  else{ 
    if(grepl('device_30', args[2])){gsub('([[:digit:]]{2}\\-[[:digit:]]{2}\\-[[:digit:]]{2}\\_FSMLab\\_)(.*)', '\\2', as.character(vitekreport$samplename_old))}
  else{
    if(grepl('device_15', args[2])){gsub('([[:digit:]]{2}\\-[[:digit:]]{2}\\-[[:digit:]]{2}\\_device_15\\_)(.*)', '\\2', as.character(vitekreport$samplename_old))}
    else{as.character(vitekreport$samplename_old)}}}}
    
#add strainnumber and species identified
sampleID<-read.csv2(args[2], sep=',', stringsAsFactors = FALSE)
sampleID['samplename_old']<-gsub(' $', '', sampleID$samplename_old)
sampleID['samplename_old']<-gsub('\\.mz.+$', '', sampleID$samplename_old)
vitekreport<-merge(vitekreport, sampleID, by='samplename_old')
vitekreport['strainnumber']<-sub('([^\\.]+)(\\.)([[:digit:]]{1,2})(.+)', '\\3', vitekreport$samplename_new)
vitekreport['strainnumber']<-if(grepl('ascii_axima_eqa', args[2])){gsub('(.*)((AXI\\-(WP1\\-)*)|(brukerextraction\\_3\\.))([[:digit:]]{1,2})(.*)', '\\6', vitekreport$samplename_new)} else {vitekreport$strainnumber}
vitekreport['strainnumber']<-ifelse(nchar(vitekreport$strainnumber) == 1, paste0('0', vitekreport$strainnumber), vitekreport$strainnumber)
vitekreport['lab']<-sub('([^\\.]+)(\\.)([[:digit:]]{1,2})(.+)', '\\1', vitekreport$samplename_new)

# correct for S- lutetiensis 
vitekreport<-data.frame(lapply(vitekreport, function(x) {
  gsub('Streptococcus infantarius ssp coli \\(Str.lutetiensis\\)', 'Streptococcus lutetiensis', x)
                }))

vitekreport['species_identified']<-sapply(strsplit(as.character(vitekreport$species_rank_1)," "), `[`, 2)
vitekreport['genus_identified']<-sapply(strsplit(as.character(vitekreport$species_rank_1)," "), `[`, 1)

# remove calibration spectra
vitekreport<-vitekreport[!grepl('Ecal.*x', vitekreport$samplename_new),]
vitekreport<-vitekreport[!vitekreport$strainnumber == '00',]

# add Species and strain by NGS
numbering<-read.csv2(args[3], sep=',')
numbering['species_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\3', numbering$Strain)
numbering['genus_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\1', numbering$Strain)
numbering['genus_NGS']<-gsub('Winkia', 'Actinomyces', numbering$genus_NGS)
numbering['Numbering_Shipment_1']<-ifelse(nchar(numbering$Numbering_Shipment_1)==1, paste0('0', numbering$Numbering_Shipment_1), numbering$Numbering_Shipment_1)

vitekreport<-merge(vitekreport, numbering, by.x ='strainnumber', by.y = 'Numbering_Shipment_1', all.x=TRUE)

# Count how many species have been identified
vitekreport['n_species']<-ifelse(!is.na(vitekreport$species_rank_4) & vitekreport$species_rank_4 !='', '4', 
                                            ifelse(!is.na(vitekreport$species_rank_3) & vitekreport$species_rank_3 !='', '3', 
                                                   ifelse(!is.na(vitekreport$species_rank_2) & vitekreport$species_rank_2 !='', '2', 
                                                          ifelse(!is.na(vitekreport$species_rank_1) & vitekreport$species_rank_1!='', '1', '0'))))
# check difference in probability between 1st and 2nd match
vitekreport['diff']<-ifelse(vitekreport$n_species > 1, abs(as.numeric(as.character(vitekreport$score_rank_1)) - as.numeric(as.character(vitekreport$score_rank_2))), NA)
vitekreport['diff_proba']<-ifelse(vitekreport$n_species > 1, abs(as.numeric(as.character(vitekreport$proba_rank_1)) - as.numeric(as.character(vitekreport$proba_rank_2))), NA)

# check if correct genus amongst identified genera
vitekreport['correct_genus_identified']<-mapply(grepl, vitekreport$genus_NGS, vitekreport$genus_identified)
vitekreport['correct_species_identified']<-mapply(grepl, vitekreport$species_NGS, vitekreport$species_identified)

# add position
vitekreport['position']<-if(grepl('device_15|device_11|device_12|device_30', args[2])){
  gsub('(.*\\_)([[:alpha:]][[:digit:]])(\\_.*)', '\\2', vitekreport$samplename_new)
} else {
  gsub('(.*(\\_|\\.))([[:digit:]][[:alpha:]][[:digit:]]{1,2})(\\.txt)*', '\\3', vitekreport$samplename_new)
}

# check which to exclude. Exclude the ones which have been assigned the wrong genus with a probability higher than 70 % with the exception of Escherichia/Shigella, Klebsiella/Raoultella and Klebsiella aerogenes, which is still called Enterobacter aerogenes in the VitekMS DB
vitekreport['include']<-ifelse(vitekreport$correct_genus_identified, TRUE, 
                           ifelse(as.numeric(as.character(vitekreport$proba_rank_1)) < 70, TRUE, #include bad spectra, only inlcude if genus from with more than 1.7 confidence
                                  ifelse(vitekreport$genus_identified == 'Escherichia' & vitekreport$genus_NGS == 'Shigella', TRUE, 
                                         ifelse(vitekreport$genus_identified == 'Raoultella' & vitekreport$genus_NGS == 'Klebsiella', TRUE, 
                                                ifelse(vitekreport$species_NGS =='aerogenes' & vitekreport$genus_identified == 'Enterobacter', TRUE, FALSE)))))

# add empty 'max_SN' and 'max_resolution' if missing
if (any(!grepl('max_SN', colnames(vitekreport)))){vitekreport['max_SN']<-NA}
if (any(!grepl('max_resolution', colnames(vitekreport)))){vitekreport['max_resolution']<-NA}

# only keep columns needed
vitekreport<-vitekreport[,c("strainnumber","samplename_new", "files","position", "nb_peaks", "max_SN","max_intensities","max_resolution","identifType", "species_rank_1",          
                            "proba_rank_1","score_rank_1","lab","Strain","n_species", "diff","diff_proba", 
                            "species_identified" ,"genus_identified","species_NGS","genus_NGS","correct_genus_identified","correct_species_identified","include")]

if (grepl('*device_35|device_36|device_37*', args[2])){
  vitekreport_device_35 <-vitekreport[grepl('device_35\\.', vitekreport$samplename_new),]
  vitekreport_device_35 <- apply(vitekreport_device_35,2,as.character)
  write.csv2(vitekreport_AXIMA_I, args[4])
  vitekreport_evice_36 <-vitekreport[grepl('evice_36\\.', vitekreport$samplename_new),]
  vitekreport_evice_36 <- apply(vitekreport_evice_36,2,as.character)
  write.csv2(vitekreport_AXIMA_II, args[5])
  vitekreport_device_37<-vitekreport[grepl('device_37\\.', vitekreport$samplename_new),]
  vitekreport_device_37 <- apply(vitekreport_device_37,2,as.character)
  write.csv2(vitekreport_AXIMA_III, args[6])
} else {
  vitekreport <- apply(vitekreport,2,as.character)
  write.csv2(vitekreport, args[4])  
}


