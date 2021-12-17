# Read out Vitek MS
# add 'include' column, add, how many species udentified, spectraname and strainnumber and merge with group
## This summarises the the VitekMS species identification (one line per spectrum). It requires the following four arguments 
#(i) input path to the VitekMS .csv file file
#(ii) Input file which shich strainnumber is which species per NGS.
#(iii) output file to the csv file
library('stringr')

args = commandArgs(trailingOnly=TRUE)

# read in VitekMS report
vitekreport<-read.csv2(args[1], sep=';')
# remove lines wchi do not include identification (empty lines)
vitekreport <- vitekreport[!(vitekreport$files == ''),]

#extract lab and strainnumber from filename
vitekreport <- cbind(vitekreport, str_split_fixed(vitekreport$files, '\\.', 6)[,c(1,3)])
colnames(vitekreport) <- c(colnames(vitekreport)[1:(length(colnames(vitekreport))-2)], c('lab', 'strainnumber'))

# correct for S- lutetiensis
vitekreport<-data.frame(lapply(vitekreport, function(x) {
  gsub('Streptococcus infantarius ssp coli \\(Str.lutetiensis\\)', 'Streptococcus lutetiensis', x)
                }))

vitekreport['species_identified']<-sapply(strsplit(as.character(vitekreport$species_rank_1)," "), `[`, 2)
vitekreport['genus_identified']<-sapply(strsplit(as.character(vitekreport$species_rank_1)," "), `[`, 1)

# remove calibration spectra
vitekreport<-vitekreport[!grepl('Ecal.*x', vitekreport$files),]
vitekreport<-vitekreport[!vitekreport$strainnumber == '00',]

# add Species and strain by NGS
numbering<-read.csv2(args[2], sep=',')
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
vitekreport['position']<-if(grepl('device\\_15|device\\_10|device\\_11|device\\_12|device\\_30', args[1])){
  gsub('(.*\\.)([[:alpha:]][[:digit:]])(\\..*)', '\\2', vitekreport$files)
} else {
  gsub('(.*(\\_|\\.))([[:digit:]][[:alpha:]][[:digit:]]{1,2})(.*\\.txt)*', '\\3', vitekreport$files)
}

# check which to exclude. Eclude the ones which habve been assigned the wrong genus with a probability higher than 70 %
vitekreport['include']<-ifelse(vitekreport$correct_genus_identified, TRUE, 
                           ifelse(as.numeric(as.character(vitekreport$proba_rank_1)) < 70, TRUE, #include bad spectra, only inlcude if genus from with more than 1.7 confidence
                                  ifelse(vitekreport$genus_identified == 'Escherichia' & vitekreport$genus_NGS == 'Shigella', TRUE, 
                                         ifelse(vitekreport$genus_identified == 'Raoultella' & vitekreport$genus_NGS == 'Klebsiella', TRUE, 
                                                ifelse(vitekreport$species_NGS =='aerogenes' & vitekreport$genus_identified == 'Enterobacter', TRUE, FALSE)))))

# add empty 'max_SN' and 'max_resolution' if missing
if (any(!grepl('max_SN', colnames(vitekreport)))){vitekreport['max_SN']<-NA}
if (any(!grepl('max_resolution', colnames(vitekreport)))){vitekreport['max_resolution']<-NA}

# only keep columns needed
vitekreport<-vitekreport[,c("strainnumber","files","position", "nb_peaks", "max_SN","max_intensities","max_resolution","identifType", "species_rank_1",          
                            "proba_rank_1","score_rank_1","lab","Strain","n_species", "diff","diff_proba", 
                            "species_identified" ,"genus_identified","species_NGS","genus_NGS","correct_genus_identified","correct_species_identified","include")]


vitekreport <- apply(vitekreport,2,as.character)
write.csv2(vitekreport, args[3])  



