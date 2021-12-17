## This file reads in a html file and exports a .csv file summarising the the bruker species identification (one line per spectrum). It requires the folloing two arguments 
## This file reads in a html file and exports a .csv file summarising the the bruker species identification (one line per spectrum). It requires the folloing three arguments 
#(i) input path to the html file
#(ii) input to the .csv file containing the brukercode translation 
#(iii) Input file which shich strainnumber is which species per NGS.
#(iii) output file to the csv file

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

rawHTML <- paste(readLines(args[1]), collapse = '\n')
rawHTML <- unlist(strsplit(rawHTML,'\\\n\\s*<'))

# add brukercode to strainnumber translation
brukerreport<-read.csv(args[2], sep=',')
# define bruker_pos variable which can be used as an anchor in the raw html
brukerreport['bruker_pos']<-if(length(unique(brukerreport$brukercode)) < length(brukerreport$brukercode)){paste0(brukerreport$brukercode, '.*',brukerreport$position_ascii)} else {brukerreport$brukercode}
brukerreport['bruker_pos']<-gsub('[[:space:]]+', '.*', brukerreport$bruker_pos)

# define empty lists for the variables to extract
Species_Score<-list()
Species_Score_all<-list()
species.names <- 1:10
col<-list()
position<-list()
run<-list()
# extract species and scores, position and run from the html files
for (i in 1:length(unique(brukerreport$bruker_pos))){
  Species<-list()
  Score <-list()
  temp<-NULL
  selected_lines<-NULL
  selected_lines<-if(any(grepl(paste0(".*MaldiBiotyperRealTimeClassification.*",unique(brukerreport$bruker_pos)[[i]], ".*"), rawHTML))){
    which(grepl(paste0(".*MaldiBiotyperRealTimeClassification.*",unique(brukerreport$bruker_pos[[i]]), ".*"), rawHTML))
  }else{
    NA}
  selected_lines2<-if(any(grepl(paste0(".*Desktop.*",unique(brukerreport$bruker_pos[[i]]), ".*"), rawHTML))){
    which(grepl(paste0(".*Desktop.*",unique(brukerreport$bruker_pos)[[i]], ".*"), rawHTML))
  }else{
    NA}
  selected_lines<-union(selected_lines, selected_lines2)
  selected_lines<-selected_lines[!is.na(selected_lines)]
  temp <- rawHTML[selected_lines]
  run<-append(run, gsub('(.*MaldiBiotyperRealTimeClassification.+)([[:alnum:]]{8}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{12})(.+)([[:alnum:]]{8}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{12})(.+)(0\\_[[:alnum:]]{2,3})(.*)', '\\2',temp))
  col <- append(col, gsub('(.*MaldiBiotyperRealTimeClassification.+)([[:alnum:]]{8}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{12})(.+)([[:alnum:]]{8}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{12})(.+)(0\\_[[:alnum:]]{2,3})(.*)', '\\4',temp))
  position <- append(position, gsub('(.*MaldiBiotyperRealTimeClassification.+)([[:alnum:]]{8}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{12})(.+)([[:alnum:]]{8}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{4}.{1}[[:alnum:]]{12})(.+)(0\\_[[:alnum:]]{2,3})(.*)', '\\6',temp))
  for (j in 1:10){
    temp2<-NULL
    temp2 <- rawHTML[selected_lines+23+((j-1)*5)]
    temp3<-NULL
    temp3 <- rawHTML[selected_lines+24+((j-1)*5)]
    Species_Score[[j]]<-list()
    Species_Score[[j]]<-append(Species_Score[[j]], paste(gsub('(.*>)(.*)(<.*|\\+|\\s$)', '\\2',temp2), gsub('(.*>)(.*)(<.*)', '\\2',temp3), sep=';'))
    
  }
  names(Species_Score) <- species.names
  Species_Score_all[[i]]<-list()
  Species_Score_all[[i]]<-append(Species_Score_all[[i]], Species_Score)
}


##merge species and scores one dataframe
Species_Score_all_df<-as.data.frame(rbindlist(Species_Score_all))
Species_Score_all_df<-Species_Score_all_df[!is.na(Species_Score_all_df$`1`),]
# convert to character
tt4<-data.frame(run = as.character(run), position = as.character(position), col = as.character(col))
tt4<-tt4[!is.na(tt4$run),]
tt4_sp_scores<-cbind(tt4, Species_Score_all_df)
# convert to long format
tt4_sp_scores_long<- tt4_sp_scores %>% gather("Rank","sp_score", -c('run', 'col', 'position'))
#split into two columns
tt4_sp_scores_long[,c('Species', 'Score')]<-str_split_fixed(tt4_sp_scores_long$sp_score, ';', 2)
tt4_sp_scores_long$sp_score<-NULL
tt4_sp_scores_long['Score']<-as.numeric(as.character(gsub('&lt; \n      0', '0', tt4_sp_scores_long$Score)))

# some have been measured double, once empty no peaks found, and then again with peaks, choose the one with peaks
tt4_sp_scores_long <- tt4_sp_scores_long %>% 
  group_by(run,position,col,Rank) %>%   
  slice(which.max(Score)) %>%
  ungroup()

tt4_all<-tt4_sp_scores_long

# set the score for 'no peaks found to 0
tt4_all['Score']<-as.numeric(ifelse(tt4_all$Species=='no peaks found', 0, as.character(tt4_all$Score)))
tt4_all<-tt4_all[!duplicated(tt4_all),]
tt4_all['col']<-ifelse(grepl('.*Desktop.*', tt4_all$col), gsub('(.*Desktop.)([^ ]+)([^\\\\]*)(\\\\)(0\\_[[:alnum:]]{2,3})(.*)','\\2', tt4_all$col), as.character(tt4_all$col))
tt4_all['run']<-ifelse(grepl('.*Desktop.*', tt4_all$run), gsub('(.*Desktop.)([^ ]+)([^\\\\]*)(\\\\)(0\\_[[:alnum:]]{2,3})(.*)','\\3', tt4_all$run), as.character(tt4_all$run))
tt4_all['position']<-ifelse(grepl('.*Desktop.*', tt4_all$position), gsub('(.*Desktop.)([^ ]+)([^\\\\]*)(\\\\)(0\\_[[:alnum:]]{2,3})(.*)','\\5', tt4_all$position), as.character(tt4_all$position))

# add strainnumber from bukerreport
tt4_all<-merge(tt4_all, brukerreport, by.x='col', by.y = 'brukercode', all.x=T)

#add 'spectra variable
tt4_all['spectra']<-paste(tt4_all$samplename, tt4_all$position, sep = '.')
tt5<-tt4_all[!duplicated(tt4_all),]

#add strainnumber
tt5['samplename']<-ifelse(grepl('([[:digit:]]{1,2})([[:punct:]])([[:digit:]]{1,2}$)', tt5$samplename), gsub('([[:digit:]]{1,2})([[:punct:]])([[:digit:]]{1,2}$)', '\\1', tt5$samplename), as.character(tt5$samplename))
tt5['samplename']<-ifelse(is.na(tt5$samplename), as.character(tt5$col), as.character(tt5$samplename))
tt5['strainnumber']<- str_extract(tt5$samplename, "[[:digit:]]{2}")
tt5['strainnumber']<-ifelse(grepl('^\\d{1}$',tt5$strainnumber), paste0('0',tt5$strainnumber),tt5$strainnumber)

# add species identified and genus identified column 
tt5['species_identified']<-sapply(strsplit(as.character(tt5$Species)," "), `[`, 2)
tt5['genus_identified']<-sapply(strsplit(as.character(tt5$Species)," "), `[`, 1)

#count how many different have been identified with a score higher than 2, and what the difference between highest and second highest species id, if any 
tt5_sum<-tt5 %>% 
  group_by(run, spectra) %>%
  mutate(n_species = n_distinct(species_identified)) %>%
  mutate(n_genera = n_distinct(genus_identified)) %>%
  mutate(n_species_over_2 = n_distinct(species_identified[Score > 2])) %>%
  mutate(n_genera_over_2 = n_distinct(genus_identified[Score > 2])) %>%
  group_by(run, spectra, species_identified) %>%
  mutate(diff = as.numeric(sort(as.numeric(Score), decreasing=T)[1] - sort(as.numeric(Score), decreasing=T)[2])) %>%
  filter(Rank == '1')

#add specied by NGS
numbering<-read.csv2(args[3], sep=',')
numbering['species_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\3', numbering$Strain)
numbering['genus_NGS']<-gsub('(.*)(\\_)(.*)(\\_)(.*)', '\\1', numbering$Strain)
numbering['genus_NGS']<-gsub('Winkia', 'Actinomyces', numbering$genus_NGS)
numbering['Numbering_Shipment_1']<-ifelse(nchar(numbering$Numbering_Shipment_1)==1, paste0('0', numbering$Numbering_Shipment_1), numbering$Numbering_Shipment_1)

#merge to df
tt5_sum<-merge(tt5_sum, numbering, by.x = 'strainnumber', by.y = 'Numbering_Shipment_1', all.x = T)

# check if correct genus amongst identified genera
tt5_sum['correct_genus_identified']<-mapply(grepl, tt5_sum$genus_NGS, tt5_sum$genus_identified)
tt5_sum['correct_species_identified']<-mapply(grepl, tt5_sum$species_NGS, tt5_sum$species_identified)

# remove spectra as contaminations if the from genus has been identified with high confidence. exceptions are Shigella/Escherichia and Klebsiella/Raoultella
tt5_sum['include']<-ifelse(tt5_sum$correct_genus_identified, TRUE, 
                           ifelse(tt5_sum$Score<1.7, TRUE, #include bad spectra, only inlcude if genus from with more than 1.7 confidence
                                  ifelse(tt5_sum$genus_identified == 'Escherichia' & tt5_sum$genus_NGS == 'Shigella', TRUE, 
                                         ifelse(tt5_sum$genus_identified == 'Raoultella' & tt5_sum$genus_NGS == 'Klebsiella', TRUE, FALSE))))
# save all as character
tt5_sum <- if(nrow(tt5_sum) > 1){apply(tt5_sum,2,as.character)} else {tt5_sum}

# export csv
write.table(tt5_sum,args[4], row.names = F, dec='.', sep=';')
