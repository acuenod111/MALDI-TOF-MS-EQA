This repository includes scripts which were used for the analysis of MALDI-TOF mass spectra, which were acquired by 36 diagnostic laboratories in the framework of a MALDI-TOF MS External Quality Assessment (EQA). 
The participating laboratories acquired mass spectra at two different timepoints: (i) for baseline quality assessment and (ii) after an intervention which consisted of detailed feedback report and instructions on how to acquire the second set of measurements. 

Scripts for each timepoint can be found in the respective directory. 
-	The scripts ‘brukercode_readout_html.R’ and  ‘sampleID_readout_html.R’ (depending on whether the mass spectra were assigned a ‘brukercode’ or were assigned their given sample name) can be used to summarise species identification retrieved by comparison of the spectra to the microflex Biotyper database (MALDI Biotyper Compass Library, Revision E (Vers. 8.0, 8468 MSP, RUO, Bruker Daltonics, Bremen, Germany).  
-	The scripts ‘readout_VitekMS_ID.R’ summarises the species identification retrieved by comparison of the spectra to the the VItekMS database (bioMérieux, Marcy-l’Étoile, France) (v3.2)
-	The scripts ‘summarise_marker_based_args.R’ summarised the species identification retrieved by comparison of the spectra to a ribosomal marker based database PAPMIDTM or PAPMIDTM subtyping modules (both Mabritec AG, Riehen Switzerland). 
-	The scripts ‘spectra_comparison_ascii_args.R’ extracts mass spectral features such as total number of peaks, total intensity, number of ribosomal subunits detected, from the peak lists (ascii format), which have previously been exported from the raw spectra using the microflex Biotyper /VitekMS / Axima Confidence softwares. 
-	The scripts ‘summarize_eval_spectra.R’ summarises all data outputted from the above listed scripts, for each timepoint separately. 
-	The script ‘sum_all_plot.R’ then merges the outputs from the two timepoints visualises the results. 
-	The script ‘questionnaire_plot.R’ combines the answers of a questionnaire the participating laboratories have filled out with the mass spectral data acquired for the baseline quality assessment and visualises observed correlations between routine laboratory practices and MALDI-TOF mass spectral quality. 
Please refer to the instructions in the scripts for more details.

