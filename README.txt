This is the GitHub repository for a project on rhythm in a 49-language sample from the DoReCo corpus.

Involved: Susanne Fuchs, Ludger Paschen, Lara S. Burchardt
Affiliations: Leibniz Center General Linguistics / Leibniz-Zentrum Allgemeine Sprachwissenschaft, Humboldt University of Berlin / Humboldt-Universität zu Berlin, IMéRa and ILCB Aix-Marseille University


Includes  six R code files:
00a_metadata: This script combines and reformats metadata retrieved from DoReCo 2.0. This allows the metadata to be easily merged with the core DoReCo data in 01_prep_doreco_data.R

00b_synthesis: This script calculates synthesis scores for all core DoReCo datasets with morphological annotations. Synthesis scores are defined as the mean number of morphs (roots, affixes, clitics) per unique word form. Disfluencies, pauses, and units adjacent to word-internal pauses are excluded from the calculations. 

01_prep_doreco_data: This script reads the word-level core datasets from DoReCo 2.0 and reformats them such that one row corresponds to one inter-onset-interval (IOI). It also merges the metadata extracted in 0a_metadata.R and the synthesis scores created in 0b_synthesis.R. Lastly, the script merges the manually coded info on tone languages in DoReCo_2_0_tone.csv and the geographical data (sourced from Glottolog) in DoReCo_2_0_geo.csv.

02_prep_ipu_format_split: This script prepares the data for the rhythm analysis. Data input includes ~100.000 iois from 49 languages and multiple speakers per language. We separate the data into uninterrupted instances of speech from one speaker, per language. We extract onsets and offsets of events. 

03_post_metadata: This script handles the post-processing of data. Rhythm analysis results from RANTO are joined back together with metadata. Additional metadata like geographical coordinates and average height per country are added. 

04_plots_analysis: This script produces and calculates all plots and statistical analyses for rhythm analysis results as reported in the manuscript. 


Includes different data files: 
Raw Input Data: DoReCo_2_0_IO_20240928.csv
Complete Rhythm Results Data: rhythm_results_doreco_ioi_meta_complete.rds

Intermediate data files used in the code 04: 
unsplit_ioi_data_for_rhythm_analysis_including_meta_data_run_Oct24.csv
iois_doreco_rerun_ir.csv
doreco_languages_metadata_2_0.csv

Additional data intermediate datafile to big for upload to GitHub: 
doreco_interger_ratios_all_pairs_all_files.csv (contact corresponding authors for data, or run analysis yourself to generate data)


