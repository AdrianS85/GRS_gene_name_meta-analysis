# required packages:: rentrez, biomaRt, tidyverse (dplyr, purrr, readr), rlist, stringr, data.table, getPass
# rm(list = ls(pattern = '^temp.*'))
# devtools::install_github("hadley/lineprof")
rm(strvec_ncbi_query_for_identifers)


######### PREPARING METADATA ######### 

source('functions_for_genename_conversions.R')

### Read in the table with description of experiments. They need to have the same identifiers as data in PRE_DATA, that is: Paper(int), Experiment(chr). Also needs Species(chr) column, if probe_id identification function or ncbi query tool is to be used
descriptions <- readr::read_tsv("descriptions.txt", col_types = "nccccccccccccccccccccccc")

######### PREPARING METADATA ######### 



##### PROBE ID - GET THE HIGHEST-HIT-RETURNING ID TYPE, THEN ANNOTATE ##### 
# Single Probe_ID may be assigned to multiple genes
DATA_FROM_HIGHEST_HIT_ANALYSIS <- get_the_highest_hit_returning_id_type(
  str_LIST_DATA = 'data_v4_Probe_ID.tsv', 
  descriptions_ = descriptions,
  int_Probe_IDs_to_test = 200,
  str_experiment_name = 'Probe_ID')

identifiers_used_for_annotation_global <- as.character(DATA_FROM_HIGHEST_HIT_ANALYSIS[[2]]$platform_to_use)
# OR
identifiers_used_for_annotation_global <- readr::read_tsv('Probe_ID/platform_to_use_for_probes_based_analysis.tsv')$platform_to_use

annotations_Probe_ID <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'data_v4_Probe_ID.tsv', str_identifier_type__ = identifiers_used_for_annotation_global, backup_experiment_name = 'Probe_ID') # 195495 vs 195495

# annotations_Probe_ID_12 <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'data_v4_Probe_ID_12.tsv', str_identifier_type__ = 'affy_mg_u74av2', backup_experiment_name = 'Probe_ID_12') # 195495 vs 195495

##### PROBE ID - GET THE HIGHEST-HIT-RETURNING ID TYPE, THEN ANNOTATE ##### 



####### GEMMA ANNOTATION (31/33) ####### 

platforms_ids_to_download <- readr::read_tsv('data_v4_gemma_platforms.tsv')

gemma_platforms <- download_platforms_from_gemma(platforms_ids_to_download$Platform_ID)

annotations_gemma_probeID <- get_gemma_annotations_for_data(str_data_file_name = 'data_v4_gemma.tsv', list_gemma_platforms = gemma_platforms) # 2002 vs 2002

dir.create('gemma_annotation')
readr::write_tsv(x = annotations_gemma_probeID, path = 'gemma_annotation/raw_gemma_annotations.tsv')

rm(platforms_ids_to_download, gemma_platforms)
####### GEMMA ANNOTATION (31/33) ####### 



####### KNOWN IDENTIFIERS ####### 
# These identifiers return unique gene name per identifier, because every id is assigned to known gene. This is different for Probe_IDs
# 'ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'refseq_mrna', 'refseq_mrna_predicted', 'embl'

#All this needs to be repeated with data redone
annotations_ensembl_gene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'data_v4_ensembl_gene_id.tsv', str_identifier_type__ = 'ensembl_gene_id') #20007 vs 20007

annotations_ensembl_transcript_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'data_v4_ensembl_transcript_id.tsv', str_identifier_type__ = 'ensembl_transcript_id') #2422 vs 2422

annotations_entrezgene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'data_v4_entrezgene_id.tsv', str_identifier_type__ = 'entrezgene_id') #29345 vs 29345 # REDO

annotations_refseq_mrna <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'data_v4_refseq_mrna.tsv', str_identifier_type__ = 'refseq_mrna') #4710 vs 4710 !!! Use this to check for repeating values in external_gene_names

annotations_refseq_mrna_predicted <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'data_v4_refseq_mrna_predicted.tsv', str_identifier_type__ = 'refseq_mrna_predicted') #464 vs 464 Doesnt work with aggregate


####### KNOWN IDENTIFIERS ####### 


### !!! add to function make_and_write_table_with_original_and_ncbi_ids chaning ncbi_response to external_gene_name
####### NCBI IDENTIFIERS ####### 
# Annotate gene names and other ids with GeneID (entrezgene_id?)
annotated_other_ids_to_geneIDs <- annotate_identifiers_to_geneID(str_filename_ = 'data_v4_other_ids.tsv', str_experiment_name = 'other_ids', descriptions_ = descriptions) #2408 vs 
# OR
# annotated_other_ids_to_geneIDs <- readr::read_tsv(file = 'other_id/data_v4_for_ncbi_output_of_ncbi_query.tsv')

annotated_names_to_geneIDs <- annotate_identifiers_to_geneID(str_filename_ = 'data_v4_names.tsv', str_experiment_name = 'names', descriptions_ = descriptions, bool_reformat_names = F, chr_gene_identifier_ = 'Gene Name') #637 vs 637

temp_annotated_geneIDs <- rbind(annotated_other_ids_to_geneIDs, annotated_names_to_geneIDs)

# Prepare input for entrezgene_id annotation
annotated_names_and_other_ids_to_geneIDs <- temp_annotated_geneIDs %>%
  dplyr::filter(.data$ncbi_response != 'none' & !is.na(.data$ncbi_response)) %>%
  dplyr::mutate(Probe_ID = ncbi_response) %>%
  dplyr::select(-ncbi_response)

# Get non-anotated values
not_annotated_data_1 <- temp_annotated_geneIDs %>%
  dplyr::filter(.data$ncbi_response == 'none' | is.na(.data$ncbi_response))

check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'ncbi_annotation_for_geneNames_and_otherIDs', df_original = temp_annotated_geneIDs, list_df_splited = list(annotated_names_and_other_ids_to_geneIDs, not_annotated_data_1))

annotations_names_and_other_ids_to_geneIDs <- master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(df_to_annotate = annotated_names_and_other_ids_to_geneIDs, str_identifier_type___ = 'entrezgene_id', str_experiment_name__ = 'ncbi_annotations_annotated_as_entrezID')

rm(temp_annotated_geneIDs, annotated_names_and_other_ids_to_geneIDs)
# annotated_short <- annotate_identifiers_to_geneID(str_filename_ = 'data_v4_names_short.tsv', str_experiment_name = 'short', descriptions_ = descriptions)

####### NCBI IDENTIFIERS ####### 

rm(data_gemma)

####### DEAL WITH NON-ANNOTATED DATA IN FINAL TABLE ####### 

# logFCs for the same gene within single Experiment are avaraged using median. No other filtering is used

# It is critical to use this workflow only after producing actual, tested proper merged dataset - we need to add all the data and check if the lenght of this set is equal to lenght of input dataset. 
final_merged_dataset <- gather_all_datasets_into_single_df(regex_pattern_to_find_datasets_with = '^annotations.*')

noNAs_final_merged_dataset <- final_merged_dataset %>%
  subset(subset = !(is.na(final_merged_dataset$external_gene_name) | final_merged_dataset$external_gene_name == ''))

not_annotated_data_2 <- final_merged_dataset %>%
  subset(subset = (is.na(final_merged_dataset$external_gene_name) | final_merged_dataset$external_gene_name == ''))

check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'final_merged_dataset', df_original = final_merged_dataset, list_df_splited = list(noNAs_final_merged_dataset, not_annotated_data_2))

not_annotated_data <- rbind(not_annotated_data_1, not_annotated_data_2)



####### DEAL WITH NON-ANNOTATED DATA IN FINAL TABLE ####### 



####### PREPARING FINAL TABLE ####### 
# Here we collapse multiple genenames to single gene name
noNAs_final_merged_dataset$final_gene_name <- as.character(lapply(
  X = noNAs_final_merged_dataset$external_gene_name,
  FUN = function(x) {
    return(select_best_geneName_wrapper_for_single_string(x)) ### !!! Gm's cna have 5 digits... damn
  }
))

readr::write_tsv(x = final_merged_dataset, path = 'final_merged_dataset.tsv')

# Uniformarize gene names by making them all all-lower case
noNAs_final_merged_dataset$lower_final_gene_name <- tolower(noNAs_final_merged_dataset$final_gene_name)

# Medianing values for the same gene within Experiment
medianed_noNAs_final_merged_dataset <- noNAs_final_merged_dataset %>%
  dplyr::group_by(Experiment, lower_final_gene_name) %>%
  dplyr::summarize(logFC_median = median(logFC))

# This is a very sparse matrix. Perhaps try working with it as such
spread_lowercase_final_merged_dataset <- tidyr::spread(data = medianed_final_merged_dataset, key = Experiment, value = logFC_median)

readr::write_tsv(x = spread_lowercase_final_merged_dataset, path = 'final_merged_dataset_for_clustering.tsv')

####### PREPARING FINAL TABLE ####### 




#######################################################









###### POST-PREPARATION OF ANNOTATED PROBE_ID TABLES ###### 




### This will be our input data for further analysis. It contains only significant columns of "Paper", "GroupID", "ensembl_gene_name", "logFC": 
SHORT_SINGLE_TEST_ANNOTATION <- SINGLE_TEST_ANNOTATION %>%
  select("Paper", "GroupID", "ensembl_gene_name", "logFC") %>%
  group_by(Paper, GroupID, ensembl_gene_name) %>%
  nest()

# I dont know why purrr::map didnt work for this. Either way, we just lapply this. Here we write UP or DOWN for every logFC for given genename in given experiment (not paper)
SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$data, FUN = function(x) { mutate(x, Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) } ) 

# Here we count how many UPs and DOWNs are in every given genename in given experiment (not paper)
SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, FUN = function (x) {table(select(x, "Symbol_direction"))} ) 

# Here we establish actual status of gene in given experiment
SHORT_SINGLE_TEST_ANNOTATION$sum_directionality <- as.character(lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, 
                                                                       FUN = function(x) {
                                                                         if (length(x) == 1 && grepl(pattern = "UP", names(x))) { "UP" }
                                                                         else if (length(x) == 1 && grepl(pattern = "DOWN", names(x))) { "DOWN" }
                                                                         else if (length(x) == 2 && grepl(pattern = "DOWN", names(x)) && grepl(pattern = "DOWN", names(x))) { "MIXED" }
                                                                         else { "ERROR" }
                                                                       }))

# Here we remove MIXED expression and multiple genenames
FILT_SHORT_SINGLE_TEST_ANNOTATION <- SHORT_SINGLE_TEST_ANNOTATION %>%
  #filter(!grepl(pattern = "(.*);(.*)", SHORT_SINGLE_TEST_ANNOTATION$ensembl_gene_name)) %>%
  filter(!sum_directionality == "MIXED") %>%
  filter(!sum_directionality == "ERROR") %>%
  filter(!is.na(ensembl_gene_name))

# Here we add mean to each !!! Perhaps this should be final input, because it has: 1) removed genes giving UP+DOWN in the same experiment, 2) removed NAs
STDINPUT_FILT_SHORT_SIN_T_ANNO <- FILT_SHORT_SINGLE_TEST_ANNOTATION %>%
  mutate(mean = as.numeric(map(data, function(x) { as.numeric(mean(x[[1]]))} ))) %>% ## This way we give actuall vector to function, not a data table(tibble)
  select("Paper", "GroupID", "ensembl_gene_name", "sum_directionality", "mean") 
### !!! THE IDEA IS THAT THE "STDINPUT_FILT_SHORT_SIN_T_ANNO" IS THE INPUT FOR ALL FURTHER INQUIRIES !!! ###

# tO raczej nie wyjdzie - nie da siÄ™ wyciÄ…gnÄ…Ä‡ wĹ‚aĹ›ciwej wartoĹ›ci ekspresji z dwĂłch rĂłznych eksperymentĂłw w tym samym paper. MoĹĽe lepiej po prostu dowiedzieÄ‡ siÄ™, ktĂłre geny sÄ… unikatowe dla danego paper i usunÄ…Ä‡ te geny z normalnie uĹĽywanego wczeĹ›niej test annotation pliku. 

###### PAPER-CENTRIC ANALYSIS ###### 


#################################################
##### ANNOTATE PROBE_IDs WITH ENSEMBL NAMES #####
#################################################




# ###### POST-PREPARATION OF ANNOTATED PROBE_ID TABLES ###### 
# 
# ## HERE WE WILL COLLAPSE duplicated probe-genename pairs WITHIN ANNOTATED GENE LIST (and also change column names to standarized ones)
# UNIQ_ANNOT_LIST_DATA <- list()
# for (n in WHICH_EXP_TO_ANAL){
#   UNIQ_ANNOT_LIST_DATA[[n]] <- unique(ANNOT_LIST_DATA[[n]])
#   colnames(UNIQ_ANNOT_LIST_DATA[[n]]) <- c("Probe_ID", "ensembl_gene_name")
# }
# 
# ## Here we will collapse all gene names for given probe
# for (n in WHICH_EXP_TO_ANAL){
#   UNIQ_ANNOT_LIST_DATA[[n]] <- aggregate(
#     ensembl_gene_name~Probe_ID, 
#     data = UNIQ_ANNOT_LIST_DATA[[n]], 
#     FUN = stringr::str_c) ## Here we aggregate the genenames into single row
#   UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- lapply(
#     X = UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name, 
#     FUN = paste, 
#     collapse = "; ") ##
#   UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- as.character(UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name) ## Here we make sure that collapsed genename column is not a list (need for writing function)
# }
# 
# # Je?eli w ensemble nie ma przypisanego genu (oznaczenie NA w przys?anym pliku "100_genow") to program pokazuje we wszystkich rubrykach oznaczenie NA chocia? te kom?rki zawiera?y wcze?niej dane (plik AAAreview11). W pe?nych danych kt?re przys?a?e? mi przed ?wi?tami (plik 1 znajduj?cy si? w skopresowanych danych "TEST_ANNOTATION") tych sond w og?le nie ma. Wygl?da na to, ?e to co przys?a?e? wczoraj odnosi si? chyba do wcze?niejszego pliku w kt?rym dane dla sond nie maj?cych przypisanego genu w ensemble zosta?y ca?kowicie usuni?te a to nie jest dobre :( - AMS W LIST_DATA S? TE DANE, ALE W TEST_ANNOTATION JUZ NIE
# 
# ## Here we merge probes annotated with ensembl with probes data that we have
# TEST_ANNOTATION <- list()
# for (n in WHICH_EXP_TO_ANAL){
#   TEST_ANNOTATION[[n]] <- merge(LIST_DATA[[n]], UNIQ_ANNOT_LIST_DATA[[n]], by = "Probe_ID", all.x = T)
# }
# 
# 
# 
# ## Here lets also remove empty lists
# ### !!! MANUAL SUPPLEMENTARY STEP - TO BE REMOVED IN FINAL VERSION, AS WE WILL NOT HAVE EMPTY LISTS !!! ###
# # TEST_ANNOTATION[[13]] <-NULL
# 
# # Here we bind all list elements into single data table
# SINGLE_TEST_ANNOTATION <- data.table::data.table(rlist::list.rbind(TEST_ANNOTATION), stringsAsFactors = F)
# readr::write_tsv(x = SINGLE_TEST_ANNOTATION, path = "SINGLE_ANNOTATION_probes.tsv")
# 
# 
# ### This will be our input data for further analysis. It contains only significant columns of "Paper", "GroupID", "ensembl_gene_name", "logFC": 
# SHORT_SINGLE_TEST_ANNOTATION <- SINGLE_TEST_ANNOTATION %>%
#   select("Paper", "GroupID", "ensembl_gene_name", "logFC") %>%
#   group_by(Paper, GroupID, ensembl_gene_name) %>%
#   nest()
# 
# # I dont know why purrr::map didnt work for this. Either way, we just lapply this. Here we write UP or DOWN for every logFC for given genename in given experiment (not paper)
# SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$data, FUN = function(x) { mutate(x, Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) } ) 
# 
# # Here we count how many UPs and DOWNs are in every given genename in given experiment (not paper)
# SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, FUN = function (x) {table(select(x, "Symbol_direction"))} ) 
# 
# # Here we establish actual status of gene in given experiment
# SHORT_SINGLE_TEST_ANNOTATION$sum_directionality <- as.character(lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, 
#                                                                        FUN = function(x) {
#                                                                          if (length(x) == 1 && grepl(pattern = "UP", names(x))) { "UP" }
#                                                                          else if (length(x) == 1 && grepl(pattern = "DOWN", names(x))) { "DOWN" }
#                                                                          else if (length(x) == 2 && grepl(pattern = "DOWN", names(x)) && grepl(pattern = "DOWN", names(x))) { "MIXED" }
#                                                                          else { "ERROR" }
#                                                                        }))
# 
# # Here we remove MIXED expression and multiple genenames
# FILT_SHORT_SINGLE_TEST_ANNOTATION <- SHORT_SINGLE_TEST_ANNOTATION %>%
#   #filter(!grepl(pattern = "(.*);(.*)", SHORT_SINGLE_TEST_ANNOTATION$ensembl_gene_name)) %>%
#   filter(!sum_directionality == "MIXED") %>%
#   filter(!sum_directionality == "ERROR") %>%
#   filter(!is.na(ensembl_gene_name))
# 
# # Here we add mean to each !!! Perhaps this should be final input, because it has: 1) removed genes giving UP+DOWN in the same experiment, 2) removed NAs
# STDINPUT_FILT_SHORT_SIN_T_ANNO <- FILT_SHORT_SINGLE_TEST_ANNOTATION %>%
#   mutate(mean = as.numeric(map(data, function(x) { as.numeric(mean(x[[1]]))} ))) %>% ## This way we give actuall vector to function, not a data table(tibble)
#   select("Paper", "GroupID", "ensembl_gene_name", "sum_directionality", "mean") 
# ### !!! THE IDEA IS THAT THE "STDINPUT_FILT_SHORT_SIN_T_ANNO" IS THE INPUT FOR ALL FURTHER INQUIRIES !!! ###
# 
# # tO raczej nie wyjdzie - nie da siÄ™ wyciÄ…gnÄ…Ä‡ wĹ‚aĹ›ciwej wartoĹ›ci ekspresji z dwĂłch rĂłznych eksperymentĂłw w tym samym paper. MoĹĽe lepiej po prostu dowiedzieÄ‡ siÄ™, ktĂłre geny sÄ… unikatowe dla danego paper i usunÄ…Ä‡ te geny z normalnie uĹĽywanego wczeĹ›niej test annotation pliku. 
# 
# ###### PAPER-CENTRIC ANALYSIS ###### 




# ###### POST-PREPARATION OF ANNOTATED PROBE_ID TABLES ###### 
# 
# ## HERE WE WILL COLLAPSE duplicated probe-genename pairs WITHIN ANNOTATED GENE LIST (and also change column names to standarized ones)
# UNIQ_ANNOT_LIST_DATA <- list()
# for (n in WHICH_EXP_TO_ANAL){
#   UNIQ_ANNOT_LIST_DATA[[n]] <- unique(ANNOT_LIST_DATA[[n]])
#   colnames(UNIQ_ANNOT_LIST_DATA[[n]]) <- c("Probe_ID", "ensembl_gene_name")
# }
# 
# ## Here we will collapse all gene names for given probe
# for (n in WHICH_EXP_TO_ANAL){
#   UNIQ_ANNOT_LIST_DATA[[n]] <- aggregate(
#     ensembl_gene_name~Probe_ID, 
#     data = UNIQ_ANNOT_LIST_DATA[[n]], 
#     FUN = stringr::str_c) ## Here we aggregate the genenames into single row
#   UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- lapply(
#     X = UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name, 
#     FUN = paste, 
#     collapse = "; ") ##
#   UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- as.character(UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name) ## Here we make sure that collapsed genename column is not a list (need for writing function)
# }
# 
# # Je?eli w ensemble nie ma przypisanego genu (oznaczenie NA w przys?anym pliku "100_genow") to program pokazuje we wszystkich rubrykach oznaczenie NA chocia? te kom?rki zawiera?y wcze?niej dane (plik AAAreview11). W pe?nych danych kt?re przys?a?e? mi przed ?wi?tami (plik 1 znajduj?cy si? w skopresowanych danych "TEST_ANNOTATION") tych sond w og?le nie ma. Wygl?da na to, ?e to co przys?a?e? wczoraj odnosi si? chyba do wcze?niejszego pliku w kt?rym dane dla sond nie maj?cych przypisanego genu w ensemble zosta?y ca?kowicie usuni?te a to nie jest dobre :( - AMS W LIST_DATA S? TE DANE, ALE W TEST_ANNOTATION JUZ NIE
# 
# ## Here we merge probes annotated with ensembl with probes data that we have
# TEST_ANNOTATION <- list()
# for (n in WHICH_EXP_TO_ANAL){
#   TEST_ANNOTATION[[n]] <- merge(LIST_DATA[[n]], UNIQ_ANNOT_LIST_DATA[[n]], by = "Probe_ID", all.x = T)
# }
# 
# 
# 
# ## Here lets also remove empty lists
# ### !!! MANUAL SUPPLEMENTARY STEP - TO BE REMOVED IN FINAL VERSION, AS WE WILL NOT HAVE EMPTY LISTS !!! ###
# # TEST_ANNOTATION[[13]] <-NULL
# 
# # Here we bind all list elements into single data table
# SINGLE_TEST_ANNOTATION <- data.table::data.table(rlist::list.rbind(TEST_ANNOTATION), stringsAsFactors = F)
# readr::write_tsv(x = SINGLE_TEST_ANNOTATION, path = "SINGLE_ANNOTATION_probes.tsv")
# 
# 
# ### This will be our input data for further analysis. It contains only significant columns of "Paper", "GroupID", "ensembl_gene_name", "logFC": 
# SHORT_SINGLE_TEST_ANNOTATION <- SINGLE_TEST_ANNOTATION %>%
#   select("Paper", "GroupID", "ensembl_gene_name", "logFC") %>%
#   group_by(Paper, GroupID, ensembl_gene_name) %>%
#   nest()
# 
# # I dont know why purrr::map didnt work for this. Either way, we just lapply this. Here we write UP or DOWN for every logFC for given genename in given experiment (not paper)
# SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$data, FUN = function(x) { mutate(x, Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) } ) 
# 
# # Here we count how many UPs and DOWNs are in every given genename in given experiment (not paper)
# SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, FUN = function (x) {table(select(x, "Symbol_direction"))} ) 
# 
# # Here we establish actual status of gene in given experiment
# SHORT_SINGLE_TEST_ANNOTATION$sum_directionality <- as.character(lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, 
#                                                                        FUN = function(x) {
#                                                                          if (length(x) == 1 && grepl(pattern = "UP", names(x))) { "UP" }
#                                                                          else if (length(x) == 1 && grepl(pattern = "DOWN", names(x))) { "DOWN" }
#                                                                          else if (length(x) == 2 && grepl(pattern = "DOWN", names(x)) && grepl(pattern = "DOWN", names(x))) { "MIXED" }
#                                                                          else { "ERROR" }
#                                                                        }))
# 
# # Here we remove MIXED expression and multiple genenames
# FILT_SHORT_SINGLE_TEST_ANNOTATION <- SHORT_SINGLE_TEST_ANNOTATION %>%
#   #filter(!grepl(pattern = "(.*);(.*)", SHORT_SINGLE_TEST_ANNOTATION$ensembl_gene_name)) %>%
#   filter(!sum_directionality == "MIXED") %>%
#   filter(!sum_directionality == "ERROR") %>%
#   filter(!is.na(ensembl_gene_name))
# 
# # Here we add mean to each !!! Perhaps this should be final input, because it has: 1) removed genes giving UP+DOWN in the same experiment, 2) removed NAs
# STDINPUT_FILT_SHORT_SIN_T_ANNO <- FILT_SHORT_SINGLE_TEST_ANNOTATION %>%
#   mutate(mean = as.numeric(map(data, function(x) { as.numeric(mean(x[[1]]))} ))) %>% ## This way we give actuall vector to function, not a data table(tibble)
#   select("Paper", "GroupID", "ensembl_gene_name", "sum_directionality", "mean") 
# ### !!! THE IDEA IS THAT THE "STDINPUT_FILT_SHORT_SIN_T_ANNO" IS THE INPUT FOR ALL FURTHER INQUIRIES !!! ###
# 
# # tO raczej nie wyjdzie - nie da siÄ™ wyciÄ…gnÄ…Ä‡ wĹ‚aĹ›ciwej wartoĹ›ci ekspresji z dwĂłch rĂłznych eksperymentĂłw w tym samym paper. MoĹĽe lepiej po prostu dowiedzieÄ‡ siÄ™, ktĂłre geny sÄ… unikatowe dla danego paper i usunÄ…Ä‡ te geny z normalnie uĹĽywanego wczeĹ›niej test annotation pliku. 
# 
# ###### PAPER-CENTRIC ANALYSIS ###### 








































##############################################################
##############################################################
##############################################################







###### WHOLE DATASET ANALYSIS ######

# Tutaj liczymy ile razy geny wystepuja w oryginalnym dataset (w eksperymentach, nie w paperach)
HOW_MANY_TIMES_EXP_STDINPUT <- STDINPUT_FILT_SHORT_SIN_T_ANNO %>%
  select(ensembl_gene_name) %>%
  group_by(ensembl_gene_name) %>%
  summarise(number = n())

# Tutaj liczymy ile razy geny wystepuja w oryginalnym dataset (w paperach)
HOW_MANY_TIMES_PAPER_STDINPUT <- STDINPUT_FILT_SHORT_SIN_T_ANNO %>%
  select(Paper, ensembl_gene_name) %>%
  unique() %>%
  select(ensembl_gene_name) %>%
  group_by(ensembl_gene_name) %>%
  summarise(number = n())

###### Tutaj robimy specjalny input do klastrowania, w ktĂłrym uĹĽywamy tylko genĂłw, ktĂłre zostaĹ‚y wykryte przynajmniej w 3 oddzielnych paperach ###### 
UP3_HOW_MANY_TIMES_PAPER_STDINPUT <- HOW_MANY_TIMES_PAPER_STDINPUT %>%
  filter(number >= 3)

UP3_PAPER_CLUSTERING_INPUT <- merge(STDINPUT_FILT_SHORT_SIN_T_ANNO, UP3_HOW_MANY_TIMES_PAPER_STDINPUT, by = "ensembl_gene_name", all.y = T)

FOR_CLUST_UP3_PAPER_CLUSTERING_INPUT <- UP3_PAPER_CLUSTERING_INPUT %>%
  FOR_CLUS()
#readr::write_tsv(x = FOR_CLUST_UP3_PAPER_CLUSTERING_INPUT, path = "FOR_CLUST_UP3_PAPER_CLUSTERING_INPUT.tsv")  









###### Tutaj robimy specjalny input do klastrowania, w ktĂłrym uĹĽywamy tylko genĂłw, ktĂłre zostaĹ‚y wykryte przynajmniej w 3 oddzielnych paperach ###### 


# Tutaj liczymy ile razy geny wyst?puj? w oryginalnym dataset, patrz?c czy s? up czy down
UorDWHOLE_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(Gene_symbol, logFC) %>%
  mutate(Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) %>%
  mutate(Symbol_direction = paste(Gene_symbol, Symbol_direction, sep = "_")) %>%
  group_by(Symbol_direction) %>%
  summarise(number = n()) %>% 
  mutate(Gene_symbol2 = str_remove(Symbol_direction, "_.*"))

###### WHOLE DATASET ANALYSIS ######



###### COMPARISONS-CENTERED ANALYSIS ######



### Here we set whether we want to analyze papers or comparisons
P_or_C = quo(Paper) #" GroupID OR Paper "



# Tutaj liczymy ile razy geny wyst?puj? W KA?DYM Z POR?WNA?, nie patrz?c czy s? up czy down
COMP_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(!!P_or_C, Gene_symbol) %>%
  group_by(!!P_or_C, Gene_symbol) %>%
  summarise(number = n())



# Tutaj liczymy ile razy geny wyst?puj? W KA?DYM Z POR?WNA?, patrz?c czy s? up czy down    
UorDCOMP_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(!!P_or_C, Gene_symbol, logFC) %>%
  mutate(Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) %>%
  mutate(Symbol_direction = paste(Gene_symbol, Symbol_direction, sep = "_")) %>%
  group_by(!!P_or_C, Symbol_direction) %>%
  summarise(Sym_dir_number = n()) %>%
  mutate(Gene_symbol2 = str_remove(Symbol_direction, "_.*"))



#Divide data into genes expressed in single direction in given comparison, vs genes expressed in different direction (bad genes)
nonUNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UorDCOMP_NO_UNIDS_ORG_DATA %>%
  group_by(!!P_or_C) %>%
  filter(duplicated(Gene_symbol2, fromLast = T) | duplicated(Gene_symbol2))

UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UorDCOMP_NO_UNIDS_ORG_DATA %>%
  group_by(!!P_or_C) %>%
  filter(!duplicated(Gene_symbol2, fromLast = T) & !duplicated(Gene_symbol2))



# Check if unique/duplicated division went well           
if (nrow(UorDCOMP_NO_UNIDS_ORG_DATA) - (nrow(nonUNIQ_UorDCOMP_NO_UNIDS_ORG_DATA) + nrow(UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA)) != 0){ 
  stop("Hey, fwend! You have some wierd values in Your counted data, buddy! Better check whats happening, or Your results will smell of moose scrotum!")
}



# Here we make a table only with genes that were replicated in few comparisons
REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA %>%
  filter(Sym_dir_number >= 3)


#Annotate base on Paper OR GroupID
ANNO_REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- merge(REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA, COMPARISONS, by = "Paper")


###### COMPARISONS-CENTERED ANALYSIS ######






















# name_global <- 'ncbi/'
# 
# PRE_DATA_global <- read_preformated_data(str_filename = 'data_v4_for_ncbi.tsv', int_numbers_are_from = 6, int_numbers_are_to = 8, col_types_ = 'nncccccc')
# 
# LIST_DATA_global <- split(PRE_DATA_global, f = PRE_DATA_global$Paper)
# 
# dir.create(name_global)
# 
# readr::write_tsv(rlist::list.rbind(LIST_DATA_global), paste0(name_global, 'input_for_ncbi.tsv'))
# 
# # List of species names based on data in description files. This is a file we will be working on. Species is in format: small letters, english plural of species
# exp_species_global <- descriptions %>%
#   select("Paper_ID", "Species") %>%
#   unique() %>%
#   filter(Paper_ID %in% unique(PRE_DATA_global$Paper))
# 
# write_lenghts_of_list_objects(LIST_DATA_global, paste0(name_global, 'list_data_lenghts_probes.tsv')) 
# readr::write_tsv(exp_species_global, paste0(name_global, 'exp_species_used_for_testing_which_platform_to_use_probes.tsv'))
# 
# short_PRE_DATA_global <- PRE_DATA_global %>%
#   filter(Paper == 44)
# 
# test_vector <- as.vector(short_PRE_DATA_global$Probe_ID)
# 
# test_ncbi_query_for_identifers <- search_for_ids_in_ncbi(test_vector)
# 
# original_and_ncbi_ids_44 <- make_and_write_table_with_original_and_ncbi_ids(df_returned_by_entrez_gene_search = test_ncbi_query_for_identifers, df_original_data = short_PRE_DATA_global, str_name_of_the_file = 'test_minus_44.tsv', experiment_directory_name = name_global)
# ### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database
# 
# 
# original_and_ncbi_ids <- rbind(original_and_ncbi_ids_minus44, original_and_ncbi_ids_44)
# original_and_ncbi_ids <- original_and_ncbi_ids[order(original_and_ncbi_ids$Nb),] 
# 
# readr::write_tsv(x = original_and_ncbi_ids, path = paste0(name_global, 'data_v4_for_ncbi_output_of_ncbi_query.tsv'))
####### OTHER (FOR_NCBI) IDENTIFIERS ####### 

# Backup
# Make sure that name format for every species is correct
# LIST_DATA <- purrr::map2(.x = LIST_DATA, .y = exp_species$Species, .f = ~ change_name_to_proper_format_given_the_species(.x, .y))

### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database
# leftovers <- readr::read_tsv('data_v4_leftovers_for_ncbi.txt')
# 
# leftovers_vector <- as.vector(leftovers$Probe_ID)
# 
# ncbi_query_for_genes <- search_for_ids_in_ncbi(leftovers_vector)
# 
# original_and_ncbi_ids <- make_and_write_table_with_original_and_ncbi_ids(ncbi_query_for_genes)
### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database


# ##### PROBE ID - DEBUG, DEPRECATED BY USING FUNCTION-BASED WORKFLOW  ##### 
# experiment_name_global <- 'Probe_ID_analysis/'
# 
# dir.create(experiment_name_global)
# 
# # Sadly, "numeric" columns may contain both . and , as decimary separator. Hence we need to load them as char, replace all , with . change into numeric and only then we have nice input. Some probes are significant despite having logFCs near 0. Also: format of input data is: single file, multiple columns, at least Paper(int), Experiment(str), Probe_ID, logFC. For now lets agree on format: Nb Paper	Experiment	Probe_ID	Symbol	adj_p	p	logFC. Note - Probe_ID will actually not always be actual probe id, because we have different ID types, but it would be a bitch to change this throughout the code
# PRE_DATA_global <- read_preformated_data(str_filename = 'data_v4_Probe_ID.tsv', int_numbers_are_from = 6, int_numbers_are_to = 8, col_types_ = 'nncccccc')
# 
# LIST_DATA_global <- split(PRE_DATA_global, f = PRE_DATA_global$Paper)
# 
# # List of species names based on data in description files. This is a file we will be working on. Species is in format: small letters, english plural of species
# exp_species_global <- descriptions %>%
#   select("Paper_ID", "Species") %>%
#   unique() %>%
#   filter(Paper_ID %in% unique(PRE_DATA_global$Paper))
# 
# # For debugging
# write_lenghts_of_list_objects(LIST_DATA_global, paste0(experiment_name_global, 'list_data_lenghts_probes.tsv')) 
# readr::write_tsv(exp_species_global, paste0(experiment_name_global, 'exp_species_used_for_testing_which_platform_to_use_probes.tsv'))
# ##### PROBE ID - DEBUG, DEPRECATED BY USING FUNCTION-BASED WORKFLOW  ##### 



# devtools::install_github('PavlidisLab/gemmaAPI.R')
# 
# ?gemmaAPI::endpointFunctions
# 
# 
# test_data <- read.csv('xxx.txt', fileEncoding="UTF-16LE")
# 
# 
# test <- biomaRt::getBM(
#   attributes = c('affy_rat230_2', "external_gene_name"),
#   filters = 'affy_rat230_2', 
#   values = test_data$Probe_ID, 
#   uniqueRows = T,
#   mart = biomaRt::useMart(
#     "ENSEMBL_MART_ENSEMBL", 
#     dataset = "rnorvegicus_gene_ensembl")
# )






# I dont know how to pass column names via string and $ operator. [''] way doesnt work. So I guess im stuck with hard-coding it for now. # this part also needs to have a capitalized first letter




# 
# 
# devtools::install_github('PavlidisLab/gemmaAPI.R')
# 
# ?gemmaAPI::endpointFunctions
# ?gemmaAPI::highLevelFunctions
# 
# test <- gemmaAPI::allPlatforms()
# 
# gemmaAPI::setGemmaUser(username = 'AdrianS85', password = '13Haslo13')
# xxx <- gemmaAPI::geneInfo('1859', request = 'evidence')
# xxx <- gemmaAPI::geneInfo('ENSRNOT00000007747', request = 'evidence')
# 
# data = 
#   gemmaAPI::datasetInfo('GSE107999',
#                         request='data', # we want this endpoint to return data. see documentation
#                         filter = FALSE, # data request accepts filter argument we want non filtered data
#                         return = TRUE, # TRUE by default, all functions have this. if false there'll be no return
#                         file = NULL # NULL by default, all functions have this. If specificed, output will be saved.
#   )
# 
# rm(data)
# 
# datasetInfo
# allPlatforms
# platformInfo
# geneInfo 'Can either be the NCBI ID (1859), Ensembl ID (ENSG00000157540) or official symbol (DYRK1A) of the gene.'
# allTaxa
# gemmaAPI::geneInfo()
# 
# 
# 
# library(httr)
# library(jsonlite)
# library(xml2)
# 
# server <- "https://rest.ensembl.org"
# ext <- "/xrefs/symbol/rattus_norvegicus/LOC497770?"
# # "/xrefs/symbol/rattus_norvegicus/MAP1B?"
# 
# 
# https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=dbfind&inputValues=A1BG,MYC,1,3&output=geneid&taxonId=9606&format=row
# 
# r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
# 
# stop_for_status(r)
# 
# head(fromJSON(toJSON(content(r))))




# ## HERE WE WILL COLLAPSE duplicated probe-genename pairs WITHIN ANNOTATED GENE LIST (and also change column names to standarized ones)
# UNIQ_ANNOT_LIST_DATA <- list()
# for (n in WHICH_EXP_TO_ANAL){
#   UNIQ_ANNOT_LIST_DATA[[n]] <- unique(ANNOT_LIST_DATA[[n]])
#   colnames(UNIQ_ANNOT_LIST_DATA[[n]]) <- c("Probe_ID", "ensembl_gene_name")
# }
# 
# ## Here we will collapse all gene names for given probe
# for (n in WHICH_EXP_TO_ANAL){
#   UNIQ_ANNOT_LIST_DATA[[n]] <- aggregate(
#     ensembl_gene_name~Probe_ID, 
#     data = UNIQ_ANNOT_LIST_DATA[[n]], 
#     FUN = stringr::str_c) ## Here we aggregate the genenames into single row
#   UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- lapply(
#     X = UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name, 
#     FUN = paste, 
#     collapse = "; ") ##
#   UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- as.character(UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name) ## Here we make sure that collapsed genename column is not a list (need for writing function)
# }
# 
# # Je?eli w ensemble nie ma przypisanego genu (oznaczenie NA w przys?anym pliku "100_genow") to program pokazuje we wszystkich rubrykach oznaczenie NA chocia? te kom?rki zawiera?y wcze?niej dane (plik AAAreview11). W pe?nych danych kt?re przys?a?e? mi przed ?wi?tami (plik 1 znajduj?cy si? w skopresowanych danych "TEST_ANNOTATION") tych sond w og?le nie ma. Wygl?da na to, ?e to co przys?a?e? wczoraj odnosi si? chyba do wcze?niejszego pliku w kt?rym dane dla sond nie maj?cych przypisanego genu w ensemble zosta?y ca?kowicie usuni?te a to nie jest dobre :( - AMS W LIST_DATA S? TE DANE, ALE W TEST_ANNOTATION JUZ NIE
# 
# ## Here we merge probes annotated with ensembl with probes data that we have
# TEST_ANNOTATION <- list()
# for (n in WHICH_EXP_TO_ANAL){
#   TEST_ANNOTATION[[n]] <- merge(LIST_DATA[[n]], UNIQ_ANNOT_LIST_DATA[[n]], by = "Probe_ID", all.x = T)
# }
# 
# 
# 
# ## Here lets also remove empty lists
# ### !!! MANUAL SUPPLEMENTARY STEP - TO BE REMOVED IN FINAL VERSION, AS WE WILL NOT HAVE EMPTY LISTS !!! ###
# # TEST_ANNOTATION[[13]] <-NULL
# 
# # Here we bind all list elements into single data table
# SINGLE_TEST_ANNOTATION <- data.table::data.table(rlist::list.rbind(TEST_ANNOTATION), stringsAsFactors = F)
# readr::write_tsv(x = SINGLE_TEST_ANNOTATION, path = "SINGLE_ANNOTATION_probes.tsv")




# temp_entrez <- rentrez::entrez_search(db="gene", term = 'ACAD9 AND "Mus musculus"[ORGANISM]')
# 
# temp_mart <- biomaRt::useMart(
#   "ENSEMBL_MART_ENSEMBL", 
#   dataset = "mmusculus_gene_ensembl")
# tempFilters <- biomaRt::listFilters(mart = temp_mart)
# temp_annotation <- biomaRt::getBM(
#   attributes = c('entrezgene_id', "external_gene_name"),
#   filters = 'entrezgene_id', 
#   values = '106264', 
#   uniqueRows = T,
#   mart = temp_mart
# )
