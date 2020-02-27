GRS_input_preparation_script.R

# parse(text = 'input_7SK')
# 
# eval(parse(text = 'input_7SK'))
# 
# ls(pattern = '^input_.*')
# 
# test2 <- lapply(X = ls(pattern = '^input_.*'), FUN = function(x){eval(parse(text = x))})
# 
# # INPUT: regex describing which objects to take and include in the list
# list_env_objects <- function(regex_)
# {
#   print(environment())
#   
#   print(environment(.GlobalEnv))
#   print(environment())
#   list_output <-
#     lapply(
#       X = ls(pattern = regex_),
#       FUN = function(x) {
#         eval(parse(text = x), envir = globalenv())
#       }
#     )
#   return(list_output)
# }
# 
# test3 <- list_env_objects(regex_ = '^input_.*')
# 
# print(globalenv())
# environment()
# .GlobalEnv
# globalenv()
# ls(pattern = '^input_.*')
# 
# test <- eval(parse(text = 'input_7SK'))
# ### !!! How to create this object by ls() query and pass a list of objects returned by enviroment query into the function
# 
# 
# write_lists(list_LIST_DATA = list(input_ProbeID, input_Gemma, input_EnsemblGeneId, input_EnsemblTranscriptId, input_EntrezGeneId, input_EntrezGeneID_LOC, input_RefSeqMRNA, input_RefSeqMRNA_corrupted, input_mgi_symbol, input_mgi_symbol_corrupted, input_accession, input_rgd, input_gm, input_RNA, input_7SK, input_bad_gene_symbol, left_to_do_16), str_experiment_name = 'setting_input')
# 
# 
# 
# eval(parse(text = ls(pattern = '^(input\\_)(.*)')))
# test_files <- as.list(parse(text = ls(pattern = '^(input\\_)(.*)')))
# 
# 
# 
# check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'raw_dataset', df_original = raw_dataset, list_df_splited = chuj)
  


####### TRASH ####### 
# inputAnalysis_list_of_input_files_names <- c('data_v4_ensembl_gene_id.tsv', 'data_v4_ensembl_transcript_id.tsv', 'data_v4_entrezgene_id.tsv', 'data_v4_gemma.tsv', 'data_v4_manual.tsv', 'data_v4_names.tsv', 'data_v4_other_ids.tsv', 'data_v4_Probe_ID.tsv', 'data_v4_refseq_mrna_predicted.tsv', 'data_v4_refseq_mrna.tsv', 'data_v4_trash.tsv')
# 
# inputAnalysis_list_of_input_files <- list()
# inputAnalysis_list_of_input_files <-
#   lapply(
#     X = inputAnalysis_list_of_input_files_names,
#     FUN = function(x) {
#       temp <- read_preformated_data(str_filename = x)
#       
#     }
#   )
# 
# names(inputAnalysis_list_of_input_files) <- inputAnalysis_list_of_input_files_names
# 
# inputAnalysis_df_of_input_files <- rlist:: list.rbind(inputAnalysis_list_of_input_files)
# 
# inputAnalysis_lenght_of_wholeDataset_vs_inputFiles <- length(raw_dataset[[1]]) - length(inputAnalysis_df_of_input_files[[1]])
# ####### CHECK IF INPUT DATA IS GOOD. This is not needed. We will be using NEW IMPROVED DATASET ####### 


















GRS.R
# test_ <- get_query_for_ncbi_geneID_annotation(char_vec_gene_id_to_query_with = 'VAMP1', char_vec_organism = 'homo')[[1]]
# 
# test_xxx <- search_for_ids_in_ncbi(test_)
# 
# test_input_gene_symbol <- input_gene_symbol %>%
#   dplyr::select(Probe_ID, Entry_number)
# 
# test <- merge(test_input_gene_symbol, y = annotated_names_to_geneIDs, by = 'Entry_number', all.x = T)
# 
# test1 <- test[duplicated(test$Entry_number),]
# 
# "VAMP1[Gene Name] AND \"Homo sapiens\"[Organism]"
####### NCBI ANNOTATION ####### 
# ####### NCBI ANNOTATION ####### 
# ####### ENSEMBL ANNOTATION OF PREVIOUS NCBI ANNOTATION ####### 
# temp_annotated_geneIDs <- rbind(annotated_other_ids_to_geneIDs, annotated_names_to_geneIDs)
# 
# 
# save(annotated_names_to_geneIDs, file = 'annotated_names_to_geneIDs')
# # Prepare input for entrezgene_id annotation
# annotated_names_and_other_ids_to_geneIDs <- temp_annotated_geneIDs %>%
#   dplyr::filter(.data$external_gene_name != 'none' & !is.na(.data$external_gene_name)) %>%
#   dplyr::mutate(Probe_ID = external_gene_name) %>%
#   dplyr::select(-external_gene_name)
# 
# # Get non-anotated values
# not_annotated_data_1 <- temp_annotated_geneIDs %>%
#   dplyr::filter(.data$external_gene_name == 'none' | is.na(.data$external_gene_name))
# 
# check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'ncbi_annotation_for_geneNames_and_otherIDs', df_original = temp_annotated_geneIDs, list_df_splited = list(annotated_names_and_other_ids_to_geneIDs, not_annotated_data_1))
# 
# annotations_names_and_other_ids_to_geneIDs <- master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(df_to_annotate = annotated_names_and_other_ids_to_geneIDs, str_identifier_type___ = 'entrezgene_id', str_experiment_name__ = 'ncbi_annotations_annotated_as_entrezID')
# 
# rm(temp_annotated_geneIDs, annotated_names_and_other_ids_to_geneIDs)
# # annotated_short <- annotate_identifiers_to_geneID(str_filename_ = 'data_v4_names_short.tsv', str_experiment_name = 'short', descriptions_ = descriptions)
# ####### ENSEMBL ANNOTATION OF PREVIOUS NCBI ANNOTATION ####### 
#############################################
####### NCBI GENE NAME PRE-ANNOTATION ####### 
#############################################






#############################################
####### NCBI GENE NAME PRE-ANNOTATION ####### 
#############################################
# annotated_names_to_geneIDs <- annotate_identifiers_to_geneID(str_filename_ = 'checked_input/input_RGD_GM_RNA_7SK.tsv', str_experiment_name = 'other_ids', descriptions_ = descriptions, chr_gene_identifier_ = 'Gene Name') #4948 vs 
# save(annotated_names_to_geneIDs, file = 'annotated_names_to_geneIDs')

# # Pre-annotate and save first part of gene names
# annotated_other_ids_to_geneIDs <- annotate_identifiers_to_geneID(str_filename_ = 'checked_input/input_RGD_GM_RNA_7SK.tsv', str_experiment_name = 'other_ids', descriptions_ = descriptions, chr_gene_identifier_ = 'Gene Name') #515 vs 515
# save(annotated_other_ids_to_geneIDs, file = 'annotated_other_ids_to_geneIDs')
# # OR
# # annotated_other_ids_to_geneIDs <- readr::read_tsv(file = 'other_id/data_v4_for_ncbi_output_of_ncbi_query.tsv')
# 
# 
# # Pre-annotate and save first part of gene names
# annotated_names_to_geneIDs <- annotate_identifiers_to_geneID(str_filename_ = 'checked_input/input_gene_symbol.tsv', str_experiment_name = 'names', descriptions_ = descriptions, bool_reformat_names = F, chr_gene_identifier_ = 'Gene Name') #4433 vs 4433
# save(annotated_names_to_geneIDs, file = 'annotated_names_to_geneIDs')

# Prepare composite table for proper annotation
# annotated_gemma_names_and_other_ids_to_geneIDs <- rbind(annotated_other_ids_to_geneIDs, annotated_names_to_geneIDs, annotated_second_stage_input_gemma_probeID_gene_names)
#############################################
####### NCBI GENE NAME PRE-ANNOTATION ####### 
#############################################










###### PAPER-CENTRIC ANALYSIS ###### 


#### IF YOU WANT TO GET THE ANNOTATION PERCENTAGES, USE THIS ####
# annotation_percentages <-
#   purrr::map2(
#     .x = list(
#       annotations_ensembl_gene_id,
#       annotations_ensembl_transcript_id,
#       annotations_entrezgene_id,
#       annotations_refseq_mrna,
#       annotations_mgi,
#       annotations_accession
#     ),
#     .y = list(
#       'annotations_ensembl_gene_id',
#       'annotations_ensembl_transcript_id',
#       'annotations_entrezgene_id',
#       'annotations_refseq_mrna',
#       'annotations_mgi',
#       'annotations_accession'
#     ),
#     .f = function(x, y) {
#       temp <- sum(!is.na(x$external_gene_name)) / length(x[[1]])
#       names(temp) <- y
#       return(temp)
#     }
#   )
#### IF YOU WANT TO GET THE ANNOTATION PERCENTAGES, USE THIS ####




####### GEMMA SECOND STAGE ANNOTATION ####### 
# save(second_stage_input_gemma_probeID_gene_names)
# save(leftovers_second_stage_input_gemma_probeID)
# 
# annotated_second_stage_input_gemma_probeID_gene_names <- subset(x = annotated_second_stage_input_gemma_probeID_gene_names, subset = !is.na(annotated_second_stage_input_gemma_probeID_gene_names$external_gene_name)) ### !!! gemma with entrez ids for ensembl
# 
# pre_finalized_gemma <- subset(x = annotated_second_stage_input_gemma_probeID_gene_names, subset = !is.na(annotated_second_stage_input_gemma_probeID_gene_names$external_gene_name)) ### !!! gemma without entrez ids for ensembl
# 
# leftovers_gemma <- rbind(leftovers_second_stage_input_gemma_probeID, leftovers_annotated_second_stage_input_gemma_probeID_gene_names)
# 
# 
# length(pre_finalized_gemma)+length(leftovers_gemma)
# 
# rm(leftovers_second_stage_input_gemma_probeID, leftovers_annotated_second_stage_input_gemma_probeID_gene_names)
# 
# 
# finalized_gemm <- master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(descriptions__ = descriptions, df_to_annotate = pre_finalized_gemma, str_identifier_type___ = 'entrezgene_id', str_experiment_name__ = 'gemma_annotation')
####### GEMMA SECOND STAGE ANNOTATION ####### 






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















functions_for_genename_conversions.R
# do_directions_for_multiple_gene_instances_within_experiment_match_old <- function(df_merged_dataset_)
# {
#   ### This will be our input data for further analysis. It contains only significant columns of "Paper", "GroupID", "ensembl_gene_name", "logFC":
#   SHORT_SINGLE_TEST_ANNOTATION <- df_merged_dataset_ %>%
#     dplyr::select("Paper", "Experiment", "external_gene_name", "logFC") %>%
#     dplyr::group_nest(Paper, Experiment, external_gene_name)
#   
#   # I dont know why purrr::map didnt work for this. Either way, we just lapply this. Here we write UP or DOWN for every logFC for given genename in given experiment (not paper)
#   SHORT_SINGLE_TEST_ANNOTATION$directionality <-
#     lapply(
#       X = SHORT_SINGLE_TEST_ANNOTATION$data,
#       FUN = function(x) {
#         dplyr::mutate(x, Symbol_direction = ifelse(logFC > 0, "UP", "DOWN"))
#       }
#     )
#   
#   # Here we count how many UPs and DOWNs are in every given genename in given experiment (not paper)
#   SHORT_SINGLE_TEST_ANNOTATION$directionality <-
#     lapply(
#       X = SHORT_SINGLE_TEST_ANNOTATION$directionality,
#       FUN = function (x) {
#         table(dplyr::select(x, "Symbol_direction"))
#       }
#     )
#   
#   # Here we establish actual status of gene in given experiment
#   SHORT_SINGLE_TEST_ANNOTATION$sum_directionality <-
#     as.character(lapply(
#       X = SHORT_SINGLE_TEST_ANNOTATION$directionality,
#       FUN = function(x) {
#         if (length(x) == 1 && grepl(pattern = "UP", names(x))) {
#           "UP"
#         }
#         else if (length(x) == 1 &&
#                  grepl(pattern = "DOWN", names(x))) {
#           "DOWN"
#         }
#         else if (length(x) == 2 &&
#                  grepl(pattern = "DOWN", names(x)) &&
#                  grepl(pattern = "DOWN", names(x))) {
#           "MIXED"
#         }
#         else {
#           "ERROR"
#         }
#       }
#     ))
#   SHORT_SINGLE_TEST_ANNOTATION$directionality <- NULL
#   
#   
#   UNNEST_SHORT_SINGLE_TEST_ANNOTATION <-
#     SHORT_SINGLE_TEST_ANNOTATION %>%
#     tidyr::unnest()
#   
#   return(UNNEST_SHORT_SINGLE_TEST_ANNOTATION)
#   
#   # # Here we remove MIXED expression and multiple genenames
#   # FILT_SHORT_SINGLE_TEST_ANNOTATION <- SHORT_SINGLE_TEST_ANNOTATION %>%
#   #   #filter(!grepl(pattern = "(.*);(.*)", SHORT_SINGLE_TEST_ANNOTATION$ensembl_gene_name)) %>%
#   #   dplyr::filter(!sum_directionality == "MIXED") %>%
#   #   dplyr::filter(!sum_directionality == "ERROR") %>%
#   #   dplyr::filter(!is.na(ensembl_gene_name))
#   # 
#   # # Here we add mean to each !!! Perhaps this should be final input, because it has: 1) removed genes giving UP+DOWN in the same experiment, 2) removed NAs
#   # STDINPUT_FILT_SHORT_SIN_T_ANNO <- FILT_SHORT_SINGLE_TEST_ANNOTATION %>%
#   #   mutate(mean = as.numeric(map(data, function(x) { as.numeric(mean(x[[1]]))} ))) %>% ## This way we give actuall vector to function, not a data table(tibble)
#   #   select("Paper", "GroupID", "ensembl_gene_name", "sum_directionality", "mean") 
# }



# # This function is not to be used, as it is not working as expected
# change_name_to_proper_format_given_the_species <- function(df_name, str_species)
# {
#   df_name$Old_Probe_ID <- df_name$Probe_ID
#   
#   if(str_species %in% c('mice', 'rats'))
#   {
#     temp1 <- subset(x = df_name, str_detect(string = df_name$Probe_ID, pattern = '(?i)(loc)|(rgd)')) #case insensitive
#     temp2 <- subset(x = df_name, str_detect(string = df_name$Probe_ID, pattern = '(?i)(rik)'))
#     temp3 <- subset(x = df_name, !str_detect(string = df_name$Probe_ID, pattern = '(?i)(rik)|(loc)|(rgd)')) 
#     
#     temp1$Probe_ID <- toupper(temp1$Probe_ID)
#     
#     temp2$Probe_ID <- toupper(temp2$Probe_ID)
#     temp2$Probe_ID <- str_replace(string = temp2$Probe_ID, pattern = '(?i)(rik)', replacement = 'Rik')
#     
#     temp3$Probe_ID <- tolower(temp3$Probe_ID)
#     substr(temp3$Probe_ID, start = 1, stop = 1) <- toupper(substr(temp3$Probe_ID, 1, 1))
#     
#     
#     temp4 <- rbind(temp1, temp2, temp3)
#     
#     return(temp4)
#   }
#   else if(str_species %in% c('humans', 'squirrelmonkeys'))
#   {
#     temp <- df_name
#     temp$Probe_ID <- toupper(temp$Probe_ID)
#     return(temp)
#   }
# }



# merge_and_remove_nas_from_list_post_annotation <- function(list_annotated_LIST_DATA, str_name_of_column_containing_annotated_gene_symbols)
# {
#   temp_df <- rlist::list.rbind(list_annotated_LIST_DATA)
#   temp_df <- subset(x = temp_df, subset = !is.na(temp_df[[str_name_of_column_containing_annotated_gene_symbols]]))
#   return(temp_df)
# }

