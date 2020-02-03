############################################
####### SET OPTIONS FOR THE ANALYSIS ####### 
############################################
# setwd('/media/adrians/USB DISK1/Projekty/GRS - GJt Review Stress/FULL_DATASET')
# source('functions_for_genename_conversions.R')
# source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# rm(list = ls(pattern = '(.*)(temp)|(test)(.*)'))
# rm(list = ls(pattern = '^(input)(.*)|^left_(.*)'))
# To see is na is introduced length( subset(x = input_EnsemblGeneId$Probe_ID, subset = is.na(input_EnsemblGeneId$Probe_ID)) )


# temp <- biomaRt::getBM(
#   attributes = c('entrezgene_id', "external_gene_name"),
#   filters = 'entrezgene_id', 
#   values = 'LOC16948', 
#   uniqueRows = T,
#   mart = biomaRt::useMart(
#     "ENSEMBL_MART_ENSEMBL", 
#     dataset = "mmusculus_gene_ensembl")
# )

# left_to_do_temp <- left_to_do_
# left_to_do_ <- left_to_do_temp
# temp <- subset(x = input_EntrezGeneId$Gene_ID, subset = stringr::str_detect(string = input_EntrezGeneId$Gene_ID, pattern = ','))



WITHIN_COL_SEPARATOR <- ','

stage = 4

if (stage == 1) 
{
  opts_everything_column_for_this_analysis <- 'everything' 
  opts_probe_IDs_already_done_columnName_list <- NA
} else if (stage == 2) 
{
  opts_everything_column_for_this_analysis <- 'everything2'
  opts_probe_IDs_already_done_columnName_list <- list('Probe_ID_left_from_first_annotating_stage', 'Probe_ID_left_from_first_ncbi_annotation_stage')
  opts_input_folder_name = 'checked_input_2_stage'
  load('leftovers')
  left_to_do_ <- leftovers # When working with leftovers (second stage) start at this etage
} else if (stage == 3) 
{
  opts_everything_column_for_this_analysis <- 'everything2'
  opts_probe_IDs_already_done_columnName_list <- list('Probe_ID_left_from_first_annotating_stage', 'Probe_ID_left_from_first_ncbi_annotation_stage', 'Probe_ID_left_from_second_annotating_stage', 'Probe_ID_left_from_second_ncbi_annotation_stage')
  opts_input_folder_name = 'checked_input_3_stage'
  load('leftovers_2')
  left_to_do_ <- leftovers
} else if (stage == 4) 
{
  opts_everything_column_for_this_analysis <- 'everything2'
  opts_probe_IDs_already_done_columnName_list <- list('Probe_ID_left_from_first_annotating_stage', 'Probe_ID_left_from_first_ncbi_annotation_stage', 'Probe_ID_left_from_second_annotating_stage', 'Probe_ID_left_from_second_ncbi_annotation_stage', 'Probe_ID_left_from_third_annotating_stage', 'Probe_ID_left_from_third_ncbi_annotation_stage')
  opts_input_folder_name = 'checked_input_4_stage'
  load('leftovers_3')
  left_to_do_ <- leftovers
}



##############################################
####### PREPARE RAW DATA FOR SUBSETING ####### 
##############################################
if (stage == 1) {
  raw_dataset <- read_preformated_data(str_filename = 'data_whole.tsv', col_types_ = 'nccccccccccccc', int_numbers_are_from = 11, int_numbers_are_to = 13) # In original file 'data_whole.tsv': 1) there are 261697 including the header; 2) there are 261475 values in first column (Paper). This is because some rows are empty (separator rows between experiments/papers). In raw_dataset there are 261696 rows.
  
  # After removing empty rows, which are defined as rows, that have no value in 'Paper' column, we are left with 261475 value-rich rows. This is in agreement with 'data_whole.tsv' file
  raw_dataset <- subset(x = raw_dataset, subset = !is.na(raw_dataset$Paper))
  
  # Prepare two columns which will include unique identfiers for given row. We need a column including all input for the entry (except entry_number, which will be mistaken for Gene_ID) for easy query of given identifer type in entire entry - column 'everything'
  raw_dataset$Entry_number <- rownames(raw_dataset)
  
  reformated_raw_dataset <- reformat_Gene_ID_column_wrapper(raw_dataset)
  
  reformated_raw_dataset <- reformat_GenBank_Accession_column_wrapper(reformated_raw_dataset)
  
  reformated_raw_dataset <- reformat_Ensembl_ID_column_wrapper(reformated_raw_dataset)
  
  reformated_raw_dataset <- reformat_Gene_symbol_column_wrapper(reformated_raw_dataset)
  
  # This column will include all the identfiers that were.
  reformated_raw_dataset$everything <- as.character(apply(X = reformated_raw_dataset, MARGIN = 1, FUN = function(x) { paste( c(x[1:9], x[11:14]), collapse = '__')  } ))
  # This column will include all the identfiers except the once that were already used
  reformated_raw_dataset$everything2 <- reformated_raw_dataset$everything
  
  # Here we need to archive original Probe_ID column, and prepare Probe_ID column to be queried against in all subsequent functions. This is because I hardcoded Probe_ID column in analysis functions. Go me.
  reformated_raw_dataset$Probe_ID_old <- reformated_raw_dataset$Probe_ID
  reformated_raw_dataset$Probe_ID <- NA

  # Count entries in raw data. These lengths were checked by GJt and are in agreement with his raw data. Thus I conclude that it is very likely that raw_dataset is good input for further analysis
  inputAnalysis_lengths_of_experiments_in_raw_dataset <- split_and_measure_length(df_ = reformated_raw_dataset, split_by_str = 'Experiment')
  ##############################################
  ####### PREPARE RAW DATA FOR SUBSETING ####### 
  ##############################################
  #############################################################
  ####### PREPARE DATA FOR MICROARRAY-CENTERED ANALYSIS ####### 
  #############################################################    
  # Subset papers to be analyzed via actual microarray Probe_ID. Make sure papers are in correct order.
  Papers_to_analyze_via_ProbeID <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 42, 43, 54, 56, 59, 60)
  
  inputAnalysis_include_in_ProbeID <- reformated_raw_dataset$Paper %in% Papers_to_analyze_via_ProbeID
  
  input_ProbeID <- subset(x = reformated_raw_dataset, subset = inputAnalysis_include_in_ProbeID)
  
  left_to_do <- subset(x = reformated_raw_dataset, subset = !inputAnalysis_include_in_ProbeID)
  
  check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_ProbeID', df_original = reformated_raw_dataset, list_df_splited = list(input_ProbeID, left_to_do))
  
  # Proper Probe_IDs are in different columns depending on the Paper. Thus we need to ascribe Probe_ID column with appropriate values
  temp_Probe_ID_old <- subset(x = input_ProbeID, subset = input_ProbeID$Paper %in% c(1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 42, 43, 54, 56, 59, 60))
  temp_Probe_ID_old$Probe_ID <- temp_Probe_ID_old$Probe_ID_old
  
  temp_GEO_ID <- subset(x = input_ProbeID, subset = input_ProbeID$Paper %in% c(2, 4))
  temp_GEO_ID$Probe_ID <- temp_GEO_ID$`GEO ID`
  
  input_ProbeID <- rbind(temp_Probe_ID_old, temp_GEO_ID)
  input_ProbeID <- remove_used_name_from_everything2_wrapper(input_ProbeID)

  input_ProbeID <- input_ProbeID[order(input_ProbeID$Paper),] 
  
  # Quality checks
  check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_ProbeID_vs_GEO', df_original = input_ProbeID, list_df_splited = list(temp_Probe_ID_old, temp_GEO_ID))
  get_all_symbols_in_chrvec(input_ProbeID$Probe_ID)
  #############################################################
  ####### PREPARE DATA FOR MICROARRAY-CENTERED ANALYSIS ####### 
  #############################################################
  ###################################################################
  ####### PREPARE DATA FOR GEMMA MICROARRAY-CENTERED ANALYSIS ####### 
  ###################################################################
  Papers_to_analyze_via_Gemma <- c(31, 33)
  
  inputAnalysis_include_in_Gemma <- left_to_do$Paper %in% Papers_to_analyze_via_Gemma
  
  input_Gemma <- subset(x = left_to_do, subset = inputAnalysis_include_in_Gemma)
  
  left_to_do_ <- subset(x = left_to_do, subset = !inputAnalysis_include_in_Gemma)
  
  input_Gemma$Probe_ID <- input_Gemma$Probe_ID_old
  
  input_Gemma <- remove_used_name_from_everything2_wrapper(input_Gemma)
  
  input_Gemma <- input_Gemma[order(input_Gemma$Paper),] 
  
  # Quality checks
  check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_Gemma', df_original = left_to_do, list_df_splited = list(input_Gemma, left_to_do_))
  get_all_symbols_in_chrvec(input_Gemma$Probe_ID)
  
  rm(inputAnalysis_include_in_ProbeID, temp_Probe_ID_old, temp_GEO_ID, Papers_to_analyze_via_ProbeID, inputAnalysis_include_in_Gemma, Papers_to_analyze_via_Gemma, left_to_do)
}
###################################################################
####### PREPARE DATA FOR GEMMA MICROARRAY-CENTERED ANALYSIS ####### 
###################################################################



######################################################
####### PREPARE DATA FOR ENSEMBL GENE ANALYSIS ####### 
######################################################
# This detects all entries in ensembl_gene_id from data_v4, Ive checked. Yet, the lenght of this input is much lower than of data_v4_ensembl_gene_id.tsv?
temp <- prepare_input(regex_ = '(ENS)(.*)G0(.*)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '(__.*)|(,.*)') ### !!! Check czy na następnych etapach skrypt bierze tutaj kolejny identyfikator


input_EnsemblGeneId <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
get_all_symbols_in_chrvec(input_EnsemblGeneId$Probe_ID)

input_EnsemblGeneId <- remove_used_name_from_everything2_wrapper(input_EnsemblGeneId)
######################################################
####### PREPARE DATA FOR ENSEMBL GENE ANALYSIS ####### 
######################################################



############################################################
####### PREPARE DATA FOR ENSEMBL TRANSCRIPT ANALYSIS ####### 
############################################################
temp <- prepare_input(regex_ = '(ENS)(.*)(T0)(.*)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '(__.*)|(,.*)') ### !!! Check czy na następnych etapach skrypt bierze tutaj kolejny identyfikator

input_EnsemblTranscriptId <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
get_all_symbols_in_chrvec(input_EnsemblTranscriptId$Probe_ID)

input_EnsemblTranscriptId <- remove_used_name_from_everything2_wrapper(input_EnsemblTranscriptId)
############################################################
####### PREPARE DATA FOR ENSEMBL TRANSCRIPT ANALYSIS ####### 
############################################################


# temp_left_to_do_ <- left_to_do_
#########################################################################
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
#########################################################################
temp <- prepare_input(regex_ = '[1-9](\\d*)', df_ = left_to_do_, col_everything = 'Gene_ID', col_get_resulting_identifer_from_here = 'Gene_ID')


input_EntrezGeneId <- temp[[1]] 
input_EntrezGeneId$Probe_ID <- get_geneNames_filtered_with_already_done_columns(input_gene_symbol_ = input_EntrezGeneId, colname_with_output_probeNames = 'Gene_ID')
### !!! testing
# test <- input_EntrezGeneId[c('7494', '49153'),]
# test$Probe_ID_left_from_first_annotating_stage <- '100036568'
# test$Probe_ID_left_from_first_ncbi_annotation_stage <- '108978'
# 
# test$Probe_ID <- get_geneNames_filtered_with_already_done_columns(input_gene_symbol_ = test, colname_with_output_probeNames = 'Gene_ID')
### !!! testing


left_to_do_ <- temp[[2]]
temp[[3]]
get_all_symbols_in_chrvec(input_EntrezGeneId$Probe_ID)

#This moves entries with empty probe_id (which means, that probe_id was removed in get_geneNames_filtered_with_already_done_columns) to left_to_do_
if(stage != 1)
{
  empty_ids <- subset(x = input_EntrezGeneId, subset = input_EntrezGeneId$Probe_ID == '')
  good_ids <- subset(x = input_EntrezGeneId, subset = input_EntrezGeneId$Probe_ID != '')
  
  check_was_the_spliting_of_df_by_filtering_ok(df_original = input_EntrezGeneId, list_df_splited = list(empty_ids, good_ids))
  
  input_EntrezGeneId <- good_ids
  left_to_do_ <- rbind(left_to_do_, empty_ids)
}

# This removes LOCs corresponding to entrezIDs that were used in this step
input_EntrezGeneId <- remove_used_entrezID_from_everything2_wrapper(input_EntrezGeneId)
####### PREPARE DATA FOR ENTREZ GENE ANALYSIS ####### 
####### PREPARE DATA FOR ENTREZ-LOC GENE ANALYSIS ####### 
# Subset entries to be analyzed via Entrez Gene from LOC numbers - LOC - genes of uncertain function. When a published symbol is not available, and orthologs have not yet been determined, Gene will provide a symbol that is constructed as 'LOC' + the GeneID. Therefore LOC is basically GeneID, and is thus unique


temp <- prepare_input(regex_ = '(LOC)(.*)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '__.*|(,.*)')

input_EntrezGeneID_LOC <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
get_all_symbols_in_chrvec(input_EntrezGeneID_LOC$Probe_ID)

input_EntrezGeneID_LOC <- remove_used_name_from_everything2_wrapper(input_EntrezGeneID_LOC)

input_EntrezGeneID_LOC$Probe_ID <- stringr::str_remove_all(string = input_EntrezGeneID_LOC$Probe_ID, pattern = 'LOC')
### !!! BEWARE! entrez loc depend on entrez gene id in such way, that if in previous stage entrez number for given LOC was used, than the number part of this LOC is also removed. This may lead to empty Probe_ID number for this LOC. Hence, we need to search for empty '' input_EntrezGeneID_LOC$Probe_ID entries and put them into left_to_do_

input_EntrezGeneId <- rbind(input_EntrezGeneId, input_EntrezGeneID_LOC)

input_EntrezGeneId <- input_EntrezGeneId[order(input_EntrezGeneId$Paper),] 

get_all_symbols_in_chrvec(input_EntrezGeneId$Probe_ID)

rm(input_EntrezGeneID_LOC, empty_ids, good_ids)
####### PREPARE DATA FOR ENTREZ-LOC GENE ANALYSIS ####### 
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
#########################################################################
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
#########################################################################



####################################################
####### PREPARE DATA FOR REFSEQMRNA ANALYSIS ####### 
####################################################
temp <- prepare_input(regex_ = '(NM_)(\\d*)(.*)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '(__.*)|(,.*)')

input_RefSeqMRNA <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
get_all_symbols_in_chrvec(input_RefSeqMRNA$Probe_ID)

input_RefSeqMRNA <- remove_used_name_from_everything2_wrapper(input_RefSeqMRNA)

# This is because biomart does not recognize transcript variants (so transcript names which have '.' after the name)
input_RefSeqMRNA$Probe_ID <- stringr::str_remove_all(string = input_RefSeqMRNA$Probe_ID, pattern = '\\..*')
####################################################
####### PREPARE DATA FOR REFSEQMRNA ANALYSIS ####### 
####################################################



###################################################
####### PREPARE DATA FOR ACCESSION ANALYSIS #######
###################################################
# Study leftover ids: XM_ - most of these names are substituted by NM_ sequences, but NCBI search does not return this new NM_ gene. Stupid.
# Study leftover ids: [letter][numbers] - a) ncbi accession nb. An accession number applies to the complete record and is usually a combination of a letter(s) and numbers, such as a single letter followed by five digits (e.g., U12345) or two letters followed by six digits (e.g., AF123456). Some accessions might be longer, depending on the type of sequence record. Accession numbers do not change, even if information in the record is changed at the author's request. Sometimes, however, an original accession number might become secondary to a newer accession number, if the authors make a new submission that combines previous sequences, or if for some reason a new submission supercedes an earlier record. These IDs are actually the same in ENA ('embl' or 'clone_based_ensembl_gene') and in ncbi nucleotide
temp <- prepare_input(regex_ = '(?<![A-Za-z]{1})[A-Z]{1,2}\\d{5,}(?=(_|,|\\.))', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '(__.*)|(,.*)')
# temp <- prepare_input(regex_ = '(?<![A-Za-z]{1})[A-Z]{1,2}\\d{5,}', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '(__.*)|(,.*)')

input_accession <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]

get_all_symbols_in_chrvec(input_accession$Probe_ID)

input_accession <- remove_used_name_from_everything2_wrapper(input_accession)
###################################################
####### PREPARE DATA FOR ACCESSION ANALYSIS #######
###################################################



##################################################################
####### PREPARE DATA FOR GOOD AND BAD GENE_SYMBOL ANALYSIS #######
##################################################################
input_bad_gene_symbol <- subset(x = left_to_do_, subset = is.na(left_to_do_$Gene_symbol))

input_gene_symbol <- subset(x = left_to_do_, subset = !is.na(left_to_do_$Gene_symbol))

# Quality checks
check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_gene_symbol', df_original = left_to_do_, list_df_splited = list(input_bad_gene_symbol, input_gene_symbol))



get_all_symbols_in_chrvec(input_gene_symbol$Probe_ID)

input_gene_symbol$Probe_ID <- input_gene_symbol$Gene_symbol
input_gene_symbol$Probe_ID <- get_geneNames_filtered_with_already_done_columns(input_gene_symbol_ = input_gene_symbol)



input_gene_symbol <- input_gene_symbol[order(input_gene_symbol$Paper),] 
input_bad_gene_symbol <- input_bad_gene_symbol[order(input_bad_gene_symbol$Paper),]


#This moves entries with empty probe_id (which means, that probe_id was removed in get_geneNames_filtered_with_already_done_columns) to left_to_do_
if(stage != 1)
{
  empty_ids <- subset(x = input_gene_symbol, subset = input_gene_symbol$Probe_ID == '')
  good_ids <- subset(x = input_gene_symbol, subset = input_gene_symbol$Probe_ID != '')
  
  check_was_the_spliting_of_df_by_filtering_ok(df_original = input_gene_symbol, list_df_splited = list(empty_ids, good_ids))
  
  input_gene_symbol <- good_ids
  input_bad_gene_symbol <- rbind(input_bad_gene_symbol, empty_ids)
}

rm(left_to_do_, empty_ids, good_ids)
##################################################################
####### PREPARE DATA FOR GOOD AND BAD GENE_SYMBOL ANALYSIS #######
##################################################################



########################################
####### CHECK AND WRITE ALL DATA #######
########################################
# Here I show that dataset composed of all the specified subsets is equal to raw_dataset
if (stage == 1)
{
  rebuild_dataset_list_first_stage <- list(input_ProbeID, input_Gemma, input_EnsemblGeneId, input_EnsemblTranscriptId, input_EntrezGeneId, input_RefSeqMRNA, input_accession, input_bad_gene_symbol, input_gene_symbol)
  rebuild_dataset_list_names_first_stage <- list('input_ProbeID', 'input_Gemma', 'input_EnsemblGeneId', 'input_EnsemblTranscriptId', 'input_EntrezGeneId', 'input_RefSeqMRNA', 'input_accession', 'input_bad_gene_symbol', 'input_gene_symbol') # FIRST STAGE
  
  rebuild_dataset <- rlist::list.rbind(rebuild_dataset_list_first_stage)
  
  write_inputs(lists_ = rebuild_dataset_list_first_stage, list_names_list_str = rebuild_dataset_list_names_first_stage, dir_name_str = 'checked_input_1_stage')

} else if (stage != 1) 
{
  rebuild_dataset_list_second_stage <- list(input_EnsemblGeneId, input_EnsemblTranscriptId, input_EntrezGeneId, input_RefSeqMRNA, input_accession, input_bad_gene_symbol, input_gene_symbol)
  rebuild_dataset_list_names_second_stage <- list('input_EnsemblGeneId', 'input_EnsemblTranscriptId', 'input_EntrezGeneId', 'input_RefSeqMRNA', 'input_accession', 'input_bad_gene_symbol', 'input_gene_symbol')
  
  rebuild_dataset <- rlist::list.rbind(rebuild_dataset_list_second_stage)
  
  write_inputs(lists_ = rebuild_dataset_list_second_stage, list_names_list_str = rebuild_dataset_list_names_second_stage, dir_name_str = opts_input_folder_name)
}

rm(list = ls(pattern = 'left_to(.*)|inputAnalysis_(.*)|pattern_for(.*)|Papers_to(.*)|rebuild_dataset_list(.*)'))
########################################
####### CHECK AND WRITE ALL DATA #######
########################################
