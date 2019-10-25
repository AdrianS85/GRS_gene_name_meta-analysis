setwd('/media/adrians/USB DISK1/Projekty/GRS - GJt Review Stress/FULL_DATASET')
source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
rm(list = ls(pattern = '(.*)(temp)|(test)(.*)'))
# To see is na is introduced length( subset(x = input_EnsemblGeneId$Probe_ID, subset = is.na(input_EnsemblGeneId$Probe_ID)) )



##############################################
####### PREPARE RAW DATA FOR SUBSETING ####### 
##############################################
raw_dataset <- read_preformated_data(str_filename = 'data_whole.tsv', col_types_ = 'nccccccccccccc', int_numbers_are_from = 11, int_numbers_are_to = 13) # In original file 'data_whole.tsv': 1) there are 261697 including the header; 2) there are 261475 values in first column (Paper). This is because some rows are empty (separator rows between experiments/papers). In raw_dataset there are 261696 rows.

# After removing empty rows, which are defined as rows, that have no value in 'Paper' column, we are left with 261475 value-rich rows. This is in agreement with 'data_whole.tsv' file
raw_dataset <- subset(x = raw_dataset, subset = !is.na(raw_dataset$Paper))

# Prepare two columns which will include unique identfiers for given row. We need a column including all input for the entry (except entry_number, which will be mistaken for Gene_ID) for easy query of given identifer type in entire entry - column 'everything'
raw_dataset$Entry_number <- rownames(raw_dataset)
raw_dataset$everything <- as.character(apply(X = raw_dataset, MARGIN = 1, FUN = function(x) { paste( c(x[1:9], x[11:14]), collapse = '__')  } ))

# Here we need to archive original Probe_ID column, and prepare Probe_ID column to be queried against in all subsequent functions. This is because I hardcoded Probe_ID column in analysis functions. Go me.
raw_dataset$Probe_ID_old <- raw_dataset$Probe_ID
raw_dataset$Probe_ID <- NA


# Count entries in raw data. These lengths were checked by GJt and are in agreement with his raw data. Thus I conclude that it is very likely that raw_dataset is good input for further analysis
inputAnalysis_lengths_of_experiments_in_raw_dataset <- split_and_measure_length(df_ = raw_dataset, split_by_str = 'Experiment')
##############################################
####### PREPARE RAW DATA FOR SUBSETING ####### 
##############################################



#############################################################
####### PREPARE DATA FOR MICROARRAY-CENTERED ANALYSIS ####### 
#############################################################    
# Subset papers to be analyzed via actual microarray Probe_ID. Make sure papers are in correct order.
Papers_to_analyze_via_ProbeID <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 42, 43, 54, 56, 59, 60, 72)

inputAnalysis_include_in_ProbeID <- raw_dataset$Paper %in% Papers_to_analyze_via_ProbeID

input_ProbeID <- subset(x = raw_dataset, subset = inputAnalysis_include_in_ProbeID)

left_to_do <- subset(x = raw_dataset, subset = !inputAnalysis_include_in_ProbeID)

inputAnalysis_was_spliting_good_1 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_ProbeID', df_original = raw_dataset, list_df_splited = list(input_ProbeID, left_to_do))

# Proper Probe_IDs are in different columns depending on the Paper. Thus we need to ascribe Probe_ID column with appropriate values
temp_Probe_ID_old <- subset(x = input_ProbeID, subset = input_ProbeID$Paper %in% c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 42, 43, 54, 56, 59, 60, 72))
temp_Probe_ID_old$Probe_ID <- temp_Probe_ID_old$Probe_ID_old

temp_GEO_ID <- subset(x = input_ProbeID, subset = input_ProbeID$Paper == 2)
temp_GEO_ID$Probe_ID <- temp_GEO_ID$`GEO ID`

inputAnalysis_was_spliting_good_ProbeID_vs_GEOID <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_ProbeID_vs_GEO', df_original = input_ProbeID, list_df_splited = list(temp_Probe_ID_old, temp_GEO_ID))

input_ProbeID <- rbind(temp_Probe_ID_old, temp_GEO_ID)
input_ProbeID <- input_ProbeID[order(input_ProbeID$Paper),] 

rm(inputAnalysis_include_in_ProbeID, temp_Probe_ID_old, temp_GEO_ID)
#############################################################
####### PREPARE DATA FOR MICROARRAY-CENTERED ANALYSIS ####### 
#############################################################



###################################################################
####### PREPARE DATA FOR GEMMA MICROARRAY-CENTERED ANALYSIS ####### 
###################################################################
Papers_to_analyze_via_Gemma <- c(31, 33)

inputAnalysis_include_in_Gemma <- left_to_do$Paper %in% Papers_to_analyze_via_Gemma

input_Gemma <- subset(x = left_to_do, subset = inputAnalysis_include_in_Gemma)

left_to_do_2 <- subset(x = left_to_do, subset = !inputAnalysis_include_in_Gemma)

inputAnalysis_was_spliting_good_2 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_Gemma', df_original = left_to_do, list_df_splited = list(input_Gemma, left_to_do_2))

rm(inputAnalysis_include_in_Gemma)

input_Gemma$Probe_ID <- input_Gemma$Probe_ID_old

input_Gemma <- input_Gemma[order(input_Gemma$Paper),] 
###################################################################
####### PREPARE DATA FOR GEMMA MICROARRAY-CENTERED ANALYSIS ####### 
###################################################################



#######################################################
####### PREPARE LEFTOVER RAW DATA FOR SUBSETING ####### 
#######################################################
# When working with leftovers, start at this etage
left_to_do <- leftovers
#######################################################
####### PREPARE LEFTOVER RAW DATA FOR SUBSETING ####### 
#######################################################



######################################################
####### PREPARE DATA FOR ENSEMBL GENE ANALYSIS ####### 
######################################################
# This detects all entries in ensembl_gene_id from data_v4, Ive checked. Yet, the lenght of this input is much lower than of data_v4_ensembl_gene_id.tsv?
temp <- prepare_input(regex_ = '(ENS)(.*)G0(.*)', df_ = left_to_do, regex_two_for_extracting_identifiers_from_string = '__.*', vec_regex_additional_pattern_removal = '( \\+ )(.*)')

input_EnsemblGeneId <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
######################################################
####### PREPARE DATA FOR ENSEMBL GENE ANALYSIS ####### 
######################################################



############################################################
####### PREPARE DATA FOR ENSEMBL TRANSCRIPT ANALYSIS ####### 
############################################################
temp <- prepare_input(regex_ = '(ENS)(.*)(T0)(.*)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '__.*')

input_EnsemblTranscriptId <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
############################################################
####### PREPARE DATA FOR ENSEMBL TRANSCRIPT ANALYSIS ####### 
############################################################



#########################################################################
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
#########################################################################
temp <- prepare_input(regex_ = '[1-9](\\d*)', df_ = left_to_do_, col_everything = 'Gene_ID', col_get_resulting_identifer_from_here = 'Gene_ID')

input_EntrezGeneId <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
####### PREPARE DATA FOR ENTREZ GENE ANALYSIS ####### 
####### PREPARE DATA FOR ENTREZ-LOC GENE ANALYSIS ####### 
# Subset entries to be analyzed via Entrez Gene from LOC numbers - LOC - genes of uncertain function. When a published symbol is not available, and orthologs have not yet been determined, Gene will provide a symbol that is constructed as 'LOC' + the GeneID. Therefore LOC is basically GeneID, and is thus unique
temp <- prepare_input(regex_ = '(LOC)(.*)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '__.*', vec_regex_additional_pattern_removal = c('(LOC)', '///(.*)'))

input_EntrezGeneID_LOC <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
####### PREPARE DATA FOR ENTREZ-LOC GENE ANALYSIS ####### 
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
input_EntrezGeneId <- rbind(input_EntrezGeneId, input_EntrezGeneID_LOC)

input_EntrezGeneId <- input_EntrezGeneId[order(input_EntrezGeneId$Paper),] 

rm(input_EntrezGeneID_LOC)
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
#########################################################################
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
#########################################################################



####################################################
####### PREPARE DATA FOR REFSEQMRNA ANALYSIS ####### 
####################################################
temp <- prepare_input(regex_ = '(NM_)(\\d*)(.*)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '__.*', vec_regex_additional_pattern_removal = c('\\.(.*)', ' /// (.*)'))

input_RefSeqMRNA <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
####### PREPARE DATA FOR REFSEQMRNA ANALYSIS ####### 
####### PREPARE DATA FOR REFSEQMRNA CORRUPTED ANALYSIS ####### 
# Subset corrupted (no _ after NM) entries to be analyzed via RefSeqMRNA
temp <- prepare_input(regex_ = '(NM )(\\d*)(.*)', df_ = left_to_do_, col_get_resulting_identifer_from_here = 'GenBank_Accession', vec_regex_additional_pattern_removal = '\\.(.*)')

input_RefSeqMRNA_corrupted <- temp[[1]] 

input_RefSeqMRNA_corrupted$Probe_ID <- stringr::str_replace(string = input_RefSeqMRNA_corrupted$Probe_ID, pattern = ' ', replacement = '_')

left_to_do_ <- temp[[2]]
temp[[3]]
####### PREPARE DATA FOR REFSEQMRNA CORRUPTED ANALYSIS ####### 
####### PREPARE DATA FOR PROPER AND CORRUPTED REFSEQMRNA ANALYSIS ####### 
input_RefSeqMRNA <- rbind(input_RefSeqMRNA, input_RefSeqMRNA_corrupted)

input_RefSeqMRNA <- input_RefSeqMRNA[order(input_RefSeqMRNA$Paper),] 

rm(input_RefSeqMRNA_corrupted)
####################################################
####### PREPARE DATA FOR REFSEQMRNA ANALYSIS ####### 
####################################################
left_to_do_8 <- left_to_do_


########################################################################
####### PREPARE DATA FOR PROPER AND CORRUPTED MGI GENES ANALYSIS ####### 
########################################################################
# Study leftover ids: Rik - MGI genes with no canonical name yet
temp <- prepare_input(regex_ = '(.*)(\\d*)(Rik)', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '(.*)__', vec_regex_additional_pattern_removal = ',(.*)')

input_mgi_symbol <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
####### PREPARE DATA FOR MGI GENES ANALYSIS ####### 
####### PREPARE DATA FOR CORRUPTED MGI GENES ANALYSIS ####### 
# Study leftover ids: Rik - MGI genes with no canonical name yet corrupted RIKs
temp <- prepare_input(regex_ = '(.*)(\\d*)(RIK)', df_ = left_to_do_, col_get_resulting_identifer_from_here = 'Gene_symbol')

input_mgi_symbol_corrupted <- temp[[1]] 

input_mgi_symbol_corrupted$Probe_ID <- stringr::str_replace(string = input_mgi_symbol_corrupted$Probe_ID, pattern = 'RIK', replacement = 'Rik')

left_to_do_ <- temp[[2]]
temp[[3]]

####### PREPARE DATA FOR CORRUPTED MGI GENES ANALYSIS ####### 
####### PREPARE DATA FOR PROPER AND CORRUPTED MGI GENES ANALYSIS ####### 
input_mgi_symbol <- rbind(input_mgi_symbol, input_mgi_symbol_corrupted)

input_mgi_symbol <- input_mgi_symbol[order(input_mgi_symbol$Paper),] 

rm(input_mgi_symbol_corrupted)
########################################################################
####### PREPARE DATA FOR PROPER AND CORRUPTED MGI GENES ANALYSIS ####### 
########################################################################



###################################################
####### PREPARE DATA FOR ACCESSION ANALYSIS #######
###################################################
# Study leftover ids: XM_ - most of these names are substituted by NM_ sequences, but NCBI search does not return this new NM_ gene. Stupid.

# Study leftover ids: [letter][numbers] - a) ncbi accession nb. An accession number applies to the complete record and is usually a combination of a letter(s) and numbers, such as a single letter followed by five digits (e.g., U12345) or two letters followed by six digits (e.g., AF123456). Some accessions might be longer, depending on the type of sequence record. Accession numbers do not change, even if information in the record is changed at the author's request. Sometimes, however, an original accession number might become secondary to a newer accession number, if the authors make a new submission that combines previous sequences, or if for some reason a new submission supercedes an earlier record. These IDs are actually the same in ENA ('embl' or 'clone_based_ensembl_gene') and in ncbi nucleotide
temp <- prepare_input(regex_ = '__[A-Z]{1,2}\\d{5,}', df_ = left_to_do_, regex_two_for_extracting_identifiers_from_string = '__')

input_accession <- temp[[1]] 
left_to_do_ <- temp[[2]]
temp[[3]]
###################################################
####### PREPARE DATA FOR ACCESSION ANALYSIS #######
###################################################



##################################################################
####### PREPARE DATA FOR GOOD AND BAD GENE_SYMBOL ANALYSIS #######
##################################################################
inputAnalysis_include_in_bad_gene_symbol <- left_to_do_$Gene_symbol %in% c('–', 'N/A', '---')

input_bad_gene_symbol <- subset(x = left_to_do_, subset = inputAnalysis_include_in_bad_gene_symbol | is.na(left_to_do_15$Gene_symbol))

input_gene_symbol <- subset(x = left_to_do_, subset = !inputAnalysis_include_in_bad_gene_symbol)
input_gene_symbol <- subset(x = input_gene_symbol, subset = !is.na(input_gene_symbol$Gene_symbol))

# 5S_rRNA.116 , 7SK.11 , A_55_P1955488, Cd247,Gm16565, , snoU2_19.1, U6.137, Y_RNA
#if include RPG in string, then remove all after _ : RGD1564074_PREDICTED
# deal_with_wierd_gene_names <- function()

check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_gene_symbol', df_original = left_to_do_, list_df_splited = list(input_bad_gene_symbol, input_gene_symbol))

input_gene_symbol$Probe_ID <- input_gene_symbol$Gene_symbol

input_gene_symbol$Probe_ID <- stringr::str_remove(string = input_gene_symbol$Probe_ID, pattern = ',(.*)') # old regex = ',(.*)|_(.*)'

# for (regex_nb in seq_along(vec_regex_additional_pattern_removal)) {
#   df_output[[col_put_resulting_identifiers_here]] <-
#     stringr::str_remove(string = df_output[[col_put_resulting_identifiers_here]], pattern = vec_regex_additional_pattern_removal[[regex_nb]])


input_gene_symbol <- input_gene_symbol[order(input_gene_symbol$Paper),] 
input_bad_gene_symbol <- input_bad_gene_symbol[order(input_bad_gene_symbol$Paper),]
##################################################################
####### PREPARE DATA FOR GOOD AND BAD GENE_SYMBOL ANALYSIS #######
##################################################################



########################################
####### CHECK AND WRITE ALL DATA #######
########################################
# Here I show that dataset composed of all the specified subsets is equal to raw_dataset
rebuild_dataset_list <- list(input_ProbeID, input_Gemma, input_EnsemblGeneId, input_EnsemblTranscriptId, input_EntrezGeneId, input_RefSeqMRNA, input_mgi_symbol, input_accession, input_RGD_GM_RNA_7SK, input_bad_gene_symbol, input_gene_symbol)

rebuild_dataset <- rlist::list.rbind(rebuild_dataset_list)

inputAnalysis_lengths_of_experiments_in_rebuild_dataset <- split_and_measure_length(df_ = rebuild_dataset, split_by_str = 'Experiment')

inputAnalysis_lengths_of_experiments_in_rebuild_dataset == inputAnalysis_lengths_of_experiments_in_raw_dataset

# Print lists. They should be good and we need not repeat this script.
rebuild_dataset_list_names <- list('input_ProbeID', 'input_Gemma', 'input_EnsemblGeneId', 'input_EnsemblTranscriptId', 'input_EntrezGeneId', 'input_RefSeqMRNA', 'input_mgi_symbol', 'input_accession', 'input_RGD_GM_RNA_7SK', 'input_bad_gene_symbol', 'input_gene_symbol')
write_inputs(lists_ = rebuild_dataset_list, list_names_list_str = rebuild_dataset_list_names)

rm(list = ls(pattern = 'left_to(.*)|inputAnalysis_(.*)|pattern_for(.*)|Papers_to(.*)|rebuild_dataset_list(.*)'))
########################################
####### CHECK AND WRITE ALL DATA #######
########################################



##################################################
####### CHECK AND WRITE ALL LEFTOOVER DATA #######
##################################################
# Here I show that dataset composed of all the specified subsets is equal to raw_dataset
rebuild_dataset_list <- list(input_EnsemblGeneId, input_EnsemblTranscriptId, input_EntrezGeneId, input_RefSeqMRNA, input_mgi_symbol, input_accession, input_RGD_GM_RNA_7SK, input_bad_gene_symbol, input_gene_symbol)

rebuild_dataset <- rlist::list.rbind(rebuild_dataset_list)



# Print lists. They should be good and we need not repeat this script.
rebuild_dataset_list_names <- list('input_EnsemblGeneId', 'input_EnsemblTranscriptId', 'input_EntrezGeneId', 'input_RefSeqMRNA', 'input_mgi_symbol', 'input_accession', 'input_RGD_GM_RNA_7SK', 'input_bad_gene_symbol', 'input_gene_symbol')
write_inputs(lists_ = rebuild_dataset_list, list_names_list_str = rebuild_dataset_list_names)

rm(list = ls(pattern = 'left_to(.*)|inputAnalysis_(.*)|pattern_for(.*)|Papers_to(.*)|rebuild_dataset_list(.*)'))
##################################################
####### CHECK AND WRITE ALL LEFTOOVER DATA #######
##################################################

# analysis workflows:
# input_ProbeID - by paper identifier
# input_Gemma - gemma annotation -> ncbi annotation -> 'entrezgene_id'
# input_EnsemblGeneId - 'ensembl_gene_id'
# input_EnsemblTranscriptId - 'ensembl_transcript_id'
# input_EntrezGeneId - 'entrezgene_id'
# input_RefSeqMRNA - 'refseq_mrna'; check the .1 notation at the end of identifiers
# input_mgi_symbol - 'mgi_symbol'
# input_accession - 'embl'
# input_RGD_GM_RNA_7SK - ncbi annotation -> 'entrezgene_id'
# input_bad_gene_symbol - manual
# input_leftover_gene_symbol - ncbi annotation -> 'entrezgene_id'

















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
