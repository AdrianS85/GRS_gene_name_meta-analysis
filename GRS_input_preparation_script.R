setwd('/media/adrians/USB DISK1/Projekty/GRS - GJt Review Stress/FULL_DATASET')
source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
rm(list = ls(pattern = '(.*)(temp)|(test)(.*)'))

# To see is na is introduced
length( subset(x = input_EnsemblGeneId$Probe_ID, subset = is.na(input_EnsemblGeneId$Probe_ID)) )



####### PREPARE RAW DATA FOR SUBSETING ####### 
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
####### PREPARE RAW DATA FOR SUBSETING ####### 



####### PREPARE DATA FOR MICROARRAY-CENTERED ANALYSIS #######     
# Subset papers to be analyzed via actual microarray Probe_ID
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

input_ProbeID <- rbind(temp_Probe_ID_old, temp_GEO_ID)
input_ProbeID <- input_ProbeID[order(input_ProbeID$Paper),] 

rm(inputAnalysis_include_in_ProbeID)
####### PREPARE DATA FOR MICROARRAY-CENTERED ANALYSIS ####### 



####### PREPARE DATA FOR GEMMA MICROARRAY-CENTERED ANALYSIS ####### 
Papers_to_analyze_via_Gemma <- c(31, 33)

inputAnalysis_include_in_Gemma <- left_to_do$Paper %in% Papers_to_analyze_via_Gemma

input_Gemma <- subset(x = left_to_do, subset = inputAnalysis_include_in_Gemma)

left_to_do_2 <- subset(x = left_to_do, subset = !inputAnalysis_include_in_Gemma)

inputAnalysis_was_spliting_good_2 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_Gemma', df_original = left_to_do, list_df_splited = list(input_Gemma, left_to_do_2))

rm(inputAnalysis_include_in_Gemma)

input_Gemma$Probe_ID <- input_Gemma$Probe_ID_old
####### PREPARE DATA FOR GEMMA MICROARRAY-CENTERED ANALYSIS ####### 



####### PREPARE DATA FOR ENSEMBL GENE ANALYSIS ####### 
# This detects all entries in ensembl_gene_id from data_v4, Ive checked. Yet, the lenght of this input is much lower than of data_v4_ensembl_gene_id.tsv?
pattern_for_EnsemblGeneID = '(ENS)(.*)G0(.*)'

inputAnalysis_include_in_EnsemblGeneID <-
  stringr::str_detect(string = left_to_do_2$everything, pattern = pattern_for_EnsemblGeneID)

input_EnsemblGeneId <- subset(x = left_to_do_2, subset = inputAnalysis_include_in_EnsemblGeneID)

left_to_do_3 <- subset(x = left_to_do_2, subset = !inputAnalysis_include_in_EnsemblGeneID)

inputAnalysis_was_spliting_good_3 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_EnsemblGeneId', df_original = left_to_do_2, list_df_splited = list(input_EnsemblGeneId, left_to_do_3))

rm(inputAnalysis_include_in_EnsemblGeneID)

input_EnsemblGeneId$Probe_ID <- extract_from_string(chr_vec = input_EnsemblGeneId$everything, regex_one = pattern_for_EnsemblGeneID, regex_two = '__.*')
####### PREPARE DATA FOR ENSEMBL GENE ANALYSIS ####### 



####### PREPARE DATA FOR ENSEMBL TRANSCRIPT ANALYSIS ####### 
pattern_for_EnsemblTranscriptID = '(ENS)(.*)(T0)(.*)'

inputAnalysis_include_in_EnsemblTranscriptID <-
  stringr::str_detect(string = left_to_do_3$everything, pattern = pattern_for_EnsemblTranscriptID)

input_EnsemblTranscriptId <- subset(x = left_to_do_3, subset = inputAnalysis_include_in_EnsemblTranscriptID)

left_to_do_4 <- subset(x = left_to_do_3, subset = !inputAnalysis_include_in_EnsemblTranscriptID)

inputAnalysis_was_spliting_good_4 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_EnsemblTranscriptId', df_original = left_to_do_3, list_df_splited = list(input_EnsemblTranscriptId, left_to_do_4))

rm(inputAnalysis_include_in_EnsemblTranscriptID)

input_EnsemblTranscriptId$Probe_ID <- extract_from_string(chr_vec = input_EnsemblTranscriptId$everything, regex_one = pattern_for_EnsemblTranscriptID, regex_two = '__.*')
####### PREPARE DATA FOR ENSEMBL TRANSCRIPT ANALYSIS ####### 



####### PREPARE DATA FOR ENTREZ GENE ANALYSIS ####### 
pattern_for_EntrezGeneID = '[1-9](\\d*)'

inputAnalysis_include_in_EntrezGeneID <- stringr::str_detect(string = left_to_do_4$Gene_ID, pattern = pattern_for_EntrezGeneID)

input_EntrezGeneId <- subset(x = left_to_do_4, subset = inputAnalysis_include_in_EntrezGeneID)

left_to_do_5 <- subset(x = left_to_do_4, subset = !inputAnalysis_include_in_EntrezGeneID | is.na(inputAnalysis_include_in_EntrezGeneID))

inputAnalysis_was_spliting_good_5 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_EntrezGeneId', df_original = left_to_do_4, list_df_splited = list(input_EntrezGeneId, left_to_do_5))

rm(inputAnalysis_include_in_EntrezGeneID)

# Here we remove paper 48 which seems to have corrupted entrez ids
input_EntrezGeneId_corrupted <- subset(x = input_EntrezGeneId, subset = input_EntrezGeneId$Paper == 48)

# Here we prepare input file
temp_input_EntrezGeneId <- subset(x = input_EntrezGeneId, subset = input_EntrezGeneId$Paper %in% c(13, 16, 28, 34, 37, 40, 47, 53, 63, 65))

inputAnalysis_was_spliting_good_entrez <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_EntrezGeneId', df_original = input_EntrezGeneId, list_df_splited = list(input_EntrezGeneId_corrupted, temp_input_EntrezGeneId))

input_EntrezGeneId <- temp_input_EntrezGeneId

input_EntrezGeneId$Probe_ID <- input_EntrezGeneId$Gene_ID
####### PREPARE DATA FOR ENTREZ GENE ANALYSIS ####### 
####### PREPARE DATA FOR ENTREZ-LOC GENE ANALYSIS ####### 
# Subset entries to be analyzed via Entrez Gene from LOC numbers - LOC - genes of uncertain function. When a published symbol is not available, and orthologs have not yet been determined, Gene will provide a symbol that is constructed as 'LOC' + the GeneID. Therefore LOC is basically GeneID, and is thus unique
pattern_for_EntrezGeneID_LOC = '(LOC)(.*)'
inputAnalysis_include_in_EntrezGeneID_LOC <- stringr::str_detect(string = left_to_do_5$everything, pattern = pattern_for_EntrezGeneID_LOC)

input_EntrezGeneID_LOC <- subset(x = left_to_do_5, subset = inputAnalysis_include_in_EntrezGeneID_LOC)

left_to_do_6 <- subset(x = left_to_do_5, subset = !inputAnalysis_include_in_EntrezGeneID_LOC)

inputAnalysis_was_spliting_good_6 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_EntrezGeneID_LOC', df_original = left_to_do_5, list_df_splited = list(input_EntrezGeneID_LOC, left_to_do_6))

rm(inputAnalysis_include_in_EntrezGeneID_LOC)

input_EntrezGeneID_LOC$Probe_ID <- extract_from_string(chr_vec = input_EntrezGeneID_LOC$everything, regex_one = pattern_for_EntrezGeneID_LOC, regex_two = '__.*')
####### PREPARE DATA FOR ENTREZ-LOC GENE ANALYSIS ####### 
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 
# Beware! input_EntrezGeneId_corrupted will not be included here, because it cannot be fixed!
input_EntrezGeneId <- rbind(input_EntrezGeneId, input_EntrezGeneID_LOC)

rm(input_EntrezGeneID_LOC)
####### PREPARE DATA FOR STANDARD AND LOC-DERIVED ENTREZ ANALYSIS ####### 



####### PREPARE DATA FOR REFSEQMRNA ANALYSIS ####### 
pattern_for_RefSeqMRNA = '(NM_)(\\d*)(.*)'
inputAnalysis_include_in_RefSeqMRNA <- stringr::str_detect(string = left_to_do_6$everything, pattern = pattern_for_RefSeqMRNA)

input_RefSeqMRNA <- subset(x = left_to_do_6, subset = inputAnalysis_include_in_RefSeqMRNA)

left_to_do_7 <- subset(x = left_to_do_6, subset = !inputAnalysis_include_in_RefSeqMRNA)

inputAnalysis_was_spliting_good_7 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_RefSeqMRNA', df_original = left_to_do_6, list_df_splited = list(input_RefSeqMRNA, left_to_do_7))

rm(inputAnalysis_include_in_RefSeqMRNA)

input_RefSeqMRNA$Probe_ID <- extract_from_string(chr_vec = input_RefSeqMRNA$everything, regex_one = pattern_for_RefSeqMRNA, regex_two = '__.*')
####### PREPARE DATA FOR REFSEQMRNA ANALYSIS ####### 
####### PREPARE DATA FOR REFSEQMRNA CORRUPTED ANALYSIS ####### 
# Subset corrupted (no _ after NM) entries to be analyzed via RefSeqMRNA
inputAnalysis_include_in_RefSeqMRNA_corrupted <- stringr::str_detect(string = left_to_do_7$everything, pattern = '(NM )(\\d*)')

input_RefSeqMRNA_corrupted <- subset(x = left_to_do_7, subset = inputAnalysis_include_in_RefSeqMRNA_corrupted)

left_to_do_8 <- subset(x = left_to_do_7, subset = !inputAnalysis_include_in_RefSeqMRNA_corrupted)

inputAnalysis_was_spliting_good_8 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_RefSeqMRNA_corrupted', df_original = left_to_do_7, list_df_splited = list(input_RefSeqMRNA_corrupted, left_to_do_8))

rm(inputAnalysis_include_in_RefSeqMRNA_corrupted)

input_RefSeqMRNA_corrupted$GenBank_Accession <- stringr::str_replace(string = input_RefSeqMRNA_corrupted$GenBank_Accession, pattern = ' ', replacement = '_')

input_RefSeqMRNA_corrupted$Probe_ID <- input_RefSeqMRNA_corrupted$GenBank_Accession
####### PREPARE DATA FOR REFSEQMRNA CORRUPTED ANALYSIS ####### 
####### PREPARE DATA FOR PROPER AND CORRUPTED REFSEQMRNA ANALYSIS ####### 
input_RefSeqMRNA <- rbind(input_RefSeqMRNA, input_RefSeqMRNA_corrupted)

rm(input_RefSeqMRNA_corrupted)
####### PREPARE DATA FOR PROPER AND CORRUPTED REFSEQMRNA ANALYSIS ####### 



####### PREPARE DATA FOR MGI GENES ANALYSIS ####### 
# Study leftover ids: Rik - MGI genes with no canonical name yet
pattern_for_mgi_symbol = '(.*)(\\d*)(Rik)'
inputAnalysis_include_in_mgi_symbol <- stringr::str_detect(string = left_to_do_8$everything, pattern = pattern_for_mgi_symbol)

input_mgi_symbol <- subset(x = left_to_do_8, subset = inputAnalysis_include_in_mgi_symbol)

left_to_do_9 <- subset(x = left_to_do_8, subset = !inputAnalysis_include_in_mgi_symbol)

inputAnalysis_was_spliting_good_9 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_mgi_symbol', df_original = left_to_do_8, list_df_splited = list(input_mgi_symbol, left_to_do_9))

rm(inputAnalysis_include_in_mgi_symbol)

input_mgi_symbol$Probe_ID <- extract_from_string(chr_vec = input_mgi_symbol$everything, regex_one = pattern_for_mgi_symbol, regex_two = '(.*)__')
####### PREPARE DATA FOR MGI GENES ANALYSIS ####### 
####### PREPARE DATA FOR CORRUPTED MGI GENES ANALYSIS ####### 
# Study leftover ids: Rik - MGI genes with no canonical name yet corrupted RIKs
inputAnalysis_include_in_mgi_symbol_corrupted <- stringr::str_detect(string = left_to_do_9$everything, pattern = '(\\d*)(RIK)')

input_mgi_symbol_corrupted <- subset(x = left_to_do_9, subset = inputAnalysis_include_in_mgi_symbol_corrupted)

left_to_do_10 <- subset(x = left_to_do_9, subset = !inputAnalysis_include_in_mgi_symbol_corrupted)

inputAnalysis_was_spliting_good_10 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_mgi_symbol_corrupted', df_original = left_to_do_9, list_df_splited = list(input_mgi_symbol_corrupted, left_to_do_10))

rm(inputAnalysis_include_in_mgi_symbol_corrupted)

input_mgi_symbol_corrupted$Gene_symbol <- stringr::str_replace(string = input_mgi_symbol_corrupted$Gene_symbol, pattern = 'RIK', replacement = 'Rik')

input_mgi_symbol_corrupted$Probe_ID <- input_mgi_symbol_corrupted$Gene_symbol
####### PREPARE DATA FOR CORRUPTED MGI GENES ANALYSIS ####### 
####### PREPARE DATA FOR PROPER AND CORRUPTED MGI GENES ANALYSIS ####### 
input_mgi_symbol <- rbind(input_mgi_symbol, input_mgi_symbol_corrupted)

rm(input_mgi_symbol_corrupted)
####### PREPARE DATA FOR PROPER AND CORRUPTED MGI GENES ANALYSIS ####### 



####### PREPARE DATA FOR ACCESSION ANALYSIS ####### 
# Study leftover ids: XM_ - most of these names are substituted by NM_ sequences, but NCBI search does not return this new NM_ gene. Stupid.

# Study leftover ids: [letter][numbers] - a) ncbi accession nb. An accession number applies to the complete record and is usually a combination of a letter(s) and numbers, such as a single letter followed by five digits (e.g., U12345) or two letters followed by six digits (e.g., AF123456). Some accessions might be longer, depending on the type of sequence record. Accession numbers do not change, even if information in the record is changed at the author's request. Sometimes, however, an original accession number might become secondary to a newer accession number, if the authors make a new submission that combines previous sequences, or if for some reason a new submission supercedes an earlier record. These IDs are actually the same in ENA ('embl' or 'clone_based_ensembl_gene') and in ncbi nucleotide
pattern_for_accession = '__[A-Z]{1,2}\\d{5,}'
inputAnalysis_include_in_accession <- stringr::str_detect(string = left_to_do_10$everything, pattern = pattern_for_accession)

input_accession <- subset(x = left_to_do_10, subset = inputAnalysis_include_in_accession)

left_to_do_11 <- subset(x = left_to_do_10, subset = !inputAnalysis_include_in_accession)

inputAnalysis_was_spliting_good_11 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_accession', df_original = left_to_do_10, list_df_splited = list(input_accession, left_to_do_11))

rm(inputAnalysis_include_in_accession)

input_accession$Probe_ID <- extract_from_string(chr_vec = input_accession$everything, regex_one = pattern_for_accession, regex_two = '__(.*)')
####### PREPARE DATA FOR ACCESSION ANALYSIS ####### 



####### PREPARE DATA FOR RGD - RAT GENOME DATABASE ANALYSIS ####### 
pattern_for_rgd = '(RGD)(.*)'
inputAnalysis_include_in_rgd <- stringr::str_detect(string = left_to_do_11$Gene_symbol, pattern = pattern_for_rgd)

input_rgd <- subset(x = left_to_do_11, subset = inputAnalysis_include_in_rgd)

left_to_do_12 <- subset(x = left_to_do_11, subset = !inputAnalysis_include_in_rgd | is.na(inputAnalysis_include_in_rgd))

inputAnalysis_was_spliting_good_12 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_rgd', df_original = left_to_do_11, list_df_splited = list(input_rgd, left_to_do_12))

rm(inputAnalysis_include_in_rgd)

input_rgd$Probe_ID <- extract_from_string(chr_vec = input_rgd$everything, regex_one = pattern_for_rgd, regex_two = '__(.*)') ### !!! Some work here
####### PREPARE DATA FOR RGD - RAT GENOME DATABASE ANALYSIS ####### 



####### PREPARE DATA FOR GM ANALYSIS ####### 
# Gm - annotated genes that do not have a canonical name (yet)
pattern_for_gm = '(Gm)(.*)'
inputAnalysis_include_in_gm <- stringr::str_detect(string = left_to_do_12$Gene_symbol, pattern = pattern_for_gm)

input_gm <- subset(x = left_to_do_12, subset = inputAnalysis_include_in_gm)

left_to_do_13 <- subset(x = left_to_do_12, subset = !inputAnalysis_include_in_gm | is.na(inputAnalysis_include_in_gm))

inputAnalysis_was_spliting_good_13 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_gm', df_original = left_to_do_12, list_df_splited = list(input_gm, left_to_do_13))

rm(inputAnalysis_include_in_gm)

input_gm$Probe_ID <- extract_from_string(chr_vec = input_gm$everything, regex_one = pattern_for_gm, regex_two = '__(.*)')

input_gm$Probe_ID <- stringr::str_remove(string = input_gm$Probe_ID, pattern = ',(.*)')
####### PREPARE DATA FOR GM ANALYSIS ####### 


####### PREPARE DATA FOR RNA ANALYSIS ####### 
inputAnalysis_include_in_RNA <- stringr::str_detect(string = left_to_do_13$Gene_symbol, pattern = '(RNA)(.*)')

input_RNA <- subset(x = left_to_do_13, subset = inputAnalysis_include_in_RNA)

left_to_do_14 <- subset(x = left_to_do_13, subset = !inputAnalysis_include_in_RNA | is.na(inputAnalysis_include_in_RNA))

inputAnalysis_was_spliting_good_14 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_5S_', df_original = left_to_do_13, list_df_splited = list(input_RNA, left_to_do_14))

rm(inputAnalysis_include_in_RNA)

input_RNA$Probe_ID <- input_RNA$Gene_symbol
####### PREPARE DATA FOR RNA ANALYSIS ####### 



####### PREPARE DATA FOR 7SK ANALYSIS ####### 
inputAnalysis_include_in_7SK <- stringr::str_detect(string = left_to_do_14$Gene_symbol, pattern = '(7SK)')

input_7SK <- subset(x = left_to_do_14, subset = inputAnalysis_include_in_7SK)

left_to_do_15 <- subset(x = left_to_do_14, subset = !inputAnalysis_include_in_7SK | is.na(inputAnalysis_include_in_7SK))

inputAnalysis_was_spliting_good_15 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_7SK', df_original = left_to_do_14, list_df_splited = list(input_7SK, left_to_do_15))

rm(inputAnalysis_include_in_7SK)

input_7SK$Probe_ID <- input_7SK$Gene_symbol
####### PREPARE DATA FOR 7SK ANALYSIS ####### 



####### PREPARE DATA FOR GOOD AND BAD GENE_SYMBOL ANALYSIS ####### 
inputAnalysis_include_in_bad_gene_symbol <- left_to_do_15$Gene_symbol %in% c('–', 'N/A')

input_bad_gene_symbol <- subset(x = left_to_do_15, subset = inputAnalysis_include_in_bad_gene_symbol | is.na(left_to_do_15$Gene_symbol))

input_leftover_gene_symbol <- subset(x = left_to_do_15, subset = !inputAnalysis_include_in_bad_gene_symbol)
input_leftover_gene_symbol <- subset(x = input_leftover_gene_symbol, subset = !is.na(input_leftover_gene_symbol$Gene_symbol))

inputAnalysis_was_spliting_good_16 <- check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'input_7SK', df_original = left_to_do_15, list_df_splited = list(input_bad_gene_symbol, input_leftover_gene_symbol))

input_leftover_gene_symbol$Probe_ID <- input_leftover_gene_symbol$Gene_symbol
####### PREPARE DATA FOR BAD GENE_SYMBOL ANALYSIS ####### 



####### PREPARE DATA FOR GENE_SYMBOL ANALYSIS ####### 
# Here I show that dataset composed of all the specified subsets is equal to raw_dataset
rebuild_dataset <- list(input_ProbeID, input_Gemma, input_EnsemblGeneId, input_EnsemblTranscriptId, input_EntrezGeneId, input_EntrezGeneId_corrupted, input_RefSeqMRNA, input_mgi_symbol, input_accession, input_rgd, input_gm, input_RNA, input_7SK, input_bad_gene_symbol, input_leftover_gene_symbol)

rebuild_dataset <- rlist::list.rbind(rebuild_dataset)

inputAnalysis_lengths_of_experiments_in_rebuild_dataset <- split_and_measure_length(df_ = rebuild_dataset, split_by_str = 'Experiment')

inputAnalysis_lengths_of_experiments_in_rebuild_dataset == inputAnalysis_lengths_of_experiments_in_raw_dataset
####### PREPARE DATA FOR GENE_SYMBOL ANALYSIS ####### 




