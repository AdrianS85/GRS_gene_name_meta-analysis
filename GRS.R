# required packages:: rentrez, biomaRt, tidyverse (dplyr, purrr, readr), rlist, stringr, data.table, getPass
# rm(list = ls(pattern = '(.*)(temp)|(test)(.*)'))
# rm(list = ls(pattern = 'rebu(.*)'))
# devtools::install_github("hadley/lineprof")



######################################
######### PREPARING METADATA ######### 
######################################
setwd('/media/adrians/USB DISK1/Projekty/GRS - GJt Review Stress/FULL_DATASET')
# setwd(getwd())
source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
  
### Read in the table with description of experiments. They need to have the same identifiers as data in PRE_DATA, that is: Paper(int), Experiment(chr). Also needs Species(chr) column, if probe_id identification function or ncbi query tool is to be used
descriptions <- readr::read_tsv("descriptions.txt", col_types = "nccccccccccccccccccccccc")
######################################
######### PREPARING METADATA ######### 
######################################



########################################################################################
################################# FIRST STAGE OF ANNOTATION ############################
########################################################################################



###########################################################################
##### PROBE ID - GET THE HIGHEST-HIT-RETURNING ID TYPE, THEN ANNOTATE ##### 
###########################################################################
# Single Probe_ID may be assigned to multiple genes
DATA_FROM_HIGHEST_HIT_ANALYSIS <- get_the_highest_hit_returning_id_type(
  str_filename_ = 'data_v4_Probe_ID.tsv', 
  descriptions_ = descriptions,
  int_Probe_IDs_to_test = 200,
  str_experiment_name = 'Probe_ID')


# identifiers_used_for_annotation_global <- as.character(DATA_FROM_HIGHEST_HIT_ANALYSIS[[2]]$platform_to_use)
# OR
identifiers_used_for_annotation_global <- readr::read_tsv('platform_to_use_for_probes_based_analysis.tsv')$platform_to_use


# Annotate and save Probe_ID file
annotations_Probe_ID <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_1_stage/input_ProbeID.tsv', str_identifier_type__ = identifiers_used_for_annotation_global, str_experiment_name = 'Probe_ID', col_types__ = 'ncccccccccccccnccc') # 195495 vs 195495
save(annotations_Probe_ID, file = 'annotations_Probe_ID')
# load(file = 'annotations_Probe_ID')


# Prepare finalized and leftover files
finalized_Probe_ID <- subset(x = annotations_Probe_ID, subset = !is.na(annotations_Probe_ID$external_gene_name))
save(finalized_Probe_ID, file = 'finalized_Probe_ID')
# load(file = 'finalized_Probe_ID')
leftover_Probe_ID <- subset(x = annotations_Probe_ID, subset = is.na(annotations_Probe_ID$external_gene_name))
save(leftover_Probe_ID, file = 'leftover_Probe_ID')
# load(file = 'leftover_Probe_ID')
check_was_the_spliting_of_df_by_filtering_ok(df_original = annotations_Probe_ID, list_df_splited = list(finalized_Probe_ID, leftover_Probe_ID))
###########################################################################
##### PROBE ID - GET THE HIGHEST-HIT-RETURNING ID TYPE, THEN ANNOTATE ##### 
###########################################################################



#################################
####### KNOWN IDENTIFIERS ####### 
#################################
####### FIRST ORDER ANNOTATION ####### 
# These identifiers return unique gene name per identifier, because every id is assigned to known gene. This is different for Probe_IDs
annotations_ensembl_gene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_1_stage/input_EnsemblGeneId.tsv', str_identifier_type__ = 'ensembl_gene_id', col_types__ = 'ncccccccccccccnccc') #15302 vs 15302
save(annotations_ensembl_gene_id, file = 'annotations_ensembl_gene_id')
# load(file = 'annotations_ensembl_gene_id')


annotations_ensembl_transcript_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_1_stage/input_EnsemblTranscriptId.tsv', str_identifier_type__ = 'ensembl_transcript_id', col_types__ = 'ncccccccccccccnccc') #2419 vs 2419
save(annotations_ensembl_transcript_id, file = 'annotations_ensembl_transcript_id')
# load(file = 'annotations_ensembl_transcript_id')


annotations_entrezgene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_1_stage/input_EntrezGeneId.tsv', str_identifier_type__ = 'entrezgene_id', col_types__ = 'ncccccccccccccnccc') #29459 vs 29459
save(annotations_entrezgene_id, file = 'annotations_entrezgene_id')
# load('annotations_entrezgene_id')


annotations_refseq_mrna <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_1_stage/input_RefSeqMRNA.tsv', str_identifier_type__ = 'refseq_mrna', col_types__ = 'ncccccccccccccnccc') #4715 vs 4715
save(annotations_refseq_mrna, file = 'annotations_refseq_mrna')
# load('annotations_refseq_mrna')


annotations_accession <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_1_stage/input_accession.tsv', str_identifier_type__ = 'embl', col_types__ = 'ncccccccccccccnccc') #3412 vs 3412
save(annotations_accession, file = 'annotations_accession')
# load('annotations_accession')


####### FIRST ORDER ANNOTATION ####### 
####### GET FINALIZED AND LEFTOVER DFS ####### 
if( min(colnames(annotations_ensembl_gene_id) == colnames(annotations_ensembl_transcript_id)) == 1 & min(colnames(annotations_ensembl_gene_id) == colnames(annotations_entrezgene_id)) == 1 & min(colnames(annotations_ensembl_gene_id) == colnames(annotations_refseq_mrna)) == 1 & min(colnames(annotations_ensembl_gene_id) == colnames(annotations_accession)) == 1 )
{
  temp_annotations_known_identifiers <- rbind(annotations_ensembl_gene_id, annotations_ensembl_transcript_id, annotations_entrezgene_id, annotations_refseq_mrna, annotations_accession)
  temp_annotations_known_identifiers <- temp_annotations_known_identifiers[order(temp_annotations_known_identifiers$Paper),] 
}

# Prepare finalized and leftover files
finalized_known_identifiers <- subset(x = temp_annotations_known_identifiers, subset = !is.na(temp_annotations_known_identifiers$external_gene_name))
save(finalized_known_identifiers, file = 'finalized_known_identifiers')
# load('finalized_known_identifiers')

leftover_known_identifiers <- subset(x = temp_annotations_known_identifiers, subset = is.na(temp_annotations_known_identifiers$external_gene_name))
save(leftover_known_identifiers, file = 'leftover_known_identifiers')
# load('leftover_known_identifiers')

check_was_the_spliting_of_df_by_filtering_ok(df_original = temp_annotations_known_identifiers, list_df_splited = list(finalized_known_identifiers, leftover_known_identifiers))

rm(temp_annotations_known_identifiers)
####### GET FINALIZED AND LEFTOVER DFS ####### 
#################################
####### KNOWN IDENTIFIERS ####### 
#################################



### !!! I have changed names here
########################################
####### GEMMA ANNOTATION (31/33) ####### 
########################################
# devtools::install_github('PavlidisLab/gemmaAPI.R')
####### GEMMA PRE-PRE-ANNOTATION ####### 
# Prepare metadata for gemma pre-annotation
platforms_ids_to_download <- readr::read_tsv('data_v4_gemma_platforms.tsv')
gemma_platforms <- download_platforms_from_gemma(platforms_ids_to_download$Platform_ID)

# Pre-annotate probe ids with gene names by gemma
pre_pre_annotations_gemmaData <- get_gemma_annotations_for_data(str_filename_ = 'checked_input_1_stage/input_Gemma.tsv', list_gemma_platforms = gemma_platforms, col_types__ = 'ncccccccccccccnccc') # 2002 vs 2002
# save(pre_pre_annotations_gemmaData, file = 'pre_pre_annotations_gemmaData')
# load(pre_pre_annotations_gemmaData)
# dir.create('gemma_pre_annotation')
# readr::write_tsv(x = pre_annotations_gemma_probeID, path = 'gemma_pre_annotation/raw_gemma_annotations.tsv')


####### GEMMA PRE-ANNOTATION ####### 
#These gene names need to be compared to ncbi gene names, so we prepare data for this:
pre_pre_annotations_gemmaData$Probe_ID <- pre_pre_annotations_gemmaData$external_gene_name
pre_pre_annotations_gemmaData$external_gene_name <- NULL

# Prepare input for annotation and leftover files. We need to select only rows that actually have gene names, and put rest where it belong
pre_annotations_gemmaData_input <- subset(x = pre_pre_annotations_gemmaData, subset = !is.na(pre_pre_annotations_gemmaData$Probe_ID)) ### !!! In this data, for every probeID gene name was produced by gemma database. Now this needs to be passed through ncbi to get entrezID, which then will be changed into final gene names by ensembl

#There are about 50 rows, where ProbeID contains multiple names. I will just leave first one, because its not worth the work
pre_annotations_gemmaData_input$Probe_ID <- stringr::str_remove(string = pre_annotations_gemmaData_input$Probe_ID, pattern = '\\|(.*)')

save(pre_annotations_gemmaData_input, file = 'pre_annotations_gemmaData_input')
# load('pre_annotations_gemmaData_input')

leftover_second_stage_input_gemma_probeID <- subset(x = pre_pre_annotations_gemmaData, subset = is.na(pre_pre_annotations_gemmaData$Probe_ID)) ### !!! This is leftovers
save(leftover_second_stage_input_gemma_probeID, file = 'leftover_second_stage_input_gemma_probeID')
# load('leftover_second_stage_input_gemma_probeID')

check_was_the_spliting_of_df_by_filtering_ok(df_original = pre_pre_annotations_gemmaData, list_df_splited = list(pre_annotations_gemmaData_input, leftover_second_stage_input_gemma_probeID))

rm(pre_pre_annotations_gemmaData, platforms_ids_to_download, gemma_platforms)

# The annotated table is further processed along with other ncbi-preannotated data tables
# pre_annotations_gemmaData <- annotate_identifiers_to_geneID_wrapper_for_using_dfs_as_input(descriptions__ = descriptions, df_to_annotate = pre_annotations_gemmaData_input, str_experiment_name__ = 'pre_annotations_gemmaData_input') # 1571 vs 1571 This file goes to ncbi pre-annotation, where it will be annotated along with other annotated in NCBI
# save(pre_annotations_gemmaData, file = 'pre_annotations_gemmaData')
########################################
####### GEMMA ANNOTATION (31/33) ####### 
########################################
#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################
### !!! add to function make_and_write_table_with_original_and_ncbi_ids chaning ncbi_response to external_gene_name
# Pre-annotate and save gene names
# Prepare df_ for preannotating - taking gemma and gene names
# input_gene_symbol <- read_preformated_data(str_filename = 'checked_input_1_stage/input_gene_symbol.tsv', col_types_ = 'ncccccccccccccnccc', int_numbers_are_from = 11, int_numbers_are_to = 13)
if( min(colnames(input_gene_symbol) == colnames(pre_annotations_gemmaData_input)) == 1 )
{
  temp_annotated_names_to_geneIDs <- rbind(input_gene_symbol, pre_annotations_gemmaData_input) # 5207/1571
  temp_annotated_names_to_geneIDs <- temp_annotated_names_to_geneIDs[order(temp_annotated_names_to_geneIDs$Paper),] 
}

          annotated_names_to_geneIDs <- annotate_identifiers_to_geneID_wrapper_for_using_dfs_as_input(df_to_annotate = temp_annotated_names_to_geneIDs, descriptions__ = descriptions, chr_gene_identifier__ = 'Gene Name', str_experiment_name__ = 'gene_names', col_types___ = 'ncccccccccccccnccc') #6778 vs 6778
          save(annotated_names_to_geneIDs, file = 'annotated_names_to_geneIDs')
# load('annotated_names_to_geneIDs')

rm(temp_annotated_names_to_geneIDs)

annotated_names_to_geneIDs <- annotated_names_to_geneIDs[order(annotated_names_to_geneIDs$Paper),]
annotated_names_to_geneIDs$Probe_ID_left_from_first_ncbi_annotation_stage <- annotated_names_to_geneIDs$Probe_ID
annotated_names_to_geneIDs$Probe_ID <- annotated_names_to_geneIDs$external_gene_name
annotated_names_to_geneIDs$Probe_ID_2 <- NULL
annotated_names_to_geneIDs$external_gene_name <- NULL


# Annotate and save gene names. Saimiri genes are not found at this step, cause ensembl generally does not link to ncbi for this species
annotations_entrezgene_id_post_ncbi_annotation <- master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(descriptions_ = descriptions, df_to_annotate = annotated_names_to_geneIDs, str_identifier_type__ = 'entrezgene_id', str_experiment_name__ = 'entrezgene_id_post_ncbi_preAnnotation', col_types___ = 'ncccccccccccccncccc', int_numbers_are_from__ = 11, int_numbers_are_to__ = 13) #6782 vs 6782
save(annotations_entrezgene_id_post_ncbi_annotation, file = 'annotations_entrezgene_id_post_ncbi_annotation')
# load('annotations_entrezgene_id_post_ncbi_annotation')


# Prepare finalized and leftover files
finalized_entrezgene_id_post_ncbi_annotation <- subset(x = annotations_entrezgene_id_post_ncbi_annotation, subset = !is.na(annotations_entrezgene_id_post_ncbi_annotation$external_gene_name))

leftover_annotations_entrezgene_id_post_ncbi_annotation <- subset(x = annotations_entrezgene_id_post_ncbi_annotation, subset = is.na(annotations_entrezgene_id_post_ncbi_annotation$external_gene_name))

check_was_the_spliting_of_df_by_filtering_ok(df_original = annotations_entrezgene_id_post_ncbi_annotation, list_df_splited = list(leftover_annotations_entrezgene_id_post_ncbi_annotation, finalized_entrezgene_id_post_ncbi_annotation))
#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################



#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################
leftover_annotations_entrezgene_id_post_ncbi_annotation$external_gene_name <- NULL
leftover_Probe_ID$external_gene_name <- NULL
leftover_known_identifiers$external_gene_name <- NULL

leftover_second_stage_input_gemma_probeID$Probe_ID_left_from_first_ncbi_annotation_stage <- NA
leftover_Probe_ID$Probe_ID_left_from_first_ncbi_annotation_stage <- NA
leftover_known_identifiers$Probe_ID_left_from_first_ncbi_annotation_stage <- NA

leftover_second_stage_input_gemma_probeID <- leftover_second_stage_input_gemma_probeID %>%
  dplyr::select(Probe_ID, dplyr::everything())

if(min(colnames(leftover_second_stage_input_gemma_probeID) == colnames(leftover_annotations_entrezgene_id_post_ncbi_annotation)) == 1 & min(colnames(leftover_second_stage_input_gemma_probeID) == colnames(leftover_Probe_ID)) == 1 & min(colnames(leftover_second_stage_input_gemma_probeID) == colnames(leftover_known_identifiers)) == 1)
{
  leftovers <-
    rbind(
      leftover_second_stage_input_gemma_probeID,
      leftover_annotations_entrezgene_id_post_ncbi_annotation,
      leftover_Probe_ID,
      leftover_known_identifiers
    )
  
  leftovers <- leftovers[order(leftovers$Paper),] 
}

leftovers$Probe_ID_left_from_first_annotating_stage <- leftovers$Probe_ID
leftovers$Probe_ID <- NA

save(leftovers, file = 'leftovers')
# load('leftovers')

# rm(list = ls(pattern = ''))
#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################



#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################
temp <- finalized_entrezgene_id_post_ncbi_annotation[['Probe_ID_left_from_first_ncbi_annotation_stage']]
finalized_entrezgene_id_post_ncbi_annotation$Probe_ID_left_from_first_ncbi_annotation_stage <- NULL
finalized_entrezgene_id_post_ncbi_annotation$Probe_ID_left_from_first_ncbi_annotation_stage <- temp

finalized_known_identifiers$Probe_ID_left_from_first_ncbi_annotation_stage <- NA
finalized_Probe_ID$Probe_ID_left_from_first_ncbi_annotation_stage <- NA


if(min(colnames(finalized_entrezgene_id_post_ncbi_annotation) == colnames(finalized_known_identifiers)) == 1 & min(colnames(finalized_entrezgene_id_post_ncbi_annotation) == colnames(finalized_Probe_ID)) == 1)
{
  finalized <-
    rbind(
      finalized_entrezgene_id_post_ncbi_annotation,
      finalized_known_identifiers,
      finalized_Probe_ID
    )
  
  finalized <- finalized[order(finalized$Paper),] 
}
save(finalized, file = 'finalized')
# load('finalized')

rm(list = ls(pattern = '(finalized_)|(annotated_)|(annotations_)|(leftover_)(.*)|temp'))
#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################

# Do number of rows in input dataset and output datasets match?
length(leftovers[[1]]) + length(read.csv(file = 'checked_input_1_stage/input_bad_gene_symbol.tsv', header = T, sep = '\t')[[1]]) + length(finalized[[1]]) == length(reformated_raw_dataset[[1]])

########################################################################################
################################# FIRST STAGE OF ANNOTATION ############################
########################################################################################






#########################################################################################
################################# FURTHER STAGE OF ANNOTATION ###########################
#########################################################################################
###########################
####### SET OPTIONS ####### 
###########################
stage_ = 4
if (stage_ == 2){
  opts_col_types___ = 'cnccccccccccccnccccc'
  opts_checked_input_stage = 'checked_input_2_stage/'
  opts_save_name_sufix = '_2'
  opts_probe_ID_left__ = 'Probe_ID_left_from_second_ncbi_annotation_stage'
  opts_leftovers_name = 'Probe_ID_left_from_second_annotating_stage'
} else if (stage_ == 3){
  opts_col_types___ = 'cnccccccccccccnccccccc'
  opts_checked_input_stage = 'checked_input_3_stage/'
  opts_save_name_sufix = '_3'
  opts_probe_ID_left__ = 'Probe_ID_left_from_third_ncbi_annotation_stage'
  opts_leftovers_name = 'Probe_ID_left_from_third_annotating_stage'
} else if (stage_ == 4){
  opts_col_types___ = 'cnccccccccccccnccccccccc'
  opts_checked_input_stage = 'checked_input_4_stage/'
  opts_save_name_sufix = '_4'
  opts_probe_ID_left__ = 'Probe_ID_left_from_fourth_ncbi_annotation_stage'
  opts_leftovers_name = 'Probe_ID_left_from_fourth_annotating_stage'
} #After 5'th stage we are left with 54 entries to be annotated, and none of these entries are later successfully annotated. Thus, it seems reasonable to stop annotating after 4th step
###########################
####### SET OPTIONS ####### 
###########################



#################################
####### KNOWN IDENTIFIERS ####### 
#################################
####### FIRST ORDER ANNOTATION ####### 
# These identifiers return unique gene name per identifier, because every id is assigned to known gene. This is different for Probe_IDs
annotations_ensembl_gene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_EnsemblGeneId.tsv'), str_identifier_type__ = 'ensembl_gene_id', col_types__ = opts_col_types___) #2-12
get_percentages_of_annotation(annotations_ensembl_gene_id)
save(annotations_ensembl_gene_id, file = paste0('annotations_ensembl_gene_id', opts_save_name_sufix))
# load(file = 'annotations_ensembl_gene_id_')


annotations_ensembl_transcript_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_EnsemblTranscriptId.tsv'), str_identifier_type__ = 'ensembl_transcript_id', col_types__ = opts_col_types___) #2-38
get_percentages_of_annotation(annotations_ensembl_transcript_id)
save(annotations_ensembl_transcript_id, file = paste0('annotations_ensembl_transcript_id', opts_save_name_sufix))
# load(file = 'annotations_ensembl_transcript_id_')

  
annotations_entrezgene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_EntrezGeneId.tsv'), str_identifier_type__ = 'entrezgene_id', col_types__ = opts_col_types___) #2-16106 3-94 4-21
get_percentages_of_annotation(annotations_entrezgene_id)
save(annotations_entrezgene_id, file = paste0('annotations_entrezgene_id', opts_save_name_sufix))
# load('annotations_entrezgene_id_')


annotations_refseq_mrna <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_RefSeqMRNA.tsv'), str_identifier_type__ = 'refseq_mrna', col_types__ = opts_col_types___) #2-293 3-144
get_percentages_of_annotation(annotations_refseq_mrna)
save(annotations_refseq_mrna, file = paste0('annotations_refseq_mrna', opts_save_name_sufix))
# load('annotations_refseq_mrna_')


annotations_accession <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_accession.tsv'), str_identifier_type__ = 'embl', col_types__ = opts_col_types___) #2-7453 3-2723 4-144
get_percentages_of_annotation(annotations_accession)
save(annotations_accession, file = paste0('annotations_accession', opts_save_name_sufix))
# load('annotations_accession_')

####### FIRST ORDER ANNOTATION ####### 
####### GET FINALIZED AND LEFTOVER DFS ####### 
# do this for all objects in the same statement using this: https://stackoverflow.com/questions/52937148/write-script-to-ignore-objects-which-can-t-be-found-in-r

set_empty_annotations_to_na()
  
if(are_vectors_the_same(chr_vec_list = list(colnames(annotations_ensembl_gene_id), colnames(annotations_ensembl_transcript_id), colnames(annotations_entrezgene_id), colnames(annotations_refseq_mrna), colnames(annotations_accession))) )
{
  temp_annotations_known_identifiers <-
    rbind(
      annotations_ensembl_gene_id,
      annotations_ensembl_transcript_id,
      annotations_entrezgene_id,
      annotations_refseq_mrna,
      annotations_accession
    )
  
  temp_annotations_known_identifiers <- subset(x = temp_annotations_known_identifiers, subset = !is.na(temp_annotations_known_identifiers$Probe_ID)) 
}

# Prepare finalized and leftover files
finalized_known_identifiers <- subset(x = temp_annotations_known_identifiers, subset = !is.na(temp_annotations_known_identifiers$external_gene_name)) #2-13103 3-70 4-
save(annotated_names_to_geneIDs, file = paste0('finalized_known_identifiers', opts_save_name_sufix))

leftover_known_identifiers <- subset(x = temp_annotations_known_identifiers, subset = is.na(temp_annotations_known_identifiers$external_gene_name)) #2-12707 3-2893 4-
save(annotated_names_to_geneIDs, file = paste0('leftover_known_identifiers', opts_save_name_sufix))

check_was_the_spliting_of_df_by_filtering_ok(df_original = temp_annotations_known_identifiers, list_df_splited = list(finalized_known_identifiers, leftover_known_identifiers))

rm(temp_annotations_known_identifiers)
####### GET FINALIZED AND LEFTOVER DFS ####### 
#################################
####### KNOWN IDENTIFIERS ####### 
#################################



#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################
annotated_names_to_geneIDs <-
  annotate_identifiers_to_geneID(
    str_filename_ = paste0(opts_checked_input_stage, 'input_gene_symbol.tsv'),
    descriptions_ = descriptions,
    chr_gene_identifier_ = 'Gene Name',
    str_experiment_name = 'gene_names_second_stage',
    col_types__ = opts_col_types___
  ) #2-2347 3-1503 4-2458 5-43
save(annotated_names_to_geneIDs, file = paste0('annotated_names_to_geneIDs', opts_save_name_sufix))
# load('annotated_names_to_geneIDs_4')


annotated_names_to_geneIDs <- post_annotate_identifiers_to_geneID(annotated_names_to_geneIDs, Probe_ID_left_ = opts_probe_ID_left__)

# Annotate and save gene names
  annotations_entrezgene_id_post_ncbi_annotation <-
    master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(
      descriptions_ = descriptions,
      df_to_annotate = annotated_names_to_geneIDs,
      str_identifier_type__ = 'entrezgene_id',
      str_experiment_name__ = 'entrezgene_id_post_ncbi_preAnnotation',
      col_types___ = paste0(opts_col_types___, 'c'),
      int_numbers_are_from__ = 11,
      int_numbers_are_to__ = 13
    ) #2-2347 3-1503 4-2458
  get_percentages_of_annotation(annotations_entrezgene_id_post_ncbi_annotation)
  save(annotations_entrezgene_id_post_ncbi_annotation, file = paste0('annotations_entrezgene_id_post_ncbi_annotation', opts_save_name_sufix))
# load('annotations_entrezgene_id_post_ncbi_annotation_')


# Prepare finalized and leftover files
finalized_entrezgene_id_post_ncbi_annotation <- subset(x = annotations_entrezgene_id_post_ncbi_annotation, subset = !is.na(annotations_entrezgene_id_post_ncbi_annotation$external_gene_name))
save(finalized_entrezgene_id_post_ncbi_annotation, file = paste0('finalized_entrezgene_id_post_ncbi_annotation', opts_save_name_sufix)) #3-107
# load('finalized_entrezgene_id_post_ncbi_annotation_')

leftover_annotations_entrezgene_id_post_ncbi_annotation <- subset(x = annotations_entrezgene_id_post_ncbi_annotation, subset = is.na(annotations_entrezgene_id_post_ncbi_annotation$external_gene_name))
save(leftover_annotations_entrezgene_id_post_ncbi_annotation, file = paste0('leftover_annotations_entrezgene_id_post_ncbi_annotation', opts_save_name_sufix)) #3-1396
# load('leftover_annotations_entrezgene_id_post_ncbi_annotation_')

check_was_the_spliting_of_df_by_filtering_ok(df_original = annotations_entrezgene_id_post_ncbi_annotation, list_df_splited = list(leftover_annotations_entrezgene_id_post_ncbi_annotation, finalized_entrezgene_id_post_ncbi_annotation))
#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################



#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################
leftover_annotations_entrezgene_id_post_ncbi_annotation$external_gene_name <- NULL
leftover_known_identifiers$external_gene_name <- NULL
leftover_known_identifiers[opts_probe_ID_left__] <- NA

if( min(colnames(leftover_annotations_entrezgene_id_post_ncbi_annotation) == colnames(leftover_known_identifiers)) == 1)
{
leftovers <-
  rbind(leftover_annotations_entrezgene_id_post_ncbi_annotation,
    leftover_known_identifiers)
}

### BEWARE OF THIS STEP. HERE WE ADD THE input_ files for which no annotation was produced, so we just need to add them to the leftovers
if (stage_ == 4) {
  annotations_accession <- read_preformated_data(str_filename = paste0(opts_checked_input_stage, 'input_accession.tsv'), col_types_ = opts_col_types___, int_numbers_are_from = 11, int_numbers_are_to = 13)
  annotations_accession[opts_probe_ID_left__] <- NA
  
  if( min(colnames(leftovers) == colnames(annotations_accession)) == 1){
    leftovers <- rbind(leftovers, annotations_accession)
  }
} 

leftovers[opts_leftovers_name] <- leftovers$Probe_ID
leftovers$Probe_ID <- NA

check_was_the_spliting_of_df_by_filtering_ok(df_original = leftovers, list_df_splited = list(leftover_annotations_entrezgene_id_post_ncbi_annotation, leftover_known_identifiers))

save(leftovers, file = paste0('leftovers', opts_save_name_sufix)) #2-12635 3-4289 4-
# load('leftovers_')

# rm(list = ls(pattern = ''))
#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################



#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################
temp <- finalized_entrezgene_id_post_ncbi_annotation[[opts_probe_ID_left__]]
finalized_entrezgene_id_post_ncbi_annotation[opts_probe_ID_left__] <- NULL
finalized_entrezgene_id_post_ncbi_annotation[opts_probe_ID_left__] <- temp

finalized_known_identifiers[opts_probe_ID_left__] <- NA

if( min(colnames(finalized_known_identifiers) == colnames(finalized_entrezgene_id_post_ncbi_annotation)) == 1)
{
  finalized <-
    rbind(
      finalized_known_identifiers,
      finalized_entrezgene_id_post_ncbi_annotation
    )
  save(finalized, file = paste0('finalized', opts_save_name_sufix)) #2-14020 3-177 4-
}
# load('finalized_')

check_was_the_spliting_of_df_by_filtering_ok(df_original = finalized, list_df_splited = list(finalized_entrezgene_id_post_ncbi_annotation, finalized_known_identifiers))

rm(list = ls(pattern = '(annotated_)|(annotations_)|(leftover_)(.*)'))
#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################

length(leftovers[[1]])
length(finalized[[1]]) + length(leftovers[[1]]) + length(read.csv(file = paste0(opts_checked_input_stage, 'input_bad_gene_symbol.tsv'), header = T, sep = '\t')[[1]])

# leftovers - 49847, leftovers_2 - 12637, leftovers_3 - , leftovers_4 - 2558



rm(leftovers_for_checking_data_integrity)

#########################################################################################
################################# FURTHER STAGE OF ANNOTATION ###########################
#########################################################################################



# RAW DATA
# entries: 261475

# LOADING STAGE 1 DATA #Check gemma annotation step #Check saimiri
# load('annotations_Probe_ID') #195022 entries #19 variables #72exp-5%annot
# load('annotations_ensembl_gene_id') #15302 #19 #16,24-0% (16,24 - deprecated ids, OK)
# load('annotations_ensembl_transcript_id') #2419 #19
# load('annotations_entrezgene_id') #29487 #19
# load('annotations_refseq_mrna') #5078 #19 #16-0% (16 - absent ids, OK)
# load('annotations_accession') #3450 #19 # low annotation, but ok i think
# load('annotations_entrezgene_id_post_ncbi_annotation') #6822 (5207 genenames+1571 gemma) #20 #18-0% (18 - absent ids, OK)
# # 431 of gemma leftovers are missing here, cause they went straight to the leftovers
# load('finalized') #208610 #20
# load('leftovers') #49401 #20
# input_bad_gene_symbol <- read.csv(file = 'checked_input_1_stage/input_bad_gene_symbol.tsv',  sep = '\t') #3464 #18
# Total input: 258011; Total output: 258011; Total previous step: 261475



# LOADING STAGE 2 DATA
# load('annotations_ensembl_gene_id_2') #12 #21
# load('annotations_ensembl_transcript_id_2') #38 #21
# load('annotations_entrezgene_id_2') #16106 #21 #24-5%,72-4% (24,72 - absent ids, OK)
# load('annotations_refseq_mrna_2') #293 #21 #13-5% (wierd... 'NM_207673 (RefSeq DNA) is used as supporting evidence for transcript ENSMUST00000169373. Likely ok')
# load('annotations_accession_2') #7453 #21
# load('annotations_entrezgene_id_post_ncbi_annotation_2') #2362 #22 #16-3%,33-7%,60-9%,63-0%,67-2%
# load('finalized_2') #13703 #22
# load('leftovers_2') #12561 #22
# input_bad_gene_symbol <- read.csv(file = 'checked_input_2_stage/input_bad_gene_symbol.tsv',  sep = '\t') #23137 #20
# Total input: 26264; Total output: 26264; Total previous step: 49401



# LOADING STAGE 3 DATA
# load('annotations_ensembl_gene_id_3')
# load('annotations_ensembl_transcript_id_3')
# load('annotations_entrezgene_id_3') #94 #23
# load('annotations_refseq_mrna_3') #144 #23 #all verylow annot. perc.
# load('annotations_accession_3') #2722 #23 #all superlow annot. perc.
# load('annotations_entrezgene_id_post_ncbi_annotation_3') #1448 #24 #all verylow annot. perc.
# load('finalized_3') #152 #24
# load('leftovers_3') #4256 #24
# input_bad_gene_symbol <- read.csv(file = 'checked_input_3_stage/input_bad_gene_symbol.tsv',  sep = '\t') #8153 #22
# Total input: 4408; Total output: 4408; Total previous step: 12561



# LOADING STAGE 4 DATA
# load('annotations_ensembl_gene_id_4')
# load('annotations_ensembl_transcript_id_4')
# load('annotations_entrezgene_id_4') #21 #25 #looked-through
# load('annotations_refseq_mrna_4')
# load('annotations_accession_4') #144 #(no output)
# load('annotations_entrezgene_id_post_ncbi_annotation_4') #2458 #26 #looked-through
# load('finalized_4') #65 #26
# load('leftovers_4') #2558 #26
# input_bad_gene_symbol <- read.csv(file = 'checked_input_4_stage/input_bad_gene_symbol.tsv',  sep = '\t') #1633 #24
# Total input: 2623; Total output: 2623; Total previous step: 4256





# 
### Check percentages of annotation for each paper in each annotation step ###
# # dir.create('percentages_of_annotation/')
# stage_ = '4'
# 
# readr::write_tsv(get_percentages_of_annotation(annotations_Probe_ID), paste0('percentages_of_annotation/annotations_Probe_ID_', stage_))
# readr::write_tsv(get_percentages_of_annotation(annotations_ensembl_gene_id), paste0('percentages_of_annotation/annotations_ensembl_gene_id_', stage_))
# readr::write_tsv(get_percentages_of_annotation(annotations_ensembl_transcript_id), paste0('percentages_of_annotation/annotations_ensembl_transcript_id_', stage_))
# readr::write_tsv(get_percentages_of_annotation(annotations_entrezgene_id), paste0('percentages_of_annotation/annotations_entrezgene_id_', stage_))
# readr::write_tsv(get_percentages_of_annotation(annotations_refseq_mrna), paste0('percentages_of_annotation/annotations_refseq_mrna_', stage_))
# readr::write_tsv(get_percentages_of_annotation(annotations_accession), paste0('percentages_of_annotation/annotations_accession_', stage_))
# readr::write_tsv(get_percentages_of_annotation(annotations_entrezgene_id_post_ncbi_annotation), paste0('percentages_of_annotation/annotations_entrezgene_id_post_ncbi_annotation_', stage_))



### Check start/end of input_ - DONE ###
# source('functions_for_genename_conversions.R')
# source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# stage_ = '4'
# input_ProbeID <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_ProbeID.tsv'), sep = '\t')
# input_EnsemblGeneId <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_EnsemblGeneId.tsv'), sep = '\t')
# input_EnsemblTranscriptId <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_EnsemblTranscriptId.tsv'), sep = '\t')
# input_EntrezGeneId <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_EntrezGeneId.tsv'), sep = '\t')
# input_RefSeqMRNA <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_RefSeqMRNA.tsv'), sep = '\t')
# input_accession <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_accession.tsv'), sep = '\t')
# input_gene_symbol <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_gene_symbol.tsv'), sep = '\t')
# input_bad_gene_symbol <- read.csv(file = paste0('checked_input_', stage_, '_stage/input_bad_gene_symbol.tsv'), sep = '\t')

### Check symbols - OK ###
# get_all_symbols_in_chrvec(input_ProbeID$Probe_ID)
# get_all_symbols_in_chrvec(input_EnsemblGeneId$Probe_ID)
# get_all_symbols_in_chrvec(input_EnsemblTranscriptId$Probe_ID)
# get_all_symbols_in_chrvec(input_EntrezGeneId$Probe_ID)
# get_all_symbols_in_chrvec(input_RefSeqMRNA$Probe_ID)
# get_all_symbols_in_chrvec(input_accession$Probe_ID)
# get_all_symbols_in_chrvec(input_gene_symbol$Probe_ID)
# get_all_symbols_in_chrvec(input_bad_gene_symbol$Probe_ID)




# temp <- subset(x = input_ProbeID, subset = is.na(input_ProbeID$Probe_ID))

### Check if ncbi annotation is done in order using annotate_identifiers_to_geneID_queries ###

# Why no XM_? give almost no results (6 per 1293) and there are relatively few of them (compared to accession)
# Why only 4 steps, if we still have 2558 leftovers after 4th step - 5th step annotated basically no of the entries (4th step annotated only 65), and it increases chaos. I could add another 3 steps to annotate 20 more, but its not worht the time (quering entrez) and effort imho
# 1) Are we using proper input for given step? - OK
# 2) Are paper column integer (for the experiments to be ordered in the same way in data and descriptions)?
# 3) Are we always merging columns of the same type?
# 4) Is exerything2 actuallized in each iterarion - OK i think
# 5) Check lengths on input_ and annotations - OK
# 5) Check number of columns annotations - OK
# Paper 72 Probe_IDs are bad - cannot find annotations for the identifiers. I dont know what they are (they are not entrez ids), thus, I will ignore them and annotate the experiment with the rest Probe-lacking experiments - OK

