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
annotations_Probe_ID <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_first_stage/input_ProbeID.tsv', str_identifier_type__ = identifiers_used_for_annotation_global, str_experiment_name = 'Probe_ID', col_types__ = 'ncccccccccccccnccc') # 195495 vs 195495
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
annotations_ensembl_gene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_first_stage/input_EnsemblGeneId.tsv', str_identifier_type__ = 'ensembl_gene_id', col_types__ = 'ncccccccccccccnccc') #15302 vs 15302
save(annotations_ensembl_gene_id, file = 'annotations_ensembl_gene_id')
# load(file = 'annotations_ensembl_gene_id')


annotations_ensembl_transcript_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_first_stage/input_EnsemblTranscriptId.tsv', str_identifier_type__ = 'ensembl_transcript_id', col_types__ = 'ncccccccccccccnccc') #2419 vs 2419
save(annotations_ensembl_transcript_id, file = 'annotations_ensembl_transcript_id')
# load(file = 'annotations_ensembl_transcript_id')


annotations_entrezgene_id <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_first_stage/input_EntrezGeneId.tsv', str_identifier_type__ = 'entrezgene_id', col_types__ = 'ncccccccccccccnccc') #29459 vs 29459
save(annotations_entrezgene_id, file = 'annotations_entrezgene_id')
# load('annotations_entrezgene_id')


annotations_refseq_mrna <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_first_stage/input_RefSeqMRNA.tsv', str_identifier_type__ = 'refseq_mrna', col_types__ = 'ncccccccccccccnccc') #4715 vs 4715
save(annotations_refseq_mrna, file = 'annotations_refseq_mrna')
# load('annotations_refseq_mrna')


annotations_accession <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_first_stage/input_accession.tsv', str_identifier_type__ = 'embl', col_types__ = 'ncccccccccccccnccc') #3412 vs 3412
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
pre_pre_annotations_gemmaData <- get_gemma_annotations_for_data(str_filename_ = 'checked_input_first_stage/input_Gemma.tsv', list_gemma_platforms = gemma_platforms, col_types__ = 'ncccccccccccccnccc') # 2002 vs 2002
save(pre_pre_annotations_gemmaData, file = 'pre_pre_annotations_gemmaData')
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
# input_gene_symbol <- read_preformated_data(str_filename = 'checked_input_first_stage/input_gene_symbol.tsv', col_types_ = 'ncccccccccccccnccc', int_numbers_are_from = 11, int_numbers_are_to = 13)
if( min(colnames(input_gene_symbol) == colnames(pre_annotations_gemmaData_input)) == 1 )
{
  temp_annotated_names_to_geneIDs <- rbind(input_gene_symbol, pre_annotations_gemmaData_input) # 5211/1571
  temp_annotated_names_to_geneIDs <- temp_annotated_names_to_geneIDs[order(temp_annotated_names_to_geneIDs$Paper),] 
}

annotated_names_to_geneIDs <- annotate_identifiers_to_geneID_wrapper_for_using_dfs_as_input(df_to_annotate = temp_annotated_names_to_geneIDs, descriptions__ = descriptions, chr_gene_identifier__ = 'Gene Name', str_experiment_name__ = 'gene_names', col_types___ = 'ncccccccccccccnccc') #6782 vs 6782
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

# # Remove identifers which already proven to be unreliable
# leftovers$everything2 <-
#   as.character(purrr::map2(
#     .x =  leftovers$everything,
#     .y = leftovers$Probe_ID,
#     .f = function(x, y) {
#       if (!is.na(y)) {
#         stringr::str_remove_all(string = x, pattern = stringr::coll(y))
#       }
#       else{
#         x
#       }
#     }
#   ))

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
length(leftovers[[1]]) + length(read.csv(file = 'checked_input_first_stage/input_bad_gene_symbol.tsv', header = T, sep = '\t')[[1]]) + length(finalized[[1]]) == length(reformated_raw_dataset[[1]])

########################################################################################
################################# FIRST STAGE OF ANNOTATION ############################
########################################################################################






########################################################################################
################################# SECOND STAGE OF ANNOTATION ###########################
########################################################################################
stage_ = 3
if (stage_ == 3){
  opts_col_types___ = 'cnccccccccccccnccccc'
  opts_checked_input_stage = 'checked_input_second_stage/'
}

#################################
####### KNOWN IDENTIFIERS ####### 
#################################
####### FIRST ORDER ANNOTATION ####### 
# These identifiers return unique gene name per identifier, because every id is assigned to known gene. This is different for Probe_IDs
annotations_ensembl_gene_id_second_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_second_stage/input_EnsemblGeneId.tsv', str_identifier_type__ = 'ensembl_gene_id', col_types__ = opts_col_types___) #12 vs 12
save(annotations_ensembl_gene_id_second_stage, file = 'annotations_ensembl_gene_id_second_stage')
# load(file = 'annotations_ensembl_gene_id_second_stage')


annotations_ensembl_transcript_id_second_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_second_stage/input_EnsemblTranscriptId.tsv', str_identifier_type__ = 'ensembl_transcript_id', col_types__ = opts_col_types___) #38 vs 38
save(annotations_ensembl_transcript_id_second_stage, file = 'annotations_ensembl_transcript_id_second_stage')
# load(file = 'annotations_ensembl_transcript_id_second_stage')


annotations_entrezgene_id_second_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_second_stage/input_EntrezGeneId.tsv', str_identifier_type__ = 'entrezgene_id', col_types__ = opts_col_types___) #16130 vs 16130
save(annotations_entrezgene_id_second_stage, file = 'annotations_entrezgene_id_second_stage')
# load('annotations_entrezgene_id_second_stage')


annotations_refseq_mrna_second_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_second_stage/input_RefSeqMRNA.tsv', str_identifier_type__ = 'refseq_mrna', col_types__ = opts_col_types___) #639 vs 639
save(annotations_refseq_mrna_second_stage, file = 'annotations_refseq_mrna_second_stage')
# load('annotations_refseq_mrna_second_stage')


annotations_accession_second_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_second_stage/input_accession.tsv', str_identifier_type__ = 'embl', col_types__ = opts_col_types___) #7489 vs 7489
save(annotations_accession_second_stage, file = 'annotations_accession_second_stage')
# load('annotations_accession_second_stage')


####### FIRST ORDER ANNOTATION ####### 
####### GET FINALIZED AND LEFTOVER DFS ####### 
if( min(colnames(annotations_ensembl_gene_id_second_stage) == colnames(annotations_ensembl_transcript_id_second_stage)) == 1 & min(colnames(annotations_ensembl_gene_id_second_stage) == colnames(annotations_entrezgene_id_second_stage)) == 1 & min(colnames(annotations_ensembl_gene_id_second_stage) == colnames(annotations_refseq_mrna_second_stage)) == 1 & min(colnames(annotations_ensembl_gene_id_second_stage) == colnames(annotations_accession_second_stage)) == 1 )
{
  temp_annotations_known_identifiers_second_stage <-
    rbind(
      annotations_ensembl_gene_id_second_stage,
      annotations_ensembl_transcript_id_second_stage,
      annotations_entrezgene_id_second_stage,
      annotations_refseq_mrna_second_stage,
      annotations_accession_second_stage
    )
}

# Prepare finalized and leftover files
finalized_known_identifiers_second_stage <- subset(x = temp_annotations_known_identifiers_second_stage, subset = !is.na(temp_annotations_known_identifiers_second_stage$external_gene_name)) #13103

leftover_known_identifiers_second_stage <- subset(x = temp_annotations_known_identifiers_second_stage, subset = is.na(temp_annotations_known_identifiers_second_stage$external_gene_name)) # 12707
# finalized_known_identifiers_second_stage + leftover_known_identifiers_second_stage should equal 25810. It does equal 25810

check_was_the_spliting_of_df_by_filtering_ok(df_original = temp_annotations_known_identifiers_second_stage, list_df_splited = list(finalized_known_identifiers_second_stage, leftover_known_identifiers_second_stage))

rm(temp_annotations_known_identifiers_second_stage)
####### GET FINALIZED AND LEFTOVER DFS ####### 
#################################
####### KNOWN IDENTIFIERS ####### 
#################################



#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################


annotated_names_to_geneIDs_second_stage <-
  annotate_identifiers_to_geneID_wrapper_for_using_dfs_as_input(
    df_to_annotate = input_gene_symbol,
    descriptions__ = descriptions,
    chr_gene_identifier__ = 'Gene Name',
    str_experiment_name__ = 'gene_names_second_stage',
    col_types___ = opts_col_types___
  ) # 2347 vs 
save(annotated_names_to_geneIDs_second_stage, file = 'annotated_names_to_geneIDs_second_stage')
# load('annotated_names_to_geneIDs_second_stage')

annotated_names_to_geneIDs_second_stage <- post_annotate_identifiers_to_geneID(annotated_names_to_geneIDs_second_stage, Probe_ID_left_ = 'Probe_ID_left_from_second_ncbi_annotation_stage')

# Annotate and save gene names
annotations_entrezgene_id_post_ncbi_annotation_second_stage <-
  master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(
    descriptions_ = descriptions,
    df_to_annotate = annotated_names_to_geneIDs_second_stage,
    str_identifier_type__ = 'entrezgene_id',
    str_experiment_name__ = 'entrezgene_id_post_ncbi_preAnnotation_second_stage',
    col_types___ = 'cnccccccccccccncccccc',
    int_numbers_are_from__ = 11,
    int_numbers_are_to__ = 13
  ) #4075 vs 4075
save(annotations_entrezgene_id_post_ncbi_annotation_second_stage, file = 'annotations_entrezgene_id_post_ncbi_annotation_second_stage')
# load('annotations_entrezgene_id_post_ncbi_annotation_second_stage')


# Prepare finalized and leftover files
finalized_entrezgene_id_post_ncbi_annotation_second_stage <- subset(x = annotations_entrezgene_id_post_ncbi_annotation_second_stage, subset = !is.na(annotations_entrezgene_id_post_ncbi_annotation_second_stage$external_gene_name))

leftover_annotations_entrezgene_id_post_ncbi_annotation_second_stage <- subset(x = annotations_entrezgene_id_post_ncbi_annotation_second_stage, subset = is.na(annotations_entrezgene_id_post_ncbi_annotation_second_stage$external_gene_name))

check_was_the_spliting_of_df_by_filtering_ok(df_original = annotations_entrezgene_id_post_ncbi_annotation_second_stage, list_df_splited = list(leftover_annotations_entrezgene_id_post_ncbi_annotation_second_stage, finalized_entrezgene_id_post_ncbi_annotation_second_stage))
#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################


#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################
leftover_annotations_entrezgene_id_post_ncbi_annotation_second_stage$external_gene_name <- NULL
leftover_known_identifiers_second_stage$external_gene_name <- NULL
leftover_known_identifiers_second_stage$Probe_ID_left_from_second_ncbi_annotation_stage <- NA



if( min(colnames(leftover_annotations_entrezgene_id_post_ncbi_annotation_second_stage) == colnames(leftover_known_identifiers_second_stage)) == 1)
{
leftovers_second_stage <-
  rbind(leftover_annotations_entrezgene_id_post_ncbi_annotation_second_stage,
    leftover_known_identifiers_second_stage)
}

# Remove identifers which already proven to be unreliable
leftovers_second_stage$everything2 <-
  as.character(purrr::map2(
    .x =  leftovers_second_stage$everything2,
    .y = leftovers_second_stage$Probe_ID,
    .f = function(x, y) {
      if (!is.na(y)) {
        stringr::str_remove_all(string = x, pattern = stringr::coll(y))
      }
      else{
        x
      }
    }
  ))

leftovers_second_stage$Probe_ID_left_from_second_annotating_stage <- leftovers_second_stage$Probe_ID
leftovers_second_stage$Probe_ID <- NA

save(leftovers_second_stage, file = 'leftovers_second_stage')
# load('leftovers_second_stage')

# rm(list = ls(pattern = ''))
#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################



#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################
temp <- finalized_entrezgene_id_post_ncbi_annotation_second_stage[['Probe_ID_left_from_second_ncbi_annotation_stage']]
finalized_entrezgene_id_post_ncbi_annotation_second_stage$Probe_ID_left_from_second_ncbi_annotation_stage <- NULL
finalized_entrezgene_id_post_ncbi_annotation_second_stage$Probe_ID_left_from_second_ncbi_annotation_stage <- temp

finalized_known_identifiers_second_stage$Probe_ID_left_from_second_ncbi_annotation_stage <- NA


if( min(colnames(finalized_known_identifiers_second_stage) == colnames(finalized_entrezgene_id_post_ncbi_annotation_second_stage)) == 1)
{
  finalized_second_stage <-
    rbind(
      finalized_known_identifiers_second_stage,
      finalized_entrezgene_id_post_ncbi_annotation_second_stage
    )
  save(finalized_second_stage, file = 'finalized_second_stage')
}
# load('finalized_second_stage')

rm(list = ls(pattern = '(annotated_)|(annotations_)|(leftover_)(.*)'))
#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################

length(leftovers[[1]]) == length(finalized_second_stage[[1]]) + length(leftovers_second_stage[[1]]) + length(read.csv(file = 'checked_input_second_stage/input_bad_gene_symbol.tsv', header = T, sep = '\t')[[1]])

rm(leftovers_for_checking_data_integrity)

########################################################################################
################################# SECOND STAGE OF ANNOTATION ###########################
########################################################################################















#########################################################################################
################################### THIRD STAGE OF ANNOTATION ###########################
#########################################################################################
stage_ = 3
if (stage_ == 3){
  opts_col_types___ = 'cnccccccccccccnccccccc'
  opts_checked_input_stage = 'checked_input_fourth_stage/'
}

#################################
####### KNOWN IDENTIFIERS ####### 
#################################
####### FIRST ORDER ANNOTATION ####### 
annotations_entrezgene_id_third_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_third_stage/input_EntrezGeneId.tsv', str_identifier_type__ = 'entrezgene_id', col_types__ = opts_col_types___) #133
save(annotations_entrezgene_id_third_stage, file = 'annotations_entrezgene_id_third_stage')
# load('annotations_entrezgene_id_third_stage')


annotations_refseq_mrna_third_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_third_stage/input_RefSeqMRNA.tsv', str_identifier_type__ = 'refseq_mrna', col_types__ = opts_col_types___) #865
save(annotations_refseq_mrna_third_stage, file = 'annotations_refseq_mrna_third_stage')
# load('annotations_refseq_mrna_third_stage')


annotations_accession_third_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = 'checked_input_third_stage/input_accession.tsv', str_identifier_type__ = 'embl', col_types__ = opts_col_types___) #4037
save(annotations_accession_third_stage, file = 'annotations_accession_third_stage')
# load('annotations_accession_third_stage')


####### FIRST ORDER ANNOTATION ####### 
####### GET FINALIZED AND LEFTOVER DFS ####### 
if( min(colnames(annotations_entrezgene_id_third_stage) == colnames(annotations_refseq_mrna_third_stage)) == 1 & min(colnames(annotations_entrezgene_id_third_stage) == colnames(annotations_accession_third_stage)) == 1 )
{
  temp_annotations_known_identifiers_third_stage <-
    rbind(
      annotations_entrezgene_id_third_stage,
      annotations_refseq_mrna_third_stage,
      annotations_accession_third_stage
    )
}

# Prepare finalized and leftover files
finalized_known_identifiers_third_stage <- subset(x = temp_annotations_known_identifiers_third_stage, subset = !is.na(temp_annotations_known_identifiers_third_stage$external_gene_name)) #361

leftover_known_identifiers_third_stage <- subset(x = temp_annotations_known_identifiers_third_stage, subset = is.na(temp_annotations_known_identifiers_third_stage$external_gene_name)) # 3031



check_was_the_spliting_of_df_by_filtering_ok(df_original = temp_annotations_known_identifiers_third_stage, list_df_splited = list(finalized_known_identifiers_third_stage, leftover_known_identifiers_third_stage))

rm(temp_annotations_known_identifiers_third_stage)
####### GET FINALIZED AND LEFTOVER DFS ####### 
#################################
####### KNOWN IDENTIFIERS ####### 
#################################
cnccccccccccccnccccccc
ccccccccccccccnccccccc
#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################
# input_gene_symbol <- read_preformated_data(str_filename = 'checked_input_third_stage/input_gene_symbol.tsv', col_types_ = 'ccccccccccccccnccccccc', int_numbers_are_from = 11, int_numbers_are_to = 13)
annotated_names_to_geneIDs_third_stage <-
  annotate_identifiers_to_geneID_wrapper_for_using_dfs_as_input(
    df_to_annotate = input_gene_symbol,
    descriptions__ = descriptions,
    chr_gene_identifier__ = 'Gene Name',
    str_experiment_name__ = 'gene_names_third_stage',
    col_types___ = opts_col_types___
  ) # 1556 vs 
save(annotated_names_to_geneIDs_third_stage, file = 'annotated_names_to_geneIDs_third_stage')
# load('annotated_names_to_geneIDs_third_stage')

annotated_names_to_geneIDs_third_stage <- post_annotate_identifiers_to_geneID(annotated_names_to_geneIDs_third_stage, Probe_ID_left_ = 'Probe_ID_left_from_third_ncbi_annotation_stage')


# Annotate and save gene names
annotations_entrezgene_id_post_ncbi_annotation_third_stage <-
  master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(
    descriptions_ = descriptions,
    df_to_annotate = annotated_names_to_geneIDs_third_stage,
    str_identifier_type__ = 'entrezgene_id',
    str_experiment_name__ = 'entrezgene_id_post_ncbi_preAnnotation_third_stage',
    col_types___ = 'cnccccccccccccncccccccc',
    int_numbers_are_from__ = 11,
    int_numbers_are_to__ = 13
  ) #4075 vs 4075
save(annotations_entrezgene_id_post_ncbi_annotation_third_stage, file = 'annotations_entrezgene_id_post_ncbi_annotation_third_stage')
# load('annotations_entrezgene_id_post_ncbi_annotation_third_stage')


# Prepare finalized and leftover files
finalized_entrezgene_id_post_ncbi_annotation_third_stage <- subset(x = annotations_entrezgene_id_post_ncbi_annotation_third_stage, subset = !is.na(annotations_entrezgene_id_post_ncbi_annotation_third_stage$external_gene_name))

leftover_annotations_entrezgene_id_post_ncbi_annotation_third_stage <- subset(x = annotations_entrezgene_id_post_ncbi_annotation_third_stage, subset = is.na(annotations_entrezgene_id_post_ncbi_annotation_third_stage$external_gene_name))

check_was_the_spliting_of_df_by_filtering_ok(df_original = annotations_entrezgene_id_post_ncbi_annotation_third_stage, list_df_splited = list(leftover_annotations_entrezgene_id_post_ncbi_annotation_third_stage, finalized_entrezgene_id_post_ncbi_annotation_third_stage))
#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################


#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################
leftover_annotations_entrezgene_id_post_ncbi_annotation_third_stage$external_gene_name <- NULL
leftover_known_identifiers_third_stage$external_gene_name <- NULL
leftover_known_identifiers_third_stage$Probe_ID_left_from_third_ncbi_annotation_stage <- NA

if( min(colnames(leftover_annotations_entrezgene_id_post_ncbi_annotation_third_stage) == colnames(leftover_known_identifiers_third_stage)) == 1)
{
leftovers_third_stage <-
  rbind(leftover_annotations_entrezgene_id_post_ncbi_annotation_third_stage,
        leftover_known_identifiers_third_stage)
}
# Remove identifers which already proven to be unreliable
leftovers_third_stage$everything2 <-
  as.character(purrr::map2(
    .x =  leftovers_third_stage$everything2,
    .y = leftovers_third_stage$Probe_ID,
    .f = function(x, y) {
      if (!is.na(y)) {
        stringr::str_remove_all(string = x, pattern = stringr::coll(y))
      }
      else{
        x
      }
    }
  ))

leftovers_third_stage$Probe_ID_left_from_third_annotating_stage <- leftovers_third_stage$Probe_ID
leftovers_third_stage$Probe_ID <- NA

save(leftovers_third_stage, file = 'leftovers_third_stage')
# load('leftovers_third_stage')

# rm(list = ls(pattern = ''))
#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################



#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################
temp <- finalized_entrezgene_id_post_ncbi_annotation_third_stage[['Probe_ID_left_from_third_ncbi_annotation_stage']]
finalized_entrezgene_id_post_ncbi_annotation_third_stage$Probe_ID_left_from_third_ncbi_annotation_stage <- NULL
finalized_entrezgene_id_post_ncbi_annotation_third_stage$Probe_ID_left_from_third_ncbi_annotation_stage <- temp

finalized_known_identifiers_third_stage$Probe_ID_left_from_third_ncbi_annotation_stage <- NA


if( min(colnames(finalized_known_identifiers_third_stage) == colnames(finalized_entrezgene_id_post_ncbi_annotation_third_stage)) == 1)
{
  finalized_third_stage <-
    rbind(
      finalized_known_identifiers_third_stage,
      finalized_entrezgene_id_post_ncbi_annotation_third_stage
    )
}
save(finalized_third_stage, file = 'finalized_third_stage')
# load('finalized_third_stage')

rm(list = ls(pattern = '(finalized_)|(annotated_)|(annotations_)|(leftover_)(.*)'))
#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################

length(leftovers_second_stage[[1]]) == length(finalized_third_stage[[1]]) + length(leftovers_third_stage[[1]]) + length(read.csv(file = 'checked_input_third_stage/input_bad_gene_symbol.tsv', header = T, sep = '\t')[[1]])


#########################################################################################
################################### THIRD STAGE OF ANNOTATION ###########################
#########################################################################################





##########################################################################################
################################### FOURTH STAGE OF ANNOTATION ###########################
##########################################################################################
stage_ = 4
if (stage_ == 4){
  opts_col_types___ = 'cnccccccccccccnccccccccc'
  opts_checked_input_stage = 'checked_input_fourth_stage/'
  
}

#################################
####### KNOWN IDENTIFIERS ####### 
#################################
####### FIRST ORDER ANNOTATION ####### 
annotations_entrezgene_id_fourth_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_EntrezGeneId.tsv'), str_identifier_type__ = 'entrezgene_id', col_types__ = opts_col_types___) #25
save(annotations_entrezgene_id_fourth_stage, file = 'annotations_entrezgene_id_fourth_stage')
# load('annotations_entrezgene_id_fourth_stage')
# Does not return any result


annotations_refseq_mrna_fourth_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_RefSeqMRNA.tsv'), str_identifier_type__ = 'refseq_mrna', col_types__ = opts_col_types___) #8
save(annotations_refseq_mrna_fourth_stage, file = 'annotations_refseq_mrna_fourth_stage')
# load('annotations_refseq_mrna_fourth_stage')


annotations_accession_fourth_stage <- master_annotator_for_known_identfiers(descriptions_ = descriptions, str_filename_ = paste0(opts_checked_input_stage, 'input_accession.tsv'), str_identifier_type__ = 'embl', col_types__ = opts_col_types___) #207
save(annotations_accession_fourth_stage, file = 'annotations_accession_fourth_stage')
# load('annotations_accession_fourth_stage')


####### FIRST ORDER ANNOTATION ####### 
####### GET FINALIZED AND LEFTOVER DFS ####### 
if( min(colnames(annotations_entrezgene_id_fourth_stage) == colnames(annotations_refseq_mrna_fourth_stage)) == 1 & min(colnames(annotations_entrezgene_id_fourth_stage) == colnames(annotations_accession_fourth_stage)) == 1 )
{
  temp_annotations_known_identifiers_fourth_stage <-
    rbind(
      annotations_refseq_mrna_fourth_stage,
      annotations_mgi_fourth_stage,
      annotations_accession_fourth_stage
  )
}
# Prepare finalized and leftover files
finalized_known_identifiers_fourth_stage <- subset(x = temp_annotations_known_identifiers_fourth_stage, subset = !is.na(temp_annotations_known_identifiers_fourth_stage$external_gene_name)) #361

leftover_known_identifiers_fourth_stage <- subset(x = temp_annotations_known_identifiers_fourth_stage, subset = is.na(temp_annotations_known_identifiers_fourth_stage$external_gene_name)) # 3031


check_was_the_spliting_of_df_by_filtering_ok(df_original = temp_annotations_known_identifiers_fourth_stage, list_df_splited = list(finalized_known_identifiers_fourth_stage, leftover_known_identifiers_fourth_stage))

rm(temp_annotations_known_identifiers_fourth_stage)
####### GET FINALIZED AND LEFTOVER DFS ####### 
#################################
####### KNOWN IDENTIFIERS ####### 
#################################


#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################
annotated_names_to_geneIDs_fourth_stage <-
  annotate_identifiers_to_geneID_wrapper_for_using_dfs_as_input(
    df_to_annotate = input_gene_symbol,
    descriptions__ = descriptions,
    chr_gene_identifier__ = 'Gene Name',
    str_experiment_name__ = 'gene_names_fourth_stage',
    col_types___ = opts_col_types___
  ) #1269
save(annotated_names_to_geneIDs_fourth_stage, file = 'annotated_names_to_geneIDs_fourth_stage')
# load('annotated_names_to_geneIDs_fourth_stage')

annotated_names_to_geneIDs_fourth_stage <- post_annotate_identifiers_to_geneID(annotated_names_to_geneIDs_fourth_stage, Probe_ID_left_ = 'Probe_ID_left_from_fourth_ncbi_annotation_stage')

# Annotate and save gene names
annotations_entrezgene_id_post_ncbi_annotation_fourth_stage <-
  master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input(
    descriptions_ = descriptions,
    df_to_annotate = annotated_names_to_geneIDs_fourth_stage,
    str_identifier_type__ = 'entrezgene_id',
    str_experiment_name__ = 'entrezgene_id_post_ncbi_preAnnotation_fourth_stage',
    col_types___ = 'cnccccccccccccncccccccccc',
    int_numbers_are_from__ = 11,
    int_numbers_are_to__ = 13
  ) #1269
save(annotations_entrezgene_id_post_ncbi_annotation_fourth_stage, file = 'annotations_entrezgene_id_post_ncbi_annotation_fourth_stage')
# load('annotations_entrezgene_id_post_ncbi_annotation_fourth_stage')


# Prepare finalized and leftover files
finalized_entrezgene_id_post_ncbi_annotation_fourth_stage <-
  subset(
    x = annotations_entrezgene_id_post_ncbi_annotation_fourth_stage,
    subset = !is.na(
      annotations_entrezgene_id_post_ncbi_annotation_fourth_stage$external_gene_name
    )
  )

leftover_annotations_entrezgene_id_post_ncbi_annotation_fourth_stage <-
  subset(
    x = annotations_entrezgene_id_post_ncbi_annotation_fourth_stage,
    subset = is.na(
      annotations_entrezgene_id_post_ncbi_annotation_fourth_stage$external_gene_name
    )
  )

check_was_the_spliting_of_df_by_filtering_ok(
  df_original = annotations_entrezgene_id_post_ncbi_annotation_fourth_stage,
  list_df_splited = list(
    leftover_annotations_entrezgene_id_post_ncbi_annotation_fourth_stage,
    finalized_entrezgene_id_post_ncbi_annotation_fourth_stage
  )
)
#########################################
####### NCBI GENE NAME ANNOTATION ####### 
#########################################


#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################
leftover_annotations_entrezgene_id_post_ncbi_annotation_fourth_stage$external_gene_name <- NULL
leftover_known_identifiers_fourth_stage$external_gene_name <- NULL
leftover_known_identifiers_fourth_stage$Probe_ID_left_from_fourth_ncbi_annotation_stage <- NA

leftovers_lenght_fourth_stage <- length(leftover_annotations_entrezgene_id_post_ncbi_annotation_fourth_stage[[1]]) + length(leftover_known_identifiers_fourth_stage[[1]]) + length(readr::read_tsv('checked_input_fourth_stage/input_bad_gene_symbol.tsv')[[1]])

colnames(leftover_annotations_entrezgene_id_post_ncbi_annotation_fourth_stage) == colnames(leftover_known_identifiers_fourth_stage)

leftovers_fourth_stage <-
  rbind(leftover_annotations_entrezgene_id_post_ncbi_annotation_fourth_stage,
        leftover_known_identifiers_fourth_stage)

# Remove identifers which already proven to be unreliable
leftovers_fourth_stage$everything2 <-
  as.character(purrr::map2(
    .x =  leftovers_fourth_stage$everything2,
    .y = leftovers_fourth_stage$Probe_ID,
    .f = function(x, y) {
      if (!is.na(y)) {
        stringr::str_remove_all(string = x, pattern = stringr::coll(y))
      }
      else{
        x
      }
    }
  ))

leftovers_fourth_stage$Probe_ID_left_from_fourth_annotating_stage <- leftovers_fourth_stage$Probe_ID
leftovers_fourth_stage$Probe_ID <- NA

save(leftovers_fourth_stage, file = 'leftovers_fourth_stage')
# load('leftovers_fourth_stage')

# rm(list = ls(pattern = ''))
#######################################
####### PREPARE LEFTOVERS TABLE ####### 
#######################################



#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################
finalized_known_identifiers_fourth_stage$Probe_ID_left_from_fourth_ncbi_annotation_stage <- NA
finalized_known_identifiers_fourth_stage <- finalized_known_identifiers_fourth_stage[, c(1:24,26,25)]

colnames(finalized_known_identifiers_fourth_stage) == colnames(finalized_entrezgene_id_post_ncbi_annotation_fourth_stage)

finalized_fourth_stage <-
  rbind(
    finalized_known_identifiers_fourth_stage,
    finalized_entrezgene_id_post_ncbi_annotation_fourth_stage
  )
save(finalized_fourth_stage, file = 'finalized_fourth_stage')
# load('finalized_fourth_stage')

rm(list = ls(pattern = '(finalized_)|(annotated_)|(annotations_)|(leftover_)(.*)'))
#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################

length(leftovers_second_stage[[1]]) == length(finalized_third_stage[[1]]) + leftovers_lenght_third_stage
##########################################################################################
################################### FOURTH STAGE OF ANNOTATION ###########################
##########################################################################################








####### GET NON-ANNOTATED DATA IN FINAL TABLE ####### 
# input_bad_gene_symbol - manual
# logFCs for the same gene within single Experiment are avaraged using median. No other filtering is used
# It is critical to use this workflow only after producing actual, tested proper merged dataset - we need to add all the data and check if the lenght of this set is equal to lenght of input dataset. 
final_merged_dataset <- gather_all_datasets_into_single_df(regex_pattern_to_find_datasets_with = '^annotations.*')

noNAs_final_merged_dataset <- final_merged_dataset %>%
  subset(subset = !(is.na(final_merged_dataset$external_gene_name) | final_merged_dataset$external_gene_name == ''))

not_annotated_data_2 <- final_merged_dataset %>%
  subset(subset = (is.na(final_merged_dataset$external_gene_name) | final_merged_dataset$external_gene_name == ''))

check_was_the_spliting_of_df_by_filtering_ok(str_what_was_splited = 'final_merged_dataset', df_original = final_merged_dataset, list_df_splited = list(noNAs_final_merged_dataset, not_annotated_data_2))

not_annotated_data <- rbind(not_annotated_data_1, not_annotated_data_2)

####### GET NON-ANNOTATED DATA IN FINAL TABLE ####### 



####### ANALYZE NON-ANNOTATED DATA ####### 
# There are some entries that have p above 0.05
not_annotated_data$composite_id <- paste(not_annotated_data$Experiment, not_annotated_data$adj_p, not_annotated_data$p, not_annotated_data$logFC, sep = "__")

whole_dataset <- read_preformated_data(str_filename = 'data_whole.tsv', col_types_ = 'nccccccccccccc', int_numbers_are_from = 11, int_numbers_are_to = 13)

whole_dataset$composite_id <- paste(whole_dataset$Experiment, whole_dataset$adj_p, whole_dataset$p, whole_dataset$logFC, sep = "__")

whole_dataset$all_names <- as.character(apply(X = whole_dataset, MARGIN = 1, FUN = function(x) { paste(x, collapse = '__')  } ))

whole_and_not_annotated <- merge(x = not_annotated_data, y = whole_dataset, by = 'composite_id', all.x = T)
# test_short <- whole_and_not_annotated[1:1000,]


is_the_ProbeID_from_notAnnotated_in_the_wholeDataset <-
  purrr::map2_lgl(
    .x = whole_and_not_annotated$Probe_ID.x,
    .y = whole_and_not_annotated$all_names,
    .f = function(x, y) {
      
      pattern_ <- paste0('(.*)\\Q', x, '\\E(.*)')
      
      stringr::str_detect(string = y,
                          pattern = pattern_)
    }
  )

proper_whole_and_not_annotated <- subset(x = whole_and_not_annotated, subset = is_the_ProbeID_from_notAnnotated_in_the_wholeDataset)



####### ANALYZE NON-ANNOTATED DATA ####### 



####### PREPARING FINAL TABLE ####### 
# Here we collapse multiple genenames to single gene name
noNAs_final_merged_dataset$final_gene_name <- as.character(lapply(
  X = noNAs_final_merged_dataset$external_gene_name,
  FUN = function(x) {
    return(select_best_geneName_wrapper_for_single_string(x)) ### !!! Gm's cna have 5 digits... damn ### !!! Ive changed \\d{5} to \\d{5,}, because the latter one matches exactly 5 and not 5 and more!
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


####### PREPARING FINAL TABLE ####### 


####### CHECK IF INPUT DATA IS GOOD ####### 




####### PREPARING NEW, IMPROVED DATASET ####### 







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




#######################################
####### PREPARE FINALIZED TABLE ####### 
#######################################

check_was_the_spliting_of_df_by_filtering_ok(df_original = raw_dataset, list_df_splited = list(finalized, leftovers_for_checking_data_integrity))

length(finalized[[1]]) + length(leftovers_for_checking_data_integrity[[1]])




















































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




test <- c(1, 10, 2)

test[order(as.character(test))]
