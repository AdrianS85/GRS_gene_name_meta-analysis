source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')

opts_ann_1_n = 'Probe_ID_left_from_first_ncbi_annotation_stage'
opts_ann_1 = 'Probe_ID_left_from_first_annotating_stage'
opts_ann_2_n = 'Probe_ID_left_from_second_ncbi_annotation_stage'
opts_ann_2 = 'Probe_ID_left_from_second_annotating_stage'
opts_ann_3_n = 'Probe_ID_left_from_third_ncbi_annotation_stage'
opts_ann_3 = 'Probe_ID_left_from_third_annotating_stage'
opts_ann_4_n = 'Probe_ID_left_from_fourth_ncbi_annotation_stage'
opts_ann_4 = 'Probe_ID_left_from_fourth_annotating_stage'
opts_nb_of_analysis_stages = 4

for (current_stage in seq(opts_nb_of_analysis_stages)) {
  if (current_stage == 1) {
    load('finalized') 
    input_bad_gene_symbol <- dplyr::select(read.csv(file = 'checked_input_1_stage/input_bad_gene_symbol.tsv',  sep = '\t'), Probe_ID, dplyr::everything())

  } else if (current_stage != 1) {
    load(paste0('finalized_', current_stage))
    input_bad_gene_symbol <- read.csv(file = paste0('checked_input_', current_stage,'_stage/input_bad_gene_symbol.tsv'),  sep = '\t')
  }

  if (current_stage == opts_nb_of_analysis_stages) {
    load(paste0('leftovers_', opts_nb_of_analysis_stages)) # Load undone leftovers
  }

  finalized <- uniformize_finalized_colnames()
  input_bad_gene_symbol <- uniformize_bad_gene_symbol_colnames()
  leftovers <- uniformize_leftovers_colnames()

  assign(x = paste0('finalized_', current_stage), value = finalized)
  assign(x = paste0('bad_gene_symbol_', current_stage), value = input_bad_gene_symbol)
  
  rm(finalized, input_bad_gene_symbol)
}

list_all_dataobjects_colnames <-
  lapply(
    ls(pattern = '^bad_gene_symbol.*|^finalized_*|^leftovers'),
    FUN = function(x) {
      colnames(eval(parse(text = x), parent.frame()))
    }
  )

if(are_vectors_the_same(chr_vec_list = list_all_dataobjects_colnames)){
  bads <- rlist::list.rbind(lapply(
    ls(pattern = '^bad_gene_symbol.*|^leftovers'),
    FUN = function(x) {
      eval(parse(text = x), parent.frame())
    }
  ))
}


bads2 <- get_usable_bad_ids(bads_ = bads)


if(are_vectors_the_same(chr_vec_list = list_all_dataobjects_colnames)){
  final_good_dataset <- rlist::list.rbind(lapply(
    ls(pattern = '^bads2|^finalized.*'),
    FUN = function(x) {
      eval(parse(text = x), parent.frame())
    }
  ))
  
  final_good_dataset <- final_good_dataset[order(final_good_dataset$Paper),]
  
  # Here we collapse multiple genenames to single gene name
  final_good_dataset$final_gene_name <- as.character(lapply(
    X = final_good_dataset$external_gene_name,
    FUN = function(x) {
      return(select_best_geneName_wrapper_for_single_string(x)) ### !!! Gm's cna have 5 digits... damn ### !!! Ive changed \\d{5} to \\d{5,}, because the latter one matches exactly 5 and not 5 and more!
    }
  ))
  
  # readr::write_tsv(x = final_merged_dataset, path = 'final_merged_dataset.tsv')
  
  # Uniformarize gene names by making them all all-lower case
  final_good_dataset$lower_final_gene_name <- tolower(final_good_dataset$final_gene_name)
  # save(final_good_dataset, file = 'final_good_dataset')
  # load('final_good_dataset')
  # readr::write_tsv(x = final_good_dataset, path = 'final_good_dataset.tsv')
}



# Medianing values for the same gene within Experiment
medianed_final_good_dataset <- final_good_dataset %>%
  dplyr::group_by(Experiment, lower_final_gene_name) %>%
  dplyr::summarize(logFC_median = median(logFC))

# This is a very sparse matrix. Perhaps try working with it as such
spread_medianed_final_good_dataset <- tidyr::spread(data = medianed_final_good_dataset, key = Experiment, value = logFC_median)
# save(spread_medianed_final_good_dataset, file = 'spread_medianed_final_good_dataset')
# load('spread_medianed_final_good_dataset')
# readr::write_tsv(x = spread_medianed_final_good_dataset, path = 'spread_medianed_final_good_dataset.tsv')



bool_entries_with_3_or_more_values <- get_subset_vector_for_entries_with_3_or_more_values_per_paper(final_good_dataset_ = final_good_dataset)

at_least_in_3_papers_spread_med_fin_g_ds <- subset(x = spread_medianed_final_good_dataset, subset = bool_entries_with_3_or_more_values)
# save(at_least_in_3_papers_spread_med_fin_g_ds, file = 'at_least_in_3_papers_spread_med_fin_g_ds')


  




####### PREPARING FINAL TABLE ####### 




  




































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
