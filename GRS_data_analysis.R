source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')


descriptions <- readr::read_tsv("descriptions.txt", col_types = "nccccccccccccccccccccccc")
load('at_least_in_3_papers_spread_med_fin_g_ds')
at_least_in_3_papers_spread_med_fin_g_ds[is.na(at_least_in_3_papers_spread_med_fin_g_ds)] <- 0
data_analysis_input <- tidyr::gather(data = at_least_in_3_papers_spread_med_fin_g_ds, key = 'Experiment', value = 'logFC', -lower_final_gene_name)



descriptions$Repetitions <- stringr::str_replace_all(string = descriptions$Repetitions, pattern = "4 weeks  3 days", replacement = '31 days')
descriptions$Repetitions_clean_days <- cleanup_differing_units(charvec_ = descriptions$Repetitions, unit_identifier_regexList = list('.*d.*', '.*w.*'), unit_multiplier_list = list(1, 7))[[1]]



descriptions$Duration <- stringr::str_replace_all(string = descriptions$Duration, pattern = ",", replacement = '.')
descriptions$Duration_clean_minutes <- cleanup_differing_units(charvec_ = descriptions$Duration, unit_identifier_regexList = list('.*m.*', '.*h.*', '.*w.*'), unit_multiplier_list = list(1, 60, 10080))[[1]]

descriptions$Measurement_latency


descriptions$Gender_clean <- stringr::str_replace(string = descriptions$Gender, pattern = 'fem.*', replacement = 'female')


descriptions$Brain_part_clean <- brain_part_cleanup_wrapper()




descriptions$Group_ID # is good as is.
descriptions$Species # is good as is.



test <- data.frame(brain_part_cleanup_wrapper())
brain_part_test <- descriptions['Brain_part']

sensitivity_cleanup_wrapper <- function(sensitivity_ = descriptions$Stress_sensitivity)
{
  sensitivity_ <- tolower(sensitivity_)
  # sensitivity_[sensitivity_ == 'social defeat stress'] <- 'social defeat'

  
  return(sensitivity_)
}
levels(as.factor(sensitivity_cleanup_wrapper()))





stress_cleanup_wrapper <- function(stress_ = descriptions$Stress)
{
  stress_ <- tolower(stress_)
  stress_[stress_ == 'social defeat stress'] <- 'social defeat'
  stress_[stress_ == 'immobilization and tail-shocks \n'] <- 'immobilization and tail shocks'
  stress_ <- stringr::str_replace(string = stress_, pattern = '^forced sw.*', replacement = 'forced swim')
  stress_[stress_ == 'restraint stress'] <- 'restraint'
  stress_[stress_ == 'fear conditioned'] <- 'fear conditioning'
  
  return(stress_)
}
levels(as.factor(stress_cleanup_wrapper()))




# TESTING DATA
load('final_good_dataset')
load('reformated_raw_dataset')


post_annotation <- subset(final_good_dataset, final_good_dataset$lower_final_gene_name == 'grm1')
pre_annotation <- subset(reformated_raw_dataset, stringr::str_detect(string = tolower(reformated_raw_dataset$everything), pattern = 'sgk1'))
# TESTING DATA


### GET NUMBERS OF GENES IN EXPS AND PAPERS ###
load('spread_medianed_final_good_dataset')

spread_medianed_final_good_dataset[is.na(spread_medianed_final_good_dataset)] <- 0

exp_number_and_percentage <- get_number_and_percentage_of_directionality_of_exp_first_column_names(spread_medianed_final_good_dataset)
save(exp_number_and_percentage, file = 'exp_number_and_percentage')
readr::write_tsv(exp_number_and_percentage, 'exp_number_and_percentage.tsv')

paper_number <- get_number_and_percentage_of_directionality_of_paper_first_column_names(spread_medianed_final_good_dataset)
save(paper_number, file = 'paper_number')
readr::write_tsv(paper_number, 'paper_number.tsv')
### GET NUMBERS OF GENES IN EXPS AND PAPERS ###



### GETING CORRELATIONS ###
gather_spread_medianed_final_good_dataset <- tidyr::gather(data = spread_medianed_final_good_dataset, key = "exp", value = "logFC", na.rm = T, -lower_final_gene_name)
save(gather_spread_medianed_final_good_dataset, file = 'gather_spread_medianed_final_good_dataset')
# gather_amygdala <- tidyr::gather(data = amygdala, key = "exp", value = "logFC", na.rm = T, -lower_final_gene_name)

descriptions_for_correlations <- descriptions %>%
  dplyr::select(Group_ID, Gender_clean, Repetitions_clean_days, Duration_clean_minutes, Brain_part_clean)

temp <- merge(x = gather_amygdala, y = descriptions_for_correlations, by.x = 'exp', by.y = 'Group_ID', all.x = T)

temp_2 <- temp %>%
  dplyr:: group_by(lower_final_gene_name) %>%
  tidyr::nest()

temp_2$cor <- as.numeric(purrr::map(.x = temp_2$data, .f = function(x){ cor(x = x$logFC, y = x$Duration_clean_minutes, use = "pairwise.complete.obs")  }))
### GETING CORRELATIONS ###



  





temp_2 <- temp %>%
  dplyr:: group_by(lower_final_gene_name) %>%
  tidyr::nest() %>%
  purrr::map(.f = function(x){ cor(x = x$logFC, y = x$Gender_clean, method = 'spearman')  })











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
