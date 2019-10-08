`%>%` <- dplyr::`%>%`


kill_corrupted_e_notation <- function(df_temp_data_, str_to_substitute_corrupted_data_with)
{
  df_temp_data_$corrupted <- stringr::str_detect(string = df_temp_data_$logFC, pattern = 'e')
  
  # Find cells with e-annotation
  df_temp_data_ <- dplyr::mutate(.data = df_temp_data_, 
                                 value = dplyr::if_else(condition = corrupted,
                                                 true = gsub(
                                                   pattern = '(.*)e\\Q+\\E', 
                                                   replacement = '',
                                                   x = logFC), 
                                                 false = '0', 
                                                 missing = '0')
  )
  
  df_temp_data_$value <- as.numeric(df_temp_data_$value)
  
  # Set all numbers with e number higher than 2 to constant value
  df_temp_data_ <- dplyr::mutate(.data = df_temp_data_, 
                                 logFC = dplyr::if_else(condition = value > 2,
                                                 true = str_to_substitute_corrupted_data_with, 
                                                 false = logFC, 
                                                 missing = logFC)
  )
  
  df_temp_data_$corrupted <- NULL
  df_temp_data_$value <- NULL
  
  return(df_temp_data_)
}



# In col types input floats that can contain either . or , as chars. Later they are converted to numeric based on int_numbers_are_ arguments. Needs at least columns: Paper, Experiment, Probe_Id logFC
read_preformated_data <- function(str_filename, col_types_ = 'nncccccc', int_numbers_are_from = 6, int_numbers_are_to = 8, str_substitute_inf_with = '15')
{
  temp_data <- readr::read_tsv(
    str_filename, 
    col_types = col_types_)
  
  # Replace , with . to make number actual R numerics
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(
    X = temp_data[int_numbers_are_from:int_numbers_are_to], 
    FUN = function(x) 
    { stringr::str_replace(
      string = x, 
      pattern = ",", 
      replacement = ".")
    }) ### COOL CONSTRUCT 
  
  # Replace inf values with constant number
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(
    X = temp_data[int_numbers_are_from:int_numbers_are_to],
    FUN = function(x) { gsub(
      pattern = "^[I|i]nf(.*)", 
      replacement = str_substitute_inf_with,
      x = x)})
  
  # Replace too large e numbers with constant number
  temp_data <- kill_corrupted_e_notation(df_temp_data_ = temp_data, str_to_substitute_corrupted_data_with = str_substitute_inf_with)
  
  # Change number columns to numeric type
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(
    temp_data[int_numbers_are_from:int_numbers_are_to], 
    as.numeric)
  
  return(temp_data)
}



# Class of input is a list. Defaults to writing lenghts of object on first level of depth of the list
write_lenghts_of_list_objects <- function(list_, string_name_of_the_file, int_length_at_this_depth = 1)
{
  temp <- lapply(
    X = list_, 
    FUN = function(x){ 
      length(x[[int_length_at_this_depth]]) })
  temp2 <- as.data.frame( rlist::list.rbind(temp) )
  write.table(temp2, string_name_of_the_file, sep = '\t')
  rm(temp, temp2)
} 



set_mart_to_be_used <- function(str_vector_of_species_names_, int_loop = 1)
{
  message( paste0('Setting mart for step ', int_loop, "...") )
  
  usedMart__ <- switch(str_vector_of_species_names_,
                       "mice" = biomaRt::useMart(
                         "ENSEMBL_MART_ENSEMBL", 
                         dataset = "mmusculus_gene_ensembl"),
                       "rats" = biomaRt::useMart(
                         "ENSEMBL_MART_ENSEMBL", 
                         dataset = "rnorvegicus_gene_ensembl"),
                       "humans" = biomaRt::useMart(
                         "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl"),
                       "squirrelmonkeys" = biomaRt::useMart(
                         "ENSEMBL_MART_ENSEMBL", 
                         dataset = "sbboliviensis_gene_ensembl"))
  
  message( paste0('Mart set as ', usedMart__@dataset, ' for step ', int_loop) )
  return(usedMart__)
}



### 
get_the_potental_identifiers <- function(usedMart___, int_loop_nb)
{
  # Here we extract all the potential gene identifiers
  ##### !!! THIS MAY NEED FURTHER WORK !!! ##### 
  message( paste0('Extracting potental_identifiers for step ', int_loop_nb) )
  
  potental_identifiers <- c(usedMart___@filters[grep(pattern = "^ensembl(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^refseq(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^affy(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^agilent(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^illumina(.*)", usedMart___@filters[[1]]) , 1])
  
  message( paste0('Potental identifiers for step ', int_loop_nb, 'extracted') )
  
  return(potental_identifiers)
}



set_0_hit_annotations_to_na <- function(list_of_dataframes) 
{ 
  if(length(list_of_dataframes) >= 1 && length(list_of_dataframes[[1]][[1]]) > 0 )
  {
    list_of_dataframes <- names(list_of_dataframes[[1]][1])
  }
  else
  {
    list_of_dataframes <- NA
  }
}



# usedMart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
# filters <- biomaRt::listFilters(mart = usedMart)
# If we have one, specific identifer, than we want to apply it to all experiments. If we have a list of identifers, than this means that each of these identifier should be applyied to corresponding experiment
set_identifiers_used_for_annotation_if_not_probeID <- function(str_identifier_type, list_LIST_DATA) 
{
  if(length(str_identifier_type) == 1)
  {
    return( rep(str_identifier_type, length(list_LIST_DATA)) )
  }
  else # This is not good, because if we have single identifier for a list of experiments, but it is not any of abovementioned identifier, than it will fuck up program. the solution should be to ask if the lenght(str_identifier_type) == 1
  {
    return(str_identifier_type)
  }
}



change_name_to_proper_format_given_the_species <- function(df_name, str_species)
{
  df_name$Old_Probe_ID <- df_name$Probe_ID
  
  if(str_species %in% c('mice', 'rats'))
  {
    temp1 <- subset(x = df_name, str_detect(string = df_name$Probe_ID, pattern = '(?i)(loc)|(rgd)')) #case insensitive
    temp2 <- subset(x = df_name, str_detect(string = df_name$Probe_ID, pattern = '(?i)(rik)'))
    temp3 <- subset(x = df_name, !str_detect(string = df_name$Probe_ID, pattern = '(?i)(rik)|(loc)|(rgd)')) 
    
    temp1$Probe_ID <- toupper(temp1$Probe_ID)
    
    temp2$Probe_ID <- toupper(temp2$Probe_ID)
    temp2$Probe_ID <- str_replace(string = temp2$Probe_ID, pattern = '(?i)(rik)', replacement = 'Rik')
    
    temp3$Probe_ID <- tolower(temp3$Probe_ID)
    substr(temp3$Probe_ID, start = 1, stop = 1) <- toupper(substr(temp3$Probe_ID, 1, 1))
    
    
    temp4 <- rbind(temp1, temp2, temp3)
    
    return(temp4)
  }
  else if(str_species %in% c('humans', 'squirrelmonkeys'))
  {
    temp <- df_name
    temp$Probe_ID <- toupper(temp$Probe_ID)
    return(temp)
  }
}



get_proper_length_vector_for_checking_annotation_percentage <- function(list_LIST_DATA__, int_Probe_IDs_to_test_)
{
  temp <- list()
  for (n in seq_along(list_LIST_DATA__))
  {
    if(length(list_LIST_DATA__[[n]][[1]]) <= int_Probe_IDs_to_test_)
    {
      temp[n] <- length(list_LIST_DATA__[[n]][[1]])
    }
    else
    {
      temp[n] <- int_Probe_IDs_to_test_
    }
  }
  return(rlist::list.rbind(temp))
}



### list_LIST_DATA_ format: just the LIST_DATA set in previous lines, str_vector_of_species_names  format: small letters, english plural of species. str_vector_of_species_names includes names for species for each experiment in list data. str_vector_of_experiment_ids includes names which identifiy any given experiment
get_the_highest_hit_returning_id_type <- function(str_LIST_DATA_name, descriptions_ = descriptions, int_Probe_IDs_to_test = 200, str_experiment_name = experiment_name)
{
  PRE_DATA <- read_preformated_data(str_filename = str_LIST_DATA_name)
  
  exp_species <- descriptions_ %>%
    dplyr::select("Paper_ID", "Species") %>%
    unique() %>%
    dplyr::filter(Paper_ID %in% unique(PRE_DATA$Paper))
  
  str_vector_of_species_names <- exp_species$Species
  
  str_vector_of_experiment_ids <- exp_species$Paper_ID
  
  
  list_LIST_DATA_ <- split(PRE_DATA, f = PRE_DATA$Paper)
  
  ### We need to first check appropriate probe ids on smaller dataset and only then do actual annotation, because its too slow otherwise. Hence this shortened list !!! ADD RANDOM SELECTION OF ROWS! !!!
  
  SHORT_LIST_DATA <- lapply(X = list_LIST_DATA_, FUN = function(x){ x[1:int_Probe_IDs_to_test,] })
  ANNOT_SHORT_LIST_DATA <- rep(list(list()), times = length(SHORT_LIST_DATA))
  all_ID_annotations <- rep(list(list()), times = length(SHORT_LIST_DATA))
  
  for(n in seq_along(SHORT_LIST_DATA))
  {
    # Here we establish which mart(species) we are using in this given dataset based on "exp_species" vector
    ##### !!! THIS NEEDS TO BE MANUALLY ESTABLISHED BASE ON exp_species !!! ##### 
    usedMart_ <- set_mart_to_be_used(str_vector_of_species_names_ = str_vector_of_species_names[n], int_loop = n)
    
    potental_identifiers <- get_the_potental_identifiers(usedMart___ = usedMart_, int_loop_nb = n)
    
    
    
    message( 'Starting to annotate the data' )
    # Here we are annotating given datasets with data from all the relevant databases
    for (m in seq_along(potental_identifiers))
    {
      message( paste0('Annotating data for step ', n, ' and ', m, '...'))
      
      ANNOT_SHORT_LIST_DATA[[n]][[m]] <- biomaRt::getBM(
        attributes = c(potental_identifiers[[m]], "external_gene_name"),
        filters = potental_identifiers[[m]], 
        values = SHORT_LIST_DATA[[n]]$Probe_ID, 
        uniqueRows = T,
        mart = usedMart_
      )
      
      message( paste0('Data for step for ', n, ' and ', m, ' annotated') )
    }
    message( 'Data annotated' )
    
    usedMart_ <- NULL
    
    
    
    message( 'Starting all_ID_annotations step' )
    
    # Here we save number of annotations from each ID
    for(k in seq_along(ANNOT_SHORT_LIST_DATA[[n]])) 
    {
      message( paste0('Starting all_ID_annotations step for ', n, ' and ', k) )
      all_ID_annotations[[n]][[k]] <- length(ANNOT_SHORT_LIST_DATA[[n]][[k]][[1]])
      message( paste0('Competed all_ID_annotations step for ', n, ' and ', k) )
    } 
    
    message( 'Competed whole all_ID_annotations step' )
  }
  
  
  
  # Lists inside main list are changed into dfs (vectors) as 'which' function demands it
  df_all_ID_annotations <- lapply(all_ID_annotations, FUN = unlist)
  
  # Here we will be returing results of appropriate microarray search
  HIGHEST_HIT_LIST <- rep(list(list()), times = length(SHORT_LIST_DATA))
  
  # Here we are getting all of the highest yielding IDs
  for(n in seq_along(list_LIST_DATA_)){
    for(m in seq_along(which(df_all_ID_annotations[[n]] == max(df_all_ID_annotations[[n]])))){
      HIGHEST_HIT_LIST[[n]][[m]] <- ANNOT_SHORT_LIST_DATA[[n]][[(which(df_all_ID_annotations[[n]] == max(df_all_ID_annotations[[n]]))[m])]]
    } }
  
  rm(df_all_ID_annotations)
  
  
  proper_length_vector_for_checking_annotation_percentage <- get_proper_length_vector_for_checking_annotation_percentage(list_LIST_DATA__ = list_LIST_DATA_, int_Probe_IDs_to_test_ = int_Probe_IDs_to_test) 
  check_annotation_percentage <- purrr::map2(
    .x = HIGHEST_HIT_LIST, 
    .y = proper_length_vector_for_checking_annotation_percentage, 
    .f = function(.x, .y) {  (length(.x[[1]][[1]]) / .y) * 100 }
  )
  check_annotation_percentage <- data.frame(str_vector_of_experiment_ids, rlist::list.rbind(check_annotation_percentage))
  colnames(check_annotation_percentage) <- c('Exp_ID', 'highest_annotated_identifier_percentages')
  readr::write_tsv(check_annotation_percentage, paste0(str_experiment_name, '/', 'highest_annotated_identifier_percentages.tsv'))
  
  # Here we simply copy/establish names of features, that were highest by themselves. The conditions ask: 1) is there at least a single hit with highest number (I dont know if there can be 0 though...) 2) Is the first (and though each) highest hit list has at least single hit?
  NAMES_HIGHEST_HIT_LIST <- lapply(X = HIGHEST_HIT_LIST, FUN = set_0_hit_annotations_to_na)
  NAMES_HIGHEST_HIT_LIST <- data.frame(str_vector_of_experiment_ids, rlist::list.rbind(NAMES_HIGHEST_HIT_LIST))
  colnames(NAMES_HIGHEST_HIT_LIST) <- c('Exp_ID', 'platform_to_use')
  readr::write_tsv(NAMES_HIGHEST_HIT_LIST, paste0(str_experiment_name, '/', 'platform_to_use_for_probes_based_analysis.tsv'))
  
  ###### Here we estalish correct lists to analyzed in further steps (currently - need to remove experiments with microarrays not captured in ensembl) ######
  ###### Yeah, i dont know what to do here
  WHICH_EXP_TO_ANAL <- seq_along(NAMES_HIGHEST_HIT_LIST[[1]])
  
  return(list(WHICH_EXP_TO_ANAL, NAMES_HIGHEST_HIT_LIST, HIGHEST_HIT_LIST, ANNOT_SHORT_LIST_DATA, all_ID_annotations))
}



### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database. We can query ncbi databases only 3 times per second - that is the reason for sys.sleep time! A
search_for_ids_in_ncbi <- function(str_vector_of_ids)
{
  ids <- list()
  counter <- 1
  for (id_ in str_vector_of_ids)
  {
    print(id_)
    print(counter)
    
    ids[[counter]] <- rentrez::entrez_search(db="gene", term = id_)
    
    counter = counter + 1
    
    Sys.sleep(0.4)
  }
  return(ids)
}



#This function should except output of search_for_ids function. BEWARE - this only returns first geneID found
make_and_write_table_with_original_and_ncbi_ids <- function(df_returned_by_entrez_gene_search, df_original_data, str_name_of_the_file = 'generic.tsv', experiment_directory_name = '.')
{
  temp_list <- list()
  for (n in seq(length(df_returned_by_entrez_gene_search)))
  {
    temp_list[[n]] <- ''
    
    if (length(df_returned_by_entrez_gene_search[[n]]$ids) != 0)
    {
      temp_list[[n]][[1]] <- df_returned_by_entrez_gene_search[[n]]$ids[[1]]
    } 
    else
    {
      temp_list[[n]][[1]] <- 'none'
    }
    temp_name <- df_returned_by_entrez_gene_search[[n]]$QueryTranslation
    
    
    temp_list[[n]][[2]] <- gsub(
      pattern = "(\\[All Fields\\])|\\)|\\(", 
      replacement = '', 
      x = df_returned_by_entrez_gene_search[[n]]$QueryTranslation)
  }
  
  
  temp_df <- as.data.frame(rlist::list.rbind(temp_list))
  colnames(temp_df) <- c('ncbi_response', 'Probe_ID')
  
  temp_df <- merge(x = df_original_data, y = temp_df, by = 'Probe_ID', all.x = T)
  temp_df <- unique(temp_df)
  
  readr::write_tsv(temp_df, paste0(experiment_directory_name, '/', str_name_of_the_file))
  
  return(temp_df)
}



# str_platforms_ids_to_download are in GEO format
download_platforms_from_gemma <- function(str_platforms_ids_to_download)
{
  username_ <- readline(prompt = "Gimme Your GEMMA username: ")
  password_ <- getPass::getPass(msg = "Gimme Your GEMMA password: ")
  
  gemmaAPI::setGemmaUser(username = username_, password = password_)
  
  temp_list <- list()
  
  temp_list <- lapply(X = str_platforms_ids_to_download, FUN = function(X)
  {
    gemmaAPI::platformInfo(platform = X, 
                           request = 'annotations')
  })
  
  gemmaAPI::setGemmaUser()
  
  return(temp_list)
}



# list_gemma_platforms is list of annotations for platforms in list_LIST_DATA downloaded from gemma: a result from download_platforms_from_gemma function
get_gemma_annotations_for_data <- function(list_LIST_DATA, list_gemma_platforms)
{
  library(dplyr)
  temp_list <-
    purrr::map2(
      .x = list_LIST_DATA,
      .y = list_gemma_platforms,
      .f = function(.x, .y)
      {
        merge(
          x = .x,
          y = .y,
          by.x = 'Probe_ID',
          by.y = 'ProbeName',
          all.x = T
        ) %>%
          dplyr::select(Probe_ID:GeneSymbols)
      }
    )
  
  return(temp_list)
}



write_lists <- function(list_LIST_DATA, str_experiment_name, str_description)
{
  temp_df <- rlist::list.rbind(list_LIST_DATA)
  readr::write_tsv(temp_df, paste0(str_experiment_name, 'table_', str_description, '.tsv'))
}



merge_and_remove_nas_from_list_post_annotation <- function(list_annotated_LIST_DATA, str_name_of_column_containing_annotated_gene_symbols)
{
  temp_df <- rlist::list.rbind(list_annotated_LIST_DATA)
  temp_df <- subset(x = temp_df, subset = !is.na(temp_df[[str_name_of_column_containing_annotated_gene_symbols]]))
  return(temp_df)
}



annotate_now <- function(list_LIST_DATA_ = LIST_DATA, str_identifier_type_, str_vector_of_species_names__, experiment_name_)
{
  WHICH_EXP_TO_ANAL <- seq(length(list_LIST_DATA_))
  ANNOT_LIST_DATA <- list()
  
  if (length(str_identifier_type_) != 1)
  {
    identifiers_used_for_annotation <- str_identifier_type_
  }
  else
  {
    identifiers_used_for_annotation <-
      set_identifiers_used_for_annotation_if_not_probeID(str_identifier_type = str_identifier_type_, list_LIST_DATA = list_LIST_DATA_)
  }
  
  for (n in WHICH_EXP_TO_ANAL)
  {
    usedMart_ = set_mart_to_be_used(str_vector_of_species_names_ = str_vector_of_species_names__[n],
                                    int_loop = n)
    
    message(paste0(
      'Annotating experiment ',
      names(list_LIST_DATA_[n]),
      ' in step ',
      n,
      '...'
    ))
    ANNOT_LIST_DATA[[n]] <- biomaRt::getBM(
      attributes = c(identifiers_used_for_annotation[[n]], "external_gene_name"),
      filters = identifiers_used_for_annotation[[n]],
      values = list_LIST_DATA_[[n]]$Probe_ID,
      uniqueRows = F,
      mart = usedMart_
    )
  }
  
  FINAL_ANNOT_LIST_DATA <-
    purrr::pmap(
      .l = list(
        list_LIST_DATA_,
        ANNOT_LIST_DATA,
        identifiers_used_for_annotation
      ),
      .f = function(.x, .y, .z)
      {
        merge(
          x = .x,
          y = .y,
          by = 'Probe_ID',
          by.y = .z,
          all.x = T
        )
      }
    ) ### !!! This does not seem to actually save all rows in original data?
  
  DF_FINAL_ANNOT_LIST_DATA <-
    rlist::list.rbind(FINAL_ANNOT_LIST_DATA)
  
  ### !!! ADD A LINE WHERE RAW ANNOTATION DATA ARE PRINTED
  
  # This is the correct way to uniqualize the resulting dataframes
  # DF_FINAL_ANNOT_LIST_DATA <- DF_FINAL_ANNOT_LIST_DATA[unique(DF_FINAL_ANNOT_LIST_DATA$Nb),]
  # Additionally Probe_ID annotation returns many duplicated rows. I am not sure why. Here we remove them
  DF_FINAL_ANNOT_LIST_DATA <- unique(DF_FINAL_ANNOT_LIST_DATA)
  
  DF_FINAL_ANNOT_LIST_DATA <-
    collapse_annotated_names_for_given_probe(DF_FINAL_ANNOT_LIST_DATA, list_LIST_DATA__ = list_LIST_DATA_)
  
  readr::write_tsv(
    DF_FINAL_ANNOT_LIST_DATA,
    paste0(
      experiment_name_,
      'annotated_data_from_',
      identifiers_used_for_annotation[1],
      '.tsv'
    )
  )
  
  return(DF_FINAL_ANNOT_LIST_DATA)
}



set_experiment_name_and_create_directory_for_output <- function(str_identifier_type___, backup_experiment_name_)
{
  if (length(str_identifier_type___) != 1)
  {
    temp_experiment_directory_name = paste0(backup_experiment_name_, '/')
    temp_str_identifier_name <- backup_experiment_name_
  }
  else
  {
    temp_experiment_directory_name <-
      paste0(str_identifier_type___, '/')
    temp_str_identifier_name <- str_identifier_type___
  }
  
  print(
    paste0(
      'Trying to create directory "',
      temp_experiment_directory_name,
      '". Ignore warning that the directory exists. It does not interfere with its creation. No need for if statements here.'
    )
  )
  dir.create(temp_experiment_directory_name, temp_str_identifier_name)
  
  temp_both_names <-
    list(temp_experiment_directory_name, temp_str_identifier_name)
  
  return(temp_both_names)
}



# We can pass two types of data into str_identifier_name: 1) one-element string vector, which is the same name as filter name in biomartr-ensembl database. 2) vector of strings with filter name for each experiment dataset in the list of experiments
master_annotator_for_known_identfiers <- function(descriptions_, str_filename_, str_identifier_type__, backup_experiment_name)
{
  # We use %>% operator somewhere in this, or downstream function - I think we do not need this, as we defined `%>%` <- dplyr::`%>%` before
  # library(dplyr)
  
  list_experiment_directory_name_and_identifier_type <-
    set_experiment_name_and_create_directory_for_output(str_identifier_type__, backup_experiment_name)
  
  
  PRE_DATA__ <-
    read_preformated_data(
      str_filename = str_filename_,
      int_numbers_are_from = 6,
      int_numbers_are_to = 8,
      col_types_ = 'nncccccc'
    )
  
  LIST_DATA__ <- split(PRE_DATA__, f = PRE_DATA__$Paper)
  
  readr::write_tsv(
    rlist::list.rbind(LIST_DATA__),
    paste0(
      list_experiment_directory_name_and_identifier_type[[1]],
      'input_for_',
      list_experiment_directory_name_and_identifier_type[[2]],
      '.tsv'
    )
  )
  
  # List of species names based on data in description files. This is a file we will be working on. Species is in format: small letters, english plural of species
  exp_species__ <- descriptions_ %>%
    dplyr::select("Paper_ID", "Species") %>%
    unique() %>%
    dplyr::filter(Paper_ID %in% unique(PRE_DATA__$Paper))
  
  write_lenghts_of_list_objects(
    LIST_DATA__,
    paste0(
      list_experiment_directory_name_and_identifier_type[[1]],
      'list_data_lenghts_probes.tsv'
    )
  )
  readr::write_tsv(
    exp_species__,
    paste0(
      list_experiment_directory_name_and_identifier_type[[1]],
      'exp_species_used_for_testing_which_platform_to_use_probes.tsv'
    )
  )
  
  annotation <-
    annotate_now(
      list_LIST_DATA_ = LIST_DATA__,
      str_identifier_type_ = str_identifier_type__,
      str_vector_of_species_names__ = exp_species__$Species,
      experiment_name_ = list_experiment_directory_name_and_identifier_type[[1]]
    )
  
  return(annotation)
}



annotate_identifiers_to_geneID <- function(str_filename_, str_experiment_name, descriptions_, bool_reformat_names = F)
{
  data <- read_preformated_data(str_filename = str_filename_, int_numbers_are_from = 6, int_numbers_are_to = 8, col_types_ = 'nncccccc')
  
  directory_name <- paste0(str_experiment_name, '/')
  
  dir.create(directory_name)
  
  strvec_ncbi_query_for_identifers <- search_for_ids_in_ncbi(as.vector(data$Probe_ID))
  
  if(bool_reformat_names == T)
  {
    print( paste0('Do not use this option. In its current implementation it does the opposite of helping. Moreover, the data I have is proper and does not need reformating') )
    # exp_species <- descriptions_ %>%
    #   select("Paper_ID", "Species") %>%
    #   unique() %>%
    #   filter(Paper_ID %in% unique(data$Paper))
    # 
    # readr::write_tsv(x = exp_species, path = paste0(directory_name, str_experiment_name, '_species_used_for_analysis.tsv'))
    # 
    # list_data <- purrr::map2(.x = split(data, f = data$Paper), .y = exp_species$Species, .f = ~ change_name_to_proper_format_given_the_species(.x, .y))
    # 
    # data <- rlist::list.rbind(list_data)
    # 
    # rm(list_data)
  }
  
  ### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database  
  annotated_data <- make_and_write_table_with_original_and_ncbi_ids(df_returned_by_entrez_gene_search = strvec_ncbi_query_for_identifers, df_original_data = data, str_name_of_the_file = paste0(str_experiment_name, '.tsv'), experiment_directory_name = directory_name)
  
  return(annotated_data)
}



gather_all_datasets_into_single_df <- function(regex_pattern_to_find_datasets_with = '^annotations.*')
{
  dataset_names <-
    as.list(parse(
      text = ls(pattern = regex_pattern_to_find_datasets_with, name = globalenv())
    ))
  
  datasets_list <- do.call(what = 'list', args = dataset_names)
  
  datasets_df <- rlist::list.rbind(datasets_list)
  
  return(datasets_df)
}



collapse_annotated_names_for_given_probe <- function(df_annotated, list_LIST_DATA__)
{
  temp_list_annotated <-
    split(x = df_annotated, f = df_annotated$Experiment)
  
  temp_bind_list_annotated <-
    aggregate(external_gene_name ~ Probe_ID,
              data = df_annotated,
              FUN = stringr::str_c)
  
  # This part collapses multiple same values to single value and saves it as temp column
  temp_bind_list_annotated$temp_gene_name <-
    as.character(purrr::map(
      .x = temp_bind_list_annotated$external_gene_name,
      .f = function(x) {
        unique(x)
      }
    ))
  
  # Check if the value for given row is proper char, or leftover list ' c("blabla") '
  was_the_value_uniqued <-
    stringr::str_detect(string = temp_bind_list_annotated$temp_gene_name, pattern = 'c\\(".*')
  
  temp_bind_list_annotated$bullshit_list_derived_strings <-
    as.character(
      purrr::map_if(
        .x = temp_bind_list_annotated$external_gene_name,
        .p = was_the_value_uniqued,
        .f = function(.x) {
          paste(.x, collapse = "; ")
        },
        .else = function(.x) {
          return('')
        }
      )
    )
  
  temp_bind_list_annotated$proper_strings <-
    as.character(
      purrr::map_if(
        .x = temp_bind_list_annotated$temp_gene_name,
        .p = !was_the_value_uniqued,
        .f = function(.x) {
          return(.x)
        },
        .else = function(.x) {
          return('')
        }
      )
    )
  
  temp_bind_list_annotated$external_gene_name <-
    paste0(
      temp_bind_list_annotated$bullshit_list_derived_strings,
      temp_bind_list_annotated$proper_strings
    )
  
  temp_bind_list_annotated <-
    dplyr::select(.data = temp_bind_list_annotated, Probe_ID, external_gene_name)
  
  bind_list_LIST_DATA__ <- rlist::list.rbind(list_LIST_DATA__)
  
  merged_bind_list_annotated <-
    merge(x = bind_list_LIST_DATA__,
          y = temp_bind_list_annotated,
          by = 'Probe_ID',
          all.x = T)
  
  merged_bind_list_annotated <- unique(merged_bind_list_annotated)
  
  return(merged_bind_list_annotated)
}



do_directions_for_multiple_gene_instances_within_experiment_match_old <- function(df_merged_dataset_)
{
  ### This will be our input data for further analysis. It contains only significant columns of "Paper", "GroupID", "ensembl_gene_name", "logFC":
  SHORT_SINGLE_TEST_ANNOTATION <- df_merged_dataset_ %>%
    dplyr::select("Paper", "Experiment", "external_gene_name", "logFC") %>%
    dplyr::group_nest(Paper, Experiment, external_gene_name)
  
  # I dont know why purrr::map didnt work for this. Either way, we just lapply this. Here we write UP or DOWN for every logFC for given genename in given experiment (not paper)
  SHORT_SINGLE_TEST_ANNOTATION$directionality <-
    lapply(
      X = SHORT_SINGLE_TEST_ANNOTATION$data,
      FUN = function(x) {
        dplyr::mutate(x, Symbol_direction = ifelse(logFC > 0, "UP", "DOWN"))
      }
    )
  
  # Here we count how many UPs and DOWNs are in every given genename in given experiment (not paper)
  SHORT_SINGLE_TEST_ANNOTATION$directionality <-
    lapply(
      X = SHORT_SINGLE_TEST_ANNOTATION$directionality,
      FUN = function (x) {
        table(dplyr::select(x, "Symbol_direction"))
      }
    )
  
  # Here we establish actual status of gene in given experiment
  SHORT_SINGLE_TEST_ANNOTATION$sum_directionality <-
    as.character(lapply(
      X = SHORT_SINGLE_TEST_ANNOTATION$directionality,
      FUN = function(x) {
        if (length(x) == 1 && grepl(pattern = "UP", names(x))) {
          "UP"
        }
        else if (length(x) == 1 &&
                 grepl(pattern = "DOWN", names(x))) {
          "DOWN"
        }
        else if (length(x) == 2 &&
                 grepl(pattern = "DOWN", names(x)) &&
                 grepl(pattern = "DOWN", names(x))) {
          "MIXED"
        }
        else {
          "ERROR"
        }
      }
    ))
  SHORT_SINGLE_TEST_ANNOTATION$directionality <- NULL
  
  
  UNNEST_SHORT_SINGLE_TEST_ANNOTATION <-
    SHORT_SINGLE_TEST_ANNOTATION %>%
    tidyr::unnest()
  
  return(UNNEST_SHORT_SINGLE_TEST_ANNOTATION)
  
  # # Here we remove MIXED expression and multiple genenames
  # FILT_SHORT_SINGLE_TEST_ANNOTATION <- SHORT_SINGLE_TEST_ANNOTATION %>%
  #   #filter(!grepl(pattern = "(.*);(.*)", SHORT_SINGLE_TEST_ANNOTATION$ensembl_gene_name)) %>%
  #   dplyr::filter(!sum_directionality == "MIXED") %>%
  #   dplyr::filter(!sum_directionality == "ERROR") %>%
  #   dplyr::filter(!is.na(ensembl_gene_name))
  # 
  # # Here we add mean to each !!! Perhaps this should be final input, because it has: 1) removed genes giving UP+DOWN in the same experiment, 2) removed NAs
  # STDINPUT_FILT_SHORT_SIN_T_ANNO <- FILT_SHORT_SINGLE_TEST_ANNOTATION %>%
  #   mutate(mean = as.numeric(map(data, function(x) { as.numeric(mean(x[[1]]))} ))) %>% ## This way we give actuall vector to function, not a data table(tibble)
  #   select("Paper", "GroupID", "ensembl_gene_name", "sum_directionality", "mean") 
}
