`%>%` <- dplyr::`%>%`
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')


kill_corrupted_e_notation <- function(df_temp_data_, str_to_substitute_corrupted_data_with)
{
  df_temp_data_$corrupted <-
    stringr::str_detect(string = df_temp_data_$logFC, pattern = 'e')
  
  # Find cells with e-annotation
  df_temp_data_ <- dplyr::mutate(
    .data = df_temp_data_,
    value = dplyr::if_else(
      condition = corrupted,
      true = gsub(
        pattern = '(.*)e\\Q+\\E',
        replacement = '',
        x = logFC
      ),
      false = '0',
      missing = '0'
    )
  )
  
  df_temp_data_$value <- as.numeric(df_temp_data_$value)
  
  # Set all numbers with e number higher than 2 to constant value
  df_temp_data_ <- dplyr::mutate(
    .data = df_temp_data_,
    logFC = dplyr::if_else(
      condition = value > 2,
      true = str_to_substitute_corrupted_data_with,
      false = logFC,
      missing = logFC
    )
  )
  
  df_temp_data_$corrupted <- NULL
  df_temp_data_$value <- NULL
  
  return(df_temp_data_)
}



# INPUT: col_types_ - add floats as 'c' here. This is because input floats that can contain either . or , as decimal. Later they are converted to numeric based on int_numbers_are_ arguments. Needs at least columns: Paper, Experiment, Probe_Id logFC
read_preformated_data <- function(str_filename, col_types_ = 'nncccccc', int_numbers_are_from = 6, int_numbers_are_to = 8, str_substitute_inf_with = '15')
{
  temp_data <- readr::read_tsv(str_filename, col_types = col_types_)

  # Replace , with . to make number actual R numerics
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(
    X = temp_data[int_numbers_are_from:int_numbers_are_to],
    FUN = function(x)
    {
      stringr::str_replace(string = x,
                           pattern = ",",
                           replacement = ".")
    }
  ) ### COOL CONSTRUCT
  
  # Replace inf values with constant number
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(
    X = temp_data[int_numbers_are_from:int_numbers_are_to],
    FUN = function(x) {
      gsub(pattern = "^[I|i]nf(.*)",
           replacement = str_substitute_inf_with,
           x = x)
    }
  )
  
  # Replace too large e numbers with constant number
  temp_data <-
    kill_corrupted_e_notation(df_temp_data_ = temp_data,
                              str_to_substitute_corrupted_data_with = str_substitute_inf_with)
  
  # Change number columns to numeric type
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(temp_data[int_numbers_are_from:int_numbers_are_to],
                                                               as.numeric)
  
  return(temp_data)
}



# INPUT: list. Defaults to writing lenghts of object on first level of depth of the list
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
get_the_highest_hit_returning_id_type <- function(descriptions_ = descriptions, int_Probe_IDs_to_test = 200, str_experiment_name = experiment_name, str_filename_, col_types__ = 'ncccccccccccccncc', int_numbers_are_from_ = 11, int_numbers_are_to_ = 13, str_substitute_inf_with_ = '15')
{
  PRE_DATA <- read_preformated_data(
    str_filename = str_filename_,
    int_numbers_are_from = int_numbers_are_from_,
    int_numbers_are_to = int_numbers_are_to_,
    col_types_ = col_types__,
    str_substitute_inf_with_ = str_substitute_inf_with__
  )
  
  
  exp_species <- descriptions_ %>%
    dplyr::select("Paper_ID", "Species") %>%
    unique() %>%
    dplyr::filter(Paper_ID %in% unique(PRE_DATA$Paper))
  
  exp_species$Paper_ID <- as.integer(exp_species$Paper_ID) ### !!! Added
  exp_species <- exp_species[order(exp_species$Paper_ID),] ### !!! Added
  
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



### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database. INPUT: delay_between_queries - 0.5 produces errors even. 3 times per second my asshole. We can query ncbi databases only 3 times per second - that is the reason for sys.sleep time! A
search_for_ids_in_ncbi <- function(str_vector_of_ids, delay_between_queries = 0.5) ### !!! We should rename this variable, because this is actually query vector, and not id vector per se. Query vector includes id, but also includes organism (and other, perhaps?)
{
  ids <- list()
  counter <- 1
  for (id_ in str_vector_of_ids)
  {
    done <- F
    
    while (done == F) {
      print(id_)
      print(counter)
      
      tryCatch( {
          ids[[counter]] <- rentrez::entrez_search(db="gene", term = id_)
          
          done <- T
          
          counter = counter + 1
          
          Sys.sleep(delay_between_queries)
        },
        error=function(e) {
          Sys.sleep(5)
          print('Failed to download data, trying again in 5 s')
        })
    }

  }
  return(ids)
}




# search_for_ids_in_ncbi <- function(str_vector_of_ids, delay_between_queries = 0.5) ### !!! We should rename this variable, because this is actually query vector, and not id vector per se. Query vector includes id, but also includes organism (and other, perhaps?)
# {
#   ids <- list()
#   counter <- 1
#   for (id_ in str_vector_of_ids)
#   {
#     print(id_)
#     print(counter)
#     
#     ids[[counter]] <- rentrez::entrez_search(db="gene", term = id_)
#     
#     counter = counter + 1
#     
#     Sys.sleep(delay_between_queries)
#   }
#   return(ids)
# }


#This function should except output of search_for_ids function. BEWARE - this only returns first geneID found INPUT: ### !!! df_returned_by_entrez_gene_search is actually a list!!!
make_and_write_table_with_original_and_ncbi_ids <- function(df_returned_by_entrez_gene_search, df_original_data, str_name_of_the_file = 'generic.tsv', experiment_directory_name = '.', str_to_cut_from_ncbi_response_to_form_back_Probe_ID, vector_of_species_names_used_, exp_species___, possible_names_for_species__)
{
  temp_list <- list()
  for (n in seq(length(df_returned_by_entrez_gene_search)))
  {
    temp_list[[n]] <- ''
    
    if (length(df_returned_by_entrez_gene_search[[n]]$ids) != 0)
    {
      temp_list[[n]][[1]] <-
        df_returned_by_entrez_gene_search[[n]]$ids[[1]]
    }
    else
    {
      temp_list[[n]][[1]] <- 'none'
    }
    temp_name <-
      df_returned_by_entrez_gene_search[[n]]$QueryTranslation
    
    
    temp_list[[n]][[2]] <- gsub(
      pattern = "(\\[All Fields\\])|\\)|\\(",
      replacement = '',
      x = df_returned_by_entrez_gene_search[[n]]$QueryTranslation
    )
  }
  
  
  temp_df <- as.data.frame(rlist::list.rbind(temp_list))
  colnames(temp_df) <- c('external_gene_name', 'Probe_ID_2')
  
  regex_str_to_cut_from_ncbi_response_to_form_back_Probe_ID <-
    paste0('(',
           Hmisc::escapeRegex(string = str_to_cut_from_ncbi_response_to_form_back_Probe_ID),
           ')(.*)')
 
  # species_for_this_entry <- vector_of_species_names_used_[stringr::str_detect(string = temp_df$Probe_ID, pattern = vector_of_species_names_used_)]
  
  temp_df$Probe_ID_2 <-
    stringr::str_remove(string = temp_df$Probe_ID_2, pattern = regex_str_to_cut_from_ncbi_response_to_form_back_Probe_ID) 
  # temp_df$Species <- stringr::str_remove_all(string = species_for_this_entry, pattern = '"')
  
  ######### HERE WE NEED TO SPLIT ORIGINAL DF INTO SPECIES, SPLIT temp_df INTO SPECIES, MERGE SPECIES LISTS AND THEM rbind the resulting lists ############# !!
  # exp_species___, possible_names_for_species__
  
  temp_df2 <- cbind(df_original_data, temp_df)
  
  # temp_df <-
  #   merge(x = df_original_data,
  #         y = temp_df,
  #         by = 'Probe_ID',
  #         all.x = T)
  # temp_df <- unique(temp_df)
  
  # readr::write_tsv(temp_df,
  #                  paste0(experiment_directory_name, '/', str_name_of_the_file))
  
  return(temp_df2)
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
get_gemma_annotations_for_data <- function(list_gemma_platforms, str_filename_, col_types__ = 'ncccccccccccccncc', int_numbers_are_from_ = 11, int_numbers_are_to_ = 13, str_substitute_inf_with_ = '15')
{
  df_data <- read_preformated_data(
    str_filename = str_filename_,
    int_numbers_are_from = int_numbers_are_from_,
    int_numbers_are_to = int_numbers_are_to_,
    col_types_ = col_types__,
    str_substitute_inf_with = str_substitute_inf_with_
  )

  df_data$Paper <- as.integer(df_data$Paper) ### !!! Added
  df_data <- df_data[order(df_data$Paper),]
  
  list_data <- split(x = df_data, f = df_data$Paper)
  
  temp_list <-
    purrr::map2(
      .x = list_data,
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
  
  df_return <- rlist::list.rbind(temp_list)
  
  df_return <- df_return %>%
    dplyr::rename(external_gene_name = GeneSymbols) %>%
    dplyr::select(Paper, Experiment, `GEO ID`, Probe_ID, dplyr::everything())
  
  return(df_return)
}



write_lists <- function(list_LIST_DATA, str_experiment_name, str_description)
{
  temp_df <- rlist::list.rbind(list_LIST_DATA)
  readr::write_tsv(temp_df, paste0(str_experiment_name, 'table_', str_description, '.tsv'))
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
    )
  
  DF_FINAL_ANNOT_LIST_DATA <-
    rlist::list.rbind(FINAL_ANNOT_LIST_DATA)
  
  readr::write_tsv(
    DF_FINAL_ANNOT_LIST_DATA,
    paste0(
      experiment_name_,
      'raw_annotated_data_from_',
      identifiers_used_for_annotation[1],
      '.tsv'
    )
  )
  
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



set_experiment_name_and_create_directory_for_output <- function(str_identifier_type___, str_experiment_name_)
{
  if (str_experiment_name_ != '')
  {
    temp_experiment_directory_name = paste0(str_experiment_name_, '/')
    temp_str_identifier_name <- str_experiment_name_
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



# set_experiment_name_and_create_directory_for_output <- function(str_identifier_type___, backup_experiment_name_)
# {
#   if (length(str_identifier_type___) != 1)
#   {
#     temp_experiment_directory_name = paste0(backup_experiment_name_, '/')
#     temp_str_identifier_name <- backup_experiment_name_
#   }
#   else
#   {
#     temp_experiment_directory_name <-
#       paste0(str_identifier_type___, '/')
#     temp_str_identifier_name <- str_identifier_type___
#   }
#   
#   print(
#     paste0(
#       'Trying to create directory "',
#       temp_experiment_directory_name,
#       '". Ignore warning that the directory exists. It does not interfere with its creation. No need for if statements here.'
#     )
#   )
#   dir.create(temp_experiment_directory_name, temp_str_identifier_name)
#   
#   temp_both_names <-
#     list(temp_experiment_directory_name, temp_str_identifier_name)
#   
#   return(temp_both_names)
# }



# We can pass two types of data into str_identifier_name: 1) one-element string vector, which is the same name as filter name in biomartr-ensembl database. 2) vector of strings with filter name for each experiment dataset in the list of experiments
master_annotator_for_known_identfiers <- function(descriptions_, str_identifier_type__, str_experiment_name = '', str_filename_, col_types__ = 'ncccccccccccccncc', int_numbers_are_from_ = 11, int_numbers_are_to_ = 13, str_substitute_inf_with_ = '15')
{
  # We use %>% operator somewhere in this, or downstream function - I think we do not need this, as we defined `%>%` <- dplyr::`%>%` before
  # library(dplyr)
  
  list_experiment_directory_name_and_identifier_type <-
    set_experiment_name_and_create_directory_for_output(str_identifier_type__, str_experiment_name)
  
  
  PRE_DATA__ <-
    read_preformated_data(
      str_filename = str_filename_,
      int_numbers_are_from = int_numbers_are_from_,
      int_numbers_are_to = int_numbers_are_to_,
      col_types_ = col_types__,
      str_substitute_inf_with = str_substitute_inf_with_
    )
  
  PRE_DATA__$Paper <- as.integer(PRE_DATA__$Paper) ### !!! Added
  PRE_DATA__ <- PRE_DATA__[order(PRE_DATA__$Paper),] ### !!! Added
  
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
  
  exp_species__$Paper_ID <- as.integer(exp_species__$Paper_ID) ### !!! Added
  exp_species__ <- exp_species__[order(exp_species__$Paper_ID),] ### !!! Added
  
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




# INPUT: two corresponding char vectors, char_vec_gene_id_to_query_with - gene names to query ncbi with, char_vec_organism - species name for every gene name. OUTPUT: List. [[1]] - actual queries to be sent to ncbi; [[2]] - vector added to Probe_ID (gene name). Without it, it is difficult for next function to return probeid/genename to its original form, which is neccesary for merging queries with original dataframe
get_query_for_ncbi_geneID_annotation <- function(char_vec_gene_id_to_query_with, char_vec_organism, chr_gene_identifier = 'Gene name') # here ive added default value to chr_gene_identifier
{
  the_species_name_vector <- c('"mus musculus"', '"rattus norvegicus"', '"homo sapiens"', 'saimiri')
  
  possible_names_for_species <- list(mouse_names = c('mus', 'mouse', 'mice'), 
                                     rat_names = c('rattus', 'rat', 'rats'), 
                                     human_names = c('homo', 'human', 'humans'), 
                                     saimiri_names = c('saimiri', 'squirrel monkey', 'squirrel monkeys', 'squirrelmonkeys'))
  
  species_vector <-
    as.character(lapply(
      X = tolower(char_vec_organism),
      FUN = function(x) {
        if (x %in% possible_names_for_species$mouse_names) {
          return(the_species_name_vector[1])
        }
        else if (x %in% possible_names_for_species$rat_names) {
          return(the_species_name_vector[2])
        }
        else if (x %in% possible_names_for_species$human_names) {
          return(the_species_name_vector[3])
        }
        else if (x %in% possible_names_for_species$saimiri_names) {
          return(the_species_name_vector[4])
        }
        else{
          stop('I did not recognize species name. You probably have to ret into this function and add Your species names manually')
        }
      }
    ))
  
  ################### WORKING ####################
  linker_string_crucial_for_returning_ProbeID_to_proper_form <- paste0('[', chr_gene_identifier, '] AND ')
  
  # linker_string_crucial_for_returning_ProbeID_to_proper_form <- '[Gene Name] AND '
  ################### WORKING ####################
  
  query_vector <- paste0(char_vec_gene_id_to_query_with, linker_string_crucial_for_returning_ProbeID_to_proper_form, species_vector, '[Organism]')

  return_ <- list(query_vector, linker_string_crucial_for_returning_ProbeID_to_proper_form, the_species_name_vector, possible_names_for_species)
  
  return(return_)
}


# INPUT: chr_gene_identifier_ - this needs to be proper ncbi search class, as can be found using advanced search in ncbi. There should be no need for it to be any different that 'Gene name'
annotate_identifiers_to_geneID <- function(str_filename_, str_experiment_name, descriptions_ = descriptions, chr_gene_identifier_ = 'Gene name', col_types__ = 'ncccccccccccccncc', int_numbers_are_from_ = 11, int_numbers_are_to_ = 13, str_substitute_inf_with_ = '15')
{
  data <- read_preformated_data(
    str_filename = str_filename_,
    int_numbers_are_from = int_numbers_are_from_,
    int_numbers_are_to = int_numbers_are_to_,
    col_types_ = col_types__,
    str_substitute_inf_with = str_substitute_inf_with_
  )

  data$Paper <- as.integer(data$Paper) ### !!! Ive added this
  data <- data[order(data$Paper),]
  
  directory_name <- paste0(str_experiment_name, '/')
  
  dir.create(directory_name)

  # Number of entires for [Organism] name: 102702 rattus, rat, rats // 273852 mouse, mus, mice // 224903 homo, humans // 224866 human // 30931 squirrel monkeys, Saimiri // 30905 saimiri boliviensis
  # List of species names based on data in description files. This is a file we will be working on. Species is in format: small letters, english plural of species
  exp_species__ <- descriptions_ %>%
    dplyr::select("Paper_ID", "Species") %>%
    unique() %>%
    dplyr::filter(Paper_ID %in% unique(data$Paper))
  
  exp_species__$Paper_ID <- as.integer(exp_species__$Paper_ID) ### !!! Ive added this
  exp_species__ <- exp_species__[order(exp_species__$Paper_ID),] ### !!! Ive added this 
  
  
  temp_data <- merge(x = data, y = exp_species__, by.x = 'Paper', by.y = 'Paper_ID')

  list_query_for_ncbi <- get_query_for_ncbi_geneID_annotation(char_vec_gene_id_to_query_with = temp_data$Probe_ID, char_vec_organism = temp_data$Species, chr_gene_identifier = chr_gene_identifier_)
  
  query_for_ncbi <- list_query_for_ncbi[[1]]
  
  write.table(x = query_for_ncbi, file = 'annotate_identifiers_to_geneID_query.txt') ### !!! this is added
  
  linker_string_crucial_for_returning_ProbeID_to_proper_form_ <- list_query_for_ncbi[[2]]
  
  vector_of_species_names_used <- list_query_for_ncbi[[3]]
  
  possible_names_for_species_ <- list_query_for_ncbi[[4]]

  strvec_ncbi_query_for_identifers <- search_for_ids_in_ncbi(query_for_ncbi)
  
  ### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database  
  annotated_data <- make_and_write_table_with_original_and_ncbi_ids(df_returned_by_entrez_gene_search = strvec_ncbi_query_for_identifers, df_original_data = data, str_name_of_the_file = paste0(str_experiment_name, '.tsv'), experiment_directory_name = directory_name, str_to_cut_from_ncbi_response_to_form_back_Probe_ID = linker_string_crucial_for_returning_ProbeID_to_proper_form_, vector_of_species_names_used_ = vector_of_species_names_used, exp_species___ = exp_species__, possible_names_for_species__ = possible_names_for_species_)
  
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
  df_annotated <-
    aggregate(external_gene_name ~ Probe_ID,
              data = df_annotated,
              FUN = stringr::str_c)
  
  # This part collapses multiple same values to single value and saves it as temp column
  df_annotated$temp_gene_name <-
    as.character(lapply(
      X = df_annotated$external_gene_name,
      FUN = function(x) {
        paste(x, collapse = "; ")
      }
    ))

  df_annotated$external_gene_name <- change_vector_of_mixed_normal_and_c_geneNames_into_unique_geneName(df_annotated$temp_gene_name)
  
  df_annotated <-
    dplyr::select(.data = df_annotated, Probe_ID, external_gene_name)
  
  bind_list_LIST_DATA__ <- rlist::list.rbind(list_LIST_DATA__)
  
  merged_df_annotated <-
    merge(x = bind_list_LIST_DATA__,
          y = df_annotated,
          by = 'Probe_ID',
          all.x = T)
  
  merged_df_annotated <- unique(merged_df_annotated)
  
  return(merged_df_annotated)
}


# Input: Vector of gene names, needs to be a proper char vector, with gene names separ
change_vector_of_mixed_normal_and_c_geneNames_into_unique_geneName <- function(char_vec, regex_pattern_to_split_with = '; ')
{
  # This part collapses multiple same values to single value and saves it as temp column
  temp_list <-
    lapply(
      X = char_vec,
      FUN = function(x) {
        temp <- as.character(stringr::str_split(
          string = x,
          pattern = regex_pattern_to_split_with,
          simplify = T
        ))
        
        temp <- unique(temp)
        
        return(temp)
      }
    )
  
  # Check if the value for given row is proper char, or leftover list ' c("blabla") '
  was_the_value_uniqued <-
    as.logical(lapply(
      X = temp_list,
      FUN = function(x) {
        if (length(x) == 1) {
          return(F)
        }
        else{
          return(T)
        }
      }
    ))

  bullshit_list_derived_strings <-
    as.character(
      purrr::map_if(
        .x = temp_list,
        .p = was_the_value_uniqued,
        .f = function(x) {
          paste(x, collapse = regex_pattern_to_split_with)
        },
        .else = function(x) {
          return('')
        }
      )
    )
  
  proper_strings <-
    as.character(
      purrr::map_if(
        .x = temp_list,
        .p = !was_the_value_uniqued,
        .f = function(x) {
          return(x)
        },
        .else = function(x) {
          return('')
        }
      )
    )
  
  resulting_char_vector <-
    paste0(bullshit_list_derived_strings,
           proper_strings)
  
  return(resulting_char_vector)
}



# INPUT: single string containing multiple gene names separeted by separator
select_best_geneName_wrapper_for_single_string <- function(string_to_be_vectorised, separator = '; ')
{
  vectorised_string <- as.character(stringr::str_split(
    string = string_to_be_vectorised,
    pattern = separator,
    simplify = T
  ))
  
  the_best_name_ <- select_best_geneName(char_vec = vectorised_string)
  
  return(the_best_name_)
}



# INPUT: char_vec - includes all gene names returned by ensembl for given gene. Hence, the function should be iterated over list of vectors, each vector/list element for single gene. regex_to_detect_bad_names_with - genes can have 4 numbers: Olr1237 OUTPUT: the_best_name - char_vec of length 1
select_best_geneName <- function(char_vec, regex_to_detect_bad_names_with = '(\\d{5})|(^Gm\\d)', regex_to_detect_less_bad_names_with = '(Rik)|(LOC)|(Gm)', ignore_case_ = T)
{
  log_is_this_name_bad <-
    stringr::str_detect(string = char_vec, pattern = regex_to_detect_bad_names_with)
  
  good_names <- subset(x = char_vec, subset = !log_is_this_name_bad)
  bad_names <- subset(x = char_vec, subset = log_is_this_name_bad)
  
  log_is_this_name_less_bad <-
    stringr::str_detect(
      string = bad_names,
      pattern = stringr::regex(regex_to_detect_less_bad_names_with, ignore_case = ignore_case_)
    )
  
  less_bad_names <- subset(x = bad_names, subset = log_is_this_name_less_bad)
  
  if(length(good_names) != 0)
  {
    the_best_name <- good_names[[1]]
  }
  else if(length(less_bad_names) != 0)
  {
    the_best_name <- less_bad_names[[1]]
  }
  else if(length(bad_names) != 0)
  {
    the_best_name <- bad_names[[1]]
  }
  else
  {
    the_best_name <- NA
  }
  
  return(the_best_name)
}




master_annotator_for_known_identfiers_wrapper_for_using_dfs_as_input <- function(descriptions__ = descriptions, df_to_annotate, str_identifier_type___, str_experiment_name__ = '', col_types___, int_numbers_are_from__, int_numbers_are_to__, str_substitute_inf_with__ = '15')
{
  temp_file <- tempfile()
  
  readr::write_tsv(x = df_to_annotate, path = temp_file)
  
  annotations_ <-
    master_annotator_for_known_identfiers(
      descriptions_ = descriptions__,
      str_filename_ = temp_file,
      str_identifier_type__ = str_identifier_type___,
      str_experiment_name = str_experiment_name__,
      col_types__ = col_types___, 
      int_numbers_are_from_ = int_numbers_are_from__, 
      int_numbers_are_to_ = int_numbers_are_to__, 
      str_substitute_inf_with_ = str_substitute_inf_with__
    )
  
  if (file.exists(temp_file))
  {
    file.remove(temp_file)
  }
  
  return(annotations_)
}



write_inputs <- function(lists_, list_names_list_str, dir_name_str = 'checked_input')
{
  dir.create(dir_name_str)
  
  purrr::walk2(
    .x = lists_,
    .y = list_names_list_str,
    .f = function(x, y) {
      file_name <- paste0(dir_name_str, '/', y, '.tsv')
      readr::write_tsv(x = x, path = file_name)
    }
  )
}



annotate_identifiers_to_geneID_wrapper_for_using_dfs_as_input <- function(descriptions__ = descriptions, df_to_annotate, str_experiment_name__ = '', chr_gene_identifier__ = 'Gene name', col_types___ = 'ncccccccccccccncc', int_numbers_are_from__ = 11, int_numbers_are_to__ = 13, str_substitute_inf_with__ = '15')
{
  temp_file <- tempfile()
  
  readr::write_tsv(x = df_to_annotate, path = temp_file)
  
  annotated_ <-
    annotate_identifiers_to_geneID(
      str_filename_ = temp_file,
      str_experiment_name = str_experiment_name__,
      descriptions_ = descriptions__,
      chr_gene_identifier_ = chr_gene_identifier__, 
      col_types__ = col_types___, 
      int_numbers_are_from_ = int_numbers_are_from__, 
      int_numbers_are_to_ = int_numbers_are_to__, 
      str_substitute_inf_with_ = str_substitute_inf_with__
    )
  
  if (file.exists(temp_file))
  {
    file.remove(temp_file)
  }
  
  return(annotated_)
}



# INPUT: regex_ - for detecting proper gene identifier, col_everything - column to search the pattern in, usually containing all the identifiers in given entry
prepare_input <- function(regex_, df_, col_everything = opts_everything_column_for_this_analysis, col_put_resulting_identifiers_here = 'Probe_ID', col_order = 'Paper', regex_two_for_extracting_identifiers_from_string = '', col_get_resulting_identifer_from_here, vec_regex_additional_pattern_removal = '')
{
  search_vector <-
    stringr::str_detect(string = df_[[col_everything]], pattern = regex_)
  
  df_output <- subset(x = df_, subset = search_vector)
  
  left_to_do <- subset(x = df_, subset = !search_vector | is.na(search_vector))
  
  filtering_goodness <-
    check_was_the_spliting_of_df_by_filtering_ok(df_original = df_,
                                                 list_df_splited = list(df_output, left_to_do))
  
  if (regex_two_for_extracting_identifiers_from_string != '') {
    df_output[[col_put_resulting_identifiers_here]] <-
      extract_from_string(
        chr_vec = df_output[[col_everything]],
        regex_one = regex_,
        regex_two = regex_two_for_extracting_identifiers_from_string
      )
  }
  else{
    df_output[[col_put_resulting_identifiers_here]] <- df_output[[col_get_resulting_identifer_from_here]]
  }

  if (vec_regex_additional_pattern_removal != '') {
    for (regex_nb in seq_along(vec_regex_additional_pattern_removal)) {
      df_output[[col_put_resulting_identifiers_here]] <-
        stringr::str_remove(string = df_output[[col_put_resulting_identifiers_here]], pattern = vec_regex_additional_pattern_removal[[regex_nb]])
    }
  }
  
  df_output <-
    df_output[order(df_output[[col_order]]), ] 
  
  output <- list(df_output, left_to_do, filtering_goodness)
  return(output)
}




capture_LOCs_orphaned_by_previous_stage_entrez_search <- function(everything_column_for_this_analysis_ = opts_everything_column_for_this_analysis, input_EntrezGeneID_LOC_, left_to_do__)
{
  if (everything_column_for_this_analysis_ == 'everything2') 
  {
    proper_rows <- subset(x = input_EntrezGeneID_LOC_, subset = input_EntrezGeneID_LOC_$Probe_ID != '')
    
    orphaned_rows <- subset(x = input_EntrezGeneID_LOC_, subset = input_EntrezGeneID_LOC_$Probe_ID == '')
    
    left_and_orphaned <- rbind(orphaned_rows, left_to_do__)
    
    spliting_good_or_not <- check_was_the_spliting_of_df_by_filtering_ok(df_original = input_EntrezGeneID_LOC_, list_df_splited = list(orphaned_rows, proper_rows))
    
    result <- list(proper_rows, left_and_orphaned, spliting_good_or_not)
    
    return(result) 
  }
}



#
check_if_geneID_was_not_used_before <- function(input_EntrezGeneId_ = input_EntrezGeneId, col_with_ID_to_check = 'Probe_ID', list_of_colnames_to_check_if_the_ID_is_present_in, left_to_do__ = left_to_do_)
{
  is_checked_value_in_any_of_columns_to_check_lgl <- list()
  for (n in seq_along(list_of_colnames_to_check_if_the_ID_is_present_in)) {
    is_checked_value_in_any_of_columns_to_check_lgl[[n]] <- input_EntrezGeneId_[[col_with_ID_to_check]] == input_EntrezGeneId_[[ list_of_colnames_to_check_if_the_ID_is_present_in[[n]] ]]
    
    is_checked_value_in_any_of_columns_to_check_lgl[[n]] <- gtools::na.replace(is_checked_value_in_any_of_columns_to_check_lgl[[n]], F)
  }
  
  

  is_checked_value_in_any_of_columns_to_check_lgl <- Reduce("|", is_checked_value_in_any_of_columns_to_check_lgl)
  new_ids <- subset(x = input_EntrezGeneId_, subset = !is_checked_value_in_any_of_columns_to_check_lgl)
  
  new_ids$Paper <- as.integer(new_ids$Paper) ### !!! Added
  new_ids <- new_ids[order(new_ids$Paper),] 
  
  old_ids <- subset(x = input_EntrezGeneId_, subset = is_checked_value_in_any_of_columns_to_check_lgl)
  
  was_splitting_ok <- check_was_the_spliting_of_df_by_filtering_ok(df_original = input_EntrezGeneId_, list_df_splited = list(new_ids, old_ids))
  
  left_to_do_return <- rlist::list.rbind(old_ids)
  left_to_do_return <- rbind(old_ids, left_to_do__)
  
  return(list(new_ids, left_to_do_return, was_splitting_ok))
}



# OUTPUT: vectorwith first gene name that was not previously analyzed based on probe_IDs_already_done_columnName_list. 
get_geneNames_filtered_with_already_done_columns <- function(input_gene_symbol_ = input_gene_symbol, colname_containing_geneSymbols_to_check = 'Gene_symbol', colname_with_output_probeNames = 'Probe_ID', probe_IDs_already_done_columnName_list_ = opts_probe_IDs_already_done_columnName_list)
{
  if (!is.na(probe_IDs_already_done_columnName_list_))
  {
    all_IDs <-
      stringr::str_split(string = input_gene_symbol_[[colname_with_output_probeNames]], pattern = ',')
    
    new_names <- list()
    for (ID_set_number in seq_along(all_IDs)) {
      ID_set <- all_IDs[[ID_set_number]]
      was_unused_geneID_found <- F
      
      for (ID_number in seq_along(ID_set)) {
        current_tested_geneID <- all_IDs[[ID_set_number]][[ID_number]]
        if(was_unused_geneID_found == F){
          new_names[[ID_set_number]] <- current_tested_geneID
        }
        was_name_found <- F
        
        # This double if is because we should detect if we have already checked entrezID for proposed tested LOC id. Only if we
        if (was_name_found == F & !stringr::str_detect(string = current_tested_geneID, pattern = 'LOC') & was_unused_geneID_found == F)
        {
          current_done_geneID_counter <- 0
          for (done_columnName in probe_IDs_already_done_columnName_list_) {
            current_done_geneID <-
              input_gene_symbol_[[done_columnName]][[ID_set_number]]

            if (current_done_geneID == current_tested_geneID && !is.na(current_done_geneID) && !is.na(current_tested_geneID)) {
              new_names[[ID_set_number]] <- ''
              was_name_found <- T
            }
            current_done_geneID_counter <- current_done_geneID_counter + 1
            if(current_done_geneID_counter == length(probe_IDs_already_done_columnName_list_) & was_name_found == F)
            {
              was_unused_geneID_found <- T
            }
          }
        } else if (was_name_found == F & stringr::str_detect(string = current_tested_geneID, pattern = 'LOC') & was_unused_geneID_found == F)
        {
          current_tested_geneID_with_LOC_removed <- stringr::str_remove_all(string = current_tested_geneID, pattern = 'LOC')
          current_done_geneID_counter <- 0
          for (done_columnName in probe_IDs_already_done_columnName_list_) {
            current_done_geneID <-
              input_gene_symbol_[[done_columnName]][[ID_set_number]]
            
            if ( (current_done_geneID == current_tested_geneID | current_done_geneID == current_tested_geneID_with_LOC_removed ) && !is.na(current_done_geneID) && !is.na(current_tested_geneID)  & was_unused_geneID_found == F) {
              new_names[[ID_set_number]] <- ''
              was_name_found <- T
            }
            current_done_geneID_counter <- current_done_geneID_counter + 1
            if(current_done_geneID_counter == length(probe_IDs_already_done_columnName_list_) & was_name_found == F)
            {
              was_unused_geneID_found <- T
            }
          }
        }
        }
    } 
    return(as.character(unlist(new_names)))
  } else if (is.na(probe_IDs_already_done_columnName_list_))
  {
    first_ID <- stringr::str_remove_all(string = input_gene_symbol_[[colname_with_output_probeNames]], pattern = ',.*')
    return(as.character(first_ID))
  }
}



# Uses columns 'Paper', 'Probe_ID', 'external_gene_name', 'Probe_ID_2'
post_annotate_identifiers_to_geneID <- function(annotated_names_to_geneIDs_, Probe_ID_left_)
{
  annotated_names_to_geneIDs_$Paper <- as.integer(annotated_names_to_geneIDs_$Paper) ### !!! added
  annotated_names_to_geneIDs_ <-
    annotated_names_to_geneIDs_[order(annotated_names_to_geneIDs_$Paper), ]
  annotated_names_to_geneIDs_[[Probe_ID_left_]] <- annotated_names_to_geneIDs_$Probe_ID
  annotated_names_to_geneIDs_$Probe_ID <-
    annotated_names_to_geneIDs_$external_gene_name
  annotated_names_to_geneIDs_$Probe_ID_2 <- NULL
  annotated_names_to_geneIDs_$external_gene_name <- NULL
  
  return(annotated_names_to_geneIDs_)
}



reformat_Gene_ID_column_wrapper <- function(df_, WITHIN_COL_SEPARATOR_ = WITHIN_COL_SEPARATOR)
{
  df_$Gene_ID <- stringr::str_remove_all(string = df_$Gene_ID, pattern = ',00')
  df_$Gene_ID[df_$Gene_ID == '?'] <- NA
  df_$Gene_ID[df_$Gene_ID == '-'] <- NA
  df_$Gene_ID[df_$Gene_ID == 'N/A'] <- NA
  df_$Gene_ID <- stringr::str_replace_all(string = df_$Gene_ID, pattern = '(///)|(;)', replacement = WITHIN_COL_SEPARATOR)
  
  return(df_)
}



# Removes all '_predicted'!
reformat_Gene_symbol_column_wrapper <- function(df_, WITHIN_COL_SEPARATOR_ = WITHIN_COL_SEPARATOR)
{
  df_$Gene_symbol[df_$Gene_symbol == 'N/A'] <- NA
  df_$Gene_symbol <- stringr::str_replace_all(string = df_$Gene_symbol, pattern = '( /// )|(///)', replacement = WITHIN_COL_SEPARATOR)
  df_$Gene_symbol <- stringr::str_remove_all(string = df_$Gene_symbol, pattern = stringr::regex('(_predicted)|(_RAT)', ignore_case = T))
  df_$Gene_symbol <- stringr::str_replace_all(string = df_$Gene_symbol, pattern = '^(//)', replacement = '')
  df_$Gene_symbol <- stringr::str_replace_all(string = df_$Gene_symbol, pattern = '^(/)', replacement = '')
  df_$Gene_symbol <- stringr::str_replace_all(string = df_$Gene_symbol, pattern = '/', replacement = WITHIN_COL_SEPARATOR)
  df_$Gene_symbol <- stringr::str_replace_all(string = df_$Gene_symbol, pattern = 'b Rpl10', replacement = 'Rpl10')
  df_$Gene_symbol <- stringr::str_replace_all(string = df_$Gene_symbol, pattern = 'Fam171 a2', replacement = 'Fam171a2')
  df_$Gene_symbol[df_$Gene_symbol == 'Image:619641'] <- NA
  df_$Gene_symbol[df_$Gene_symbol == '–'] <- NA # This is long hypen ubuntu utf
  df_$Gene_symbol[df_$Gene_symbol == '?'] <- NA # This is long hypen win10 utf
  df_$Gene_symbol[df_$Gene_symbol == '-'] <- NA ### !!! This was added
  df_$Gene_symbol[df_$Gene_symbol == '---'] <- NA
  df_$Gene_symbol <- stringr::str_replace_all(string = df_$Gene_symbol, pattern = stringr::regex("(?<=\\d{2})rik", ignore_case = TRUE), replacement = 'Rik') # This removes problems with different spelling of Rik symbols
  
  
  return(df_)
}



reformat_GenBank_Accession_column_wrapper <- function(df_, WITHIN_COL_SEPARATOR_ = WITHIN_COL_SEPARATOR)
{
  df_$GenBank_Accession <- stringr::str_replace_all(string = df_$GenBank_Accession, pattern = '( /// )|(///)|;', replacement = WITHIN_COL_SEPARATOR)
  df_$GenBank_Accession <- stringr::str_replace_all(string = df_$GenBank_Accession, pattern = 'NM ', replacement = 'NM_')
  df_$GenBank_Accession[df_$GenBank_Accession == '---'] <- NA
  df_$GenBank_Accession <- stringr::str_remove_all(string = df_$GenBank_Accession, pattern = '\\?')
  df_$GenBank_Accession <- stringr::str_replace_all(string = df_$GenBank_Accession, pattern = '^,', replacement = '')
  df_$GenBank_Accession[df_$GenBank_Accession == ''] <- NA
  
  
  return(df_)
}


reformat_Ensembl_ID_column_wrapper <- function(df_, WITHIN_COL_SEPARATOR_ = WITHIN_COL_SEPARATOR)
{
  df_$Ensembl_ID <- stringr::str_replace_all(string = df_$Ensembl_ID, pattern = '( \\+ )|(///)|;', replacement = WITHIN_COL_SEPARATOR)
  
  return(df_)
}


# Take current probeID and remove it from everthing2 column
remove_used_name_from_everything2_wrapper <- function(input_, col_identifer_to_remove = 'Probe_ID')
{
  input_$everything2 <- as.character(
    purrr::map2(
      .x =  input_$everything2,
      .y = input_[[col_identifer_to_remove]],
      .f = function(x, y) {
        if (!is.na(y)) {
          stringr::str_remove_all(string = x, pattern = stringr::coll(y))
        }
        else{
          x
        }
      }
    )
  )
  return(input_)
}



# Take current probeID and remove it from everthing2 column
remove_used_entrezID_from_everything2_wrapper <- function(input__, col_identifer_to_remove_ = 'Probe_ID')
{

  input__$temp_to_remove <- paste0('LOC', input__[[col_identifer_to_remove_]])
  
  
  input__ <- remove_used_name_from_everything2_wrapper(input_ = input__, col_identifer_to_remove = 'temp_to_remove')
  
  input__$temp_to_remove <- NULL
  
  return(input__)
}



get_percentages_of_annotation <- function(df_)
{
  found <- df_ %>%
    dplyr::filter(!is.na(external_gene_name)) %>%
    dplyr::group_by(Paper) %>%
    dplyr::summarise(dplyr::n())
  
  queried <- df_%>%
    dplyr::filter(!is.na(Probe_ID)) %>%
    dplyr::group_by(Paper) %>%
    dplyr::summarise(dplyr::n())
  
  percentages <- merge(x = found, y = queried, by = 'Paper', all.y = T)
  percentages[is.na(percentages)] <- 0
  percentages$percentage <- percentages[[2]]/percentages[[3]]
  
  return(percentages)
}



set_empty_annotations_to_na <- function()
{
  sapply(c("annotations_ensembl_gene_id", "annotations_ensembl_transcript_id", "annotations_entrezgene_id", "annotations_refseq_mrna", "annotations_accession"), function(x) 
    if(!exists(x)) { assign(x, NA, envir=.GlobalEnv) } )
}



uniformize_finalized_colnames <- function(finalized_ = finalized, opts_nb_of_analysis_stages_ = opts_nb_of_analysis_stages, current_stage_ = current_stage)
{
  finalized_ <- dplyr::select(finalized_, external_gene_name, dplyr::everything())
  
  if (current_stage_ < opts_nb_of_analysis_stages_) { # Exclude final finalized here
    finalized_[eval(parse(text = paste0('opts_ann_', current_stage_)), parent.frame())] <- NA
    if (current_stage_+1  < opts_nb_of_analysis_stages_) { # Exclude secons-to-last finalized
      for(stages_above_1 in seq(from = current_stage_+1, to = opts_nb_of_analysis_stages_-1)){
        finalized_[eval(parse(text = paste0('opts_ann_', stages_above_1, '_n')), parent.frame())] <- NA
        finalized_[eval(parse(text = paste0('opts_ann_', stages_above_1)), parent.frame())] <- NA
      }
    }
    finalized_[eval(parse(text = paste0('opts_ann_', opts_nb_of_analysis_stages_, '_n')), parent.frame())] <- NA
  }
  
  finalized_[eval(parse(text = paste0('opts_ann_', opts_nb_of_analysis_stages_)), parent.frame())] <- NA
  
  return(finalized_)
}



uniformize_bad_gene_symbol_colnames <- function(bad = input_bad_gene_symbol, opts_nb_of_analysis_stages_ = opts_nb_of_analysis_stages, current_stage_ = current_stage)
{
  bad$external_gene_name <- NA
  bad <- dplyr::select(bad, external_gene_name, dplyr::everything())
  for(stages_above_1 in seq(from = current_stage_, to = opts_nb_of_analysis_stages_)){
    bad[eval(parse(text = paste0('opts_ann_', stages_above_1, '_n')), parent.frame())] <- NA
    bad[eval(parse(text = paste0('opts_ann_', stages_above_1)), parent.frame())] <- NA
  }
  
  colnames(bad)[5] <- 'GEO ID'
  
  return(bad)
}



uniformize_leftovers_colnames <- function(leftovers_ = leftovers){
  leftovers_$external_gene_name <- NA
  leftovers_ <- dplyr::select(leftovers_, external_gene_name, dplyr::everything())
  return(leftovers_)
}



get_usable_bad_ids <- function(bads_ = bads){
  bads_$external_gene_name <- stringr::str_remove(string = bads_$Gene_symbol, pattern = ',.*')
  
  bads_ <- subset(x = bads_, subset = !is.na(bads_$external_gene_name))
  
  return(bads_)
}



get_subset_vector_for_entries_with_3_or_more_values_per_paper <- function(final_good_dataset_){
  temp <- final_good_dataset_ %>%
    dplyr::select(Paper, lower_final_gene_name) %>%
    unique()
  
  temp$present <- T
  
  spread_temp <- tidyr::spread(data = temp, key = Paper, value = present)
  
  entries_with_3_or_more_values <- purrrlyr::by_row(.d = spread_temp[,-1], .collate = "rows", ..f = function(x) {sum(!is.na(x)) >= 3})$.out
  
  return(entries_with_3_or_more_values)
}



create_subset_of_exps_with_at_least_3_papers_wrapper <- function(final_good_dataset__ = final_good_dataset, experiments_to_include_ = experiments_to_include, save_as_chr)
{
  subset_matrix <-subset(final_good_dataset__, subset = final_good_dataset__$Experiment %in% experiments_to_include_)
  
  medianed_subset_matrix <- subset_matrix %>%
    dplyr::group_by(Experiment, lower_final_gene_name) %>%
    dplyr::summarize(logFC_median = median(logFC))
  
  spread_medianed_subset_matrix <- tidyr::spread(data = medianed_subset_matrix, key = Experiment, value = logFC_median)
  
  bool_entries_with_3_or_more_values_subset_matrix <- get_subset_vector_for_entries_with_3_or_more_values_per_paper(final_good_dataset_ = subset_matrix)
  
  at_least_in_3_papers_spread_med_sub_mat <- subset(x = spread_medianed_subset_matrix, subset = bool_entries_with_3_or_more_values_subset_matrix)
  
  # save(at_least_in_3_papers_spread_med_sub_mat, file = save_as_chr)
  
  return(at_least_in_3_papers_spread_med_sub_mat)
}



# OUTPUT: df
get_number_and_percentage_of_directionality_of_exp_first_column_names <- function(spread_dataset_){
  
  if (!is.character(spread_dataset_[[1]])) {
    stop('First column needs to include gene names')
  }
  
  no_of_exps <- purrrlyr::by_row(.d = spread_dataset_[,-1], .collate = "rows", ..f = function(x) {sum(x != 0)})$.out
  
  temp_2 <- purrrlyr::by_row(.d = spread_dataset_[,-1], .collate = "rows", ..f = function(x) {sum(x > 0)})$.out
  
  perc_of_upregulated <- temp_2/no_of_exps
  
  spread_dataset_ <- cbind(spread_dataset_[1], no_of_exps, perc_of_upregulated)
  
  return(spread_dataset_)
}



get_number_and_percentage_of_directionality_of_paper_first_column_names <- function(spread_data_)
{
  spread_data_[spread_data_ == 0] <- NA
  
  temp_ <- tidyr::gather(data = spread_data_, key = "exp", value = "logFC", na.rm = T, -lower_final_gene_name)
  temp_ <- temp_[,-3]
  
  temp_$paper <- stringr::str_remove(string = temp_$exp, pattern = '_.*')
  temp_$exp <- NULL
  
  temp_ <- temp_ %>%
    unique() %>%
    dplyr::select(lower_final_gene_name) %>%
    dplyr::group_by(lower_final_gene_name) %>%
    dplyr::summarise(number = dplyr::n())
  
  return(temp_)
}



# OUTPUT: [[1]] - actual prepared vector, [[2]] - dataframe for validation
cleanup_differing_units <- function(charvec_, unit_identifier_regexList, unit_multiplier_list)
{
  temp <- data.frame(charvec_)
  
  temp$order <- c(1:length(charvec_))
  
  temp$number <- stringr::str_extract(string = temp$charvec_, pattern = '[0-9]{1,}')
  
  temp$number <- sapply(X = temp$number, FUN = function(x) {paste0(x, collapse = '') } )
  
  for (unit_nb in seq_along(unit_identifier_regexList)) {
    temp2 <- subset(x = temp , subset = stringr::str_detect(string = charvec_, pattern = unit_identifier_regexList[[unit_nb]]))
    
    temp2$temp <- as.numeric(temp2$number) * unit_multiplier_list[[unit_nb]]
    
    if (unit_nb == 1) {
      ret_urn <- temp2
    } else if (unit_nb > 1) {
      ret_urn <- rbind(ret_urn, temp2)
    }
  }
  
  ret_urn <- merge(temp, ret_urn, by = 'order', all.x = T)
  
  if(length(charvec_) != length(ret_urn[[1]])){
    stop('Output differs from input. Probably You have an entry which matches to two or more of the unit_identifier_regexList')
  }
  
  ret_urn <- dplyr::select(.data = ret_urn, order, charvec_.x, temp)
  
  ret_urn <- ret_urn[order(ret_urn$order),]
  
  return(list(ret_urn$temp, ret_urn))
}



brain_part_cleanup_wrapper <- function(brain_parts_ = descriptions$Brain_part)
{
  brain_parts_ <- tolower(brain_parts_)
  brain_parts_[brain_parts_ == 'ventral tegmental areas'] <- 'ventral tegmental area'
  brain_parts_[brain_parts_ == 'hypothalamic paraventricular\nnucleus'] <- 'hypothalamic paraventricular nucleus'
  brain_parts_[brain_parts_ == 'paraventricular nucleus of the hypothalamus'] <- 'hypothalamic paraventricular nucleus'
  brain_parts_ <- stringr::str_replace_all(string = brain_parts_, pattern = '.*myg.*', replacement = 'amygdala')
  brain_parts_ <- stringr::str_replace_all(string = brain_parts_, pattern = '.*hipp.*', replacement = 'hippocampus')
  brain_parts_ <- stringr::str_replace_all(string = brain_parts_, pattern = '.*cumbens.*', replacement = 'nucleus_accumbens')
  
  return(brain_parts_)
  
}
