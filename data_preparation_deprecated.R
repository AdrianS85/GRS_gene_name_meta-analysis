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


##### FUNCTION FOR FORMATING DATA FOR CLUSTERING ##### 

read_preformated_data <- function(str_filename, int_numbers_are_from = 8, int_numbers_are_to = 8, col_types_ = 'nncccccc')
{
  temp_data <- readr::read_tsv(
    str_filename, 
    col_types = col_types_)
  
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(
    X = temp_data[int_numbers_are_from:int_numbers_are_to], 
    FUN = function(x) { stringr::str_replace(
      string = x, 
      pattern = ",", 
      replacement = ".")}) ### !!! THIS IS VERY COOL CONSTRUCT THAT I NEED TO NOTE!!
  
  temp_data[int_numbers_are_from:int_numbers_are_to] <- lapply(
    temp_data[int_numbers_are_from:int_numbers_are_to], 
    as.numeric)
  
  return(temp_data)
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
                         dataset = "sboliviensis_gene_ensembl"))
  
  message( paste0('Mart set as ', usedMart__@dataset, ' for step ', int_loop) )
  return(usedMart__)
}



### 
get_the_potental_identifiers <- function(usedMart___)
{
  # Here we extract all the potential gene identifiers
  ##### !!! THIS MAY NEED FURTHER WORK !!! ##### 
  message( paste0('Extracting potental_identifiers for step ', n) )
  
  potental_identifiers <- c(usedMart___@filters[grep(pattern = "^ensembl(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^refseq(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^affy(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^agilent(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^illumina(.*)", usedMart___@filters[[1]]) , 1])
  
  message( paste0('Potental identifiers for step ', n, 'extracted') )
  
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


# testusemart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
# testfilters <- biomaRt::listFilters(mart = testusemart)

set_identifiers_used_for_annotation_if_not_probeID <- function(str_identifier_type, list_LIST_DATA = LIST_DATA, identifiers_used_for_annotation_if_probeID_is_used = DATA_FROM_HIGHEST_HIT_ANALYSIS[[2]])
{
  if(str_identifier_type == 'Probe_ID')
  {
    return(identifiers_used_for_annotation_if_probeID_is_used)
  }
  else
  {
    
    switch(str_identifier_type,
           "ensembl_g" = return( rep(str_identifier_type, length(LIST_DATA)) ),
           "genebank_ID" = return( rep(str_identifier_type, length(LIST_DATA)) ),
           "Gene_ID" = return( rep(str_identifier_type, length(LIST_DATA)) ),
           "external_gene_name" = return( rep('external_gene_name', length(LIST_DATA)))
    )
  }
  
}



change_name_to_proper_format_given_the_species <- function(df_name, str_species)
{
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


#This function is used to query ncbi-gene db and return first (best, most current?) gene-id identifier
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
    
    Sys.sleep(0.5)
  }
  return(ids)
}


#This function should except output of search_for_ids function. BEWARE - this only returns first geneID found
make_table_with_original_and_ncbi_ids <- function(list_returned_by_entrez_gene_search)
{
  temp_list <- list()
  for (n in seq(length(list_returned_by_entrez_gene_search)))
  {
    temp_list[[n]] <- ''
    
    if (length(list_returned_by_entrez_gene_search[[n]]$ids) != 0)
    {
      temp_list[[n]][[1]] <- list_returned_by_entrez_gene_search[[n]]$ids[[1]]
    } 
    else
    {
      temp_list[[n]][[1]] <- 'none'
    }
    temp_name <- list_returned_by_entrez_gene_search[[n]]$QueryTranslation
    
    
    temp_list[[n]][[2]] <- gsub(
      pattern = "(\\[All Fields\\])|\\)|\\(", 
      replacement = '', 
      x = list_returned_by_entrez_gene_search[[n]]$QueryTranslation)
  }
  
  temp_df <- as.data.frame(rlist::list.rbind(temp_list))
  colnames(temp_df) <- c('ncbi_response', 'query')
  
  return(temp_df)
}

get_the_highest_hit_returning_id_type <- function(list_LIST_DATA_, str_vector_of_species_names, int_Probe_IDs_to_test = 200)
{
  ### We need to first check appropriate probe ids on smaller dataset and only then do actual annotation, because its too slow otherwise. Hence this shortened list !!! ADD RANDOM SELECTION OF ROWS! !!!
  
  SHORT_LIST_DATA <- lapply(X = list_LIST_DATA_, FUN = function(x){ x[1:int_Probe_IDs_to_test,] })
  ANNOT_SHORT_LIST_DATA <- rep(list(list()), times = length(SHORT_LIST_DATA))
  all_ID_annotations <- rep(list(list()), times = length(SHORT_LIST_DATA))
  
  for(n in seq_along(SHORT_LIST_DATA))
  {
    # Here we establish which mart(species) we are using in this given dataset based on "exp_species" vector
    ##### !!! THIS NEEDS TO BE MANUALLY ESTABLISHED BASE ON exp_species !!! ##### 
    usedMart_ <- set_mart_to_be_used(str_vector_of_species_names_ = str_vector_of_species_names[n], int_loop = n)
    
    potental_identifiers <- get_the_potental_identifiers(usedMart___ = usedMart_)
    
    
    
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
  for(n in seq_along(LIST_DATA)){
    for(m in seq_along(which(df_all_ID_annotations[[n]] == max(df_all_ID_annotations[[n]])))){
      HIGHEST_HIT_LIST[[n]][[m]] <- ANNOT_SHORT_LIST_DATA[[n]][[(which(df_all_ID_annotations[[n]] == max(df_all_ID_annotations[[n]]))[m])]]
    } }
  
  rm(df_all_ID_annotations)
  
  # This is list allowing us to check how many hits did we get for highest-yielding feature for given experiment
  check_annotation_percentage <- lapply(HIGHEST_HIT_LIST, FUN = function(x) {  length(x[[1]][[1]]) / int_Probe_IDs_to_test })
  
  # Here we simply copy/establish names of features, that were highest by themselves. The conditions ask: 1) is there at least a single hit with highest number (I dont know if there can be 0 though...) 2) Is the first (and though each) highest hit list has at least single hit?
  NAMES_HIGHEST_HIT_LIST <- lapply(X = HIGHEST_HIT_LIST, FUN = set_0_hit_annotations_to_na)
  
  readr::write_tsv(as.data.frame(NAMES_HIGHEST_HIT_LIST), 'data_on_which_platform_to_use_for_probes_based_probes.tsv')
  
  ###### Here we estalish correct lists to analyzed in further steps (currently - need to remove experiments with microarrays not captured in ensembl) ######
  ###### Yeah, i dont know what to do here
  WHICH_EXP_TO_ANAL <- seq_along(NAMES_HIGHEST_HIT_LIST)
  
  return(list(WHICH_EXP_TO_ANAL, NAMES_HIGHEST_HIT_LIST, HIGHEST_HIT_LIST, ANNOT_SHORT_LIST_DATA, all_ID_annotations))
}



PRE_DATA <- read_preformated_data("zbitka_from_prepared_data.txt")

#################################################
##### ANNOTATE PROBE_IDs WITH ENSEMBL NAMES #####
#################################################

library(tidyverse)

##### PREPARING LISTS TO BE USED ##### 

### Read in the table with description of experiments. They need to have the same identifiers as data in PRE_DATA, that is: Paper(int), Experiment(str)
descriptions <- readr::read_tsv("descriptions.txt", col_types = "nccccccccccccccccccccccc")
### List of species names based on data in description files. This is a file we will be working on. Species is in format: small letters, english plural of species
exp_species <- descriptions %>%
  select("Paper_ID", "Species") %>%
  unique() %>%
  filter(Paper_ID %in% unique(PRE_DATA$Paper))
  # filter(Paper_ID %in% c(40, 53))
  # filter(Paper_ID %in% c(1, 2))


readr::write_tsv(exp_species, 'exp_species_used_for_testing_which_platform_to_use_probes.tsv')

# The list object is needed for annotation of probes with ensembl names
LIST_DATA <- split(PRE_DATA, f = PRE_DATA$Paper) ###!!! <---
write_lenghts_of_list_objects(LIST_DATA, 'list_data_lenghts_probes.tsv')

write_lenghts_of_list_objects(list_ = LIST_DATA, string_name_of_the_file = 'zbitka_from_prepared_data_lenght.txt')

# How many Probe_IDs should we test to established appropriate microarray for given paper

### Here we make list with number of lists equal to number of experiments
###!!! <---

### This is just a help-list, it should probably be lower, when it is actually used

##### PREPARING LISTS TO BE USED ##### 



##### GET THE HIGHEST-HIT-RETURNING ID TYPE ##### 
### list_LIST_DATA_ format: just the LIST_DATA set in previous lines, vector_of_species_names  format: small letters, english plural of species

DATA_FROM_HIGHEST_HIT_ANALYSIS <- get_the_highest_hit_returning_id_type(
  list_LIST_DATA_ = LIST_DATA, 
  str_vector_of_species_names = exp_species$Species)

WHICH_EXP_TO_ANAL <- DATA_FROM_HIGHEST_HIT_ANALYSIS[[1]]
# identifiers_used_for_annotation <- DATA_FROM_HIGHEST_HIT_ANALYSIS[[2]]

temp2 <- as.data.frame( rlist::list.rbind(ANNOT_SHORT_LIST_DATA) )


##### GET THE HIGHEST-HIT-RETURNING ID TYPE #####

# I dont know how to pass column names via string and $ operator. [''] way doesnt work. So I guess im stuck with hard-coding it for now.
# this part also needs to have a capitalized first letter





LIST_DATA <- split(PRE_DATA, f = PRE_DATA$Paper) ###!!! <---

LIST_DATA <- purrr::map2(.x = LIST_DATA, .y = exp_species$Species, .f = ~ change_name_to_proper_format_given_the_species(.x, .y))

temp <- rlist::list.rbind(LIST_DATA)
readr::write_tsv(temp, 'data_v3_names_rn_better_names.txt')


###### PROBE_ID ANNOTATION WITH ENSEMBL GENE NAMES ###### 
WHICH_EXP_TO_ANAL <- seq( length(LIST_DATA) )
ANNOT_LIST_DATA <- list()
identifiers_used_for_annotation <- set_identifiers_used_for_annotation_if_not_probeID('external_gene_name')
# identifiers_used_for_annotation <- rep('external_gene_name', length(LIST_DATA))
try(
  {
    for (n in WHICH_EXP_TO_ANAL){
      
      usedMart_ = set_mart_to_be_used(str_vector_of_species_names_ = exp_species$Species[n], int_loop = n)
      
      message( paste0('Starting ANNOT_LIST_DATA step for ', n) )    
      ANNOT_LIST_DATA[[n]] <- biomaRt::getBM(
        attributes = c(identifiers_used_for_annotation[[n]], "external_gene_name"),
        filters = identifiers_used_for_annotation[[n]],
        values = LIST_DATA[[n]]$Probe_ID,
        uniqueRows = F,
        mart = usedMart_)
      
      message( paste0('Competed usedMart step for ', n) )
    }
  }
)




### !!! INTER-SPECIES GENE-NAME STANDARDIZATION? !!! ###

###### PROBE_ID ANNOTATION WITH ENSEMBL GENE NAMES ###### 



###### POST-PREPARATION OF ANNOTATED PROBE_ID TABLES ###### 

## HERE WE WILL COLLAPSE duplicated probe-genename pairs WITHIN ANNOTATED GENE LIST (and also change column names to standarized ones)
UNIQ_ANNOT_LIST_DATA <- list()
for (n in WHICH_EXP_TO_ANAL){
  UNIQ_ANNOT_LIST_DATA[[n]] <- unique(ANNOT_LIST_DATA[[n]])
  colnames(UNIQ_ANNOT_LIST_DATA[[n]]) <- c("Probe_ID", "ensembl_gene_name")
}

## Here we will collapse all gene names for given probe
for (n in WHICH_EXP_TO_ANAL){
  UNIQ_ANNOT_LIST_DATA[[n]] <- aggregate(
    ensembl_gene_name~Probe_ID, 
    data = UNIQ_ANNOT_LIST_DATA[[n]], 
    FUN = stringr::str_c) ## Here we aggregate the genenames into single row
  UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- lapply(
    X = UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name, 
    FUN = paste, 
    collapse = "; ") ##
  UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- as.character(UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name) ## Here we make sure that collapsed genename column is not a list (need for writing function)
}

# Je?eli w ensemble nie ma przypisanego genu (oznaczenie NA w przys?anym pliku "100_genow") to program pokazuje we wszystkich rubrykach oznaczenie NA chocia? te kom?rki zawiera?y wcze?niej dane (plik AAAreview11). W pe?nych danych kt?re przys?a?e? mi przed ?wi?tami (plik 1 znajduj?cy si? w skopresowanych danych "TEST_ANNOTATION") tych sond w og?le nie ma. Wygl?da na to, ?e to co przys?a?e? wczoraj odnosi si? chyba do wcze?niejszego pliku w kt?rym dane dla sond nie maj?cych przypisanego genu w ensemble zosta?y ca?kowicie usuni?te a to nie jest dobre :( - AMS W LIST_DATA S? TE DANE, ALE W TEST_ANNOTATION JUZ NIE

## Here we merge probes annotated with ensembl with probes data that we have
TEST_ANNOTATION <- list()
for (n in WHICH_EXP_TO_ANAL){
  TEST_ANNOTATION[[n]] <- merge(LIST_DATA[[n]], UNIQ_ANNOT_LIST_DATA[[n]], by = "Probe_ID", all.x = T)
}


## Here lets also remove empty lists
### !!! MANUAL SUPPLEMENTARY STEP - TO BE REMOVED IN FINAL VERSION, AS WE WILL NOT HAVE EMPTY LISTS !!! ###
# TEST_ANNOTATION[[13]] <-NULL

# Here we bind all list elements into single data table
SINGLE_TEST_ANNOTATION <- data.table::data.table(rlist::list.rbind(TEST_ANNOTATION), stringsAsFactors = F)
readr::write_tsv(x = SINGLE_TEST_ANNOTATION, path = "SINGLE_ANNOTATION_probes.tsv")


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



