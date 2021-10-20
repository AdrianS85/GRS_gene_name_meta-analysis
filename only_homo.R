library(magrittr)
library(Hmisc)

save(homo, file = "homo.save")

homo <- list( "input" = openxlsx::read.xlsx(xlsxFile = "Human PTSD datasets.xlsx") )

homo$input_ensembl <- subset(homo$input, !is.na(homo$input$Ensembl_ID))

homo$input_no_ensembl <- subset(homo$input, is.na(homo$input$Ensembl_ID))
homo$input_no_ensembl$external_gene_name <- NA

homo$input_t2 <- openxlsx::read.xlsx(xlsxFile = "Supplementary Table 2 v4.xlsx")

homo$mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

homo$mart_mus <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")






###############
### ENSEMBL ###

homo$ens$getBM <- biomaRt::getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = homo$input_ensembl$Ensembl_ID,
  mart = homo$mart)

homo$ens$input_getBM <- merge(
  x = homo$input_ensembl, 
  y = homo$ens$getBM, 
  by.x = "Ensembl_ID", 
  by.y = "ensembl_gene_id", 
  all.x = T)

homo$ens$input_getBM$external_gene_name[homo$ens$input_getBM$external_gene_name == ""] <- NA

# here we get standardized homo names in "external_gene_name" column
homo$ens$done <- subset(x = homo$ens$input_getBM, !is.na(homo$ens$input_getBM$external_gene_name))

homo$ens$failed <- subset(x = homo$ens$input_getBM, is.na(homo$ens$input_getBM$external_gene_name))


### ENSEMBL ###
###############






##############
### ENTREZ ###

homo$entrez$input <- rbind(homo[["input_no_ensembl"]], homo$ens$failed)

homo$entrez$failed_1 <- subset(homo$entrez$input, is.na(homo$entrez$input$Gene_symbol))
homo$entrez$failed_1$external_gene_name <- tolower(homo$entrez$failed_1$Ensembl_ID)

homo$entrez$for_entrez <- subset(homo$entrez$input, !is.na(homo$entrez$input$Gene_symbol))


# Getting entrez ids for each gene name that exists
homo$entrez$query <- get_query_for_ncbi_geneID_annotation(
  char_vec_gene_id_to_query_with = homo$entrez$for_entrez$Gene_symbol, 
  char_vec_organism = rep('homo sapiens', length(homo$entrez$for_entrez$external_gene_name)), 
  chr_gene_identifier = 'Gene name')

homo$entrez$search <- search_for_ids_in_ncbi(homo$entrez$query)

homo$entrez$output <- make_and_write_table_with_original_and_ncbi_ids(
  entrez_gene_search_output = homo$entrez$search, 
  df_original_data = homo$entrez$for_entrez)

homo$entrez$output <- unique(homo$entrez$output)


# Annotate only the entrez ids that were found with standardized names.
### !!! What about the 7 gene names for which entrez id was not found?
### !!! What about many entrez id for which no standardized name was found?
homo$entrez$annot <- biomaRt::getBM(
  attributes = c("external_gene_name", "entrezgene_id"),
  filters = "entrezgene_id",
  values = homo$entrez$output$Gene_ID,
  mart = homo$mart)

homo$entrez$annot <- subset(homo$entrez$annot, homo$entrez$annot$external_gene_name != "")

homo$entrez$annot_output <- merge(
  x = homo$entrez$output, 
  y = homo$entrez$annot, 
  by.x = "Gene_ID", 
  by.y = "entrezgene_id", 
  all.x = T)

homo$entrez$for_entrez_annot_output <- merge(
  x = homo$entrez$for_entrez, 
  y = homo$entrez$annot_output, 
  by.x = "Gene_symbol", 
  by.y = "input_id", 
  all.x = T)

homo$entrez$for_entrez_annot_output <- dplyr::select(homo$entrez$for_entrez_annot_output, -c(26:29))

colnames(homo$entrez$for_entrez_annot_output)[26] <- "external_gene_name"


homo$entrez$done <- subset(x = homo$entrez$for_entrez_annot_output, !is.na(homo$entrez$for_entrez_annot_output$external_gene_name))

homo$entrez$failed_2 <- subset(x = homo$entrez$for_entrez_annot_output, is.na(homo$entrez$for_entrez_annot_output$external_gene_name))
homo$entrez$failed_2$external_gene_name <- tolower(homo$entrez$failed_2$Gene_symbol)

homo$entrez$failed <- rbind(homo$entrez$failed_1, homo$entrez$failed_2 )
homo$entrez$failed$external_gene_name[homo$entrez$failed$external_gene_name == "loc153328"] <- "slc25a48"

### ENTREZ ###
##############






###################
### PTSD TABLES ###

homo$output$homo_names <- rbind(homo$ens$done, homo$entrez$done, homo$entrez$failed)

homo$output$homo_names$external_gene_name_homo <- tolower(homo$output$homo_names$external_gene_name)
homo$output$homo_names$external_gene_name <- NULL

homo$output$medianed_final_good_dataset <- dplyr::summarize(
  dplyr::group_by(homo$output$homo_names, `Paper/Group.code`, external_gene_name_homo), 
  logFC_median = median(logFC))

homo$output$spread_medianed_final_good_dataset <- tidyr::spread(
  data = homo$output$medianed_final_good_dataset, 
  key = `Paper/Group.code`, 
  value = logFC_median)

homo$output$spread_medianed_final_good_dataset <- dplyr::rename(homo$output$spread_medianed_final_good_dataset, Standarised_gene_symbol = external_gene_name_homo)

homo$output$spread_medianed_final_good_dataset[is.na(homo$output$spread_medianed_final_good_dataset)] <- 0

homo$output$pap_exp <- get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(homo$output$spread_medianed_final_good_dataset, '')

homo$output$pap_exp$merge <- merge(
  x = homo$output$pap_exp$exps, 
  y = homo$output$pap_exp$pap, 
  by = "Standarised_gene_symbol")



colnames(homo$output$homo_names) <- stringr::str_replace_all(colnames(homo$output$homo_names), "\\.", " ")
homo$output$homo_names$Standarised_gene_symbol <- homo$output$homo_names$external_gene_name_homo
openxlsx::write.xlsx(homo$output$homo_names, file = "Supplementary Dataset 1 ptsd.xlsx")

colnames(homo$output$pap_exp$merge)[2:4] <- c("Total number of transcriptomic comparisons with altered expression", "Fraction of comparisons with up-regulated expression", "Number of reporting papers")
openxlsx::write.xlsx(homo$output$pap_exp$merge , file = "Supplementary Table 2 ptsd.xlsx")

### PTSD TABLES ###
###################






####################
### MUS HOMOLOGS ###

homo$mouse_homs$getBM <- biomaRt::getBM(
  attributes = c("external_gene_name", "mmusculus_homolog_ensembl_gene"),
  filters = "external_gene_name",
  values = unique(homo$output$homo_names$external_gene_name),
  mart = homo$mart)

homo$mouse_homs$getBM_proper <- subset(homo$mouse_homs$getBM, homo$mouse_homs$getBM$mmusculus_homolog_ensembl_gene != "")


homo$mouse_homs$homologs_and_symbols <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = homo$mouse_homs$getBM_proper$mmusculus_homolog_ensembl_gene,
  mart = homo$mart_mus)

homo$mouse_homs$homologs_and_symbols_proper <- subset(homo$mouse_homs$homologs_and_symbols, homo$mouse_homs$homologs_and_symbols$external_gene_name != "")


homo$mouse_homs$getBM_homologs_and_symbols <- merge(
  x = homo$mouse_homs$getBM_proper, 
  y = homo$mouse_homs$homologs_and_symbols_proper, 
  by.x = "mmusculus_homolog_ensembl_gene", 
  by.y = "ensembl_gene_id")
colnames(homo$mouse_homs$getBM_homologs_and_symbols) <- c("mmusculus_homolog_ensembl_gene", "external_gene_name_homo", "external_gene_name_mus")

homo$mouse_homs$getBM_homologs_and_symbols$external_gene_name_homo <- tolower(homo$mouse_homs$getBM_homologs_and_symbols$external_gene_name_homo)

homo$mouse_homs$getBM_homologs_and_symbols$external_gene_name_mus <- tolower(homo$mouse_homs$getBM_homologs_and_symbols$external_gene_name_mus)

homo$mouse_homs$getBM_homologs_and_symbols$mmusculus_homolog_ensembl_gene <- NULL

homo$mouse_homs$with_main_data <- merge(
  x = homo$mouse_homs$getBM_homologs_and_symbols, 
  y = homo$input_t2, 
  by.x = "external_gene_name_mus", 
  by.y = "Gene.symbol")

colnames(homo$mouse_homs$with_main_data)[3:36] <- stringr::str_c(colnames(homo$mouse_homs$with_main_data)[3:36], " main")

homo$mouse_homs$with_main_data_and_seq <- merge(
  x = homo$output$pap_exp$merge, 
  y = homo$mouse_homs$with_main_data, 
  by.x = "Standarised_gene_symbol", 
  by.y = "external_gene_name_homo")

openxlsx::write.xlsx(homo$mouse_homs$with_main_data_and_seq , file = "ptsd vs main dataset.xlsx")




temp <- homo$mouse_homs$with_main_data_and_seq[homo$mouse_homs$with_main_data_and_seq$Standarised_gene_symbol != homo$mouse_homs$with_main_data_and_seq$external_gene_name_mus,]









# homo$mouse_homs$getBM_homologs_and_symbols_nest <- tidyr::nest( dplyr::group_by(homo$mouse_homs$getBM_homologs_and_symbols, external_gene_name_homo) )
# 
# colnames(homo$mouse_homs$getBM_homologs_and_symbols_nest)[1] <- "external_gene_name_homo"
# 
# homo$mouse_homs$getBM_homologs_and_symbols_nest$bestName <- NA
# 
# for (rowNb in seq_along(homo$mouse_homs$getBM_homologs_and_symbols_nest[[1]]) ) {
# 
#   homo$mouse_homs$getBM_homologs_and_symbols_nest$bestName[[rowNb]] <- select_best_geneName(homo$mouse_homs$getBM_homologs_and_symbols_nest[[2]][[rowNb]][[1]])
# }

# test <- subset(homo$mouse_homs$getBM_homologs_and_symbols_nest, homo$mouse_homs$getBM_homologs_and_symbols_nest$external_gene_name == "b4gat1")

# homo$mouse_homs$getBM_homologs_and_symbols$external_gene_name[duplicated(homo$mouse_homs$getBM_homologs_and_symbols$external_gene_name)]

homo$mouse_homs$getBM_homologs_and_symbols_nest$data <- NULL


homo$output$three <- merge(
  x = homo$output$two, 
  y = homo$mouse_homs$getBM_homologs_and_symbols_nest, 
  by = "external_gene_name",
  all.x = T)





















#################
### FUNCTIONs ###

# get_query_for_ncbi_geneID_annotation <- function(
#   char_vec_gene_id_to_query_with, 
#   char_vec_organism, 
#   chr_gene_identifier)# 'Gene name'
# {
#   linker_string_crucial_for_returning_ProbeID_to_proper_form <- paste0('[', chr_gene_identifier, '] AND ') 
#   
#   query_vector <- paste0(char_vec_gene_id_to_query_with, linker_string_crucial_for_returning_ProbeID_to_proper_form, char_vec_organism, '[Organism]') 
#   
#   return(query_vector)
# }



get_query_for_ncbi_geneID_annotation <- function(
  char_vec_gene_id_to_query_with, 
  char_vec_organism, 
  chr_gene_identifier)# 'Gene name'
{
  linker_string_crucial_for_returning_ProbeID_to_proper_form <- paste0(' AND ') 
  
  query_vector <- paste0(char_vec_gene_id_to_query_with, linker_string_crucial_for_returning_ProbeID_to_proper_form, char_vec_organism, '[Organism]') 
  
  return(query_vector)
}




search_for_ids_in_ncbi <- function(query_vector, entrez_db = 'gene', delay_between_queries = 0.5)
{
  ids <- list()
  counter <- 1
  for (id_ in query_vector)
  {
    done <- F
    
    while (done == F) {
      message(id_)
      message(counter)
      
      tryCatch( {
        ids[[counter]] <- rentrez::entrez_search(db = entrez_db, term = id_)
        
        done <- T
        
        counter = counter + 1
        
        Sys.sleep(delay_between_queries)
      },
      error=function(e) {
        Sys.sleep(5)
        warning('Failed to download data, trying again in 5 s')
      })
    }
  }
  return(ids)
}










make_and_write_table_with_original_and_ncbi_ids <- function(
  entrez_gene_search_output, 
  df_original_data)
{
  
  extracted_data <- purrr::map(.x = entrez_gene_search_output, .f = function(x){
    if (length(x$ids) != 0) {
      entrez_id <- x$ids[[1]]
    } else { entrez_id <- NA }
    
    
    query <- stringr::str_split(string = x$QueryTranslation, pattern = ' AND ', simplify = T)
    
    input_id <- stringr::str_remove(string = query[1], pattern = '\\[.*')
    
    organism <- stringr::str_remove(string = tolower(query[2]), pattern = '\\[.*')
    organism <- stringr::str_remove_all(string = organism, pattern = '"')
    
    return(c(entrez_id, input_id, organism))
    
  })
  
  extracted_data_df <- as.data.frame(t(as.data.frame(extracted_data)))
  
  colnames(extracted_data_df) <- c('Gene_ID', 'input_id', 'organism')
  rownames(extracted_data_df) <- seq_along(extracted_data_df[[1]])
  
  extracted_data_df$dummy_exp <- as.numeric(as.factor(extracted_data_df$organism))
  
  extracted_data_df$Gene_ID <- as.character(extracted_data_df$Gene_ID)
  extracted_data_df$input_id <- as.character(extracted_data_df$input_id)
  extracted_data_df$organism <- as.character(extracted_data_df$organism)
  
  return(extracted_data_df)
}












