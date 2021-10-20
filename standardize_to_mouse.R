# The basic method for determining transcript orthologs across species (identical gene data) misses many orthologs, and may create bias in the results against detecting genes with differing datas across species. For example, within the Mouse Genome Database (https://eu-west-1.protection.sophos.com?d=jax.org&u=aHR0cDovL3d3dy5pbmZvcm1hdGljcy5qYXgub3Jn&i=NjA4MmJmMDQzNTk5ZWQ0NzczNjRkZGVh&t=N256bkVTMjFqd2FocHZIVmVPOTgzdnVJZGNYR0RHYktuQlk4Y2JLLzBOQT0=&h=0dfc16b099bc44c0b010704a8e9d74bb), around 75% of gene symbols are shared by mice and rats (15,741/20,638 genes with official symbols in each species), and the symbols that are different are often genes that are generally less well-studied. Re-running the analysis to consider orthologs with different datas in different species would make this paper stronger, and could potentially be a revision to the code completed at the same time as updating the p-values (above). However, the benefit would be incremental (at most, increasing the # of results by 25%) and may not be worth the investment of effort if there is a time crunch. A less time intensive approach would be to quantify the impact of the bias in the analysis by determining what percentage of the transcripts that are "underrepresented" in the differential expression findings are genes that have different gene symbols in different species (e.g., using the ortholog database above) and include that in the discussion section.
library(Hmisc)
library(magrittr)
options(scipen = 4)
source("https://github.com/AdrianS85/helper_R_functions/blob/master/little_helpers.R")
source("no_fomo_func.R")
rm(list = ls(pattern = '(.*)(temp)|(test)(.*)'))
# data <- NULL
# load("medianed_final_good_dataset")
# load("final_good_dataset_1_and_2")

save(data, file = "data.save")
# load("data.save")


# I load GJ-given Supplementary Dataset 1.xlsx, Supplementary Table 1.xlsx, Supplementary Table 2.xlsx
data <- list(
  "mart" = list(
    "mus" = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"),
    "rat" = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl"),
    "saimiri" = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "sbboliviensis_gene_ensembl")),
  "inputs" = list(
    "DS1" = openxlsx::read.xlsx(xlsxFile = "exel/Supplementary Dataset 1.xlsx"),
    "T1" = openxlsx::read.xlsx(xlsxFile = "exel/Supplementary Table 1.xlsx"),
    "T2" = openxlsx::read.xlsx(xlsxFile = "exel/Supplementary Table 2.xlsx") ### !!! not used
  ))






###########
### ST1 ###

# I remove experiment 27.1 from `Supplementary Table 1` and save it as "Supplementary Table 1 v2.xlsx"
data$output$T1 <- subset(
  data$input$T1, 
  subset = data$input$T1[["Paper/Group.code"]] != "27.1")

colnames(data$output$T1) <- stringr::str_replace_all(colnames(data$output$T1), pattern = "\\.", replacement = " ")

openxlsx::write.xlsx(data$output$T1, file = "Supplementary Table 1 v2.xlsx")

### ST1 ###
###########






#######################
### DS1 PREPARATION ###

# I remove experiment 27_1 from `Supplementary Dataset 1`
data$hs_removed$DS1 <- subset(
  data$input$DS1, 
  subset = data$input$DS1[["Paper/Group.code"]] != "27_1")



### !!! Add gene names that were not detected via annotation (Standarised_gene_symbol), but are present in Gene_symbol column Check: LOC102547879,LOC100910801,RGD1561143 gene_symbol
# I add non-filtered gene symbols to standrdized list
data$hsRm_SymbolToStd$DS1 <- data$hs_removed$DS1

data$hsRm_SymbolToStd$DS1$old_Standarised_gene_symbol <- data$hsRm_SymbolToStd$DS1$Standarised_gene_symbol

for (rowNb in seq_along(data$hsRm_SymbolToStd$DS1[[1]]) ) {
  
  if( is.na(data$hsRm_SymbolToStd$DS1$Standarised_gene_symbol)[[rowNb]] ){
    
    data$hsRm_SymbolToStd$DS1$Standarised_gene_symbol[[rowNb]] <- tolower( select_best_geneName_wrapper_for_single_string(data$hsRm_SymbolToStd$DS1$Gene_symbol[[rowNb]] ) )
    
  } else if ( !is.na(data$hsRm_SymbolToStd$DS1$Standarised_gene_symbol)[[rowNb]] ) {
    
    data$hsRm_SymbolToStd$DS1$Standarised_gene_symbol[[rowNb]] <- data$hsRm_SymbolToStd$DS1$old_Standarised_gene_symbol[[rowNb]]
  }
}



# Got all orthologs (in ensembl id form) from ensembl for both mus vs rat/saimiri and rat vs mus and saimiri vs mus
data$orthologs$all_ortholog_ids <- purrr::map2(
  .x = list(
  "mus" = c("rnorvegicus_homolog_ensembl_gene", "sbboliviensis_homolog_ensembl_gene"), 
  "rat" = "mmusculus_homolog_ensembl_gene", 
  "saimiri" = "mmusculus_homolog_ensembl_gene"), 
  .y = data$mart,
  .f = function(orth, mart_){
    
    biomaRt::getBM(
      attributes = c("external_gene_name", "ensembl_gene_id", orth),
      mart = mart_)
  }) 



# Kept only those that had gene symbol for the first species (cause if transcript doesnt have mouse name, than we cannot normalize to this name that doesnt exist)
data$orthologs$all_ortholog_ids_clean <- purrr::map2(
  .x = data$orthologs$all_ortholog_ids,
  .y = names(data$orthologs$all_ortholog_ids),
  .f = function(orth, orth_name){
    
    colnames(orth) <- stringr::str_c(colnames(orth), "_", orth_name)

    orth <- dplyr::filter(.data = orth, orth[[1]] != "")
    
    orth[[1]] <- tolower(orth[[1]])
    
    return(orth)
  })



# We extract rows specific to given species
data$orthologs$speciesSpecificSymbols <- purrr::map(
  .x = list("mice" = "mice", "rats" = "rats", "squirrelmonkeys" = "squirrelmonkeys"), 
  .f = function(species){
    data.frame("Standarised_gene_symbol" = unique(
        dplyr::filter(data$hsRm_SymbolToStd$DS1, data$hsRm_SymbolToStd$DS1$Species %in% species)$Standarised_gene_symbol )
    )
  })



# We subset only such ortholog data, which have gene names in our data for appropriate species - for all rat genes and their orthologs we keep only such genes and their orthologs that are present in rat-derived data in our dataset
data$orthologs$speciesSpecificSymbols_andOrths <- purrr::map2(
  .x = data$orthologs$speciesSpecificSymbols, 
  .y = data$orthologs$all_ortholog_ids_clean, 
  .f = function(specificSymbols, orths){
    
    merge(
      specificSymbols, 
      orths,
      by.x = "Standarised_gene_symbol",
      by.y = colnames(orths)[[1]])
  })



# We merge ensembl ids for mouse and rat/saimiri in two directions - mouse id vs it rat/saimiri homologues and rat/saimir ids vs its mouse homogoues
data$orthologs$merge$musRat1 <- merge(
  data$orthologs$speciesSpecificSymbols_andOrths$mice, 
  data$orthologs$speciesSpecificSymbols_andOrths$rats,
  by.x = "ensembl_gene_id_mus",
  by.y = "mmusculus_homolog_ensembl_gene_rat")
data$orthologs$merge$musRat1$sbboliviensis_homolog_ensembl_gene_mus <- NULL
colnames(data$orthologs$merge$musRat1)[c(2,4)] <- c("Standarised_gene_symbol_mus", "Standarised_gene_symbol_rat")

data$orthologs$merge$musRat2 <- merge(
  data$orthologs$speciesSpecificSymbols_andOrths$mice, 
  data$orthologs$speciesSpecificSymbols_andOrths$rats,
  by.x = "rnorvegicus_homolog_ensembl_gene_mus",
  by.y = "ensembl_gene_id_rat")
data$orthologs$merge$musRat2$sbboliviensis_homolog_ensembl_gene_mus <- NULL
colnames(data$orthologs$merge$musRat2)[c(2,4)] <- c("Standarised_gene_symbol_mus", "Standarised_gene_symbol_rat")

data$orthologs$merge$musSai1 <- merge(
  data$orthologs$speciesSpecificSymbols_andOrths$mice, 
  data$orthologs$speciesSpecificSymbols_andOrths$squirrelmonkeys,
  by.x = "ensembl_gene_id_mus",
  by.y = "mmusculus_homolog_ensembl_gene_saimiri")
data$orthologs$merge$musSai1$rnorvegicus_homolog_ensembl_gene_mus <- NULL
colnames(data$orthologs$merge$musSai1)[c(2,4)] <- c("Standarised_gene_symbol_mus", "Standarised_gene_symbol_sai")

data$orthologs$merge$musSai2 <- merge(
  data$orthologs$speciesSpecificSymbols_andOrths$mice, 
  data$orthologs$speciesSpecificSymbols_andOrths$squirrelmonkeys,
  by.x = "sbboliviensis_homolog_ensembl_gene_mus",
  by.y = "ensembl_gene_id_saimiri")
data$orthologs$merge$musSai2$rnorvegicus_homolog_ensembl_gene_mus <- NULL
colnames(data$orthologs$merge$musSai2)[c(2,4)] <- c("Standarised_gene_symbol_mus", "Standarised_gene_symbol_sai")



# Then we remove rows that are duplicates in these 3 columns: mouse id, symbol and rat/saimir symbol
data$orthologs$mergeClean <- purrr::map(
  .x = data$orthologs$merge, 
  .f = function(ds){
    
    ds[[3]] <- NULL
    ds[[4]] <- NULL
    ds <- unique(ds)
    
    return(ds)
  })

data$orthologs$mergeCleanFull$musRat <- unique( rbind(data$orthologs$mergeClean$musRat1[,c(2,3)], data$orthologs$mergeClean$musRat2[,c(2,3)]) )

data$orthologs$mergeCleanFull$musSai <- data$orthologs$mergeClean$musSai1[,c(2,3)]



# We keep only rows for which mouse and rat/saimir name are different
data$orthologs$mergeCleanFullDiff <- purrr::map(
  .x = data$orthologs$mergeCleanFull, 
  .f = function(ds){
    
    ds$areDifferent <- ifelse(
      test = ds[[1]] == ds [[2]], 
      yes = F, 
      no = T)

    return(ds)
  })

data$orthologs$mergeCleanFullDiffClean <- purrr::map(
  .x = data$orthologs$mergeCleanFullDiff, 
  .f = function(ds){
    
    dplyr::filter(ds, ds$areDifferent == T)[,-3]
  })



# We keep only rat and saimiri ortholog gene symbols that ARE NOT THE SAME as gene symbols already present in our mouse-derived data. That is because this rat/saimir gene was already merged with a mouse gene 
data$orthologs$mergeCleanFullDiffCleanConfusingNames <- purrr::map(
  .x = data$orthologs$mergeCleanFullDiffClean,
  .f = function(ds){
    
    ds$alreadyInMus <- ifelse(
      test = ds[[2]] %in% data$orthologs$speciesSpecificSymbols$mice[[1]], 
      yes = T, 
      no = F)
    
    return(ds)
  })

data$orthologs$mergeCleanFullDiffCleanConfusingNamesClean <- purrr::map(
  .x = data$orthologs$mergeCleanFullDiffCleanConfusingNames,
  .f = function(ds){
    
    dplyr::filter(ds, ds$alreadyInMus == F)[,-3]
  })



# Some rat/saimir gene names have 2 or more homologoues in mus. Here we select best homologoues gene name to represent given rat/saimiri gene
data$orthologs$bestNames <- purrr::map(
  .x = data$orthologs$mergeCleanFullDiffCleanConfusingNamesClean,
  .f = function(ds){
    
    colname_ds <- colnames(ds)[2]
    
    ds <- tidyr::nest( dplyr::group_by(ds, eval(parse(text = colnames(ds)[2])) ) )
    
    colnames(ds)[1] <- colname_ds
    
    ds$bestName <- NA
    
    for (rowNb in seq_along(ds[[1]]) ) {
      
      ds$bestName[[rowNb]] <- select_best_geneName(ds[[2]][[rowNb]][[1]])
    }
    
    return(ds)
  })

data$orthologs$bestNamesClean <- purrr::map(
  .x = data$orthologs$bestNames,
  .f = function(ds){
    
    return(ds[,c(1,3)])
  })




data$hsRm_namesChanged$DS1 <- data$hsRm_SymbolToStd$DS1

data$hsRm_namesChanged$mice <- subset(data$hsRm_namesChanged$DS1, data$hsRm_namesChanged$DS1$Species == "mice")
data$hsRm_namesChanged$mice$old_Standarised_gene_symbol_2 <- data$hsRm_namesChanged$mice$Standarised_gene_symbol

data$hsRm_namesChanged$rats <- subset(data$hsRm_namesChanged$DS1, data$hsRm_namesChanged$DS1$Species == "rats")

data$hsRm_namesChanged$squirrelmonkeys <- subset(data$hsRm_namesChanged$DS1, data$hsRm_namesChanged$DS1$Species == "squirrelmonkeys")




data$hsRm_namesChanged$namechange <- purrr::map2(
  .x = list("rats" = data$hsRm_namesChanged$rats, "squirrelmonkeys" = data$hsRm_namesChanged$squirrelmonkeys),
  .y =  list("rats" = data$orthologs$bestNamesClean$musRat, "squirrelmonkeys" = data$orthologs$bestNamesClean$musSai),
  .f = function(ds, namesExchange){
    
    ds$old_Standarised_gene_symbol_2 <- ds$Standarised_gene_symbol
    ds$Standarised_gene_symbol <- NA
    
    for (rowNb in seq_along(ds[[1]])) {
      
      ds$Standarised_gene_symbol[[rowNb]] <- ifelse(
        test = ds$old_Standarised_gene_symbol_2[[rowNb]] %in% namesExchange[[1]], 
        yes = recode_values_based_on_key(
          to_recode_chrvec = ds$old_Standarised_gene_symbol_2[[rowNb]], 
          replace_this_chrvec = namesExchange[[1]], 
          with_this_chrvec = namesExchange$bestName), 
        no = ds$old_Standarised_gene_symbol_2[[rowNb]])
    }
    
    return(ds)
  })

data$hsRm_namesChanged$DS1_updated <- rbind(data$hsRm_namesChanged$mice, data$hsRm_namesChanged$namechange$rats, data$hsRm_namesChanged$namechange$squirrelmonkeys)

openxlsx::write.xlsx(data$hsRm_namesChanged$DS1_updated, file = "Supplementary Dataset 1 v3.xlsx")

### DS1 PREPARATION ###
#######################







data$datasets$medianed_final_good_dataset <- dplyr::summarize(
  dplyr::group_by(data$hsRm_namesChanged$DS1_updated, `Paper/Group.code`, Standarised_gene_symbol), 
  logFC_median = median(logFC))


data$datasets$spread_medianed_final_good_dataset <- tidyr::spread(
  data = data$datasets$medianed_final_good_dataset, 
  key = `Paper/Group.code`, 
  value = logFC_median)


data$datasets$bool_entries_with_3_or_more_values <- get_subset_vector_for_entries_with_3_or_more_values_per_paper(final_good_dataset_ = data$hsRm_namesChanged$DS1_updated)


data$datasets$at_least_in_3_papers_spread_med_fin_g_ds <- subset(
  x = data$datasets$spread_medianed_final_good_dataset, 
  subset = data$datasets$bool_entries_with_3_or_more_values)


data$datasets$at_least_in_3_papers_spread_med_fin_g_ds <- subset(
  data$datasets$at_least_in_3_papers_spread_med_fin_g_ds, 
  subset = !is.na(data$datasets$at_least_in_3_papers_spread_med_fin_g_ds$Standarised_gene_symbol))


# Are any columns empty?
test <- data$datasets$at_least_in_3_papers_spread_med_fin_g_ds
test[is.na(test)] <- 0
test <- t(purrr::map_df(.x = test[,2:ncol(test)], .f = sum))
test2 <- subset(test, test == 0)
# Are any columns empty?
### Remove 




################################
### SAVE DATA FOR CLUSTERING ###

data$datasets$for_clustering <- as.data.frame(data$datasets$at_least_in_3_papers_spread_med_fin_g_ds)

data$datasets$for_clustering[is.na(data$datasets$for_clustering)] <- 0

readr::write_tsv(x = data$datasets$for_clustering, 'for_clustering_v5.tsv')

### SAVE DATA FOR CLUSTERING ###
################################








data$datasets$spread_medianed_final_good_dataset <- subset(
  data$datasets$spread_medianed_final_good_dataset, 
  subset = !is.na(data$datasets$spread_medianed_final_good_dataset$Standarised_gene_symbol))

data$datasets$spread_medianed_final_good_dataset[is.na(data$datasets$spread_medianed_final_good_dataset)] <- 0

get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(data$datasets$spread_medianed_final_good_dataset, '')






### SD1 ###
###########




##########
### S2 ###

data$datasets[["ST2 temp"]] <- merge(
  x = exp_number_and_percentage_, 
  y = paper_number_, 
  by = "Standarised_gene_symbol")

data$datasets[["ST2"]] <- merge(
  x = data$datasets[["ST2 temp"]], 
  y = data$input$T2, 
  by.x = "Standarised_gene_symbol",
  by.y = "Gene.symbol",
  all.x = T)

data$datasets[["ST2"]] <- dplyr::select(
  data$datasets[["ST2"]], 
  "Gene symbol" = "Standarised_gene_symbol",
  "Number of reporting papers" = "number",
  "Total number of transcriptomic comparisons with altered expression" = "no_of_exps",
  "Fraction of comparisons with up-regulated expression" = "perc_of_upregulated", "Glucocorticoid-responsive.from.core.list", "Glucocorticoid-responsive.from.extended.list", "Hemoglobin.cluster", "Meningeal.cluster", "Choroid.cluster" 
)

data$datasets[["ST2"]][["Fraction of comparisons with down-regulated expression"]] <- 1 - data$datasets[["ST2"]]$`Fraction of comparisons with up-regulated expression`

data$datasets[["ST2"]][["Number of papers x fraction"]] <- ifelse(
  test = data$datasets[["ST2"]]$`Fraction of comparisons with up-regulated expression` >= data$datasets[["ST2"]]$`Fraction of comparisons with down-regulated expression`,
  yes = data$datasets[["ST2"]]$`Number of reporting papers` * data$datasets[["ST2"]]$`Fraction of comparisons with up-regulated expression`,
  no = data$datasets[["ST2"]]$`Number of reporting papers` * data$datasets[["ST2"]]$`Fraction of comparisons with down-regulated expression`)

colnames(data$datasets[["ST2"]]) <- stringr::str_replace_all(string = colnames(data$datasets[["ST2"]]), pattern = "\\.", replacement = " " )

data$datasets[["ST2"]] <- dplyr::select(
  data$datasets[["ST2"]],
  1:4, 10, 11, 5:9)



data$datasets[["ST2"]]$freq_signif <- 1 - pbinom(
  data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] - 1, 
  size = 223, 
  prob = 0.05)

data$datasets[["ST2"]]$freq_signif_left <- pbinom(
  data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
  size = 223, 
  prob = 0.05)

data$datasets[["ST2"]]$freq_signif_non_cumul <- dbinom(
  data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
  size = 223, 
  prob = 0.05)

data$datasets[["ST2"]]$freq_up <- 1 - pbinom(
  as.integer( data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] * data$datasets[["ST2"]][["Fraction of comparisons with up-regulated expression"]]) - 1, 
  size = data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
  prob = 0.5)

data$datasets[["ST2"]]$freq_up_non_cumul <- dbinom(
  as.integer( data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] * data$datasets[["ST2"]][["Fraction of comparisons with up-regulated expression"]]), 
  size = data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
  prob = 0.5)

data$datasets[["ST2"]]$freq_down <- pbinom(
  as.integer( data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] * data$datasets[["ST2"]][["Fraction of comparisons with up-regulated expression"]] ), 
  size = data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
  prob = 0.5) 






# for (rowNb in seq_along(data$datasets[["ST2"]][[1]])) {
#   data$datasets[["ST2"]]$freq_signif_non_cumul[[rowNb]] <- binom.test(
#     x = data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]][[rowNb]], 
#     n = 223,
#     p = 0.05, 
#     alternative = "two.sided")$p.value
#   
#   data$datasets[["ST2"]]$freq_up_non_cumul[[rowNb]] <- binom.test(
#     x = as.integer( data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]][[rowNb]] * data$datasets[["ST2"]][["Fraction of comparisons with up-regulated expression"]][[rowNb]]), 
#     n = data$datasets[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]][[rowNb]],
#     p = 0.5, 
#     alternative = "greater")$p.value
# }
# 
# test <- binom.test(
#   x = 8, 
#   n = 223,
#   p = 0.05, 
#   alternative = "two.sided")
# 
# dbinom(8, 223, 0.05)

# 1)  w takim bądź razie policz lewostronny i prawostronny skumulowany test dwumianowy dla częstości wykrywania genów przyjmując 223 eksperymenty z których mamy dane po standaryzacji oraz FDR ale tylko dla wyników testu prawostronnego. 
# 
# 3) Czy jesteś w stanie policzyć test dwumianowy nieskumulowany dla tych danych? Jeśli tak to dobrze, jeśli nie to sam to sam uzupełnie. 


temp <- purrr::map2_dfc(
  .x = list(
    "signif" =  data$datasets[["ST2"]]$freq_signif, 
    "up" =  data$datasets[["ST2"]]$freq_up, 
    "down" =  data$datasets[["ST2"]]$freq_down),
  .y =c("signif", "up", "down"),
  .f = function(column, colname_){
    
    colnames_ <- stringr::str_c(c("bon_", "fdr_", "q_", "dq_", "rob_"), colname_)
    colnames_ <- c(colname_, colnames_)
    
    temp <- as.data.frame(column)
    
    temp$bon <- p.adjust(temp[[1]], method = "bonferroni", )
    temp$fdr <- p.adjust(temp[[1]], method = "BH")
    temp$q <- qvalue::qvalue(p = temp[[1]], pfdr = T)$qvalues
    temp$dq <- DiscreteQvalue::DQ(temp[[1]])$q.values
    temp$rob <- robust.fdr(p = temp[[1]], sides = 1, discrete = T)$q
    
    colnames(temp) <- colnames_
    
    return(temp)
  })

data$datasets[["ST2"]] <- cbind(data$datasets[["ST2"]], temp)

openxlsx::write.xlsx(data$datasets$`ST2`, file = "Supplementary Table 2 v4.xlsx")

### S2 ###
##########









##########################################
### GET NUMBERS OF GENES IN ALL PAPERS ###

paper_number_numbers <- paper_number_ %>%
  dplyr::group_by(number) %>%
  dplyr::mutate(nb_of_genes_detected_in_this_nb_of_papers = dplyr::n()) %>%
  dplyr::mutate(percent_of_genes_detected_in_this_nb_of_papers = (dplyr::n()/length(paper_number_[[1]])) ) %>%
  dplyr::select(-Standarised_gene_symbol) %>%
  unique()
readr::write_tsv(paper_number_numbers, file = 'exp_and_paper_numbers/paper_number_numbers.tsv')
### GET NUMBERS OF GENES IN ALL PAPERS ###
### GET NUMBERS OF GENES IN ALL PAPERS DRAWING ###
# library(ggplot2)
paper_number_numbers <- paper_number_numbers[order(paper_number_numbers$number),]

paper_number_numbers_for_plot <- paper_number_numbers
paper_number_numbers_for_plot$nb_of_genes_detected_in_this_nb_of_papers[19] <-sum(paper_number_numbers$nb_of_genes_detected_in_this_nb_of_papers[19:25])
paper_number_numbers_for_plot$percent_of_genes_detected_in_this_nb_of_papers[19] <-sum(paper_number_numbers$percent_of_genes_detected_in_this_nb_of_papers[19:25])
paper_number_numbers_for_plot <- paper_number_numbers_for_plot[-c(20:25),]

plot <- paper_number_numbers_for_plot %>%
  ggplot(aes(x = number, y = nb_of_genes_detected_in_this_nb_of_papers)) +
  geom_col() +
  geom_text(data = paper_number_numbers_for_plot,
            aes(
              label = paste0(round(
                100 * paper_number_numbers_for_plot$percent_of_genes_detected_in_this_nb_of_papers,
                digits = 2
              ), ' %'),
              y = 1000,
              angle = 90
            ),
            size = 4)+
  scale_x_continuous(breaks = seq(1, 19, 1), labels = c(as.character(seq(1, 18, 1)), '19 <'))+
  labs(title = 'number of papers in which gene was detected')+
  xlab('number of papers')+
  ylab('number of genes')

ggplot2::ggsave("plot.png", plot = plot, device = "png")



exp_number_numbers <- exp_number_and_percentage_ %>%
  dplyr::select(Standarised_gene_symbol, no_of_exps)  %>%
  dplyr::group_by(no_of_exps) %>%
  dplyr::mutate(nb_of_genes_detected_in_this_nb_of_exps = dplyr::n()) %>%
  dplyr::mutate(percent_of_genes_detected_in_this_nb_of_exps = (dplyr::n()/length(exp_number_and_percentage_[[1]])) ) %>%
  dplyr::select(-Standarised_gene_symbol) %>%
  unique()
readr::write_tsv(exp_number_numbers, file = 'exp_and_paper_numbers/exp_number_numbers.tsv')

### GET NUMBERS OF GENES IN ALL PAPERS ###
##########################################



#########################
### GET NEW SUBGROUPS ###

load('descriptions_1_and_2')

descriptions_1_and_2 <- subset(descriptions_1_and_2, descriptions_1_and_2$Group_ID != "27_1")


descriptions_1_and_2$Stress_sensitivity_clean[descriptions_1_and_2$Group_ID == "9_1"] <- as.factor("vulnerable")

descriptions_1_and_2 <- descriptions_1_and_2[descriptions_1_and_2$Group_ID %nin% c("61_1", "62_1", "32_1", "32_2"),]


data$subsets$search_for <- list(
  'acute_all' = list('acute', "Stress_duration"), 
  'medium_all' = list('medium', "Stress_duration"), 
  'prolonged_all' = list('prolonged', "Stress_duration"),
  'vulnerable_all' = list('vulnerable', "Stress_sensitivity_clean"),
  'resilient_all' = list('resilient', "Stress_sensitivity_clean"),
  'male_all' = list('male', "Gender_clean"),
  'female_all' = list('female', "Gender_clean")
)


data$subsets$experiments_to_subset <- lapply(
  X = data$subsets$search_for, 
  function(x)
    {
    descriptions_1_and_2 %>%
      dplyr::select(Group_ID, x[[2]]) %>%
      dplyr::filter(descriptions_1_and_2[[ x[[2]] ]] == x[[1]])
})






data$subsets$sets <- purrr::map2(
  .x = data$subsets$experiments_to_subset,
  .y = names(data$subsets$search_for),
  .f = function(x, y) {
    create_subset_of_exps_with_all_papers_wrapper_wrapper(
      final_good_dataset___ = data$hsRm_namesChanged$DS1_updated,
      experiments_to_include_df_with_group_id_col = x,
      save_as_chr_ = y,
      group_id_col_str_ = 'Group_ID'
    ) } )


#16_1, 27_1 - remove these publication from female_all, cause their female status is unsure
data$subsets$sets$female_all <- dplyr::select(.data = data$subsets$sets$female_all, -`16_1`)


#load new, upgraded data for 16_1 experiment, and merge it to female_all as new column. This file there are some NA name values, it should not matter for further analysis, as we 
data$subsets$sets$female_all <- merge(
  x = data$subsets$sets$female_all, 
  y = readr::read_tsv(file = 'roszkowski wydzielone.txt'), 
  by.x = 'Standarised_gene_symbol', 
  by.y = "lower_final_gene_name",
  all = T)


data$subsets$sets$female_all$`16_1`[is.na(data$subsets$sets$female_all$`16_1`)] <- 0


data$subsets$percent <- purrr::map2(
  .x = data$subsets$sets,
  .y = names(data$subsets$sets),
  .f = function(x, y) {
    get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(x, y)
  } )


data$subsets$st1 <- purrr::pmap(
  .l = list(data$subsets$percent, data$subsets$experiments_to_subset, names(data$subsets$percent)),
  .f = function(ds, experiments_for_nb, name){
    
    make_dataset_wrapper(
      exp_number_and_percentage__ = ds$exps, 
      paper_number__ = ds$pap,
      nb_of_exps = length(experiments_for_nb[[1]]),
      ds_name = name)
  })

### GET NEW SUBGROUPS ###
#########################






################
### DURATION ###

data$subsets$descs$t4 <- data.frame(c("", "T1", "T2", "T3"), c("", "Genes with altered expression after acute stress defined as procedures performed during one day", "Genes with altered expression after medium stress that included procedures repeated for at least 2 days but not longer than one week", "Genes with altered expression after prolonged stress that included procedures repeated for at least 8 days."))
colnames(data$subsets$descs$t4) <- c("Ranking lists of genes displaying significant changes in expression after stress", "")


openxlsx::write.xlsx(
  x = list(
    "Contents" = data$subsets$descs$t4,
    "T1 (Acute)" = data$subsets$st1$acute_all,
    "T2 (Medium)" = data$subsets$st1$medium_all,
    "T3 (Prolonged)" = data$subsets$st1$prolonged_all),
  file = "Supplementary Table 4 v3.xlsx")

### DURATION ###
################


###########################
### VUL VS RES ANALYSIS ###

data$subsets$descs$t5 <- data.frame(c("", "T1", "T2", "T3", "T4", "T5"), c("", "Genes with altered expression found only in animals that are vulnerable to stress.", "Genes with altered expression found only in animals that are resistant to stress.", "Genes displaying opposite direction of altered expression vulnerable and resistant animals.", "Ranking list of all genes displaying significant change in expression in vulnerable animals.", "Ranking list of all genes displaying significant change in expression in resistant animals."))
colnames(data$subsets$descs$t5) <- c("Ranking lists of genes displaying significant changes in expression after stress", "")



data$subsets$v_vs_r$merged <- merge(
  x = data$subsets$st1$resilient_all, 
  y = data$subsets$st1$vulnerable_all, 
  by = 'Gene symbol', 
  all = T)

colnames(data$subsets$v_vs_r$merged) <- stringr::str_replace(
  string = colnames(data$subsets$v_vs_r$merged), 
  pattern = "\\.x",
  replacement = " res")

colnames(data$subsets$v_vs_r$merged) <- stringr::str_replace(
  string = colnames(data$subsets$v_vs_r$merged), 
  pattern = "\\.y",
  replacement = " vul")



data$subsets$v_vs_r$vul_only_genes_subset <- ifelse(
  test = !is.na(data$subsets$v_vs_r$merged$`Total number of transcriptomic comparisons with altered expression vul`) & is.na(data$subsets$v_vs_r$merged$`Total number of transcriptomic comparisons with altered expression res`),
  yes = T,
  no = F)

data$subsets$v_vs_r$vul_only_genes <- subset(x = data$subsets$v_vs_r$merged, subset = data$subsets$v_vs_r$vul_only_genes_subset)

data$subsets$v_vs_r$vul_only_genes <- dplyr::select(data$subsets$v_vs_r$vul_only_genes, -c(2:35))



data$subsets$v_vs_r$res_only_genes_subset <- ifelse(
 test = !is.na(data$subsets$v_vs_r$merged$`Total number of transcriptomic comparisons with altered expression res`) & is.na(data$subsets$v_vs_r$merged$`Total number of transcriptomic comparisons with altered expression vul`),
 yes = T,
 no = F)

data$subsets$v_vs_r$res_only_genes <- subset(x = data$subsets$v_vs_r$merged, subset = data$subsets$v_vs_r$res_only_genes_subset)

data$subsets$v_vs_r$res_only_genes <- dplyr::select(data$subsets$v_vs_r$res_only_genes, c(1:35))



data$subsets$v_vs_r$opposite_direction_genes_subset <- ifelse(
  test = (data$subsets$v_vs_r$merged$`Fraction of comparisons with up-regulated expression res` > 0.5 & data$subsets$v_vs_r$merged$`Fraction of comparisons with up-regulated expression vul` < 0.5) | (data$subsets$v_vs_r$merged$`Fraction of comparisons with up-regulated expression res` < 0.5 & data$subsets$v_vs_r$merged$`Fraction of comparisons with up-regulated expression vul` > 0.5),
  yes = T,
  no = F)

data$subsets$v_vs_r$opposite_direction_genes <- subset(x = data$subsets$v_vs_r$merged, subset = data$subsets$v_vs_r$opposite_direction_genes_subset)



openxlsx::write.xlsx(
  x = list(
    "Contents" = data$subsets$descs$t5,
    "T1 (vul.)" = data$subsets$v_vs_r$vul_only_genes,
    "T2 (res.)" = data$subsets$v_vs_r$res_only_genes,
    "T3 (opp.)" = data$subsets$v_vs_r$opposite_direction_genes,
    "T4 (vul. all)" = data$subsets$st1$vulnerable_all,
    "T5 (res. all)" = data$subsets$st1$resilient_all),
  file = "Supplementary Table 5 v3.xlsx")

### VUL VS RES ANALYSIS ###
###########################




###############################
### SEX COMPARISON ANALYSIS ###

# Get only those genes that were present in 3 papers or more
data$subsets$sex$input$female <- subset(x = data$subsets$st1$female_all, subset = data$subsets$st1$female_all$`Number of reporting papers` >= 3) 

colnames(data$subsets$sex$input$female) <- stringr::str_c(colnames(data$subsets$sex$input$female), " f")



data$subsets$sex$input$male <- subset(x = data$subsets$st1$male_all, subset = data$subsets$st1$male_all$`Number of reporting papers` >= 3) 

colnames(data$subsets$sex$input$male) <- stringr::str_c(colnames(data$subsets$sex$input$male), " m")



data$subsets$sex$merged <- merge(
  x = data$subsets$sex$input$male, 
  y = data$subsets$sex$input$female, 
  by.x = 'Gene symbol m',
  by.y = 'Gene symbol f')



data$subsets$sex$opposite_direction_genes_subset <- ifelse(
  test = (data$subsets$sex$merged$`Fraction of comparisons with up-regulated expression f` > 0.5 & data$subsets$sex$merged$`Fraction of comparisons with up-regulated expression m` < 0.5) | (data$subsets$sex$merged$`Fraction of comparisons with up-regulated expression f` < 0.5 & data$subsets$sex$merged$`Fraction of comparisons with up-regulated expression m` > 0.5), 
  yes = T, 
  no = F)

data$subsets$sex$opposite_direction_genes <- subset(x = data$subsets$sex$merged, subset = data$subsets$sex$opposite_direction_genes_subset)

openxlsx::write.xlsx(x = list("opposite_direction_genes values" = data$subsets$sex$opposite_direction_genes), file = "Supplementary Table 6 v3.xlsx")

# data$subsets$sex$stable_direction_genes <- subset(x = data$subsets$sex$merged, subset = !data$subsets$sex$opposite_direction_genes_subset)
# 
# data$subsets$sex$merge_all_f <- merge(x = data$subsets$sex$input$male, y = data$subsets$sex$input$female, by = 'lower_final_gene_name', all.y = T)
# 
# data$subsets$sex$f_only_genes_subset <- ifelse(
#   test = !is.na(data$subsets$sex$merge_all_f$no_of_exps_f) & is.na(data$subsets$sex$merge_all_f$no_of_exps_m), 
#   yes = T, 
#   no = F)
# 
# data$subsets$sex$f_only_genes <- subset(x = data$subsets$sex$merge_all_f, subset = data$subsets$sex$f_only_genes_subset)

### SEX COMPARISON ANALYSIS ###
###############################






####################################
### WHICH GENE NAMES DISAPPEARED ###

data$absent$new <- data$datasets$`ST2`

data$absent$old <- data$inputs$T2

data$absent$old_but_not_new <- data$absent$old[data$absent$old$Gene.symbol %nin% data$absent$new$`Gene symbol`,]

openxlsx::write.xlsx(x = list("opposite_direction_genes values" = data$absent$old_but_not_new), file = "Table 2 old_but_not_new v2.xlsx")

### WHICH GENE NAMES DISAPPEARED ###
####################################






############################
### PROTEIN CODING GENES ###

data$prot_coding$temp_bm_results <- biomaRt::getBM(
  attributes = c('ensembl_peptide_id', "external_gene_name"),
  uniqueRows = T,
  mart = data$mart$mus)

data$prot_coding$temp_bm_results_subset <- subset(
  data$prot_coding$temp_bm_results,
  data$prot_coding$temp_bm_results$ensembl_peptide_id != '')

data$prot_coding$temp_bm_results_subset <- data.frame(
  'external_gene_name' = tolower(unique(data$prot_coding$temp_bm_results_subset$external_gene_name)), 
  'protein' = T)

data$prot_coding$our_genes <- data.frame(
  'external_gene_name' = data[["datasets"]]$ST2$Standarised_gene_symbol, 
  'in_our_data' = T)

data$prot_coding$mergement <- merge(
  x = data$prot_coding$our_genes, 
  y = data$prot_coding$temp_bm_results_subset, 
  by = 'external_gene_name')
length(data$prot_coding$mergement[[1]])/length(data$prot_coding$temp_bm_results_subset[[1]])

### PROTEIN CODING GENES ###
############################






# at least 8 comparisons
#################################
### ADDITIONAL FOR CLUSTERING ###

data$subsets2$search_for <- list(
  'mice' = list('mice', "Species"), 
  'rats' = list('rats', "Species"), 
  'amy' = list('amygdala', "Brain_part_clean"),
  'fc' = list('frontal cortex', "Brain_part_clean"),
  'hp' = list('hippocampus', "Brain_part_clean"),
  'ht' = list('hypothalamus', "Brain_part_clean"),
  'nac' = list('nucleus_accumbens', "Brain_part_clean"),
  'pfc' = list('prefrontal cortex', "Brain_part_clean"),
  'striatum' = list('striatum', "Brain_part_clean")
)


data$subsets2$experiments_to_subset <- lapply(
  X = data$subsets2$search_for, 
  function(x)
  {
    descriptions_1_and_2 %>%
      dplyr::select(Group_ID, x[[2]]) %>%
      dplyr::filter(descriptions_1_and_2[[ x[[2]] ]] == x[[1]])
  })






data$subsets2$sets <- purrr::map2(
  .x = data$subsets2$experiments_to_subset,
  .y = names(data$subsets2$search_for),
  .f = function(x, y) {
    create_subset_of_exps_with_all_papers_wrapper_wrapper(
      final_good_dataset___ = data$hsRm_namesChanged$DS1_updated,
      experiments_to_include_df_with_group_id_col = x,
      save_as_chr_ = y,
      group_id_col_str_ = 'Group_ID'
    ) } )


data$subsets2$percent <- purrr::map2(
  .x = data$subsets2$sets,
  .y = names(data$subsets2$sets),
  .f = function(x, y) {
    get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(x, y)
  } )


data$subsets2$st1 <- purrr::pmap(
  .l = list(data$subsets2$percent, data$subsets2$experiments_to_subset, names(data$subsets2$percent)),
  .f = function(ds, experiments_for_nb, name){
    
    make_dataset_wrapper(
      exp_number_and_percentage__ = ds$exps, 
      paper_number__ = ds$pap,
      nb_of_exps = length(experiments_for_nb[[1]]),
      ds_name = name)
  })

### ADDITIONAL FOR CLUSTERING ###
#################################

















#################
### FUNCTIONS ###


recode_values_based_on_key <- function(to_recode_chrvec, replace_this_chrvec, with_this_chrvec)
{
  assertthat::assert_that(length(replace_this_chrvec) == length(with_this_chrvec), msg = 'replace_this_chrvec and with_this_chrvec arguments need to be the same lenght. This is because every element in first vector will be recoded as corresponding element in the second vector')
  assertthat::assert_that(!any(duplicated(replace_this_chrvec)), msg = 'replace_this_chrvec argument includes duplicated values. This cannot be, because we use one specific value to be changed into another specific value')
  
  replace_this_chrvec <- as.character(replace_this_chrvec)
  with_this_chrvec <- as.character(with_this_chrvec)
  to_recode_chrvec <- as.character(to_recode_chrvec)
  
  key_value_pairs <-
    data.frame(replace_this_chrvec, with_this_chrvec)
  
  to_recode_chrvec <- data.frame(to_recode_chrvec)
  
  to_recode_chrvec$order <- as.integer(rownames(to_recode_chrvec))
  
  result_chrvec <-
    merge(
      x = to_recode_chrvec,
      y = key_value_pairs,
      by.x = 'to_recode_chrvec',
      by.y = 'replace_this_chrvec',
      all.x = T
    )
  
  result_chrvec <- result_chrvec[order(result_chrvec$order),]
  
  return(as.character(result_chrvec$with_this_chrvec)) ### !!! ADDED LATER - as.character()
}



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



robust.fdr<-function(p,sides=1,p2=1-p,discrete=F,use8=T)
{
  library(MASS)
  
  m<-length(p)
  ord<-order(p)
  pcumtab<-cumsum(table(p))/m
  F05<-mean(p<=0.5)	
  edf<-approx(as.numeric(names(pcumtab)),pcumtab,xout=p,rule=2)$y
  if (sides==2)
  {
    pi<-min(1,2*mean(p))
    loc.fdr<-pi*p/edf
  }
  else
  {
    p.new<-2*(p*(p<=p2)+p2*(p>p2))
    pi<-min(1,2*mean(p.new))
    if (discrete) 
    {
      if (use8) pi<-min(1,8*mean(p.new))
      else
      {
        lam<-max(p[p<=0.5])
        k<-1/(lam^2+(0.5-lam)^2)
        pi<-min(k*mean(p.new),1)
        
      }
    }
    loc.fdr<-pi*p/edf
    loc.fdr[p>0.5]<-(0.5*pi+edf[p>0.5]-F05)/edf[p>0.5]
  }
  up<-unique(p)
  ufdr<-approx(p[ord],loc.fdr[ord],xout=up)$y
  res<-lts.rank.regr(up,ufdr)
  fdr<-approx(up,res$yhat,xout=p)$y
  fdr[fdr>1]<-1
  fp<-m*edf*fdr
  fn<-m*(F05-edf-pi*(0.5-p))*(p<0.5)*(sides==1)+(sides==2)*m*(1-edf-pi*(1-p))
  fn[fn<0]<-0
  te<-fp+fn
  q<-fdr
  q[rev(ord)]<-cummin(fdr[rev(ord)])
  
  return(list(p=p,fdr=fdr,q=q,cdf=edf,loc.fdr=loc.fdr,fp=fp,fn=fn,te=te,pi=pi,ord=ord))
}

lts.rank.regr<-function(x,y)
  
{
  rx<-rank(x)
  ry<-rank(y)
  rfit<-ltsreg(rx,ry)
  ryhat<-rfit$fitted.values
  yhat<-approx(ry,y,xout=ryhat,rule=2)$y
  return(list(x=x,y=y,yhat=yhat))
}








# create_subset_of_exps_with_at_least_3_papers_wrapper_wrapper <- function(final_good_dataset___, experiments_to_include_df_with_group_id_col, save_as_chr_, group_id_col_str_ = 'Group_ID')
# {
#   include_these_exps <- experiments_to_include_df_with_group_id_col[[group_id_col_str_]]
#   
#   temp_ <- create_subset_of_exps_with_at_least_3_papers_wrapper(final_good_dataset__ = final_good_dataset___, experiments_to_include_ = include_these_exps, save_as_chr = save_as_chr_)
#   temp_[is.na(temp_)] <- 0
#   
#   readr::write_tsv(x = temp_, path = , paste0('for_clustering_', save_as_chr_, '.tsv'))
#   save(temp_, file = paste0('for_clustering_', save_as_chr_))
#   assign(x = save_as_chr_, value = temp_, envir = globalenv())
# }
# 
# 
# 
# 
# 
# create_subset_of_exps_with_at_least_3_papers_wrapper <- function(final_good_dataset__ = final_good_dataset, experiments_to_include_ = experiments_to_include, save_as_chr)
# {
#   subset_matrix <-subset(final_good_dataset__, subset = final_good_dataset__$Experiment %in% experiments_to_include_)
#   
#   subset_matrix <- dplyr::group_by(subset_matrix, Experiment, lower_final_gene_name)
#   
#   medianed_subset_matrix <-  dplyr::summarize(logFC_median = median(logFC))
#   
#   spread_medianed_subset_matrix <- tidyr::spread(data = medianed_subset_matrix, key = Experiment, value = logFC_median)
#   
#   bool_entries_with_3_or_more_values_subset_matrix <- get_subset_vector_for_entries_with_3_or_more_values_per_paper(final_good_dataset_ = subset_matrix)
#   
#   at_least_in_3_papers_spread_med_sub_mat <- subset(x = spread_medianed_subset_matrix, subset = bool_entries_with_3_or_more_values_subset_matrix)
#   
#   return(at_least_in_3_papers_spread_med_sub_mat)
# }









create_subset_of_exps_with_all_papers_wrapper_wrapper <- function(final_good_dataset___, experiments_to_include_df_with_group_id_col, save_as_chr_, group_id_col_str_ = 'Group_ID')
{
  include_these_exps <- experiments_to_include_df_with_group_id_col[[group_id_col_str_]]
  
  temp_ <- create_subset_of_exps_with_all_papers_wrapper(final_good_dataset__ = final_good_dataset___, experiments_to_include_ = include_these_exps, save_as_chr = save_as_chr_)

  temp_[is.na(temp_)] <- 0
  
  readr::write_tsv(x = temp_, path = , paste0('for_clustering_', save_as_chr_, '.tsv'))
  
  return(temp_)
}









create_subset_of_exps_with_all_papers_wrapper <- function(final_good_dataset__ = final_good_dataset, experiments_to_include_ = experiments_to_include, save_as_chr)
{
  subset_matrix <-subset(final_good_dataset__, subset = final_good_dataset__$`Paper/Group.code` %in% experiments_to_include_)
  
  medianed_subset_matrix <- subset_matrix %>%
    dplyr::group_by(`Paper/Group.code`, Standarised_gene_symbol) %>%
    dplyr::summarize(logFC_median = median(logFC))
  
  spread_medianed_subset_matrix <- tidyr::spread(data = medianed_subset_matrix, key = `Paper/Group.code`, value = logFC_median)
  
  spread_medianed_subset_matrix <- dplyr::filter(spread_medianed_subset_matrix, !is.na(Standarised_gene_symbol))
  
  return(spread_medianed_subset_matrix)
}









get_subset_vector_for_entries_with_3_or_more_values_per_paper <- function(final_good_dataset_){
  temp <- final_good_dataset_ %>%
    dplyr::select(`Paper.code`, Standarised_gene_symbol) %>%
    unique()
  
  temp$present <- T
  
  spread_temp <- tidyr::spread(data = temp, key = `Paper.code`, value = present)
  
  entries_with_3_or_more_values <- purrrlyr::by_row(.d = spread_temp[,-1], .collate = "rows", ..f = function(x) {sum(!is.na(x)) >= 3})$.out
  
  return(entries_with_3_or_more_values)
}




get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper <- function(spread_df_, name_of_df_)
{
  if (!dir.exists('exp_and_paper_numbers')) {
    dir.create('exp_and_paper_numbers')
  }
  
  temp_e <- get_number_and_percentage_of_directionality_of_exp_first_column_names(spread_df_)
  save(temp_e, file = paste0('exp_number_and_percentage_', name_of_df_))
  
  temp_e <- subset(temp_e, temp_e$no_of_exps != 0)
  readr::write_tsv(temp_e,
                   paste0('exp_and_paper_numbers/exp_number_and_percentage_', name_of_df_, '.tsv'))
  
  temp_p <- get_number_and_percentage_of_directionality_of_paper_first_column_names(spread_df_)
  readr::write_tsv(temp_p, paste0('exp_and_paper_numbers/paper_number_', name_of_df_, '.tsv'))
  
  return(list("exps" = temp_e, "pap" = temp_p))
  
}




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
  
  temp_ <- tidyr::gather(data = spread_data_, key = "exp", value = "logFC", na.rm = T, -Standarised_gene_symbol)
  temp_ <- temp_[,-3]
  
  temp_$paper <- stringr::str_remove(string = temp_$exp, pattern = '_.*')
  temp_$exp <- NULL
  
  temp_ <- temp_ %>%
    unique() %>%
    dplyr::select(Standarised_gene_symbol) %>%
    dplyr::group_by(Standarised_gene_symbol) %>%
    dplyr::summarise(number = dplyr::n())
  
  return(temp_)
}












make_dataset_wrapper <- function(exp_number_and_percentage__, paper_number__, nb_of_exps, ds_name)
{
  ds_ <- list()
  
  ds_[["ST2 temp"]] <- merge(
    x = exp_number_and_percentage__, 
    y = paper_number__, 
    by = "Standarised_gene_symbol")
  
  ds_[["ST2"]] <- merge(
    x = ds_[["ST2 temp"]], 
    y = data$input$T2, 
    by.x = "Standarised_gene_symbol",
    by.y = "Gene.symbol",
    all.x = T)
  
  ds_[["ST2"]] <- dplyr::select(
    ds_[["ST2"]], 
    "Gene symbol" = "Standarised_gene_symbol",
    "Number of reporting papers" = "number",
    "Total number of transcriptomic comparisons with altered expression" = "no_of_exps",
    "Fraction of comparisons with up-regulated expression" = "perc_of_upregulated", "Glucocorticoid-responsive.from.core.list", "Glucocorticoid-responsive.from.extended.list", "Hemoglobin.cluster", "Meningeal.cluster", "Choroid.cluster" 
  )
  
  ds_[["ST2"]][["Fraction of comparisons with down-regulated expression"]] <- 1 - ds_[["ST2"]]$`Fraction of comparisons with up-regulated expression`
  
  ds_[["ST2"]][["Number of papers x fraction"]] <- ifelse(
    test = ds_[["ST2"]]$`Fraction of comparisons with up-regulated expression` >= ds_[["ST2"]]$`Fraction of comparisons with down-regulated expression`,
    yes = ds_[["ST2"]]$`Number of reporting papers` * ds_[["ST2"]]$`Fraction of comparisons with up-regulated expression`,
    no = ds_[["ST2"]]$`Number of reporting papers` * ds_[["ST2"]]$`Fraction of comparisons with down-regulated expression`)
  
  colnames(ds_[["ST2"]]) <- stringr::str_replace_all(string = colnames(ds_[["ST2"]]), pattern = "\\.", replacement = " " )
  
  ds_[["ST2"]] <- dplyr::select(
    ds_[["ST2"]],
    1:4, 10, 11, 5:9)
  
  
  
  ds_[["ST2"]]$freq_signif <- 1 - pbinom(
    ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] - 1, 
    size = nb_of_exps, 
    prob = 0.05)
  
  ds_[["ST2"]]$freq_signif_left <- pbinom(
    ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
    size = nb_of_exps, 
    prob = 0.05)
  
  ds_[["ST2"]]$freq_signif_non_cumul <- dbinom(
    ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
    size = nb_of_exps, 
    prob = 0.05)
  
  ds_[["ST2"]]$freq_up <- 1 - pbinom(
    as.integer( ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] * ds_[["ST2"]][["Fraction of comparisons with up-regulated expression"]]) - 1, 
    size = ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
    prob = 0.5)
  
  ds_[["ST2"]]$freq_up_non_cumul <- dbinom(
    as.integer( ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] * ds_[["ST2"]][["Fraction of comparisons with up-regulated expression"]]), 
    size = ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
    prob = 0.5)
  
  ds_[["ST2"]]$freq_down <- pbinom(
    as.integer( ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]] * ds_[["ST2"]][["Fraction of comparisons with up-regulated expression"]] ), 
    size = ds_[["ST2"]][["Total number of transcriptomic comparisons with altered expression"]], 
    prob = 0.5) 
  
  
  temp <- purrr::map2_dfc(
    .x = list(
      "signif" =  ds_[["ST2"]]$freq_signif, 
      "up" =  ds_[["ST2"]]$freq_up, 
      "down" =  ds_[["ST2"]]$freq_down),
    .y =c("signif", "up", "down"),
    .f = function(column, colname_){
      
      colnames_ <- stringr::str_c(c("bon_", "fdr_", "q_", "dq_", "rob_"), colname_)
      colnames_ <- c(colname_, colnames_)
      
      temp <- as.data.frame(column)
      
      temp$bon <- p.adjust(temp[[1]], method = "bonferroni", )
      temp$fdr <- p.adjust(temp[[1]], method = "BH")
      
      temp$q <- tryCatch(
        {
          qvalue::qvalue(p = temp[[1]], pfdr = T)$qvalues
        },
        error=function(cond) {
          qvalue::qvalue(p = temp[[1]], pfdr = T, lambda = 0)$qvalues
        } )
      
      temp$dq <- DiscreteQvalue::DQ( temp[[1]] )$q.values
      temp$rob <- robust.fdr(p = temp[[1]], sides = 1, discrete = T)$q
      
      colnames(temp) <- colnames_
      
      return(temp)
    })
  
  ds_[["ST2"]] <- cbind(ds_[["ST2"]], temp)
  
  return(ds_[["ST2"]])
}
