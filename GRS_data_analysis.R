source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# rm(list = ls(pattern = 'temp.*|test.*'))
load('final_good_dataset_1_and_2') #No genes for experiments 32_1, 32_2, 61_1, 61_2!!
load('reformated_raw_dataset_1_and_2')
load('descriptions_1_and_2') #Includes descriptions for experiments 32_1, 32_2, 61_1. Also, later columns from experiments from second batch are moved one column to the right. This means that 1 resilient and 1 vulnerable group are not included in the analysis. Ive repaired this, but the stress/sensitivity analysis that I made before are obsolete, I removed them and did not made them back

####### SET NAMES OF CLUSTERS OF INTEREST ####### 
groups_of_genes <- list('hemoglobin_cluster' = c('alas2', 'hbb-bt', 'hba-a1', 'hbb-bs', 'hba-a2', 'ch25h', 'lcn2', 'lrg1', 's100a8', 's100a9'), 
                        'choroid_cluster' = c('cldn1',	'epn3',	'msx1',	'col8a2',	'lbp',	'ace',	'clic6',	'mfrp',	'krt8',	'drc7',	'ecrg4',	'prr32',	'aqp1',	'col8a1',	'steap1',	'f5',	'enpp2',	'trpv4',	'otx2',	'folr1',	'cldn2',	'kcne2',	'tmem72',	'slc4a5',	'kl',	'sostdc1',	'ttr',	'col9a3',	'slc39a4',	'sema3b',	'prlr',	'cox8b',	'oca2',	'slc2a12',	'igf2',	'igfbp2',	'slc13a4',	'pcolce',	'wdr86'), 
                        'early_genes_cluster' = c('cdkn1a',	'dusp1',	'egr1',	'egr2',	'fosb',	'fos',	'arc',	'irs2',	'junb',	'midn',	'b4galt1',	'fam107a',	'coq10b',	'gch1',	'csrnp1',	'per1',	'sgk1',	'ddit4',	'tsc22d3',	'pnpla2'), 
                        'additional_1' = c('col1a1', 'lyz2', 'cytl1', 'alx4', 'efemp1', 'col9a2', 'car13', 'h2-q1', 'wnt16', 'bmp4', 'aebp1', 'aldh1a2', 'fam180a', 'slc6a13', 'ptgds', 'slc22a6', 'slc13a3', 'fgfbp1', 'h2-aa', 'mpzl2', 'fmod', 'mrgprf', 'cdh1', 'prdm6', 'itih2', 'asgr1', 'wnt6', 'slc22a2', 'crabp2', 'slc47a1', 'slc6a12', 'mgp'), 
                        'additional_2' = c('cd4', 'drd2', 'gpr88', 'gm14703', 'gm15852', 'xlr3b', 'hist1h2af', 'i830134h01rik', 'mt-nd3', 'smco3'), 
                        'additional_4' = c('btg2', 'ier2', 'atf3', 'apold1', 'ccn1', 'gimap6', 'tagap', 'tmem252', 'lrrc32'), 
                        'additional_5' = c('atp5l', 'fbl', 'h2afz', 'a530032d15rik', 'gapdh', 'rps13', 'gm6180', 'fau', 'mrpl23', 'btg3', 'rpl23a', 'rpl17', 'rps23', 'rps27', 'ppia', 'rpl28', 'sec61g', 'rpl27', 'lsm7', 'pin4', 'rpl9', 'prl8a1', 'rps16', 'hnrnpa3', 'gm14401', 'rpl10a', 'rpl27a', 'rpl31', 'rpl5', 'rps27a', 'rps12', 'olfr613', 'sp140', 'rpl11', 'rpl37', 'rpl38'),
                        '88_best' = c('ddit4','errfi1','klf9','bcl6','fkbp5','mt2','mt2a','nfkbia','pdk4','adcy9','cxxc5','dusp1','eva1a','litaf','nedd9','rhob','sgk1','sult1a1','tiparp','aldoc','arhgef3','arl4d','bcl6b','cables1','calm2','ccnd1','cdkn1a','cdo1','chst1','cyp7b1','ehd3','fzd1','gab1','gap43','gjb6','hepacam','id1','il6r','il6ra','irf1','jun','klf15','lhfp','lyve1','mertk','mgst1','mical2','myh2','ndrg2','npy1r','nr3c1','nudt9','osbpl3','pim3','plscr1','prr5','rasl11b','rdx','rhou','sall2','scamp2','sdc4','sesn1','slc25a33','sox2','sox4','sox9','spsb1','svil','tgfbr1','thra','tle4','tmem109','tob2','tsc22d3','vps37b','wipf3','wnt16','wnt7a','gpd1','ctgf','plekhf','dgkz','mtmr2','zfp36l1','azin1','cklf','ppp5c','sema6d','tle3'),
                        'ribosomal' = c('4831440e17rik', 'a530032d15rik', 'atp5l', 'b130024g19rik', 'btg3', 'fau', 'fbl', 'gapdh', 'gm14401', 'gm14409', 'gm6180', 'h2afz', 'hnrnpa3', 'il1f5', 'lsm7', 'mrpl23', 'olfr613', 'pin4', 'ppia', 'prl8a1', 'rpl10a', 'rpl11', 'rpl17', 'rpl23a', 'rpl27', 'rpl27a', 'rpl28', 'rpl31', 'rpl37', 'rpl38', 'rpl5', 'rpl9', 'rps12', 'rps13', 'rps16', 'rps23', 'rps27', 'rps27a', 'sec61g', 'sp140', 'vwa7', 'zkscan4'),
                        'early_genes_2' = c('btg2', 'ier2', 'atf3', 'apold1', 'ccn1', 'gimap6', 'tagap', 'tmem252', 'lrrc32'),
                        'early_genes_2_no_artifacts' = c('btg2', 'ier2', 'atf3', 'apold1', 'ccn1', 'gimap6', 'tmem252'))


####### SET NAMES OF CLUSTERS OF INTEREST ####### 





####### PREPARING SUBSET TABLES FOR SUBSTRUCTURES ####### 
search_for <- list(
  'hippocampus' = stringr::str_detect(
    string = tolower(descriptions_1_and_2$Brain_part_clean),
    pattern = '.*hipp*'
  ),
  'nucleus_accumbens' = stringr::str_detect(
    string = tolower(descriptions_1_and_2$Brain_part_clean),
    pattern = '.*cumbens*'
  ),
  'amygdala' = stringr::str_detect(
    string = tolower(descriptions_1_and_2$Brain_part_clean),
    pattern = '.*amyg*'
  ),
  'prefrontal_cortex' = stringr::str_detect(
    string = tolower(descriptions_1_and_2$Brain_part_clean),
    pattern = '.*fro*')
)

experiments_to_subset <- lapply(X = search_for, function(x){
  descriptions_1_and_2 %>%
    dplyr::select(Group_ID, Brain_part_clean) %>%
    dplyr::filter(x)
})

#This prints files into working folder
purrr::walk2(
  .x = experiments_to_subset,
  .y = names(search_for),
  .f = function(x, y) {
    create_subset_of_exps_with_at_least_3_papers_wrapper_wrapper(
      final_good_dataset___ = final_good_dataset_1_and_2,
      experiments_to_include_df_with_group_id_col = x,
      save_as_chr_ = y,
      group_id_col_str_ = 'Group_ID'
    )
  }
)
####### PREPARING SUBSET TABLES FOR SUBSTRUCTURES ####### 
####### PREPARING SUBSET TABLES FOR SPECIES ####### 
search_for <- list('mice' = 'mice', 'rats' = 'rats')

experiments_to_subset <- lapply(X = search_for, function(x){
  descriptions_1_and_2 %>%
    dplyr::select(Group_ID, Species) %>%
    dplyr::filter(descriptions_1_and_2$Species == x)
})

purrr::walk2(
  .x = experiments_to_subset,
  .y = names(search_for),
  .f = function(x, y) {
    create_subset_of_exps_with_at_least_3_papers_wrapper_wrapper(
      final_good_dataset___ = final_good_dataset_1_and_2,
      experiments_to_include_df_with_group_id_col = x,
      save_as_chr_ = y,
      group_id_col_str_ = 'Group_ID'
    )
  }
)
####### PREPARING SUBSET TABLES FOR SPECIES ####### 
####### PREPARING SUBSET TABLES FOR STRESS SENSITIVITY ####### 
search_for <- list('vulnerable' = 'vulnerable', 'resilient' = 'resilient')

experiments_to_subset <- lapply(X = search_for, function(x){
  descriptions_1_and_2 %>%
    dplyr::select(Group_ID, Stress_sensitivity_clean) %>%
    dplyr::filter(descriptions_1_and_2$Stress_sensitivity_clean == x)
})

purrr::walk2(
  .x = experiments_to_subset,
  .y = names(search_for),
  .f = function(x, y) {
    create_subset_of_exps_with_at_least_3_papers_wrapper_wrapper(
      final_good_dataset___ = final_good_dataset_1_and_2,
      experiments_to_include_df_with_group_id_col = x,
      save_as_chr_ = y,
      group_id_col_str_ = 'Group_ID'
    )
  }
)
####### PREPARING SUBSET TABLES FOR STRESS SENSITIVITY ####### 
####### PREPARING SUBSET TABLES FOR STRESS ####### 
search_for <- list('chronic_unpredictable_stress' = 'chronic unpredictable stress', 'fear_conditioning' = 'fear conditioning', 'forced_swimming' = 'forced swimming', 'immobilization_stress' = 'immobilization stress', 'social_stress' = 'social stress')

experiments_to_subset <- lapply(X = search_for, function(x){
  descriptions_1_and_2 %>%
    dplyr::select(Group_ID, Stress_clean) %>%
    dplyr::filter(descriptions_1_and_2$Stress_clean == x)
})

purrr::walk2(
  .x = experiments_to_subset,
  .y = names(search_for),
  .f = function(x, y) {
    create_subset_of_exps_with_at_least_3_papers_wrapper_wrapper(
      final_good_dataset___ = final_good_dataset_1_and_2,
      experiments_to_include_df_with_group_id_col = x,
      save_as_chr_ = y,
      group_id_col_str_ = 'Group_ID'
    )
  }
)
####### PREPARING SUBSET TABLES FOR STRESS ####### 
####### PREPARING SUBSET TABLES FOR ACUTE-CHRONIC ####### 
### Where to put cutoffs
table(descriptions_1_and_2$Repetitions_clean_days)

search_for <- list('acute' = 'acute', 'medium' = 'medium', 'prolonged' = 'prolonged')

experiments_to_subset <- lapply(X = search_for, function(x){
  descriptions_1_and_2 %>%
    dplyr::select(Group_ID, Stress_duration) %>%
    dplyr::filter(descriptions_1_and_2$Stress_duration == x)
})

purrr::walk2(
  .x = experiments_to_subset,
  .y = names(search_for),
  .f = function(x, y) {
    create_subset_of_exps_with_at_least_3_papers_wrapper_wrapper(
      final_good_dataset___ = final_good_dataset_1_and_2,
      experiments_to_include_df_with_group_id_col = x,
      save_as_chr_ = y,
      group_id_col_str_ = 'Group_ID'
    )
  }
)

####### PREPARING SUBSET TABLES FOR ACUTE-CHRONIC ####### 





### ARE ANY COLMUN IN SUBTABLES EMPTY? ###
temp_data <- list(hippocampus, nucleus_accumbens, amygdala, prefrontal_cortex, mice, rats, vulnerable, resilient, chronic_unpredictable_stress, fear_conditioning, forced_swimming, immobilization_stress, social_stress, acute, medium, prolonged)

temp_data_2 <- lapply(temp_data, function(x) { t(purrr::map_df(.x = x[,-1], .f = sum)) } )
empty_tables <- lapply(temp_data_2, function(x) {subset(x, x == 0)})
rm(temp_data, temp_data_2)
# nucleus_accumbens - 1 kolumna, rats - 2 kolumny
### ARE ANY COLMUN IN SUBTABLES EMPTY? ###





### GET NUMBERS OF GENES IN EXPS AND PAPERS ###
load('spread_medianed_final_good_dataset')

spread_medianed_final_good_dataset[is.na(spread_medianed_final_good_dataset)] <- 0

get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(spread_medianed_final_good_dataset, '')

### GET NUMBERS OF GENES IN EXPS AND PAPERS IN SUBGROUPS ###
objects_ <- list('hippocampus' = hippocampus, 'nucleus_accumbens' = nucleus_accumbens, 'amygdala' = amygdala, 'prefrontal_cortex' = prefrontal_cortex, 'mice' = mice, 'rats' = rats, 'chronic_unpredictable_stress' = chronic_unpredictable_stress, 'fear_conditioning' = fear_conditioning, 'forced_swimming' = forced_swimming, 'immobilization_stress' = immobilization_stress, 'social_stress' = social_stress, 'vulnerable' = vulnerable, 'resilient' = resilient, 'acute' = acute, 'medium' = medium, 'prolonged' = prolonged)

purrr::walk2(
  .x = objects_,
  .y = names(objects_),
  .f = function(x, y) {
    get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(x, y)
  }
)
### GET NUMBERS OF GENES IN EXPS AND PAPERS IN SUBGROUPS ###
### GET NUMBERS OF GENES IN ALL PAPERS ###
paper_number_numbers <- paper_number_ %>%
  dplyr::group_by(number) %>%
  dplyr::mutate(nb_of_genes_detected_in_this_nb_of_papers = dplyr::n()) %>%
  dplyr::mutate(percent_of_genes_detected_in_this_nb_of_papers = (dplyr::n()/length(paper_number_[[1]])) ) %>%
  dplyr::select(-lower_final_gene_name) %>%
  unique()
readr::write_tsv(paper_number_numbers, path = 'exp_and_paper_numbers/paper_number_numbers.tsv')
### GET NUMBERS OF GENES IN ALL PAPERS ###
### GET NUMBERS OF GENES IN ALL PAPERS DRAWING ###
# library(ggplot2)
paper_number_numbers <- paper_number_numbers[order(paper_number_numbers$number),]

paper_number_numbers_for_plot <- paper_number_numbers
paper_number_numbers_for_plot$nb_of_genes_detected_in_this_nb_of_papers[19] <-sum(paper_number_numbers$nb_of_genes_detected_in_this_nb_of_papers[19:25])
paper_number_numbers_for_plot$percent_of_genes_detected_in_this_nb_of_papers[19] <-sum(paper_number_numbers$percent_of_genes_detected_in_this_nb_of_papers[19:25])
paper_number_numbers_for_plot <- paper_number_numbers_for_plot[-c(20:25),]

paper_number_numbers_for_plot %>%
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
  scale_x_continuous(breaks = seq(1, 19, 1), labels = c(as.character(seq(1, 18, 1)), '18 <'))+
  labs(title = 'number of papers in which gene was detected')+
  xlab('number of papers')+
  ylab('number of genes')
### GET NUMBERS OF GENES IN ALL PAPERS DRAWING ###
### GET SUBSET OF EXP NUMBERS ###
temp <- readr::read_tsv(file = 'exp_and_paper_numbers/exp_number_and_percentage_hippocampus.tsv')
temp_cluster_of_interest <- groups_of_genes$ribosomal

temp_2 <- subset(x = temp, subset = temp$lower_final_gene_name %in% temp_cluster_of_interest)
### GET SUBSET OF EXP NUMBERS ###

  
### COMPARE NUMBERS OF GENES IN ALL PAPERS AND EXPERIMENTS TO GC DATA ###
gc_data <- readr::read_tsv(file = 'GC_publication/gc.tsv')
gc_data$in_gc_article <- T
gc_data$in_gc_88_best_genes <- ifelse(test = gc_data$Gene %in% groups_of_genes$`88_best`, yes = T, no = F)

temp_biomart <-biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
temp_bm_results <- biomaRt::getBM(
  attributes = "external_gene_name",
  filters = "external_gene_name",
  values = gc_data$Gene,
  uniqueRows = T,
  mart = temp_biomart
)
temp_bm_results$Gene <- tolower(temp_bm_results$external_gene_name)
temp_bm_and_gc <- merge(gc_data, temp_bm_results, by = 'Gene', all.x = T)

temp_bm_and_gc$Gene[temp_bm_and_gc$Gene == 'ctgf'] <- 'cnn2'
temp_bm_and_gc$Gene[temp_bm_and_gc$Gene == 'fam188a'] <- 'mindy3'
temp_bm_and_gc$Gene[temp_bm_and_gc$Gene == 'gyg1'] <- 'gyg'
temp_bm_and_gc$Gene[temp_bm_and_gc$Gene == 'mb21d1'] <- 'cgas'
temp_bm_and_gc$external_gene_name <- NULL

gc_and_experiements <- merge(x = temp_bm_and_gc, y = exp_number_and_percentage_, by.x = 'Gene', by.y = 'lower_final_gene_name', all.x = T)
gc_and_papers <- merge(x = temp_bm_and_gc, y = paper_number_, by.x = 'Gene', by.y = 'lower_final_gene_name', all.x = T)

readr::write_tsv(gc_and_experiements, 'GC_publication/gc_and_experiements.tsv')
readr::write_tsv(gc_and_papers, 'GC_publication/gc_and_papers.tsv')
### COMPARE NUMBERS OF GENES IN ALL PAPERS AND EXPERIMENTS TO GC DATA ###



### GET NEW NAMES FOR WHOLE GC DATA ###
opts_whole_gc <- list()
opts_whole_gc$input <- readr::read_tsv(file = 'whole_gc_data/gc_genes_all_input.txt')

opts_whole_gc <- getting_all_gc_from_input_with_cols_nb_and_gene(opts_whole_gc)

opts_whole_gc$output$biomart_detected <- tolower(opts_whole_gc$output$biomart_detected)

readr::write_tsv(opts_whole_gc$output, path = 'whole_gc_data/gc_genes_all_output.tsv')
### GET NEW NAMES FOR WHOLE GC DATA ###



### COMPARE NUMBERS OF GENES IN ALL PAPERS AND EXPERIMENTS TO EXTENDED GC DATA ###

opts_extended_gc <- list()

opts_extended_gc$input <- readr::read_tsv(file = 'whole_gc_data/extended_gc_input.txt')
opts_extended_gc$input$lower_final_gene_name <- tolower(opts_extended_gc$input$lower_final_gene_name)

opts_extended_gc$pn <- readr::read_tsv('exp_and_paper_numbers/paper_number_.tsv')
opts_extended_gc$exp <- readr::read_tsv('exp_and_paper_numbers/exp_number_and_percentage_.tsv')

opts_extended_gc$merge <- Reduce(function(x,y){merge(x, y, by = 'lower_final_gene_name', all.x = T)}, opts_extended_gc)

readr::write_tsv(x = opts_extended_gc$merge, path = 'whole_gc_data/extended_gc_output.tsv')
### COMPARE NUMBERS OF GENES IN ALL PAPERS AND EXPERIMENTS TO EXTENDED GC DATA ###



### ANOVA  ###
gather_spread_medianed_final_good_dataset <- tidyr::gather(data = spread_medianed_final_good_dataset, key = "exp", value = "logFC", na.rm = T, -lower_final_gene_name)
save(gather_spread_medianed_final_good_dataset, file = 'gather_spread_medianed_final_good_dataset')
# load('gather_spread_medianed_final_good_dataset')

descriptions_1_and_2_for_correlations <- descriptions_1_and_2 %>%
  dplyr::select(Group_ID, Species, Gender_clean, Repetitions_clean_days, Duration_clean_minutes, Brain_part_clean, Stress_sensitivity_clean, Stress_clean, Stress_duration, Measurement_latency_clean)

### In spread_medianed_final_good_datasetnie ma zduplikowanych eksperymentów. W correlations są. precorrelations rozwala liczbę eksperymentów!!!
  
pre_correlations <- merge(x = gather_spread_medianed_final_good_dataset, y = descriptions_1_and_2_for_correlations, by.x = 'exp', by.y = 'Group_ID', all.x = T)
pre_correlations$amg_vs_rest <- as.factor(ifelse(pre_correlations$Brain_part_clean == 'amygdala', 'amygdala', 'rest'))
pre_correlations$hp_vs_rest <- as.factor(ifelse(pre_correlations$Brain_part_clean == 'hippocampus', 'hippocampus', 'rest'))
pre_correlations$nac_vs_rest <- as.factor(ifelse(pre_correlations$Brain_part_clean == 'nucleus accumbens', 'nucleus accumbens', 'rest'))
pre_correlations$fcrx_vs_rest <- as.factor(ifelse(pre_correlations$Brain_part_clean %in% c('prefrontal cortex', 'frontal cortex'), 'frontal cortex', 'rest'))
pre_correlations$cus_vs_rest <- as.factor(ifelse(pre_correlations$Stress_clean == 'chronic unpredictable stress', 'chronic unpredictable stress', 'rest'))
pre_correlations$fc_vs_rest <- as.factor(ifelse(pre_correlations$Stress_clean == 'fear conditioning', 'fear conditioning', 'rest'))
pre_correlations$fs_vs_rest <- as.factor(ifelse(pre_correlations$Stress_clean == 'forced swimming', 'forced swimming', 'rest'))
pre_correlations$is_vs_rest <- as.factor(ifelse(pre_correlations$Stress_clean == 'immobilization stress', 'immobilization stress', 'rest'))
pre_correlations$ss_vs_rest <- as.factor(ifelse(pre_correlations$Stress_clean == 'social stress', 'social stress', 'rest'))

correlations_ <- pre_correlations %>%
  dplyr:: group_by(lower_final_gene_name) %>%
  tidyr::nest()
save(correlations_, file = 'correlations_')
# load('correlations_')

analyses_names <- c('Gender_clean', 'Species', 'Stress_sensitivity_clean', 'amg_vs_rest' , 'hp_vs_rest', 'nac_vs_rest', 'fcrx_vs_rest', 'cus_vs_rest', 'fc_vs_rest', 'fs_vs_rest' , 'is_vs_rest', 'ss_vs_rest', 'Stress_duration', 'Measurement_latency_clean')
# analyses_names <- c('amg_vs_rest' , 'hp_vs_rest', 'nac_vs_rest', 'fcrx_vs_rest')

# anova_test_name = 'anova'
# post_hoc_test_name = 'tukeyhsd'
# anova_test_name = 'welch'
# post_hoc_test_name = 'tukeyhsd'
anova_test_name = 'k-w'
post_hoc_test_name = 'conover'

analyses <-
  purrr::map(
    .x = analyses_names,
    .f = function(x) {
      anova_on_nested_df(
        df_ = correlations_,
        value_col_name = 'logFC',
        trait_col_name = x,
        main_test = anova_test_name,
        post_hoc = post_hoc_test_name
      )
    }
  ) # http://www.biostathandbook.com - they suggest, that one-way anova is very resistant to non-normal distributions (http://www.biostathandbook.com/kruskalwallis.html)

# broom::tidy(oneway.test())
# onewaytests::aov.test(logFC ~ ss_vs_rest, data = correlations_$data[[1]])
# test <- 
  
  dir.create(anova_test_name)
  purrr::walk2(.x = analyses, .y = analyses_names, .f = function(x,y){  jsonlite::write_json(x = x, path = paste0(anova_test_name, '/', anova_test_name, '_', y, '.json')) })
  
  purrr::walk2(
    .x = analyses,
    .y = analyses_names,
    .f = function(x, y) {
      x <- dplyr::select(x, -data)
      readr::write_tsv(x = x, path = paste0(anova_test_name, '/', anova_test_name, '_', y, '.tsv'))
    }
  )

  ### PRINT GENERAL FIGURES FOR ANOVA ANALYSIS. ALSO CREATES OPTS_GGPLOT ###
  
  for(n in seq(length(analyses)))
  {
    opts_ggplot <- list('analysis_number' = n)
    opts_ggplot <- rlist::list.append(opts_ggplot, 'full_dataset' = analyses[[opts_ggplot$analysis_number]], 'groupname' = analyses_names[[opts_ggplot$analysis_number]], 'y' = 'logFC')
    opts_ggplot <- rlist::list.append(opts_ggplot, 'nested_dataset' = opts_ggplot$full_dataset$data, 'names_for_nested_dataset' = opts_ggplot$full_dataset$lower_final_gene_name)
    
    opts_ggplot$dir_name <-paste0(anova_test_name, '/stat_analyses_', opts_ggplot$groupname)
    
    library(ggplot2)
    dir.create(opts_ggplot$dir_name)
    purrr::walk2(
      .x = opts_ggplot$nested_dataset,
      .y = opts_ggplot$names_for_nested_dataset,
      .f = function(x, y) {
        x %>%
          ggplot(aes(x = eval(parse(text = opts_ggplot$groupname)), y = eval(parse(text = opts_ggplot$y)))) +
          geom_jitter() +
          labs(title = y) +
          xlab(opts_ggplot$groupname) +
          ylab(opts_ggplot$y)


        ggsave(paste0(opts_ggplot$dir_name, '/', y, '.jpg'))
      }
    )
  }
### PRINT GENERAL FIGURES FOR ANOVA ANALYSIS ###
opts_individual <- list('logFC_col' = 'log')
# rm(opts_ggplot)
### PRINT INDIVIDUAL GENE FIGURES FOR ANOVA ANALYSIS ###
# opts_individual$individual_genes_for_filtering_only <- c('apold1', 'egr1', 'egr2', 'fosb', 'fos', 'dusp1', 'arc', 'epyc', 'gng8', 'ankk1', 'fezf1', 'sh2d1a', 'slc22a3', 'dio2', 'klf2')
# opts_individual$individual_genes_groupname <- 'Brain_part_clean'

# opts_individual$individual_genes_for_filtering_only <- c('slc8b1', 'ss18', 'ch25h', 'lcn2', 'lrg1', 'egr1', 'egr2', 'fos', 'fosb', 'fosl2', 'htra1', 'ilf2', 'kcnq2', 'proz', 'slc9a3r1', 'socs5', 'vmn2r1')
# opts_individual$individual_genes_groupname <- 'Stress_clean'


# opts_individual$individual_genes_for_filtering_only <- c('apold1', 'arc', 'arntl2', 'blacf1', 'btg2', 'ccn1', 'crkl', 'dusp1', 'egr1', 'egr2', 'fos', 'fosb', 'gcnt', 'gtf2a1', 'ier2', 'klf2', 'npas4', 'nr4a1', 'polq', 'sgk1')
# opts_individual$individual_genes_groupname <- 'Repetitions_clean_days'
# opts_individual$individual_genes_groupname <- 'Duration_clean_minutes'
# opts_individual$individual_genes_groupname <- 'Measurement_latency_clean'

# opts_individual$individual_genes_for_filtering_only <- c('fos', 'sgk1', 'ccl5', 'hba-a1', 'aqp1', 'fkbp5', 'fam180a', 'fmc1', 'ptgds', 'grm1')
# opts_individual$individual_genes_groupname <- 'exp'

# opts_individual$individual_genes_for_filtering_only <- groups_of_genes$hemoglobin_cluster
# opts_individual$individual_genes_groupname <- 'exp'

# opts_individual$individual_genes_for_filtering_only <- groups_of_genes$ribosomal
# opts_individual$individual_genes_groupname <- 'Stress_sensitivity_clean'

opts_individual$individual_genes_for_filtering_only <- groups_of_genes$early_genes_2_no_artifacts
opts_individual$individual_genes_groupname <- 'Brain_part_clean'
opts_individual$analysis_name <- 'early_genes_2'
# opts_individual$individual_genes_for_filtering_only <- groups_of_genes$early_genes_2_no_artifacts
# opts_individual$individual_genes_groupname <- 'Brain_part_clean'

#To remove
# temp <- readr::read_tsv(file = 'exp_and_paper_numbers/exp_number_and_percentage_.tsv')
# temp_2 <- subset(x = temp, subset = temp$lower_final_gene_name %in% opts_individual$individual_genes_for_filtering_only)
#To remove

opts_individual$individual_genes_table <- subset(correlations_, correlations_$lower_final_gene_name %in% opts_individual$individual_genes_for_filtering_only)
opts_individual$individual_genes_for_filtering_only <- NULL

 
dir.create(opts_individual$analysis_name)
dir.create(paste0(opts_individual$analysis_name, '/individual_genes'))

# library(ggplot2)
purrr::walk2(
  .x = opts_individual$individual_genes_table$data,
  .y = opts_individual$individual_genes_table$lower_final_gene_name,
  .f = function(x, y) {
      x %>%
        ggplot(aes(x = as.factor(eval(parse(text = opts_individual$individual_genes_groupname))), y = eval(parse(text = opts_individual$logFC_col)))) +
        geom_jitter(size = 0.3, width = 0.25) +
        labs(title = y) +
        xlab(opts_individual$individual_genes_groupname) +
        ylab(opts_individual$logFC_col) +
        theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1))

    dir.create(paste0(opts_individual$analysis_name, '/individual_genes/', opts_individual$individual_genes_groupname))
    
    if (opts_individual$individual_genes_groupname == 'Stress_clean') {
      ggsave(paste0(opts_individual$analysis_name, '/individual_genes/', opts_individual$individual_genes_groupname, '/', y, '_', opts_individual$individual_genes_groupname, '.jpg'), width = 3, height = 4)
    } else {
      ggsave(paste0(opts_individual$analysis_name, '/individual_genes/', opts_individual$individual_genes_groupname, '/', y, '_', opts_individual$individual_genes_groupname, '.jpg'))
    }
  }
)


### PRINT INDIVIDUAL GENE FIGURES FOR ANOVA ANALYSIS ###
### ADD INFO ON ENRICHMENT OF CLUSTER IN GIVEN EXPERIMENT CHARACTERISTIC ###
# Needs subset of opts_individual$individual_genes_table .
opts_individual <- is_characteristic_enriched_wrapper(opts_ggplot_ = opts_individual, value_of_interest_ = 'hippocampus', table_of_interest_ = 'Brain_part_clean')

opts_individual$how_many_times_gene_detected_in_exp_per_paper <- purrr::map(
    .x = opts_individual$individual_genes_table$data,
    .f = function(x) {
      temp_ <- subset(x = x, subset = x$logFC != 0)
      temp_2 <- stringr::str_remove(temp_$exp, pattern = '_.*') 
      return(as.data.frame(table(temp_2)))
    }
  )

opts_individual <- get_number_of_exps_for_paper_per_gene_wrapper(opts_ggplot_ = opts_individual, papers_to_check_number_of_exprs_in = c('7', '44', '52'))

opts_individual$nb_of_exps_for_value <- get_exps_in_which_gene_was_present_for_given_value_of_interest_wrapper( opts_individual_ = opts_individual, logFC_col = 'logFC', exp_col = 'exp', descriptions_ = descriptions_1_and_2_for_correlations, exp_col_in_descr_ = 'Group_ID', list_for_subsetting_name_is_colname_value_is_value_of_interest = list('Brain_part_clean' = 'hippocampus', 'Stress_duration' = 'acute') )
### ADD INFO ON ENRICHMENT OF CLUSTER IN GIVEN EXPERIMENT CHARACTERISTIC ###
### ANOVA  ###

opts_individual$names_for_nested_dataset



### GET CLUSTERS NB_1 - HUSHED PART ON THE BOTTOM ###
opts_cluster <- list('folder' = 'clust_visualizations')
dir.create(opts_cluster$folder)

opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$additional_1)
opts_cluster$cluster_to_vis_name <- 'additional_1'
opts_cluster$exp_to_remove <- NA

opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$additional_2)
opts_cluster$cluster_to_vis_name <- 'additional_2'
opts_cluster$exp_to_remove <- NA



opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$additional_4)
opts_cluster$cluster_to_vis_name <- 'additional_4'
opts_cluster$exp_to_remove <- NA

opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$additional_5)
opts_cluster$cluster_to_vis_name <- 'additional_5'
opts_cluster$exp_to_remove <- c('6_1', '2_2', '2_1', '15_2', '15_5', '25_2', '15_1', '13_1', '13_2', '1_3', '12_4', '12_2', '44_3', '78_2', '3_1', '67_3', '1_1', '78_4', '14_1', '2_4', '45_5')

opts_cluster$cluster_to_vis_prepared <- prepare_for_clustering_wrapper(spread_df = opts_cluster$cluster_to_vis, remove_exp_with_this_many_hits = 1)
      
opts_cluster$cluster_to_vis_prepared <- opts_cluster$cluster_to_vis_prepared[, !(colnames(opts_cluster$cluster_to_vis_prepared) %in% opts_cluster$exp_to_remove)]

dev.off()
tiff(paste0(opts_cluster$folder, '/', opts_cluster$cluster_to_vis_name, '_only_exps_with_at_least_2_non0_values_removed_emptish_columns.tiff'), width = 1920,  height = 1080)
gplots::heatmap.2(x = opts_cluster$cluster_to_vis_prepared, trace="none", dendrogram = 'column', cexCol = 0.75, lwid=c(0.1,4), col = colorRamps::matlab.like, breaks = 200)
dev.off()

# Table for clustering
write.table(
  x = as.data.frame(opts_cluster$cluster_to_vis_prepared), 
  file = paste0(opts_cluster$folder, '/', opts_cluster$cluster_to_vis_name, '_only_exps_with_at_least_2_non0_values_removed_emptish_columns.tsv'),
  sep = '\t',
  row.names = T,
  col.names = NA,
  dec = '.'
)  
### GET CLUSTERS NB_1 - HUSHED PART ON THE BOTTOM ###




### GET EXPRESSION STABILITY FOR GIVEN CLUSTER PER EXPERIMENT ###
temp <- as.data.frame(t(opts_cluster$cluster_to_vis_prepared))
temp_names <- as.data.frame(rownames(temp), stringsAsFactors = F)
temp_2 <- cbind(temp_names, temp)
colnames(temp_2) <- stringr::str_replace(string = colnames(temp_2), pattern = 'rownames\\(temp\\)', replacement = 'lower_final_gene_name')

get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(temp_2, opts_cluster$cluster_to_vis_name)
### GET EXPRESSION STABILITY FOR GIVEN CLUSTER PER EXPERIMENT ###




### GET EXPRESSION STABILITY FOR GIVEN CLUSTER FROM ALL ANALYSES ###
opts_expr_stability <- list()
opts_expr_stability$file_names <- as.list(list.files(path = 'exp_and_paper_numbers', pattern = 'exp_number_and_percentage.*'))
opts_expr_stability$genes_to_subset <- groups_of_genes$early_genes_2_no_artifacts
opts_expr_stability$name_of_col_with_gene_names <- 'lower_final_gene_name'
opts_expr_stability$name_of_col_with_percentage <- 'perc_of_upregulated' 

opts_expr_stability$files <-
  lapply(opts_expr_stability$file_names, function(x) {
    temp <- readr::read_tsv(file = paste0('exp_and_paper_numbers/', x))
    temp <- subset(x = temp, subset = temp[[opts_expr_stability$name_of_col_with_gene_names]] %in% opts_expr_stability$genes_to_subset)  })

opts_expr_stability$files_mean_SD <-
  lapply(opts_expr_stability$files, function(x) {
    temp <- mean(x[[opts_expr_stability$name_of_col_with_percentage]])
    temp2 <- sd(x[[opts_expr_stability$name_of_col_with_percentage]])
    return(paste0('Mean: ', temp, ', SD: ', temp2)) })
### GET EXPRESSION STABILITY FOR GIVEN CLUSTER FROM ALL ANALYSES ###




### GET PROTEIN CODING GENES  ###
temp_biomart <-biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

temp_bm_results <- biomaRt::getBM(
  attributes = c('ensembl_peptide_id', "external_gene_name"),
  uniqueRows = T,
  mart = temp_biomart)

temp_bm_results_subset <- subset(temp_bm_results, temp_bm_results$ensembl_peptide_id != '')
temp_bm_results_subset <- data.frame('external_gene_name' = tolower(unique(temp_bm_results_subset$external_gene_name)), 'protein' = T)

our_genes <- data.frame('external_gene_name' = spread_medianed_final_good_dataset$lower_final_gene_name, 'in_our_data' = T)

mergement <- merge(x = our_genes, y = temp_bm_results_subset, by = 'external_gene_name')
length(mergement[[1]])/length(spread_medianed_final_good_dataset[[1]])

rm(our_genes, mergement)
### GET PROTEIN CODING GENES  ###



### GET HUMAN GENES  ###
opts_human <- list()
opts_human$genes <- subset(final_good_dataset_1_and_2, subset = final_good_dataset_1_and_2$Paper == 27) %>%
  dplyr::select(lower_final_gene_name)

opts_human$pn <- readr::read_tsv('exp_and_paper_numbers/paper_number_.tsv')
opts_human$exp <- readr::read_tsv('exp_and_paper_numbers/exp_number_and_percentage_.tsv')
opts_human$pn_pfc <- readr::read_tsv('exp_and_paper_numbers/paper_number_prefrontal_cortex.tsv')
opts_human$exp_pfc <- readr::read_tsv('exp_and_paper_numbers/exp_number_and_percentage_prefrontal_cortex.tsv')

opts_human$merge <- Reduce(function(x,y){merge(x, y, by = 'lower_final_gene_name', all.x = T)}, list(opts_human$genes, opts_human$pn, opts_human$exp, opts_human$pn_pfc, opts_human$exp_pfc))

colnames(opts_human$merge) <- stringr::str_replace(string = colnames(opts_human$merge), pattern = '.x', replacement = '_all_data')
colnames(opts_human$merge) <- stringr::str_replace(string = colnames(opts_human$merge), pattern = '.y', replacement = '_pfc')

dir.create('human_genes')
readr::write_tsv(x = opts_human$merge, path = 'human_genes/human_genes.tsv')

### GET HUMAN GENES  ###



### GET INTEGRATED DATA  ###
library(Hmisc)
load('final_good_dataset_1_and_2') #232 146
load('reformated_raw_dataset_1_and_2') #268 669
load('descriptions_1_and_2')

temp_final_good_dataset_1_and_2 <- final_good_dataset_1_and_2
temp_reformated_raw_dataset_1_and_2 <- reformated_raw_dataset_1_and_2
temp_descriptions_1_and_2 <- descriptions_1_and_2

temp_final_good_dataset_1_and_2$for_merge <- paste0(temp_final_good_dataset_1_and_2$Paper, temp_final_good_dataset_1_and_2$Experiment, temp_final_good_dataset_1_and_2$logFC, temp_final_good_dataset_1_and_2$adj_p, temp_final_good_dataset_1_and_2$Gene_symbol, temp_final_good_dataset_1_and_2$GenBank_Accession, temp_final_good_dataset_1_and_2$Gene_title, temp_final_good_dataset_1_and_2$`GEO ID`, temp_final_good_dataset_1_and_2$Ensembl_ID, temp_final_good_dataset_1_and_2$p)

temp_reformated_raw_dataset_1_and_2$for_merge <- paste0(temp_reformated_raw_dataset_1_and_2$Paper, temp_reformated_raw_dataset_1_and_2$Experiment, temp_reformated_raw_dataset_1_and_2$logFC, temp_reformated_raw_dataset_1_and_2$adj_p, temp_reformated_raw_dataset_1_and_2$Gene_symbol, temp_reformated_raw_dataset_1_and_2$GenBank_Accession, temp_reformated_raw_dataset_1_and_2$Gene_title, temp_reformated_raw_dataset_1_and_2$`GEO ID`, temp_reformated_raw_dataset_1_and_2$Ensembl_ID, temp_reformated_raw_dataset_1_and_2$p)

temp_super <- merge(x = temp_reformated_raw_dataset_1_and_2, y = temp_final_good_dataset_1_and_2, by = 'for_merge', all.x = T) #268 669 vs 268 671

temp_super_2 <- temp_super[duplicated(temp_super$for_merge),]

temp_superv2 <- temp_super[-c(63312, 63313),] %>% # 2 duplicated entries!!
  dplyr::select(dplyr::ends_with(".x"), lower_final_gene_name) %>%
  dplyr::select(-c(everything.x, Entry_number.x))

colnames(temp_superv2) <- stringr:: str_remove(string = colnames(temp_superv2), pattern = '\\.x')

temp_superv2$for_merge_2 <- paste0(temp_superv2$Paper, '_', temp_superv2$Experiment)

temp_descriptions_1_and_2$for_merge_2 <- paste0(temp_descriptions_1_and_2$Paper_ID, '_', temp_descriptions_1_and_2$Group_ID)

temp <- integrated_data

integrated_data <- merge(x = temp_superv2, y = temp_descriptions_1_and_2, by = 'for_merge_2')
integrated_data$for_merge_2 <- NULL
integrated_data$Measurment_latency <- stringr::str_replace_all(string = integrated_data$Measurment_latency, pattern = '\t', replacement = ' ')
integrated_data$Measurment_latency <- stringr::str_replace_all(string = integrated_data$Measurment_latency, pattern = '\n', replacement = ' ')
integrated_data$Repetitions <- stringr::str_replace_all(string = integrated_data$Repetitions, pattern = '\n', replacement = ' ')
integrated_data$reference_1 <- stringr::str_replace_all(string = integrated_data$reference_1, pattern = '\n', replacement = ' ')
integrated_data$Stress <- stringr::str_replace_all(string = integrated_data$Stress, pattern = '\n', replacement = ' ')
integrated_data$Brain_part <- stringr::str_replace_all(string = integrated_data$Brain_part, pattern = '\n', replacement = ' ')
integrated_data$Paper_link_2 <- stringr::str_replace_all(string = integrated_data$Paper_link_2, pattern = '\n', replacement = ' ')

save(integrated_data, file = 'integrated_data')

dir.create('integrated')

write.table(x = integrated_data, file = 'integrated/integrated_data_v4.tsv', sep = '\t', dec = ',', row.names = F, quote = F)

rm(temp_final_good_dataset_1_and_2, temp_reformated_raw_dataset_1_and_2, temp_descriptions_1_and_2)
### GET INTEGRATED DATA  ###



### TESTING DATA - GENES ###
load('final_good_dataset_1_and_2') #232146
# load('reformated_raw_dataset')
# reformated_raw_dataset_1 <- reformated_raw_dataset
# library(Hmisc)
# reformated_raw_dataset_1_no_10_49 <- subset(x = reformated_raw_dataset_1, subset = reformated_raw_dataset_1$Paper %nin% c(10, 49))
# load('reformated_raw_dataset_2')
# reformated_raw_dataset_2 <- reformated_raw_dataset
# reformated_raw_dataset_2$everything2 <- NULL
# reformated_raw_dataset_2$Probe_ID_old <- NULL
# reformated_raw_dataset_1_and_2 <- rbind(reformated_raw_dataset_1_no_10_49, reformated_raw_dataset_2)
# save(reformated_raw_dataset_1_and_2, file = 'reformated_raw_dataset_1_and_2')
load('reformated_raw_dataset_1_and_2') #268669

test <- subset(x = reformated_raw_dataset_1_and_2, abs(reformated_raw_dataset_1_and_2$logFC) < 0.05)

post_annotation <- subset(final_good_dataset_1_and_2, final_good_dataset_1_and_2$lower_final_gene_name == 'slc22a3', select = c(Experiment, logFC, everything))
pre_annotation <- subset(reformated_raw_dataset_1_and_2, stringr::str_detect(string = tolower(reformated_raw_dataset_1_and_2$everything), pattern = 'slc22a3'), select = c(Experiment, logFC, everything))
post_annotation_nb_of_exps <- data.frame(unique(post_annotation$Experiment))
### TESTING DATA - GENES ###
### TESTING DATA - EXPERIMENTS ###
opts_general <- list('column' = 'Stress_duration', 'value' = 'acute')

test_desc <- subset(descriptions_1_and_2, descriptions_1_and_2[[opts_general$column]] == opts_general$value, select = c(Group_ID, eval(parse(text = opts_general$column))))
test_desc$xx <- unique(test_desc$Group_ID)

test_desc2 <- data.frame('Group_ID' = colnames(eval(parse(text = opts_general$value))))

test_difference <- dplyr::anti_join(x = test_desc, y = test_desc2, by = 'Group_ID')
### TESTING DATA - EXPERIMENTS ###
### TESTING DATA - WHICH EXPREIMENTS WERE REMOVED ###
test_corr <- subset(x = correlations_, subset = correlations_$lower_final_gene_name == 'sox2')

test_gather <- subset(gather_spread_medianed_final_good_dataset, gather_spread_medianed_final_good_dataset$lower_final_gene_name == 'fezf1')

test_pre_corr <- subset(pre_correlations, pre_correlations$lower_final_gene_name == 'sox2')

test_gene <- subset(x = amygdala, subset = amygdala$lower_final_gene_name == 'slc22a3')

test_desc <- subset(x = descriptions_1_and_2_for_correlations, descriptions_1_and_2_for_correlations$Group_ID == '10_1')
### TESTING DATA - WHICH EXPREIMENTS WERE REMOVED ###
### TESTING DATA - DESCRIPTIONS ###



### TESTING DATA - DESCRIPTIONS ###



load('final_good_dataset_1_and_2')
load('reformated_raw_dataset') #261475
load('reformated_raw_dataset_2') #9173
load('reformated_raw_dataset_1_and_2') #268669
#reformated_raw_dataset - c(10, 49) # 259496
#reformated_raw_dataset + reformated_raw_dataset_2 = 268669

test <- subset(x = reformated_raw_dataset, subset = reformated_raw_dataset$Paper %nin% c(10, 49))


library(Hmisc)
test <- subset(x = descriptions_1_and_2_for_correlations, subset = (descriptions_1_and_2_for_correlations$Stress_duration == 'acute') & (descriptions_1_and_2_for_correlations$Group_ID %nin% c('32_1', '32_2', '61_1', '61_2')))
test <- subset(x = descriptions_1_and_2_for_correlations, subset = (descriptions_1_and_2_for_correlations$Brain_part_clean == 'hippocampus') & (descriptions_1_and_2_for_correlations$Group_ID %nin% c('32_1', '32_2', '61_1', '61_2')))

test <- subset(x = descriptions_1_and_2_for_correlations, subset = (descriptions_1_and_2_for_correlations$Stress_duration == 'acute') & (descriptions_1_and_2_for_correlations$Group_ID %nin% c('32_1', '32_2', '61_1', '61_2')) & !(stringr::str_detect(string = descriptions_1_and_2_for_correlations$Group_ID, pattern = '^52.*') ))
test <- subset(x = descriptions_1_and_2_for_correlations, subset = (descriptions_1_and_2_for_correlations$Brain_part_clean == 'hippocampus') & (descriptions_1_and_2_for_correlations$Group_ID %nin% c('32_1', '32_2', '61_1', '61_2')) & !(stringr::str_detect(string = descriptions_1_and_2_for_correlations$Group_ID, pattern = '^52.*') ))


test <- subset(x = descriptions_1_and_2_for_correlations, subset = (descriptions_1_and_2_for_correlations$Stress_duration == 'acute') & (descriptions_1_and_2_for_correlations$Group_ID %nin% c('32_1', '32_2', '61_1', '61_2')) & descriptions_1_and_2_for_correlations$Brain_part_clean == 'hippocampus')

test2 <- 










test <- unique(final_good_dataset_1_and_2$lower_final_gene_name)




### GET CLUSTERS NB_2 - HUSHED PART SOMEWHERE UP ###
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$hemoglobin_cluster)
# opts_cluster$cluster_to_vis_name <- 'hemoglobin_cluster'
# opts_cluster$exp_to_remove <- c('45_24', '15_2', '46_1', '2_4', '78_2', '78_4', '16_1')
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$choroid_cluster)
# opts_cluster$cluster_to_vis_name <- 'choroid_cluster'
# opts_cluster$exp_to_remove <- c('7_21', '53_1', '5_1', '45_5', '44_4', '1_1', '1_3', '14_1', '12_2', '12_4', '15_1', '13_1', '78_2', '13_3', '25_2', '45_3')
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$early_genes_cluster)
# opts_cluster$cluster_to_vis_name <- 'early_genes_cluster'
# opts_cluster$exp_to_remove <- c('67_1', '25_4', '7_17', '19_1', '19_2', '19_3', '19_4', '14_1', '13_3', '2_4')

# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('alas2', 'hbb-bt', 'hba-a1', 'hbb-bs', 'hba-a2', 'ch25h', 'lcn2', 'lrg1', 's100a8', 's100a9', 'cd68', 'cd206'))
# opts_cluster$cluster_to_vis_name <- 'hemoglobin_cluster_and_markers'
# opts_cluster$exp_to_remove <- c('45_24', '15_2', '46_1', '2_4', '78_2', '78_4', '16_1')

# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('slamf1', 'itga2b', 'gata1'))
# opts_cluster$cluster_to_vis_name <- 'platelets'
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('aif1' , 'ccl22' , 'cd14' , 'fcgr3a' , 'fcgr2a' , 'cd33' , 'cd40' , 'ptprc' , 'fcgr1a' , 'cd68' , 'cd163' , 'ptgs1' , 'cx3cr1' , 'ftl', 'slc2a5' ,'hla-dra' , 'hla-drb1' , 'p2ry12' , 'spi1' , 'tlr2' , 'tmem119' , 'trem2' ))
# opts_cluster$cluster_to_vis_name <- 'microglia'
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('mme','csf3r','itgam','itgax','il3ra','anpep','cd14','fut4','fcgr3','pecam1','fcgr2b','cd33','itga4','sell','fcgr1a','ceacam8','c5ar1','cxcr1','cxcr2','fpr1','cd68','gr1','hla-dr','jaml','lcn2','prtn3','tlr2'))
# opts_cluster$cluster_to_vis_name <- 'neutrophile'
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('cd4','il2ra','cd45ra','cd8a','sell','cd27','il7r','foxp3','ccr7','ptprc','cd8a','cd8b','ccr6','itgam','tnfrsf8','ptprc','cd6','ctla4','il2ra','ablim1','actn1','bcl2','c1orf162','snhg32','ccr10','ccr4','cd101','itgb2','cd28','pecam1','entpd1','cd3d','itga4','cd55','ncam1','b3gat1','fas','cxcr3','eif3l','fam117b','fcrl3','lrrc32','gpr183','ifng','il4','il7r','ldlrap1','snhg29','lrrn3','mal','myc','nell2','nosip','prf1','serinc5','tcf7','tmem204','trabd2a','txk'))
# opts_cluster$cluster_to_vis_name <- 't_cell_homo'
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('tcf7','trbc2','cd4','cd8a','cd3e','cd8a','cd8b1','izumo1r','cd3e','ccl5','cd4','igfbp4','nrp1'))
# opts_cluster$cluster_to_vis_name <- 't_cell_mus'
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('cd68', 'mrc1'))
# opts_cluster$cluster_to_vis_name <- 'macrophages'
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('arg2','bhlhe40','ccl22','cd14','cd163','siglec1','cd19','cd1a','cd1c','cd200r1','cd68','cd74','cd80','cd86','cd93','cib1','ciita','crem','cx3cr1','cxcl10','cybb','cyth1','dock2','dok3','dse','emb','cyria','fgr','flt3','fosl2','fpr3','ftl','fxyd5','gpr132','gpr65','hla-dmb','hla-dqa1','hla-dra','hla-drb1','hla-drb5','ido1','ifitm2','il10','il1rn','iqgap1','itga4','itgam','kit','kynu','lyve1','lyz','metrnl','mrc1','ms4a6a','ms4a7','mxd1','myb','nfil3','nos2','pde4b','pim1','plac8','plbd1','pltp','slc66a3','ptpn7','runx1','s100a11','samhd1','sh3bgrl','socs3','spint2','syngr2','tgfbi','tgm2','thbd','timd4','tlr1','tlr2','tlr4','tlr8','tmem119','tmem123','tnfsf13','trem1','vopp1'))
# opts_cluster$cluster_to_vis_name <- 'macrophages_many'
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('gfap','aldh1l1','atp13a4','cbs','sox9','slc1a3','slc1a2'))
# opts_cluster$cluster_to_vis_name <- 'astrocytes'

# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% groups_of_genes$ribosomal)
# opts_cluster$cluster_to_vis_name <- 'ribosomal'
# opts_cluster$exp_to_remove <- NA
### GET CLUSTERS NB_2 - HUSHED PART SOMEWHERE UP ###