source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# rm(list = ls(pattern = 'temp.*|test.*'))
load('final_good_dataset_1_and_2') 
load('reformated_raw_dataset_1_and_2')
load('descriptions_1_and_2')



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



### COMPARE NUMBERS OF GENES IN ALL PAPERS AND EXPERIMENTS TO GC DATA ###
`88_best` <- c('ddit4','errfi1','klf9','bcl6','fkbp5','mt2','mt2a','nfkbia','pdk4','adcy9','cxxc5','dusp1','eva1a','litaf','nedd9','rhob','sgk1','sult1a1','tiparp','aldoc','arhgef3','arl4d','bcl6b','cables1','calm2','ccnd1','cdkn1a','cdo1','chst1','cyp7b1','ehd3','fzd1','gab1','gap43','gjb6','hepacam','id1','il6r','il6ra','irf1','jun','klf15','lhfp','lyve1','mertk','mgst1','mical2','myh2','ndrg2','npy1r','nr3c1','nudt9','osbpl3','pim3','plscr1','prr5','rasl11b','rdx','rhou','sall2','scamp2','sdc4','sesn1','slc25a33','sox2','sox4','sox9','spsb1','svil','tgfbr1','thra','tle4','tmem109','tob2','tsc22d3','vps37b','wipf3','wnt16','wnt7a','gpd1','ctgf','plekhf','dgkz','mtmr2','zfp36l1','azin1','cklf','ppp5c','sema6d','tle3')
gc_data <- readr::read_tsv(file = 'GC_publication/gc.tsv')
gc_data$in_gc_article <- T
gc_data$in_gc_88_best_genes <- ifelse(test = gc_data$Gene %in% `88_best`, yes = T, no = F)

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
  
  
  
    
  ### PRINT GENERAL FIGURES FOR ANOVA ANALYSIS ###
  
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


### PRINT INDIVIDUAL GENE FIGURES FOR ANOVA ANALYSIS ###
opts_ggplot$individual_genes_for_filtering_only <- c('apold1', 'egr1', 'egr2', 'fosb', 'fos', 'dusp1', 'arc', 'epyc', 'gng8', 'ankk1', 'fezf1', 'sh2d1a', 'slc22a3', 'dio2', 'klf2')
opts_ggplot$individual_genes_groupname <- 'Brain_part_clean'

# opts_ggplot$individual_genes_for_filtering_only <- c('slc8b1', 'ss18', 'ch25h', 'lcn2', 'lrg1', 'egr1', 'egr2', 'fos', 'fosb', 'fosl2', 'htra1', 'ilf2', 'kcnq2', 'proz', 'slc9a3r1', 'socs5', 'vmn2r1')
# opts_ggplot$individual_genes_groupname <- 'Stress_clean'


# opts_ggplot$individual_genes_for_filtering_only <- c('apold1', 'arc', 'arntl2', 'blacf1', 'btg2', 'ccn1', 'crkl', 'dusp1', 'egr1', 'egr2', 'fos', 'fosb', 'gcnt', 'gtf2a1', 'ier2', 'klf2', 'npas4', 'nr4a1', 'polq', 'sgk1')
# opts_ggplot$individual_genes_groupname <- 'Repetitions_clean_days'
# opts_ggplot$individual_genes_groupname <- 'Duration_clean_minutes'
# opts_ggplot$individual_genes_groupname <- 'Measurement_latency_clean'

# opts_ggplot$individual_genes_for_filtering_only <- c('fos', 'sgk1', 'ccl5', 'hba-a1', 'aqp1', 'fkbp5', 'fam180a', 'fmc1', 'ptgds', 'grm1')
# opts_ggplot$individual_genes_groupname <- 'exp'

opts_ggplot$individual_genes_table <- subset(correlations_, correlations_$lower_final_gene_name %in% opts_ggplot$individual_genes_for_filtering_only)

dir.create(paste0(anova_test_name, '/individual_genes'))

purrr::walk2(
  .x = opts_ggplot$individual_genes_table$data,
  .y = opts_ggplot$individual_genes_table$lower_final_gene_name,
  .f = function(x, y) {
      x %>%
        ggplot(aes(x = as.factor(eval(parse(text = opts_ggplot$individual_genes_groupname))), y = eval(parse(text = opts_ggplot$y)))) +
        geom_jitter(size = 0.3, width = 0.25) +
        labs(title = y) +
        xlab(opts_ggplot$individual_genes_groupname) +
        ylab(opts_ggplot$y) +
        theme(axis.text.x = element_text(angle = 45, size = 5, hjust = 1))

    dir.create(paste0(anova_test_name, '/individual_genes/', opts_ggplot$individual_genes_groupname))
    
    if (opts_ggplot$individual_genes_groupname == 'Stress_clean') {
      ggsave(paste0(anova_test_name, '/individual_genes/', opts_ggplot$individual_genes_groupname, '/', y, '_', opts_ggplot$individual_genes_groupname, '.jpg'), width = 3, height = 4)
    } else {
      ggsave(paste0(anova_test_name, '/individual_genes/', opts_ggplot$individual_genes_groupname, '/', y, '_', opts_ggplot$individual_genes_groupname, '.jpg'))
    }
  }
)
### PRINT INDIVIDUAL GENE FIGURES FOR ANOVA ANALYSIS ###
### ANOVA  ###



### GET CLUSTERS  ###
opts_cluster <- list('folder' = 'clust_visualizations')
dir.create(opts_cluster$folder)

# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('alas2', 'hbb-bt', 'hba-a1', 'hbb-bs', 'hba-a2', 'ch25h', 'lcn2', 'lrg1', 's100a8', 's100a9'))
# opts_cluster$cluster_to_vis_name <- 'hemoglobin_cluster'
# opts_cluster$exp_to_remove <- c('45_24', '15_2', '46_1', '2_4', '78_2', '78_4', '16_1')
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('cldn1',	'epn3',	'msx1',	'col8a2',	'lbp',	'ace',	'clic6',	'mfrp',	'krt8',	'drc7',	'ecrg4',	'prr32',	'aqp1',	'col8a1',	'steap1',	'f5',	'enpp2',	'trpv4',	'otx2',	'folr1',	'cldn2',	'kcne2',	'tmem72',	'slc4a5',	'kl',	'sostdc1',	'ttr',	'col9a3',	'slc39a4',	'sema3b',	'prlr',	'cox8b',	'oca2',	'slc2a12',	'igf2',	'igfbp2',	'slc13a4',	'pcolce',	'wdr86'))
# opts_cluster$cluster_to_vis_name <- 'choroid_cluster'
# opts_cluster$exp_to_remove <- c('7_21', '53_1', '5_1', '45_5', '44_4', '1_1', '1_3', '14_1', '12_2', '12_4', '15_1', '13_1', '78_2', '13_3', '25_2', '45_3')
# 
# opts_cluster$cluster_to_vis <- subset(x = spread_medianed_final_good_dataset, spread_medianed_final_good_dataset$lower_final_gene_name %in% c('cdkn1a',	'dusp1',	'egr1',	'egr2',	'fosb',	'fos',	'arc',	'irs2',	'junb',	'midn',	'b4galt1',	'fam107a',	'coq10b',	'gch1',	'csrnp1',	'per1',	'sgk1',	'ddit4',	'tsc22d3',	'pnpla2'))
# opts_cluster$cluster_to_vis_name <- 'early_genes_cluster'
# opts_cluster$exp_to_remove <- c('67_1', '25_4', '7_17', '19_1', '19_2', '19_3', '19_4', '14_1', '13_3', '2_4')

unique(our_genes$external_gene_name)

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



### GET CLUSTERS  ###



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






















# correlations_gender <- anova_on_nested_df(df_ = correlations_, value_col_name = 'logFC', trait_col_name = 'Gender_clean', main_test = 'anova', post_hoc = 'tukeyhsd') #4021 vs 2000 vs 1991
# 
# correlations_species <- anova_on_nested_df(df_ = correlations_, value_col_name = 'logFC', trait_col_name = 'Species', main_test = 'anova', post_hoc = 'tukeyhsd') #568 vs 240 vs 230
# 
# correlations_stress_sens <- anova_on_nested_df(df_ = correlations_, value_col_name = 'logFC', trait_col_name = 'Stress_sensitivity_clean', main_test = 'anova', post_hoc = 'tukeyhsd') # 34 vs 0 
# 
# correlations_stress <- anova_on_nested_df(df_ = correlations_, value_col_name = 'logFC', trait_col_name = 'Stress_clean') # 4091 vs 2687
# 
# correlations_brain_part <- anova_on_nested_df(df_ = correlations_, value_col_name = 'logFC', trait_col_name = 'Brain_part_clean') #16213 vs 15183
# correlations_gender <- correlations_[1:100,]
# 
# correlations_gender$k_w_pval <- as.numeric(purrr::map(.x = correlations_gender$data, ~ { broom::tidy(kruskal.test(logFC ~ Gender_clean, data = .x))$p.value[[1]] }))
# 
# correlations_gender$k_w_should_i_post <- ifelse(test = correlations_gender$k_w_pval < 0.05, T, F)
# 
# correlations_gender <- subset(x = correlations_gender, subset = correlations_gender$k_w_should_i_post == T) %>%
#   dplyr::select(-k_w_should_i_post)
# 
# correlations_gender$conover <-
#   purrr::map(.x = correlations_gender$data, ~ {
#     tryCatch(
#       {conover.test::conover.test(
#         x = .x$logFC,
#         g = .x$Gender_clean,
#         method = 'bh'
#       ) },
#       error = function(c){return(NA)})})
# correlations_$stress_cor <- as.numeric(purrr::map(.x = correlations_$data, .f = function(x){ cor(x = x$logFC, y = as.integer(x$Stress_clean), use = "pairwise.complete.obs")  }))
# correlations_$brain_part_cor <- as.numeric(purrr::map(.x = correlations_$data, .f = function(x){ cor(x = x$logFC, y = as.integer(x$Brain_part_clean), use = "pairwise.complete.obs")  }))
# 
# correlations_$species_cor <- as.numeric(purrr::map(.x = correlations_$data, .f = function(x){ cor(x = x$logFC, y = as.integer(x$Species), use = "pairwise.complete.obs")  }))
# 
# correlations_$stress_sens_Pear_cor_pw_cm_obs <- as.numeric(purrr::map(.x = correlations_$data, .f = function(x){ cor(x = x$logFC, y = as.integer(x$Stress_sensitivity_clean), use = "pairwise.complete.obs")  }))
# correlations_$stress_sens_Spear_cor <- as.numeric(purrr::map(.x = correlations_$data, .f = function(x){ cor(x = x$logFC, y = as.integer(x$Stress_sensitivity_clean), method = "spearman")  }))

# correlations_$gender_Pear_cor_pw_cm_obs <- as.numeric(purrr::map(.x = correlations_$data, .f = function(x){ cor(x = x$logFC, y = as.integer(x$Gender_clean), use = "pairwise.complete.obs")  }))
# correlations_$gender_Spear_cor <- as.numeric(purrr::map(.x = correlations_$data, .f = function(x){ cor(x = x$logFC, y = as.integer(x$Gender_clean), method = "spearman")  }))
# t-test/non-parametric/normality
# 
# test_1 <- purrr::map(.x = correlations_$data, ~ { broom::tidy(aov(logFC ~ Gender_clean, data = .x))$p.value[1] })







  

### GETING CORRELATIONS ###








# temp <- subset(exp_number_and_percentage_, exp_number_and_percentage_$no_of_exps == 0)
# temp_outer_join <- dplyr::anti_join(x = exp_number_and_percentage_, y = paper_number_, "lower_final_gene_name")
# 
# temp_final_dataset <- subset(x = final_good_dataset_1_and_2, subset = final_good_dataset_1_and_2$lower_final_gene_name %in% temp_outer_join$lower_final_gene_name)
# 
# readr::write_tsv(x = temp_outer_join, path = 'genes_present_only_in_experiments.tsv')
# readr::write_tsv(x = temp_final_dataset, path = 'genes_present_only_in_experiments_dataset_used_as_input_for_data_analysis.tsv')
#     
# 
# length(paper_number_[[1]]) + length(temp_outer_join[[1]]) == length(exp_number_and_percentage_[[1]])
# 
# 
# temp_final_dataset <- subset(x = final_good_dataset_1_and_2, subset = final_good_dataset_1_and_2$Paper == 45)






# ### This will be our input data for further analysis. It contains only significant columns of "Paper", "GroupID", "ensembl_gene_name", "logFC": 
# SHORT_SINGLE_TEST_ANNOTATION <- SINGLE_TEST_ANNOTATION %>%
#   select("Paper", "GroupID", "ensembl_gene_name", "logFC") %>%
#   group_by(Paper, GroupID, ensembl_gene_name) %>%
#   nest()
# 
# # I dont know why purrr::map didnt work for this. Either way, we just lapply this. Here we write UP or DOWN for every logFC for given genename in given experiment (not paper)
# SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$data, FUN = function(x) { mutate(x, Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) } ) 
# 
# # Here we count how many UPs and DOWNs are in every given genename in given experiment (not paper)
# SHORT_SINGLE_TEST_ANNOTATION$directionality <- lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, FUN = function (x) {table(select(x, "Symbol_direction"))} ) 
# 
# # Here we establish actual status of gene in given experiment
# SHORT_SINGLE_TEST_ANNOTATION$sum_directionality <- as.character(lapply(X = SHORT_SINGLE_TEST_ANNOTATION$directionality, 
#                                                                        FUN = function(x) {
#                                                                          if (length(x) == 1 && grepl(pattern = "UP", names(x))) { "UP" }
#                                                                          else if (length(x) == 1 && grepl(pattern = "DOWN", names(x))) { "DOWN" }
#                                                                          else if (length(x) == 2 && grepl(pattern = "DOWN", names(x)) && grepl(pattern = "DOWN", names(x))) { "MIXED" }
#                                                                          else { "ERROR" }
#                                                                        }))
# 
# # Here we remove MIXED expression and multiple genenames
# FILT_SHORT_SINGLE_TEST_ANNOTATION <- SHORT_SINGLE_TEST_ANNOTATION %>%
#   #filter(!grepl(pattern = "(.*);(.*)", SHORT_SINGLE_TEST_ANNOTATION$ensembl_gene_name)) %>%
#   filter(!sum_directionality == "MIXED") %>%
#   filter(!sum_directionality == "ERROR") %>%
#   filter(!is.na(ensembl_gene_name))
# 
# # Here we add mean to each !!! Perhaps this should be final input, because it has: 1) removed genes giving UP+DOWN in the same experiment, 2) removed NAs
# STDINPUT_FILT_SHORT_SIN_T_ANNO <- FILT_SHORT_SINGLE_TEST_ANNOTATION %>%
#   mutate(mean = as.numeric(map(data, function(x) { as.numeric(mean(x[[1]]))} ))) %>% ## This way we give actuall vector to function, not a data table(tibble)
#   select("Paper", "GroupID", "ensembl_gene_name", "sum_directionality", "mean") 
# ### !!! THE IDEA IS THAT THE "STDINPUT_FILT_SHORT_SIN_T_ANNO" IS THE INPUT FOR ALL FURTHER INQUIRIES !!! ###
# 
# # tO raczej nie wyjdzie - nie da siÄ™ wyciÄ…gnÄ…Ä‡ wĹ‚aĹ›ciwej wartoĹ›ci ekspresji z dwĂłch rĂłznych eksperymentĂłw w tym samym paper. MoĹĽe lepiej po prostu dowiedzieÄ‡ siÄ™, ktĂłre geny sÄ… unikatowe dla danego paper i usunÄ…Ä‡ te geny z normalnie uĹĽywanego wczeĹ›niej test annotation pliku. 
# 
# 
# 
# 
# #######################################
# ####### PREPARE FINALIZED TABLE ####### 
# #######################################
# 
# check_was_the_spliting_of_df_by_filtering_ok(df_original = raw_dataset, list_df_splited = list(finalized, leftovers_for_checking_data_integrity))
# 
# length(finalized[[1]]) + length(leftovers_for_checking_data_integrity[[1]])




















































##############################################################
##############################################################
##############################################################







# ###### WHOLE DATASET ANALYSIS ######
# 
# # Tutaj liczymy ile razy geny wystepuja w oryginalnym dataset (w eksperymentach, nie w paperach)
# HOW_MANY_TIMES_EXP_STDINPUT <- STDINPUT_FILT_SHORT_SIN_T_ANNO %>%
#   select(ensembl_gene_name) %>%
#   group_by(ensembl_gene_name) %>%
#   summarise(number = n())
# 
# # Tutaj liczymy ile razy geny wystepuja w oryginalnym dataset (w paperach)
# HOW_MANY_TIMES_PAPER_STDINPUT <- STDINPUT_FILT_SHORT_SIN_T_ANNO %>%
#   select(Paper, ensembl_gene_name) %>%
#   unique() %>%
#   select(ensembl_gene_name) %>%
#   group_by(ensembl_gene_name) %>%
#   summarise(number = n())
# 
# ###### Tutaj robimy specjalny input do klastrowania, w ktĂłrym uĹĽywamy tylko genĂłw, ktĂłre zostaĹ‚y wykryte przynajmniej w 3 oddzielnych paperach ###### 
# UP3_HOW_MANY_TIMES_PAPER_STDINPUT <- HOW_MANY_TIMES_PAPER_STDINPUT %>%
#   filter(number >= 3)
# 
# UP3_PAPER_CLUSTERING_INPUT <- merge(STDINPUT_FILT_SHORT_SIN_T_ANNO, UP3_HOW_MANY_TIMES_PAPER_STDINPUT, by = "ensembl_gene_name", all.y = T)
# 
# FOR_CLUST_UP3_PAPER_CLUSTERING_INPUT <- UP3_PAPER_CLUSTERING_INPUT %>%
#   FOR_CLUS()
# #readr::write_tsv(x = FOR_CLUST_UP3_PAPER_CLUSTERING_INPUT, path = "FOR_CLUST_UP3_PAPER_CLUSTERING_INPUT.tsv")  
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###### Tutaj robimy specjalny input do klastrowania, w ktĂłrym uĹĽywamy tylko genĂłw, ktĂłre zostaĹ‚y wykryte przynajmniej w 3 oddzielnych paperach ###### 
# 
# 
# # Tutaj liczymy ile razy geny wyst?puj? w oryginalnym dataset, patrz?c czy s? up czy down
# UorDWHOLE_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
#   select(Gene_symbol, logFC) %>%
#   mutate(Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) %>%
#   mutate(Symbol_direction = paste(Gene_symbol, Symbol_direction, sep = "_")) %>%
#   group_by(Symbol_direction) %>%
#   summarise(number = n()) %>% 
#   mutate(Gene_symbol2 = str_remove(Symbol_direction, "_.*"))
# 
# ###### WHOLE DATASET ANALYSIS ######
# 
# 
# 
# ###### COMPARISONS-CENTERED ANALYSIS ######
# 
# 
# 
# ### Here we set whether we want to analyze papers or comparisons
# P_or_C = quo(Paper) #" GroupID OR Paper "
# 
# 
# 
# # Tutaj liczymy ile razy geny wyst?puj? W KA?DYM Z POR?WNA?, nie patrz?c czy s? up czy down
# COMP_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
#   select(!!P_or_C, Gene_symbol) %>%
#   group_by(!!P_or_C, Gene_symbol) %>%
#   summarise(number = n())
# 
# 
# 
# # Tutaj liczymy ile razy geny wyst?puj? W KA?DYM Z POR?WNA?, patrz?c czy s? up czy down    
# UorDCOMP_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
#   select(!!P_or_C, Gene_symbol, logFC) %>%
#   mutate(Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) %>%
#   mutate(Symbol_direction = paste(Gene_symbol, Symbol_direction, sep = "_")) %>%
#   group_by(!!P_or_C, Symbol_direction) %>%
#   summarise(Sym_dir_number = n()) %>%
#   mutate(Gene_symbol2 = str_remove(Symbol_direction, "_.*"))
# 
# 
# 
# #Divide data into genes expressed in single direction in given comparison, vs genes expressed in different direction (bad genes)
# nonUNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UorDCOMP_NO_UNIDS_ORG_DATA %>%
#   group_by(!!P_or_C) %>%
#   filter(duplicated(Gene_symbol2, fromLast = T) | duplicated(Gene_symbol2))
# 
# UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UorDCOMP_NO_UNIDS_ORG_DATA %>%
#   group_by(!!P_or_C) %>%
#   filter(!duplicated(Gene_symbol2, fromLast = T) & !duplicated(Gene_symbol2))
# 
# 
# 
# # Check if unique/duplicated division went well           
# if (nrow(UorDCOMP_NO_UNIDS_ORG_DATA) - (nrow(nonUNIQ_UorDCOMP_NO_UNIDS_ORG_DATA) + nrow(UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA)) != 0){ 
#   stop("Hey, fwend! You have some wierd values in Your counted data, buddy! Better check whats happening, or Your results will smell of moose scrotum!")
# }
# 
# 
# 
# # Here we make a table only with genes that were replicated in few comparisons
# REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA %>%
#   filter(Sym_dir_number >= 3)
# 
# 
# #Annotate base on Paper OR GroupID
# ANNO_REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- merge(REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA, COMPARISONS, by = "Paper")
# 
# 
# ###### COMPARISONS-CENTERED ANALYSIS ######