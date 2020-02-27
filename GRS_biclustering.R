# Required packages: biclust, tidyverse, isa2, runibic

#############################################
############## FUNCTIONS ####################
#############################################

## Watch out for slash/backslash in path name
get_the_biclusters <- function(named_matrix_input_matrix, BiClust_resulting_biclustering)
{
  dir_name <- paste0('clusters_', substitute(BiClust_resulting_biclustering))
  
  dir.create(dir_name)
  
  biclusters <- list()

  for (n in seq(BiClust_resulting_biclustering@Number)) 
  {
    biclusters[[n]] <- bicluster(
      named_matrix_input_matrix, 
      BiClust_resulting_biclustering, 
      number = n)
    
      write.table(
      x = as.data.frame(biclusters[[n]]), 
      file = paste0(dir_name, '/', substitute(BiClust_resulting_biclustering), '_', n,  '.tsv'), 
      sep = '\t',
      row.names = T,
      col.names = NA,
      dec = ','
    )
  }
  return(biclusters)
}



dev_off_looped <- function()
{
  if (!is.null(dev.list() )) 
  {
    dev.off()
  }
}



draw_bicluster_with_some_function <- function(function_to_use, named_matrix_input_matrix, BiClust_resulting_biclustering, str_name = '')
{
  dev_off_looped()
  
  function_name <- deparse(substitute(function_to_use))
  function_name <- gsub('(.)*::', '', function_name)
  
  dir_name <- paste0( 
    function_name, 
    '_', 
    substitute(named_matrix_input_matrix),
    '_', 
    substitute(BiClust_resulting_biclustering))
  
  dir.create(dir_name)
  
  
  for (cluster in seq(BiClust_resulting_biclustering@Number))
  {
    pdf_name <- paste0( dir_name, '/', substitute(BiClust_resulting_biclustering), '_', cluster, ".pdf" )
    pdf(pdf_name)
    
    ## Putting try here allows walthrough for this problem: ISA2 SPITS OUT CLUSTERS OF LENGHT 1, WHICH DO NOT HAVE ROW OR COL NAMES. THESE CUNTS CAUSE PROBLEMS.
    try(
      switch (function_name,
              'function_name' = function_to_use(
                x = named_matrix_input_matrix,
                bicResult = BiClust_resulting_biclustering,
                number = cluster),
              'quheatmap' = function_to_use(
                x = named_matrix_input_matrix,
                bicResult = BiClust_resulting_biclustering,
                number = cluster,
                showlabel = T),
              'qunetwork' = qunetwork_loop_wrapper(
                x = named_matrix_input_matrix, 
                BicRes = BiClust_resulting_biclustering, 
                number = cluster,
                pdf_name_ = pdf_name
              )
      ))
    
    dev_off_looped()
  }
}



qunetwork_loop_wrapper <- function(x, BicRes, number, pdf_name_)
{
  net <- QUBIC::qunetwork(x = x, BicRes = BicRes, number = number)
  
  dev_off_looped()
  pdf(pdf_name_)
  
  qgraph::qgraph(net[[1]], 
                 groups = net[[2]], 
                 layout = 'spring', 
                 minimum = 0.6,
                 color = cbind(rainbow(length(net[[2]]) - 1), 'gray'), 
                 edge.label = FALSE)
  
  dev_off_looped()
}



# Extract a dataframe from list, add column with list element name and write/return it as a single dataframe
print_diagnostic_table <- function(diagnosticTest_result, str_name = '')
{
  temp_list <- list()
  for (bicluster in seq_along(diagnosticTest_result))
  {
    temp_df <- diagnosticTest_result[[bicluster]]$table
    temp_df <- dplyr::mutate(.data = temp_df, 
                             Bicluster = names(diagnostic_testisa2[bicluster]))
    temp_list[[bicluster]] <- temp_df
  }
  temp_list_df <- rlist::list.rbind(temp_list)
  readr::write_tsv(temp_list_df, paste0('diagnosticTest_results_', str_name, '.tsv'))
  return(temp_list_df)
}



# TS bicluster stratification score - Highly negative TS1(b) (< -1) shows B-type co-expression, highly positive TS1(b) (> 1) shows T-type co-expression and TS1(b) close to zero (-1< TS1(b) <1) indicates ?-type co-expression
# SB differential co-expression score - A bicluster with no co-expression of any type for the bicluster genes in the non-bicluster conditions is the true bicluster. Comparable co-expression in the non-bicluster conditions means the conditions in the bicluster are not distinctive enough from the remaining conditions and hence do not qualify to be a bicluster. In such a case, the bicluster genes with all conditions in the study can be considered as a gene cluster with a strong co-expression across all conditions. Strong positive SB(b) indicates strong co-expression in G1 and weaker or no co-expression in G2 vice versa.
get_ChiaKaruturi_results_for_all_clusters <- function(named_matrix_input_matrix, BiClust_resulting_biclustering, str_name = '')
{
  temp_list <- list()
  
  for (bicluster in seq(BiClust_resulting_biclustering@Number)) 
  {
    temp_df <- biclust::ChiaKaruturi(x = named_matrix_input_matrix, 
                                     bicResult = BiClust_resulting_biclustering, 
                                     number = bicluster)
    temp_df <- dplyr::mutate(.data = temp_df, Bicluster = bicluster)
    temp_list[[bicluster]] <- temp_df
  }
  
  temp_list_df <- rlist::list.rbind(temp_list)
  readr::write_tsv(temp_list_df, paste0('ChiaKaruturi_results_', str_name, '.tsv'))
  return(temp_list_df)
}



# int_cutoof_for_nb_of_genes_int_term - if enricher returns term with number of annotated genes lower than this number, these terms are not returned
enrich_biclusters <- function(list_get_the_biclusters_results, int_cutoof_for_nb_of_genes_int_term = 1)
{
  dbs_Onto_Path <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "MGI_Mammalian_Phenotype_2017", "Human_Phenotype_Ontology", "KEGG_2016", "WikiPathways_2016", "Panther_2016", "Reactome_2016", "BioCarta_2016", "NCI-Nature_2016", "ARCHS4_Kinases_Coexp", "HumanCyc_2016", "BioPlex_2017", "SILAC_Phosphoproteomics")
  dbs_Regul <- c("Genome_Browser_PWMs", "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "ChEA_2016",  "TF-LOF_Expression_from_GEO", "PPI_Hub_Proteins", "ENCODE_TF_ChIP-seq_2015", "ENCODE_Histone_Modifications_2015", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "CORUM", "Pfam_InterPro_Domains", "Phosphatase_Substrates_from_DEPOD", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp", "miRTarBase_2017", "TargetScan_microRNA_2017", "Enrichr_Submissions_TF-Gene_Coocurrence", "Epigenomics_Roadmap_HM_ChIP-seq")
  dbs_Drug_Tissue_Other <- c("Jensen_TISSUES", "ARCHS4_IDG_Coexp", "DrugMatrix", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "OMIM_Disease", "Jensen_DISEASES", "DSigDB",  "Jensen_COMPARTMENTS", "ARCHS4_Tissues", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB", "Mouse_Gene_Atlas", "ESCAPE", "Chromosome_Location", "MSigDB_Computational", "dbGaP", "Genes_Associated_with_NIH_Grants", "GeneSigDB")
  dbs <- c(dbs_Onto_Path, dbs_Regul, dbs_Drug_Tissue_Other)
  
  temp_list <- list()
  
  temp_list <- lapply(X = list_get_the_biclusters_results, FUN = function(X)
  {
    ALL_ENRICHR <- rlist::list.rbind(enrichR::enrichr(rownames(X[[1]]), dbs))
    ALL_ENRICHR <- subset(x = ALL_ENRICHR, subset = as.integer(sub('/(.)*', '', ALL_ENRICHR$Overlap)) > int_cutoof_for_nb_of_genes_int_term)
  })
  
  return(temp_list)
}


#############################################
############## FUNCTIONS ####################
#############################################

library(tidyverse)


#### LOAD DATA #### 
# data(SyntrenEcoli)
# 
# test <- SyntrenEcoli
# 
# test <- readr::read_tsv('data_v3_names_mm.txt')
# test2 <- test %>%
#   select(Experiment, Probe_ID, logFC) %>%
#   mutate(logFC = as.numeric(stringr::str_replace(logFC, ',', '.'))) %>%
#   mutate(filt = paste0(Experiment, Probe_ID)) %>%
#   mutate(remov = duplicated(filt)) %>%
#   filter(remov == F) %>%
#   select(Experiment, Probe_ID, logFC)
# test3 <- test2 %>%
#   tidyr::spread(key = Experiment, value = logFC, fill = 0)
# test4 <- as.matrix(test3[2:24])
# rownames(test4) <- test3$Probe_ID
# 

# MOVED TO GRS_post_annotation.R
# for_clustering <- as.data.frame(at_least_in_3_papers_spread_med_fin_g_ds)
# rownames(for_clustering) <- for_clustering$lower_final_gene_name
# for_clustering$lower_final_gene_name <- NULL
# for_clustering <- as.matrix(for_clustering)
# for_clustering[is.na(for_clustering)] <- 0
# 
# 
# # For cluster-based hierachical clustering
# write.table(
#   x = as.data.frame(for_clustering), 
#   file = 'for_clustering.tsv',
#   sep = '\t',
#   row.names = T,
#   col.names = NA,
#   dec = '.'
# )  

  # ~/cluster-1.57/src/cluster -f for_clustering.tsv -m a -g 2 #http://bonsai.hgc.jp/~mdehoon/software/cluster/command.txt 

### Clustering subsets of data
search_for = stringr::str_detect(string = tolower(descriptions$Brain_part), pattern = '.*cumbens*')

brain_areas <- descriptions %>%
  dplyr::select(Group_ID, Brain_part) %>%
  dplyr::filter(search_for)

temp_other <- descriptions %>%
  dplyr::select(Group_ID, Brain_part) %>%
  dplyr::filter(!search_for)

filter_out_missing_exps <- sapply(X = brain_areas$Group_ID, function(x) {x %in% colnames(for_clustering)})
experiments_to_include <- brain_areas$Group_ID[filter_out_missing_exps]

subset_matrix <- for_clustering[, experiments_to_include]

write.table(
  x = as.data.frame(subset_matrix), 
  file = 'for_clustering_hippocampus.tsv',
  sep = '\t',
  row.names = T,
  col.names = NA,
  dec = '.'
)





#### LOAD DATA #### 


#### LEARN ABOUT THE DATA? #### 

pca_between_individuals <- prcomp(x = for_clustering, scale = FALSE)
pca_between_individuals <- prcomp(x = t(test), scale = FALSE)
pca_between_variables <- princomp(x = for_clustering, cor = FALSE, scores = TRUE)

# Show the percentage of variances explained by each principal component
factoextra::fviz_eig(pca_between_individuals)
factoextra::fviz_eig(pca_between_variables)

# Individuals with a similar profile are grouped together
factoextra::fviz_pca_ind(pca_between_individuals,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

factoextra::fviz_pca_var(pca_between_individuals,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


#### LEARN ABOUT THE DATA? #### 


### ACTUALL ANALYSIS ###

# This algorithm is chosen because it acted generally well (A systematic comparative evaluation of biclustering techniques.). It gave reasonably small amount of clusters. It detects Upregulated biclusters well
# isa_for_clustering <- isa2::isa(for_clustering) # perform clustering on a named matrix using ISA
# isa_for_clustering <- isa2::isa.biclust(isa_for_clustering) # convert results to biclust results



# This algorithm is chosen because it acted generally well (A systematic comparative evaluation of biclustering techniques.). It gave reasonably small amount of clusters. It detects Shift (and Upregulated) biclusters well
# testplaid <- biclust::biclust(x = for_clustering, method = biclust::BCPlaid(), 
#                               cluster="b", fit.model = y ~ m + a + b, background = TRUE, background.layer = NA, 
#                               background.df = 1, row.release = 0.7,col.release = 0.7, shuffle = 3, 
#                               back.fit = 0, max.layers = 20, iter.startup = 5,iter.layer = 10, verbose = TRUE)

# This algorithm is chosen because it should deal well with constant biclusters (A systematic comparative evaluation of biclustering techniques.): biclusters with a constant expression level close to dataset mean. The constant expression values of the biclusters were chosen to be 0; background values were independent, identically-distributed (i.i.d.) draws from the standard normal: N(0, 1). I dont think it works that way! But it does detect constant clusters well, while isa and plaid do not
# testcc <- biclust::biclust(x = test, method = biclust::BCCC(),
#                            delta = 1.0, alpha=1.5, number=100)



    #t - consistency level of the block (0.5-1.0] (the lower the t, the fewer clusters are found?). q - a double value for quantile discretization f - filtering overlapping blocks (default 1 do not remove any blocks) nbic - maximum number of biclusters in output div - number of ranks for up(down)-regulated genes: default: 0==ncol(x)
# Is 'block' the same as 'bicluster'? What is 'consistency level'?
unibic_for_clustering_t05_f0 <- biclust::biclust(x = for_clustering, method = runibic::BCUnibic(),
                               t = 0.5, q = 0, f = 0, nbic = 100, div = 0, useLegacy = F)
save(unibic_for_clustering_t05_f0, file = 'unibic_for_clustering_t05_f0')
biclusters_unibic_for_clustering <- get_the_biclusters(named_matrix_input_matrix = for_clustering, BiClust_resulting_biclustering = unibic_for_clustering_t05_f0)

unibic_for_clustering_t08_f0 <- biclust::biclust(x = for_clustering, method = runibic::BCUnibic(),
                                                 t = 0.8, q = 0, f = 0, nbic = 100, div = 0, useLegacy = F)
save(unibic_for_clustering_t08_f0, file = 'unibic_for_clustering_t08_f0')
biclusters_unibic_for_clustering <- get_the_biclusters(named_matrix_input_matrix = for_clustering, BiClust_resulting_biclustering = unibic_for_clustering_t08_f0)


# https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/34/24/10.1093_bioinformatics_bty512/1/bty512_supplementary_material.pdf?Expires=1567872409&Signature=1utz6dd-Vf9AQJNpsNkAuGaRkZqa1LO72cWhkiyu5u-CiD5xA5mcFBso128bdfXawVD2Wkyryd~goOc3JYHzgRMV2N4AMBD1j6rIXoph1MlrtaZFC-ubbUDQoZpiLuIXG9tp4uTI4AMTaKeIbkuyXFdNp554ihZ7uw-HopQRJKeBVLBUBzxu6KILolJPwBEkcb3Yioera6zgYkWelE-0D-yktilWL1rY4ZMq~2SWISvsIdAMgb489DBqOndjihTbkOHdO8FRaUceEw2E~P~Vur1NKxFALjbyZSyasRpzFHFaLu4IwmL5k9CEuOMVG7LBSc12aHQD~mlc4RTtO35CnA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA



draw_bicluster_with_some_function(function_to_use = biclust::drawHeatmap, 
                                  named_matrix_input_matrix = test, 
                                  BiClust_resulting_biclustering = testisa2)
draw_bicluster_with_some_function(function_to_use = QUBIC::quheatmap, 
                                  named_matrix_input_matrix = test, 
                                  BiClust_resulting_biclustering = testisa2)
draw_bicluster_with_some_function(function_to_use = QUBIC::qunetwork, 
                                  named_matrix_input_matrix = test, 
                                  BiClust_resulting_biclustering = testisa2)

biclust::heatmapBC(x = test, bicResult = testisa2, number = 0)

biclust::parallelCoordinates(x = test, bicResult = testisa2, number = 2)



### BICLUSTERING METRICS ###

# Which samples compose which cluster?
biclust::biclustbarchart(x = test, Bicres = testunibic)
biclust::biclustmember(bicResult = testunibic, x = test)

# Which clusters are similar to eachother (also between clusterings)? Color represents the bi-cluster set to which bicluster pertains (up to three bicluster sets can be represented simultaneously).Brightness represents the bicluster homogeneity (darker, less homogeneous).  Size represents thesize of the bicluster, as (number of genes)x(number of conditions).  Location is a 2D-projection ofgene and condition profiles.
biclust::bubbleplot(x = test, 
           bicResult1 = testunibic, 
           bicResult2 = testisa2, 
           bicResult3 = testplaid, 
           projection = "mean", 
           showLabels = T)


# How similar are different clusterings? Jaccard Index (measure of similarity between two sets) Perhaps similarity between https://rdrr.io/rforge/biclust/man/jaccard.html
biclust::jaccardind(bicres1 = testisa2, bicres2 = testunibic)

jaccard2


### BICLUSTERING METRICS ###




### BICLUSTERING DIAGNOSTICS ###
diagnostic_testisa2 <- biclust::diagnosticTest(BCresult = testisa2, 
                                       data = test, 
                                       statistics = c("F", "Tukey"),
                                       nSim = 200, 
                                       alpha = 0.05, 
                                       save_F = T)
# biclust::diagnosticPlot2(diagnosticTest = diagnostic_testisa2, 
#                 number = 1, 
#                 StatVal = T,
#                 binwidth = NULL)


results_diagnostic_testisa2 <- print_diagnostic_table(diagnostic_testisa2, 'testisa2')




diagnostic_testunibic <- biclust::diagnosticTest(BCresult = testunibic, 
                                                 data = test, 
                                                 statistics = c("F", "Tukey"),
                                                 nSim = 200, 
                                                 alpha = 0.05, 
                                                 save_F = T)

diagnostic_CK_testunibic <- get_ChiaKaruturi_results_for_all_clusters(test, testunibic)

# diagnostic_CR_testunibic <- diagnoseColRow(x = test, 
#                                            bicResult = testunibic, 
#                                            number = 1, 
#                                            nResamplings = 200, 
#                                            replace = TRUE)
# diagnosticPlot(bootstrapOutput = diagnostic_CR_testunibic)
# diagnostic_OF_testunibic <- computeObservedFstat(x = test, 
#                                                  bicResult = testunibic, 
#                                                  number = 1)



http://bioconductor.org/packages/release/bioc/html/GOstats.html
### BICLUSTERING DIAGNOSTICS ###






QUBIC::showinfo(matrix = test, bic = c(testunibic, testisa2))





# biclust::writeBiclusterResults(fileName = 'testisa2.txt', bicResult = testisa2, bicName = 'ass', geneNames = rownames(test4), arrayNames = colnames(test4), append=FALSE, delimiter="\t")
# 
# biclust::writeBiclusterResults(fileName = 'testunibic.txt', bicResult = testunibic, bicName = 'ass', geneNames = rownames(test4), arrayNames = colnames(test4), append=FALSE, delimiter="\t")





# ISA has two parameters, tg and tc, to prune genes and
# experimental conditions out of the bicluster being formed.
# For the experiments with the first synthetic data collec-
#   tion, the values suggested by [24] were too restrictive and
# did not allow the algorithm to detect any bicluster model.
# So, in all of our experiments, we also tested the default
# values of the isa2 package [41], which considers all the
# possible tg and tc combinations in the interval [ 1.0, 3.0]
# with steps of 0.5.








xxx <- enrich_biclusters(testisa2_biclutrs)



#Filter results based on this This works, but does it work in loop?
xxx2 <- filter(x = xxx[[1]], filter = as.integer(sub('/(.)*', '', xxx[[1]]$Overlap)) > 1)









KAJA:
  
  setwd("E:/Projekty/Kaja Review LDH")

library(enrichR)

dbs_Onto_Path <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "MGI_Mammalian_Phenotype_2017", "Human_Phenotype_Ontology", "KEGG_2016", "WikiPathways_2016", "Panther_2016", "Reactome_2016", "BioCarta_2016", "NCI-Nature_2016", "ARCHS4_Kinases_Coexp", "HumanCyc_2016", "BioPlex_2017", "SILAC_Phosphoproteomics")
dbs_Regul <- c("Genome_Browser_PWMs", "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "ChEA_2016",  "TF-LOF_Expression_from_GEO", "PPI_Hub_Proteins", "ENCODE_TF_ChIP-seq_2015", "ENCODE_Histone_Modifications_2015", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "CORUM", "Pfam_InterPro_Domains", "Phosphatase_Substrates_from_DEPOD", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp", "miRTarBase_2017", "TargetScan_microRNA_2017", "Enrichr_Submissions_TF-Gene_Coocurrence", "Epigenomics_Roadmap_HM_ChIP-seq")
dbs_Drug_Tissue_Other <- c("Jensen_TISSUES", "ARCHS4_IDG_Coexp", "DrugMatrix", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "OMIM_Disease", "Jensen_DISEASES", "DSigDB",  "Jensen_COMPARTMENTS", "ARCHS4_Tissues", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB", "Mouse_Gene_Atlas", "ESCAPE", "Chromosome_Location", "MSigDB_Computational", "dbGaP", "Genes_Associated_with_NIH_Grants", "GeneSigDB")
dbs <- c(dbs_Onto_Path, dbs_Regul, dbs_Drug_Tissue_Other)

genes <- c("ldha", "ldhb")

ALL_ENRICHR <- enrichR::enrichr(genes, dbs)
ALL_ENRICHR_DATA <- rlist::list.rbind(ALL_ENRICHR)
write.table(ALL_ENRICHR_DATA, "ALL_ENRICHR_DATA.txt", sep="\t")

getwd()

tp <- read.table("TF_PPI_FROM_ENRICHR.txt", header = T, stringsAsFactors = F)

TP_ENRICHR <- enrichR::enrichr(tp$Enrichr_TFs_PPIs_unique, dbs_Onto_Path)
ALL_TP_ENRICHR_DATA <- rlist::list.rbind(TP_ENRICHR)
write.table(ALL_TP_ENRICHR_DATA, "ALL_TP_ENRICHR_DATA.txt", sep="\t")























######################
##### CLUSTERING #####
######################



FOR_CLUS <- function(data__) 
{
  CLU_data__ <- data__ %>%
    select("GroupID", "ensembl_gene_name", "mean") %>% 
    spread(key = GroupID, value = mean)  %>%
    as.data.frame()
  
  # We change NAs to 0, because we assume that all genes that we do not have in our data set are not differentially expressed, so logFC should be equal to 0
  ZEROED_CLU_data__ <- CLU_data__
  ZEROED_CLU_data__[is.na(ZEROED_CLU_data__)] <- 0 ###!!!
  
  # Nie rozumiem zupeĹ‚nie czemu, ale matrixy wariujÄ… jak siÄ™ applyuje im po rowach, data framey dziaĹ‚ajÄ… perfekcyjnie. Anyway, tutaj wyciÄ…gamy do klastrowania tylko geny ktĂłre majÄ… przynajmniej 3 wartoĹ›ci (pomijajÄ…c entrez_gene column)
  ZEROED_CLU_data__$filter <- apply(X = ZEROED_CLU_data__[,-1], MARGIN = 1, FUN = function (x) 
  { sum(abs(x) > 0.5) })  ###!!!
  
  UP3_ZEROED_CLU_data__ <- dplyr::filter(ZEROED_CLU_data__, filter >= 3)
  UP3_ZEROED_CLU_data__$filter <- NULL
  
  return(UP3_ZEROED_CLU_data__)
}




UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO <- FOR_CLUS(STDINPUT_FILT_SHORT_SIN_T_ANNO)
#readr::write_tsv(x = UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, path = "MULTIPLE_NAMES_INCLUDED_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO.tsv")

#Mamy tutaj duĹĽo sond, ktĂłre nie majÄ… wartoĹ›ci w ĹĽadnym z eksperymentĂłw. Czy to sondy natywne dla usuniÄ™tych dwĂłch eksperymentĂłw, czy coĹ› dziwnego siÄ™ dzieje? - dzieje siÄ™ to, ĹĽe w oryginalnych danych sÄ… sondy, ktĂłre majÄ… log FC bliski zeru do dwĂłch miejsc po przecinku - 




### runibic ###

# Change the prepared data frame to matrix for clustering purposes
rownames(UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO) <- UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO$ensembl_gene_name
UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO$ensembl_gene_name <- NULL
MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO <- as.matrix(UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO)


BiocManager::install(c("runibic", "biclust", "gplots"))
library(runibic)

runibicRES2 <- biclust::biclust(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, method = runibic::BCUnibic())

dir.create("ClusterSummaryRunibic")

# Return figure representing which groups are included in which clusters
png(filename = paste0("ClusterSummaryRunibic/", names(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO), "_Cluster_Summary.png"), width = 1920, height = 1080, units = "px")
par(cex.axis = 1, oma=c(0,4,0,8))
biclust::biclustmember(runibicRES2, MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, main = names(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO))
dev.off()
# png(filename = paste0("ClusterSummaryRunibic/", the_names[n], "_Cluster_Summary.png"), width = 1920, height = 1080, units = "px")
# par(cex.axis = 1, oma=c(0,4,0,8))
# biclustmember(runibicRES2[[n]], matrixed_annot_the_list[[n]], main = the_names[n])
# dev.off()


# Shows which clusters are similar to other clusters (multidimensional scaling - similar to PCA)
png(filename = paste0("ClusterSummaryRunibic/", names(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO), "_Cluster_Relatedness.png"), width = 3840, height = 2160, units = "px")
biclust::bubbleplot(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, runibicRES2, showLabels = T)
dev.off()

# This shows entire heatmap
par(cex.axis = 0.5)
biclust::heatmapBC(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, runibicRES2)
axis(1, at=1:dim(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO)[2], labels = colnames(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO), las=2)
axis(2, at=1:dim(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO)[1], labels = rownames(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO), las=2)


biclust::drawHeatmap(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, runibicRES2, number = 1, beamercolor = F, paleta = col_for_heatmaps(1000))



biclust::writeBiclusterResults(fileName = paste0("ClusterSummaryRunibic/", names(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO), "_ClustersRunibic.txt"), 
                               bicResult = runibicRES2, 
                               bicName = paste0(names(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO)),
                               geneNames = rownames(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO), 
                               arrayNames = colnames(MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO),
                               append =FALSE, 
                               delimiter="\t")

### runibic ###



### hierachical ###

dist1 <- dist(x = MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, method = "euclidean")

hclust1 <- hclust(dist1, method="complete")



SHORT_MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO <- MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO[1:100,]


# For server
# UP3 <- read.delim("UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO.tsv", header = T, sep = "\t")
# rownames(UP3) <- UP3$ensembl_gene_name
# UP3$ensembl_gene_name <- NULL
# UP3 <- as.matrix(UP3)
# ../cluster -f UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO.tsv -m a -g 2 

pdf("heatmap.pdf", pointsize = 4)
gplots::heatmap.2(SHORT_MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO, 
                  distfun = function (x) { dist(x, ) }, 
                  hclustfun = function (x) { hclust(x, method = "average") },
                  trace = "none"
)
dev.off()

png(filename = "heatmap.png", width = 3840, height = 2160, units = "px")
heatmap(UP3)
dev.off()

png(filename = "plot.png", width = 3840, height = 2160, units = "px")
plot(hclust1)
dev.off()

install.packages("proxy")

simil1 <- proxy::simil(UP3)
png(filename = "_heatmap.png", width = 3840, height = 2160, units = "px")
heatmap(simil1)
dev.off()

# For server




plot(hclust1)
heatmap(SHORT_MATRIX_UP3_ZEROED_CLU_STDINPUT_FILT_SHORT_SIN_T_ANNO)


TOTAL_GENE_COUNT_PER_EXPERIMENT <- as.tibble(SINGLE_TEST_ANNOTATION) %>%
  group_by(ensembl_gene_name) %>%
  summarise(n())

TOTAL_GENE_COUNT_PER_PUBLICATION <- as.tibble(SINGLE_TEST_ANNOTATION) %>%
  group_by(Paper) %>%
  group_by(ensembl_gene_name) %>%
  summarise(n())


######################
##### CLUSTERING #####
######################





# Here we put info on what to do provided each visualizing function. I dont know why quheatmap doesnt work. It should.
# dplyr::case_when(
#   function_name %in% c('drawHeatmap', 'quheatmap') ~ function_to_use(
#     x = named_matrix_input_matrix,
#     bicResult = BiClust_resulting_biclustering,
#     number = cluster),
#   function_name %in% c('qunetwork') ~ function_to_use(
#     x = named_matrix_input_matrix,
#     BicRes = BiClust_resulting_biclustering,
#     number = cluster)
# )

