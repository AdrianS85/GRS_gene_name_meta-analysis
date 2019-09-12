# Required packages: biclust, tidyverse, isa2, runibic

#############################################
############## FUNCTIONS ####################
#############################################

## Watch out for slash/backslash in path name
get_the_biclusters <- function(named_matrix_input_matrix, BiClust_resulting_biclustering, str_name)
{
  biclusters <- list()
  dir.create(str_name)
  
  for (n in seq(BiClust_resulting_biclustering@Number)) 
  {
    biclusters[[n]] <- bicluster(
      named_matrix_input_matrix, 
      BiClust_resulting_biclustering, 
      number = n)
    
    write.table(
      x = as.data.frame(biclusters[[n]]), 
      file = paste0(str_name, '/', str_name, '_', n,  '.tsv'), 
      sep = '\t',
      row.names = T,
      col.names = T,
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

#############################################
############## FUNCTIONS ####################
#############################################

library(tidyverse)


#### LOAD DATA #### 
data(SyntrenEcoli)

test <- SyntrenEcoli

test <- readr::read_tsv('data_v3_names_mm.txt')

test2 <- test %>%
  select(Experiment, Probe_ID, logFC) %>%
  mutate(logFC = as.numeric(stringr::str_replace(logFC, ',', '.'))) %>%
  mutate(filt = paste0(Experiment, Probe_ID)) %>%
  mutate(remov = duplicated(filt)) %>%
  filter(remov == F) %>%
  select(Experiment, Probe_ID, logFC)


  
test3 <- test2 %>%
  tidyr::spread(key = Experiment, value = logFC, fill = 0)


test4 <- as.matrix(test3[2:24])
rownames(test4) <- test3$Probe_ID

colnames(test4)
library(isa2)

class(test4)
#### LOAD DATA #### 





### ACTUALL ANALYSIS ###

testisa <- isa2::isa(test) # perform clustering on a named matrix using ISA
testisa2 <- isa2::isa.biclust(testisa) # convert results to biclust results

testunibic <- biclust::biclust(x = test, method = runibic::BCUnibic())



testisa2_biclutrs <- get_the_biclusters(named_matrix_input_matrix = test, BiClust_resulting_biclustering = testisa2, str_name = 'testisa2')

testunibic_biclutrs <- get_the_biclusters(named_matrix_input_matrix = test, BiClust_resulting_biclustering = testunibic, str_name = 'testunibic')

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
