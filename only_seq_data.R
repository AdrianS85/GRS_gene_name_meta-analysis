library(Hmisc)
library(magrittr)

# save(seq_data, file = "seq_data.save")

seq_data <- list()

seq_data$input <- purrr::map(
  .x = list(
    "Supplementary Dataset 1" = "Supplementary Dataset 1 v3.xlsx", 
    "Supplementary Table 1" = "Supplementary Table 1 v2 GJ.xlsx", 
    "Supplementary Table 2" = "Supplementary Table 2 v4.xlsx"), 
  .f = function(file){
    
    openxlsx::read.xlsx(xlsxFile = file)
  })



###########
### ST1 ###

seq_data$seq_only$t1 <- subset(
  seq_data$input$`Supplementary Table 1`,
  subset = seq_data$input$`Supplementary Table 1`$Transcriptomics.method == "RNA-Seq")

seq_data$output$t1 <- seq_data$seq_only$t1

colnames(seq_data$seq_only$t1) <- stringr::str_replace_all(string = colnames(seq_data$seq_only$t1), pattern = "\\.", replacement = " ")

### ST1 ###
###########



###########
### SD1 ###

temp <- unique( stringr::str_replace(string = seq_data$seq_only$t1$`Paper/Group code`, pattern = "\\.", replacement = "_") )

seq_data$seq_only$dt1 <- subset(
  seq_data$input$`Supplementary Dataset 1`,
  subset = seq_data$input$`Supplementary Dataset 1`$`Paper/Group.code` %in% temp)

seq_data$output$dt1 <- seq_data$seq_only$dt1

colnames(seq_data$seq_only$dt1) <- stringr::str_replace_all(string = colnames(seq_data$seq_only$dt1), pattern = "\\.", replacement = " ")

### SD1 ###
###########






seq_data$dataset1$medianed_final_good_dataset <- dplyr::summarize(
  dplyr::group_by(seq_data$output$dt1, `Paper/Group.code`, Standarised_gene_symbol), 
  logFC_median = median(logFC))


seq_data$dataset1$spread_medianed_final_good_dataset <- tidyr::spread(
  data = seq_data$dataset1$medianed_final_good_dataset, 
  key = `Paper/Group.code`, 
  value = logFC_median)


seq_data$dataset1$bool_entries_with_3_or_more_values <- get_subset_vector_for_entries_with_3_or_more_values_per_paper(final_good_dataset_ = seq_data[["output"]]$dt1)


seq_data$dataset1$at_least_in_3_papers_spread_med_fin_g_ds <- subset(
  x = seq_data$dataset1$spread_medianed_final_good_dataset, 
  subset = seq_data$dataset1$bool_entries_with_3_or_more_values)


seq_data$dataset1$at_least_in_3_papers_spread_med_fin_g_ds <- subset(
  seq_data$dataset1$at_least_in_3_papers_spread_med_fin_g_ds, 
  subset = !is.na(seq_data$dataset1$at_least_in_3_papers_spread_med_fin_g_ds$Standarised_gene_symbol))


# Are any columns empty?
test <- seq_data$dataset1$at_least_in_3_papers_spread_med_fin_g_ds
test[is.na(test)] <- 0
test <- t(purrr::map_df(.x = test[,2:ncol(test)], .f = sum))
test2 <- subset(test, test == 0)
# Are any columns empty?
### Remove 




################################
### SAVE DATA FOR CLUSTERING ###

# seq_data$dataset1$for_clustering <- as.data.frame(seq_data$dataset1$at_least_in_3_papers_spread_med_fin_g_ds)
# 
# 
# seq_data$dataset1$for_clustering <- dplyr::select(seq_data$dataset1$for_clustering, !c("32_1", "32_2", "61_1", "62_1"))
# 
# seq_data$dataset1$for_clustering[is.na(seq_data$dataset1$for_clustering)] <- 0
# 
# readr::write_tsv(x = seq_data$dataset1$for_clustering, 'for_clustering_v3.tsv')

### SAVE DATA FOR CLUSTERING ###
################################








seq_data$dataset1$spread_medianed_final_good_dataset <- subset(
  seq_data$dataset1$spread_medianed_final_good_dataset, 
  subset = !is.na(seq_data$dataset1$spread_medianed_final_good_dataset$Standarised_gene_symbol))

seq_data$dataset1$spread_medianed_final_good_dataset[is.na(seq_data$dataset1$spread_medianed_final_good_dataset)] <- 0

seq_data$exp_number_and_percentage <- get_number_and_percentage_of_directionality_of_paper_and_exp_first_column_names_wrapper(seq_data$dataset1$spread_medianed_final_good_dataset, '')



##########
### S2 ###

seq_data$dataset1[["Supplementary Table 2 temp"]] <- merge(
  x = seq_data$exp_number_and_percentage$exps, 
  y = seq_data$exp_number_and_percentage$pap, 
  by = "Standarised_gene_symbol")

seq_data$dataset1[["Supplementary Table 2"]] <- merge(
  x = seq_data$dataset1[["Supplementary Table 2 temp"]], 
  y = seq_data$input$`Supplementary Table 2`, 
  by.x = "Standarised_gene_symbol",
  by.y = "Gene.symbol",
  all.x = T)

seq_data$output[["Supplementary Table 2"]] <- dplyr::select(
  seq_data$dataset1[["Supplementary Table 2"]], 
  "Gene symbol" = "Standarised_gene_symbol",
  "Number of reporting papers" = "number",
  "Total number of transcriptomic comparisons with altered expression" = "no_of_exps",
  "Fraction of comparisons with up-regulated expression" = "perc_of_upregulated", "Glucocorticoid-responsive.from.core.list", "Glucocorticoid-responsive.from.extended.list", "Hemoglobin.cluster", "Meningeal.cluster", "Choroid.cluster" 
  )

seq_data$output[["Supplementary Table 2"]][["Fraction.of.comparisons.with.down-regulated.expression"]] <- 1 - seq_data$output[["Supplementary Table 2"]]$`Fraction of comparisons with up-regulated expression`

seq_data$output[["Supplementary Table 2"]][["Number of papers x fraction"]] <- ifelse(
  test = seq_data$output[["Supplementary Table 2"]]$`Fraction of comparisons with up-regulated expression` >= seq_data$output[["Supplementary Table 2"]]$`Fraction.of.comparisons.with.down-regulated.expression`,
  yes = seq_data$output[["Supplementary Table 2"]]$`Number of reporting papers` * seq_data$output[["Supplementary Table 2"]]$`Fraction of comparisons with up-regulated expression`,
  no = seq_data$output[["Supplementary Table 2"]]$`Number of reporting papers` * seq_data$output[["Supplementary Table 2"]]$`Fraction.of.comparisons.with.down-regulated.expression`)

colnames(seq_data$output[["Supplementary Table 2"]]) <- stringr::str_replace_all(string = colnames(seq_data$output[["Supplementary Table 2"]]), pattern = "\\.", replacement = " " )

seq_data$output[["Supplementary Table 2"]] <- dplyr::select(
  seq_data$output[["Supplementary Table 2"]],
  1:4, 10, 11, 5:9)



seq_data$output[["Supplementary Table 2"]]$freq_signif <- 1 - pbinom(
  seq_data$output[["Supplementary Table 2"]][["Total number of transcriptomic comparisons with altered expression"]] - 1, 
  size = 89, 
  prob = 0.05) 

seq_data$output[["Supplementary Table 2"]]$freq_up <- 1 - pbinom(
  as.integer( seq_data$output[["Supplementary Table 2"]][["Total number of transcriptomic comparisons with altered expression"]] * seq_data$output[["Supplementary Table 2"]][["Fraction of comparisons with up-regulated expression"]]) - 1, 
  size = seq_data$output[["Supplementary Table 2"]][["Total number of transcriptomic comparisons with altered expression"]], 
  prob = 0.5) 

seq_data$output[["Supplementary Table 2"]]$freq_down <- pbinom(
  as.integer( seq_data$output[["Supplementary Table 2"]][["Total number of transcriptomic comparisons with altered expression"]] * seq_data$output[["Supplementary Table 2"]][["Fraction of comparisons with up-regulated expression"]] ), 
  size = seq_data$output[["Supplementary Table 2"]][["Total number of transcriptomic comparisons with altered expression"]], 
  prob = 0.5) 

temp <- purrr::map2_dfc(
  .x = list(
    "signif" = seq_data$output[["Supplementary Table 2"]]$freq_signif, 
    "up" = seq_data$output[["Supplementary Table 2"]]$freq_up, 
    "down" = seq_data$output[["Supplementary Table 2"]]$freq_down),
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

seq_data$output[["Supplementary Table 2"]] <- cbind(seq_data$output[["Supplementary Table 2"]], temp)

openxlsx::write.xlsx(seq_data$output$`Supplementary Table 2`, file = "Supplementary Table 2 v3 seq only.xlsx")

### S2 ###
##########


seq_data$merge$seq <- dplyr::select(seq_data$output[["Supplementary Table 2"]], 1:4)

colnames(seq_data$merge$seq)[2:4] <- stringr::str_c(colnames(seq_data$merge$seq)[2:4], "_seq")

seq_data$merge$all <- dplyr::select(seq_data$input$`Supplementary Table 2`, 1:4, 7:11)

colnames(seq_data$merge$all) <- stringr::str_replace_all(colnames(seq_data$merge$all), pattern = "\\.", replacement = " ")

colnames(seq_data$merge$all)[2:4] <- stringr::str_c(colnames(seq_data$merge$all)[2:4], "_all")

seq_data$merge$merged <- merge(x = seq_data$merge$seq, y = seq_data$merge$all, by = "Gene symbol")

seq_data$merge$merged$trans_comps_ratio <- seq_data$merge$merged$`Total number of transcriptomic comparisons with altered expression_all` / seq_data$merge$merged$`Total number of transcriptomic comparisons with altered expression_seq`

seq_data$merge$merged$upreg_ratio <- seq_data$merge$merged$`Fraction of comparisons with up-regulated expression_all` / seq_data$merge$merged$`Fraction of comparisons with up-regulated expression_seq`

openxlsx::write.xlsx(seq_data$merge$merged, file = "Supplementary Table 2 v3 seq vs all data.xlsx")







temp <- subset(seq_data[["dataset1"]][["spread_medianed_final_good_dataset"]], subset = seq_data[["dataset1"]][["spread_medianed_final_good_dataset"]]$Standarised_gene_symbol == "hydin")

#
















##########################################
### GET NUMBERS OF GENES IN ALL PAPERS ###

paper_number_numbers <- seq_data[["exp_number_and_percentage"]][["pap"]] %>%
  dplyr::group_by(number) %>%
  dplyr::mutate(nb_of_genes_detected_in_this_nb_of_papers = dplyr::n()) %>%
  dplyr::mutate(percent_of_genes_detected_in_this_nb_of_papers = (dplyr::n()/length(seq_data[["exp_number_and_percentage"]][["pap"]][[1]])) ) %>%
  dplyr::select(-Standarised_gene_symbol) %>%
  unique()
readr::write_tsv(paper_number_numbers, file = 'exp_and_paper_numbers/paper_number_numbers.tsv')
### GET NUMBERS OF GENES IN ALL PAPERS ###
### GET NUMBERS OF GENES IN ALL PAPERS DRAWING ###
# library(ggplot2)
paper_number_numbers <- paper_number_numbers[order(paper_number_numbers$number),]

paper_number_numbers_for_plot <- paper_number_numbers
paper_number_numbers_for_plot$nb_of_genes_detected_in_this_nb_of_papers[8] <-sum(paper_number_numbers$nb_of_genes_detected_in_this_nb_of_papers[8:12])
paper_number_numbers_for_plot$percent_of_genes_detected_in_this_nb_of_papers[8] <-sum(paper_number_numbers$percent_of_genes_detected_in_this_nb_of_papers[8:12])
paper_number_numbers_for_plot <- paper_number_numbers_for_plot[-c(9:12),]

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
  scale_x_continuous(breaks = seq(1, 8, 1), labels = c(as.character(seq(1, 7, 1)), '8 <'))+
  labs(title = 'number of papers in which gene was detected')+
  xlab('number of papers')+
  ylab('number of genes')

ggplot2::ggsave("plot.png", plot = plot, device = "png")

### GET NUMBERS OF GENES IN ALL PAPERS ###
##########################################

























#################
### FUNCTIONS ###


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
