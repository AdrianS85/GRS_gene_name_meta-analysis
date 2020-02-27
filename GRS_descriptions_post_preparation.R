source('functions_for_genename_conversions.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')



####### PREPARING SUBSET TABLES FOR SUBSTRUCTURES ####### 
descriptions_1 <- readr::read_tsv("descriptions.txt", col_types = "nccccccccccccccccccccccc")

### PREPARE DATA WITH 10 AND 49 REMOVED AND descriptions_2 added ###
descriptions_2 <- readr::read_tsv("descriptions_2.txt", col_types = "nccccccccccccccccccccccc")

if( min(colnames(descriptions_1) == colnames(descriptions_2)) == 1)
{
  library(Hmisc)
  descriptions_1 <- subset(x = descriptions_1, subset = descriptions_1$Paper_ID %nin% c(10,49))
  
  descriptions_1_and_2 <- rbind(descriptions_1, descriptions_2)
  descriptions_1_and_2 <- descriptions_1_and_2[order(descriptions_1_and_2$Paper_ID),]
  readr::write_tsv(x = descriptions_1_and_2, path = 'descriptions_1_and_2.tsv')
}
### PREPARE DATA WITH 10 AND 49 REMOVED AND descriptions_2 added ###


descriptions_1_and_2 <- readr::read_tsv("descriptions_1_and_2.tsv", col_types = "nccccccccccccccccccccccc")
# load('at_least_in_3_papers_spread_med_fin_g_ds')
# at_least_in_3_papers_spread_med_fin_g_ds[is.na(at_least_in_3_papers_spread_med_fin_g_ds)] <- 0
# data_analysis_input <- tidyr::gather(data = at_least_in_3_papers_spread_med_fin_g_ds, key = 'Experiment', value = 'logFC', -lower_final_gene_name)



### Move badly placed descriptions ### 
descriptions_1_and_2 <- descriptions_1_and_2 %>% dplyr::mutate(Additional_group_info = ifelse(Group_ID == '23_2', "mice selectivly bred for high stress-induced analgesia", Additional_group_info))
descriptions_1_and_2 <- descriptions_1_and_2 %>% dplyr::mutate(Additional_group_info = ifelse(Group_ID == '23_1', "mice selectivly bred for low stress-induced analgesia", Additional_group_info))
descriptions_1_and_2$Stress_sensitivity[descriptions_1_and_2$Stress_sensitivity == 'Mice selectivly bred for high stress-induced analgesia'] <- NA
descriptions_1_and_2$Stress_sensitivity[descriptions_1_and_2$Stress_sensitivity == 'Mice selectivly bred for low stress-induced analgesia'] <- NA
### Move badly placed descriptions ### 



descriptions_1_and_2$Repetitions <- stringr::str_replace_all(string = descriptions_1_and_2$Repetitions, pattern = "4 weeks  3 days", replacement = '31 days')
descriptions_1_and_2$Repetitions_clean_days <- as.factor(cleanup_differing_units(charvec_ = descriptions_1_and_2$Repetitions, unit_identifier_regexList = list('.*d.*', '.*w.*'), unit_multiplier_list = list(1, 7))[[1]])
levels(descriptions_1_and_2$Repetitions_clean_days)
test <- data.frame(descriptions_1_and_2$Repetitions, descriptions_1_and_2$Repetitions_clean_days)



descriptions_1_and_2$Duration <- stringr::str_replace_all(string = descriptions_1_and_2$Duration, pattern = ",", replacement = '.')
descriptions_1_and_2$Duration_clean_minutes <- as.factor(cleanup_differing_units(charvec_ = descriptions_1_and_2$Duration, unit_identifier_regexList = list('.*m.*', '.*h.*', '.*w.*'), unit_multiplier_list = list(1, 60, 10080))[[1]])
descriptions_1_and_2$Duration_clean_minutes[as.character(descriptions_1_and_2$Duration) == '3 weeks / 9 weeks'] <- NA #66
levels(descriptions_1_and_2$Duration_clean_minutes)
test <- data.frame(descriptions_1_and_2$Duration, descriptions_1_and_2$Duration_clean_minutes)



descriptions_1_and_2$Gender_clean <- stringr::str_replace(string = descriptions_1_and_2$Gender, pattern = 'fem.*', replacement = 'female')
descriptions_1_and_2$Gender_clean <- stringr::str_replace(string = descriptions_1_and_2$Gender_clean, pattern = 'Male', replacement = 'male')
descriptions_1_and_2$Gender_clean <- stringr::str_replace(string = descriptions_1_and_2$Gender_clean, pattern = 'maleandfemale', replacement = 'male_and_female')
descriptions_1_and_2$Gender_clean <- as.factor(descriptions_1_and_2$Gender_clean)
levels(descriptions_1_and_2$Gender_clean)
test <- data.frame(descriptions_1_and_2$Gender, descriptions_1_and_2$Gender_clean)



descriptions_1_and_2$Brain_part_clean <- as.factor(brain_part_cleanup_wrapper(brain_parts_ = descriptions_1_and_2$Brain_part))
levels(descriptions_1_and_2$Brain_part_clean)
test <- data.frame(descriptions_1_and_2$Brain_part, descriptions_1_and_2$Brain_part_clean)


descriptions_1_and_2$Species <- as.factor(descriptions_1_and_2$Species)
levels(descriptions_1_and_2$Species)



descriptions_1_and_2$Stress_sensitivity_clean <-  as.factor(sensitivity_cleanup_wrapper(sensitivity_ = descriptions_1_and_2$Stress_sensitivity))
levels(descriptions_1_and_2$Stress_sensitivity_clean)
test <- data.frame(descriptions_1_and_2$Stress_sensitivity, descriptions_1_and_2$Stress_sensitivity_clean)



descriptions_1_and_2$Stress_clean <-  as.factor(stress_cleanup_wrapper(stress_ = descriptions_1_and_2$Stress))
levels(descriptions_1_and_2$Stress_clean)
test <- data.frame(descriptions_1_and_2$Stress, descriptions_1_and_2$Stress_clean)



descriptions_1_and_2$Stress_duration <-
  as.factor(dplyr::case_when(
    is.na(descriptions_1_and_2$Repetitions_clean_days) &
      !is.na(descriptions_1_and_2$Duration_clean_minutes) ~ 'acute',
    as.numeric(as.character(descriptions_1_and_2$Repetitions_clean_days)) <= 7 ~ 'medium',
    as.numeric(as.character(descriptions_1_and_2$Repetitions_clean_days)) > 7 ~ 'prolonged'
  ))
levels(descriptions_1_and_2$Stress_duration)
test <- data.frame(descriptions_1_and_2$Repetitions_clean_days, descriptions_1_and_2$Duration_clean_minutes,descriptions_1_and_2$Stress_duration)




descriptions_1_and_2$Measurement_latency_clean <- as.factor(latency_cleanup_wrapper(descriptions_1_and_2$Measurment_latency))
levels(descriptions_1_and_2$Measurement_latency_clean)
test <- data.frame(descriptions_1_and_2$Measurment_latency, descriptions_1_and_2$Measurement_latency_clean)



descriptions_1_and_2 <- unique(descriptions_1_and_2) ### For some reason some experiments are duplicated



save(descriptions_1_and_2, file = 'descriptions_1_and_2')