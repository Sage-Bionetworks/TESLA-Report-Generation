library(tidyverse)
library(magrittr)
library(synapser)
library(bigrquery)

source("../markdown_generation/query_functions.R")

synapser::synLogin()



#### input

args <- commandArgs(trailingOnly=TRUE)
version <- args[[1]]

if(version == "test"){
    id_column <- "Test_report_project_id"
} else {
    id_column <- "Report_project_id"
}

rounds <- c("1", "2")
sources <- c("fastq", "vcf")

project_df <- "syn11612493" %>% 
    synapser::synGet() %>% 
    magrittr::use_series("path") %>% 
    readr::read_csv() %>% 
    dplyr::select("TEAM" = "Bird_alias", "parent" = id_column)

prediction_dbi <- 
    make_combined_prediction_dbi(rounds, sources, ranked = F) %>% 
    dplyr::select(-c(PEP_LEN, HLA_REF_BINDING))

lji_binding_dbi <- 
    make_lji_binding_dbi() %>% 
    dplyr::inner_join(prediction_dbi) %>% 
    dplyr::rename(PREDICTED_BINDING = HLA_ALT_BINDING) %>% 
    dplyr::rename(MEASURED_BINDING = LJI_BINDING) %>% 
    dplyr::rename(EPITOPE = ALT_EPI_SEQ) 


validated_bindings_dbi <-
    make_binding_assay_dbi() %>% 
    dplyr::inner_join(
        prediction_dbi,
        by = c("PATIENT_ID", "HLA_ALLELE", "ALT_EPI_SEQ")
    ) 
    
validated_epitopes_dbi <-
    make_epitope_assay_dbi() %>% 
    dplyr::inner_join(
        prediction_dbi,
        by = c("PATIENT_ID", "ALT_EPI_SEQ")
    ) %>% 
    dplyr::select(-HLA_ALLELE) %>% 
    dplyr::distinct()

assay_dbi <-
    dplyr::union_all(validated_bindings_dbi, validated_epitopes_dbi) %>% 
    dplyr::rename(EPITOPE = ALT_EPI_SEQ) 
    
    
### data processing   

assay_param_df <- assay_dbi %>% 
    dplyr::as_tibble() %>% 
    dplyr::group_by(SOURCE, PATIENT_ID, HLA_ALLELE, EPITOPE, ASSAY) %>% 
    dplyr::mutate(NUM_PREDICTED_EPITOPES = dplyr::n()) %>% 
    dplyr::mutate(MEDIAN_RANK = median(RANK, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(
        TEAM,
        PATIENT_ID,
        HLA_ALLELE,
        EPITOPE,
        SOURCE,
        ASSAY,
        RESULT,
        RANK,
        NUM_PREDICTED_EPITOPES,
        MEDIAN_RANK) %>% 
    tidyr::nest(-TEAM, .key = values) %>% 
    dplyr::mutate(name = "Validation_assays")
    


binding_df_all <- lji_binding_dbi %>% 
    dplyr::group_by(SOURCE, PATIENT_ID, HLA_ALLELE, EPITOPE) %>% 
    dplyr::mutate(NUM_PREDICTED_EPITOPES = dplyr::n()) %>% 
    dplyr::mutate(STDEV_PREDICTED_BINDING = sd(PREDICTED_BINDING, na.rm = T)) %>% 
    dplyr::mutate(MEAN_PREDICTED_BINDING = mean(PREDICTED_BINDING, na.rm = T)) %>% 
    as_tibble() %>% 
    dplyr::mutate(STDEV_PREDICTED_BINDING = round(STDEV_PREDICTED_BINDING, 2)) %>% 
    dplyr::mutate(MEAN_PREDICTED_BINDING = round(MEAN_PREDICTED_BINDING, 2)) 

binding_df_top_20 <- lji_binding_dbi %>% 
    dplyr::filter(RANK <= 20) %>% 
    dplyr::group_by(PATIENT_ID, HLA_ALLELE, EPITOPE) %>% 
    dplyr::mutate(TOP_20_NUM_PREDICTED_EPITOPES = dplyr::n()) %>% 
    dplyr::mutate(TOP_20_STDEV_PREDICTED_BINDING = sd(PREDICTED_BINDING, na.rm = T)) %>% 
    dplyr::mutate(TOP_20_MEAN_PREDICTED_BINDING = mean(PREDICTED_BINDING, na.rm = T)) %>% 
    as_tibble() %>% 
    dplyr::mutate(TOP_20_STDEV_PREDICTED_BINDING = round(TOP_20_STDEV_PREDICTED_BINDING, 2)) %>% 
    dplyr::mutate(TOP_20_MEAN_PREDICTED_BINDING = round(TOP_20_MEAN_PREDICTED_BINDING, 2)) 

binding_param_df <- 
    dplyr::full_join(binding_df_all, binding_df_top_20) %>% 
    dplyr::select(
        TEAM, 
        PATIENT_ID, 
        EPITOPE, 
        SOURCE,
        HLA_ALLELE,
        RANK,
        PREDICTED_BINDING,
        MEASURED_BINDING,
        NUM_PREDICTED_EPITOPES,
        STDEV_PREDICTED_BINDING,
        MEAN_PREDICTED_BINDING,
        TOP_20_NUM_PREDICTED_EPITOPES,
        TOP_20_STDEV_PREDICTED_BINDING,
        TOP_20_MEAN_PREDICTED_BINDING
    ) %>% 
    tidyr::nest(-TEAM, .key = values) %>% 
    dplyr::mutate(name = "MHC_binding_assay")
   
param_df <- 
    dplyr::bind_rows(assay_param_df, binding_param_df) %>% 
    dplyr::left_join(project_df) %>% 
    dplyr::select(-TEAM) 
    
    
### output
synapse_table_objects <- purrr::pmap(param_df, synapser::synBuildTable)
purrr::map(synapse_table_objects, synapser::synStore)
