library(tidyverse)
library(magrittr)
library(synapser)
library(bigrquery)

source("../markdown_generation/query_functions.R")

synapser::synLogin()



#### input

# args <- commandArgs(trailingOnly=TRUE)
# version <- args[[1]]

version <- "test"

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
    make_submission_plot_prediction_dbi(rounds, sources, ranked = F) %>% 
    dplyr::select(-c(PEP_LEN, HLA_REF_BINDING))

validated_bindings_dbi <- BQ_DBI %>% 
    dplyr::tbl("Validated_Bindings") %>% 
    dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ, TCR_NANOPARTICLE, TCR_FLOW_I, TCR_FLOW_II, LJI_BINDING) %>% 
    dplyr::inner_join(prediction_dbi, by = c("PATIENT_ID", "HLA_ALLELE", "ALT_EPI_SEQ")) %>% 
    dplyr::rename(EPITOPE = ALT_EPI_SEQ) 

validated_epitopes_dbi <- BQ_DBI %>% 
    dplyr::tbl("Validated_Epitopes") %>% 
    dplyr::select(PATIENT_ID, ALT_EPI_SEQ, TCELL_REACTIVITY) %>% 
    dplyr::inner_join(prediction_dbi, by = c("PATIENT_ID", "ALT_EPI_SEQ")) %>% 
    dplyr::select(-HLA_ALLELE) %>% 
    dplyr::distinct()

binding_dbi <- validated_bindings_dbi %>% 
    dplyr::select(-c(TCR_NANOPARTICLE, TCR_FLOW_I, TCR_FLOW_II)) %>% 
    dplyr::filter(!is.na(LJI_BINDING)) %>% 
    dplyr::distinct() %>% 
    dplyr::rename(PREDICTED_BINDING = HLA_ALT_BINDING) %>% 
    dplyr::rename(MEASURED_BINDING = LJI_BINDING)
    

### data processing   

assay_df <- 
    dplyr::bind_rows(
        dplyr::as_tibble(validated_bindings_dbi),
        dplyr::as_tibble(validated_epitopes_dbi)) %>% 
    dplyr::select(
        PATIENT_ID, 
        HLA_ALLELE, 
        EPITOPE, 
        TCR_NANOPARTICLE, 
        TCR_FLOW_I, 
        TCR_FLOW_II,
        TCELL_REACTIVITY, 
        TEAM, 
        RANK) %>% 
    tidyr::gather(
        key = "ASSAY", 
        value = "RESULT", 
        TCR_NANOPARTICLE, 
        TCR_FLOW_I, 
        TCR_FLOW_II, 
        TCELL_REACTIVITY) %>% 
    dplyr::filter(!is.na(RESULT)) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(PATIENT_ID, HLA_ALLELE, EPITOPE, ASSAY) %>% 
    dplyr::mutate(NUM_PREDICTED_EPITOPES = dplyr::n()) %>% 
    dplyr::mutate(MEDIAN_RANK = median(RANK, na.rm = T)) %>% 
    dplyr::ungroup()

assay_param_df <- assay_df %>% 
    tidyr::nest(-TEAM, .key = values) %>% 
    dplyr::mutate(name = "Validation_assays")
    


binding_df_all <- binding_dbi %>% 
    dplyr::group_by(PATIENT_ID, HLA_ALLELE, EPITOPE) %>% 
    dplyr::mutate(NUM_PREDICTED_EPITOPES = dplyr::n()) %>% 
    dplyr::mutate(STDEV_PREDICTED_BINDING = sd(PREDICTED_BINDING, na.rm = T)) %>% 
    dplyr::mutate(MEAN_PREDICTED_BINDING = mean(PREDICTED_BINDING, na.rm = T)) %>% 
    as_tibble() %>% 
    dplyr::mutate(STDEV_PREDICTED_BINDING = round(STDEV_PREDICTED_BINDING, 2)) %>% 
    dplyr::mutate(MEAN_PREDICTED_BINDING = round(MEAN_PREDICTED_BINDING, 2)) 

binding_df_top_20 <- binding_dbi %>% 
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
