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
    dplyr::select("TEAM" = "Bird_alias", "Project_ID" = id_column)

prediction_dbi <- 
    make_submission_plot_prediction_dbi(rounds, sources, ranked = F) %>% 
    dplyr::select(-PEP_LEN) 

validated_bindings_dbi <- BQ_DBI %>% 
    dplyr::tbl("Validated_Bindings") %>% 
    dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ, TCR_NANOPARTICLE, TCR_FLOW_I, TCR_FLOW_II, LJI_BINDING) %>% 
    dplyr::inner_join(prediction_dbi, by = c("PATIENT_ID", "HLA_ALLELE", "ALT_EPI_SEQ"))

validated_epitopes_dbi <- BQ_DBI %>% 
    dplyr::tbl("Validated_Epitopes") %>% 
    dplyr::select(PATIENT_ID, ALT_EPI_SEQ, TCELL_REACTIVITY) %>% 
    dplyr::inner_join(prediction_dbi, by = c("PATIENT_ID", "ALT_EPI_SEQ")) %>% 
    dplyr::select(-HLA_ALLELE) %>% 
    dplyr::distinct()

binding_dbi <- validated_bindings_dbi %>% 
    dplyr::select(-c(TCR_NANOPARTICLE, TCR_FLOW_I, TCR_FLOW_II)) %>% 
    dplyr::filter(!is.na(LJI_BINDING)) %>% 
    dplyr::distinct()

### data processing   

assay_df <- 
    dplyr::bind_rows(
        dplyr::as_tibble(validated_bindings_dbi),
        dplyr::as_tibble(validated_epitopes_dbi)) %>% 
    dplyr::select(
        PATIENT_ID, 
        HLA_ALLELE, 
        ALT_EPI_SEQ, 
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
    dplyr::distinct()

x <- assay_df %>% 
    dplyr::group_by(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ) %>% 
    dplyr::mutate(N_PREDICTED_EPITOPES = dplyr::n()) %>% 
    dplyr::mutate(MEDIAN_RANK = median(RANK))

binding_df <- binding_dbi %>% 
    as_tibble() %>% 
    dplyr::group_by(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ) %>% 
    dplyr::mutate(N_PREDICTED_EPITOPES = dplyr::n())

binding_score_df 
   
    
    
### output


