library(tidyverse)
library(magrittr)
library(synapser)
library(bigrquery)

source("../markdown_generation/query_functions.R")

synapser::synLogin()



#### input

# args <- commandArgs(trailingOnly=TRUE)
# version <- args[[1]]
# 
# if(version == "test"){
#     id_column <- "Test_report_project_id"
# } else {
#     id_column <- "Report_project_id"
# }

id_column <- "Report_project_id"

project_df <- "syn11612493" %>% 
    synapser::synGet() %>% 
    magrittr::use_series("path") %>% 
    readr::read_csv() %>% 
    dplyr::select("TEAM" = "Bird_alias", "parent" = id_column)

param_tbl <- BQ_DBI %>% 
    dplyr::tbl("Submissions") %>% 
    dplyr::filter(ROUND == "x") %>%
    dplyr::select(TEAM, SUBMISSION_ID, PATIENT_ID, AUPRC) %>%
    dplyr::arrange(dplyr::desc(SUBMISSION_ID)) %>% 
    dplyr::as_tibble() %>% 
    dplyr::left_join(project_df) %>% 
    dplyr::select(-TEAM) %>% 
    tidyr::nest(values = -parent) %>% 
    dplyr::mutate(name = "Round X results")


### output
synapse_table_objects <- purrr::pmap(param_tbl, synapser::synBuildTable)
purrr::map(synapse_table_objects, synapser::synStore)
