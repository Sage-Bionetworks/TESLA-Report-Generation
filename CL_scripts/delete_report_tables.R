library(tidyverse)
library(magrittr)
library(synapser)

synapser::synLogin()

#### input

args <- commandArgs(trailingOnly=TRUE)
version <- args[[1]]

if(version == "test"){
    id_column <- "Test_report_project_id"
} else {
    id_column <- "Report_project_id"
}

table_generators <- "syn11612493" %>% 
    synapser::synGet() %>% 
    magrittr::use_series("path") %>% 
    readr::read_csv() %>% 
    magrittr::extract2(id_column) %>% 
    purrr::map(synapser::synGetChildren, includeTypes = list("table"))
    

### data processing    
    
table_df <- table_generators %>% 
    purrr::map(synapser::as.list) %>% 
    purrr::flatten() %>% 
    purrr::map(as.data.frame) %>% 
    purrr::map(dplyr::as_tibble) %>% 
    dplyr::bind_rows()

table_ids <- table_df$id


### output

purrr::walk(table_ids, synapser::synDelete)
