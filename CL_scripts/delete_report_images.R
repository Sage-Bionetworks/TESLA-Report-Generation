library(tidyverse)
library(magrittr)
library(synapser)

args <- commandArgs(trailingOnly=TRUE)
version <- args[[1]]

if(version == "test"){
    id_column <- "Test_report_project_id"
} else {
    id_column <- "Report_project_id"
}


synapser::synLogin()

project_map_id <- "syn11612493"



project_ids <- project_map_id %>% 
    synapser::synGet() %>% 
    magrittr::use_series("path") %>% 
    readr::read_csv() %>% 
    magrittr::extract2(id_column)


get_file_ids_from_project <- function(id){
    id %>%  
        synapser::synGetChildren(includeTypes=list("file")) %>% 
        synapser::as.list() %>% 
        purrr::map_chr("id")
}

project_ids %>%
    purrr::map(get_file_ids_from_project) %>% 
    unlist() %>% 
    purrr::walk(synapser::synDelete)
    

