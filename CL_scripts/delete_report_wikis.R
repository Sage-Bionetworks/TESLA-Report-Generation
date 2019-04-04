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

delete_wikis_from_project <- function(id){
    wiki_ids <- try({
        id %>% 
            synapser::synGetWikiHeaders() %>% 
            purrr::map("id") 
    })
    if(!inherits(wiki_ids, "try-error")){
        wiki_ids %>% 
            stringr::str_c("/entity/", id, "/wiki/", .) %>% 
            purrr::map(synapser::synRestDELETE)
    }
}

purrr::map(project_ids, delete_wikis_from_project) 
    

