library(synapser)
library(tidyverse)
library(magrittr)
library(knit2synapse)
library(bigrquery)

synapser::synLogin()

source("upload_markdown_functions.R")

args <- commandArgs(trailingOnly=TRUE)
version <- args[[1]]

r2_markdown <- c(
    "round2_submissions_from_fastq.Rmd",
    "round2_validation_from_fastq1.Rmd",
    "round2_validation_from_fastq3.Rmd",  
    "round2_validation_from_vcf2.Rmd", 
    "round2_variant_counts.Rmd",
    "round2_submissions_from_vcf.Rmd",
    "round2_validation_from_fastq2.Rmd",  
    "round2_validation_from_vcf1.Rmd",
    "round2_validation_from_vcf3.Rmd",
    "round2_variant_overlap.Rmd"
)

if(version == "live") {
    id_column <- "Report_project_id"
} 
if(version == "test") {
    id_column <- "Test_report_project_id"
}

project_df <- "syn11612493" %>% 
    synapser::synGet() %>% 
    magrittr::use_series("path") %>% 
    readr::read_csv() %>% 
    dplyr::select(team = "Bird_alias", id_column) %>% 
    magrittr::set_colnames(c("team", "owner"))

submission_df <- 
    DBI::dbConnect(bigrquery::bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Submissions") %>% 
    dplyr::select(TEAM, ROUND) %>% 
    dplyr::distinct() %>% 
    dplyr::as_tibble() 

r2_teams <- submission_df %>% 
    dplyr::filter(ROUND == "2") %>% 
    magrittr::use_series(TEAM)


if(!all(r2_teams %in% project_df$team)) {
    stop("Not all teams have existing projects")
}

param_dfs <- project_df %>% 
    filter(team %in% r2_teams) %>% 
    merge(r2_markdown) %>% 
    dplyr::as_tibble() %>% 
    dplyr::rename(file = y) %>% 
    dplyr::mutate(file = as.character(file)) %>% 
    dplyr::mutate(wikiName = stringr::str_remove_all(file, ".Rmd")) %>% 
    dplyr::mutate(wikiName = basename(wikiName)) %>% 
    dplyr::mutate(wikiName = stringr::str_replace_all(wikiName, "_", " ")) %>% 
    dplyr::mutate(wikiName = stringr::str_to_title(wikiName)) %>% 
    dplyr::group_split(team) 

purrr::map(param_dfs, knit_markdown_by_team)
                     




    
