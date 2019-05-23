library(synapser)
library(tidyverse)
library(magrittr)
library(knit2synapse)
library(bigrquery)

synapser::synLogin()

args <- commandArgs(trailingOnly=TRUE)
version <- args[[1]]

r1_markdown <- c(
    "round1_submissions_from_fastq.Rmd",
    "round1_validation1.Rmd", 
    "round1_validation2.Rmd",  
    "round1_validation3.Rmd",
    "round1_variant_counts.Rmd", 
    "round1_variant_overlap.Rmd"
)

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

survey_markdown <- c(
    "survey_results.Rmd"
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

r1_teams <- submission_df %>% 
    dplyr::filter(ROUND == "1") %>% 
    magrittr::use_series(TEAM)

r2_teams <- submission_df %>% 
    dplyr::filter(ROUND == "2") %>% 
    magrittr::use_series(TEAM)

project_teams <- 
    c(r1_teams, r2_teams) %>% 
    unique

if(!all(project_teams %in% project_df$team)) {
    stop("Not all teams have existing projects")
}

survey_teams <- 
    DBI::dbConnect(bigrquery::bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Survey_Answers") %>% 
    dplyr::select(TEAM) %>% 
    dplyr::distinct() %>% 
    dplyr::as_tibble() %>% 
    magrittr::use_series(TEAM) %>% 
    purrr::discard(., ! . %in% project_teams)

r1_df <- project_df %>% 
    filter(team %in% r1_teams) %>% 
    merge(r1_markdown)

r2_df <- project_df %>% 
    filter(team %in% r2_teams) %>% 
    merge(r2_markdown)

survey_df <- project_df %>% 
    filter(team %in% survey_teams) %>% 
    merge(survey_markdown)

get_or_create_root_wiki_id <- function(project_id){
    wiki_id <- try(synGetWiki(project_id)$id)
    if(class(wiki_id) == "try-error"){
        wiki <- synapser::Wiki(title = "Report", owner = project_id, markdown = "")
        wiki <- synapser::synStore(wiki)
        wiki_id <- wiki$id
    }
    return(wiki_id)
}

create_config_yaml <- function(team_name){
    config_string <- stringr::str_c('team: "', team_name, '"')
    writeLines(config_string, "config.yaml")
}

clean_up_directory <- function(){
    file.remove("config.yaml")
    dir(pattern = "*_cache$") %>% 
        str_c("rm -rf ", .) %>% 
        purrr::map(system)
}

knit_markdown_by_team <- function(df){
    synapse_project_id <- df$owner[[1]]
    team_name <- df$team[[1]]
    wiki_id <- get_or_create_root_wiki_id(synapse_project_id)
    create_config_yaml(team_name)

    df %>% 
        dplyr::select(-team) %>% 
        dplyr::mutate(parentWikiId = wiki_id) %>% 
        purrr::pmap(knit2synapse::knitfile2synapse)
    
    clean_up_directory()
}

param_dfs <- 
    dplyr::bind_rows(r1_df, r2_df, survey_df) %>% 
    dplyr::as_tibble() %>% 
    dplyr::rename(file = y) %>% 
    dplyr::mutate(wikiName = stringr::str_remove_all(file, ".Rmd")) %>% 
    dplyr::mutate(wikiName = stringr::str_replace_all(wikiName, "_", " ")) %>% 
    dplyr::mutate(wikiName = stringr::str_to_title(wikiName)) %>% 
    dplyr::group_split(team) 

purrr::map(param_dfs[1], knit_markdown_by_team)
                     




    
