# library(argparse)
# 
# parser = ArgumentParser()
# 
# parser$add_argument(
#     "--version",
#     type = "character",
#     default = "test")
# 
# parser$add_argument(
#     "--markdown_types",
#     type = "character",
#     default = "all")
#
# parser$add_argument(
#     "--teams",
#     type = "character",
#     default = "all")
#
# parser$add_argument(
#     "replace_behavior",
#     type = "character",
#     default = "replace")
# 
# args = parser$parse_args()

args = list(
    "version" = "test", 
    "markdown_types" = "all",
    "teams" = "merlin",
    "replace_behavior" = "add")


library(synapser)
library(tidyverse)
library(magrittr)
library(knit2synapse)
library(bigrquery)
library(yaml)

synapser::synLogin()

source("upload_markdown_functions.R")

wiki_df <- readr::read_tsv("../markdown_generation/markdown_file_list.tsv")

if(args$version == "live") {
    id_column <- "Report_project_id"
} 
if(args$version == "test") {
    id_column <- "Test_report_project_id"
}

if(args$markdown_types == "all"){
    markdown_types <- c("round1", "round2", "survey")
} else {
    markdown_types <- args$markdown_types
}

project_df <- "syn11612493" %>% 
    synapser::synGet() %>% 
    magrittr::use_series("path") %>% 
    readr::read_csv() %>% 
    dplyr::select(team = "Bird_alias", id_column) %>% 
    magrittr::set_colnames(c("team", "owner")) %>% 
    dplyr::mutate(root_wiki = purrr::map(owner, get_or_create_root_wiki_id))

submission_df <- 
    DBI::dbConnect(bigrquery::bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Submissions") %>% 
    dplyr::select(TEAM, ROUND) %>% 
    dplyr::distinct() %>% 
    dplyr::as_tibble() 

if(args$teams != "all"){
    project_df = dplyr::filter(project_df, team == args$teams)
    submission_df = dplyr::filter(submission_df, TEAM == args$teams)
}

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
    merge(wiki_df) %>% 
    dplyr::as_tibble() %>% 
    dplyr::filter(type == "round1") 

r2_df <- project_df %>% 
    filter(team %in% r2_teams) %>% 
    merge(wiki_df) %>% 
    dplyr::as_tibble() %>% 
    dplyr::filter(type == "round2") 

survey_df <- project_df %>% 
    filter(team %in% survey_teams) %>% 
    merge(wiki_df) %>% 
    dplyr::as_tibble() %>% 
    dplyr::filter(type == "survey")

param_df <- 
    dplyr::bind_rows(r1_df, r2_df, survey_df) %>% 
    dplyr::filter(type %in% markdown_types) %>% 
    dplyr::select(-type) %>%
    tidyr::nest(-c(team, round, source), .key = df)

# purrr::pmap(param_df, knit_markdown_by_group)
                     




    
