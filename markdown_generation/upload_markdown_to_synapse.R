library(argparse)

parser = ArgumentParser()

parser$add_argument(
    "--version",
    type = "character",
    default = "test")

parser$add_argument(
    "--markdown_types",
    type = "character",
    default = "all")

parser$add_argument(
    "--teams",
    type = "character",
    default = "all")

parser$add_argument(
    "--replace_behavior",
    type = "character",
    default = "replace")

args = parser$parse_args()

# testing code
# args = list(
#     "version" = "test",
#     "markdown_types" = "all",
#     "teams" = "merlin",
#     "replace_behavior" = "replace")


library(synapser)
library(tidyverse)
library(magrittr)
library(knit2synapse)
library(bigrquery)
library(yaml)

synapser::synLogin()


source("upload_markdown_functions.R")

if(args$version == "live") {
    id_column <- "Report_project_id"
} 
if(args$version == "test") {
    id_column <- "Test_report_project_id"
}

if(args$markdown_types == "all"){
    markdown_types <- c("round1", "round2", "survey", "root")
} else {
    markdown_types <- args$markdown_types
}

# input

markdown_df <- readr::read_tsv("../markdown_generation/markdown_file_list.tsv")

project_df <- "syn11612493" %>% 
    synapser::synGet() %>% 
    magrittr::use_series("path") %>% 
    readr::read_csv() %>% 
    dplyr::select(team = "Bird_alias", id_column) %>% 
    magrittr::set_colnames(c("team", "owner")) %>% 
    dplyr::mutate(parentWikiId = purrr::map_chr(owner, get_or_create_root_wiki_id))

if(args$teams != "all"){
    project_df = dplyr::filter(project_df, team == args$teams)
}

submission_dbi <- 
    DBI::dbConnect(bigrquery::bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Submissions") %>% 
    dplyr::select(TEAM, ROUND) %>% 
    dplyr::distinct() 

survey_dbi <- 
    DBI::dbConnect(bigrquery::bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Survey_Answers") %>% 
    dplyr::select(TEAM) %>% 
    dplyr::distinct()


if(args$replace_behavior != "add"){
    wiki_df <- project_df %>% 
        magrittr::use_series(owner) %>% 
        purrr::map(get_wiki_df_by_project) %>% 
        dplyr::bind_rows()
}

# data processing ----

r1_teams <- get_teams_from_submissions_dbi(submission_dbi, "1")
r2_teams <- get_teams_from_submissions_dbi(submission_dbi, "2")
project_teams <- get_teams_from_submissions_dbi(submission_dbi)
survey_teams <- get_survey_teams(project_teams, survey_dbi)

if(!all(project_teams %in% project_df$team)) {
    stop("Not all teams have existing projects")
}

team_lists <- list(
    r1_teams,
    r2_teams,
    survey_teams,
    project_teams
)

types <- c("round1", "round2", "survey", "root")

create_wiki_param_df <- 
    purrr::map2(
        team_lists, 
        types, 
        ~join_markdown_to_project_df(project_df, markdown_df, .x, .y)
    ) %>% 
    dplyr::bind_rows() %>%
    dplyr::filter(type %in% markdown_types) %>%
    dplyr::select(-type)

if(args$replace_behavior == "replace"){
    delete_wiki_param_df <-  wiki_df %>% 
        dplyr::filter(!is.na(parent_id)) %>% 
        dplyr::left_join(
            create_wiki_param_df, 
            by = c("project_id" = "owner", "parent_id" = "parentWikiId", "title" = "wikiName")) %>% 
        dplyr::select(project_id, id) %>% 
        dplyr::rename(wiki_id = id)
}

if(args$replace_behavior == "keep"){
    create_wiki_param_df <- dplyr::anti_join(
        create_wiki_param_df,
        wiki_df,
        by = c("owner" = "project_id", "parentWikiId" = "parent_id", "wikiName" = "title"))
}



nested_create_wiki_param_df <- create_wiki_param_df %>% 
    dplyr::filter(wikiName != "Report") %>%
    tidyr::nest(
        -c(team, round, source), 
        .key = df)
    
nested_create_root_wiki_param_df <- create_wiki_param_df %>% 
    dplyr::filter(wikiName == "Report") %>% 
    dplyr::select(-parentWikiId) %>% 
    tidyr::nest(
        -c(team, round, source), 
        .key = df)

# output

if(args$replace_behavior == "replace" && nrow(delete_wiki_param_df) > 0){
    purrr::pmap(delete_wiki_param_df, delete_project_wiki)
}
purrr::pmap(nested_create_wiki_param_df, knit_markdown_by_group)
purrr::pmap(nested_create_root_wiki_param_df, knit_markdown_by_group)





    
