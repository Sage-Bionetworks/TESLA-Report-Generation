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

#testing code
# args = list(
#     "version" = "test",
#     "markdown_types" = "round2",
#     "teams" = "all",
#     "replace_behavior" = "replace")


library(synapser)
library(tidyverse)
library(magrittr)
library(knit2synapse)
library(bigrquery)
library(yaml)

synapser::synLogin()


source("upload_markdown_functions.R")
devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")

if(args$version == "live") {
    id_column <- "Report_project_id"
} 
if(args$version == "test") {
    id_column <- "Test_report_project_id"
}

if(args$markdown_types == "all"){
    markdown_types <- c("round1", "round2", "round3", "survey", "root")
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

insect_tbl <- query_synapse_table("select * from syn20782002")
bird_tbl   <- query_synapse_table("select * from syn8220615") 
translation_tbl <-
    dplyr::inner_join(insect_tbl, bird_tbl, by = "realTeam") %>% 
    dplyr::select(Insect_name = alias.x, Bird_name = alias.y)

r1_teams <- get_teams_from_submissions_dbi(submission_dbi, "1")
r2_teams <- get_teams_from_submissions_dbi(submission_dbi, "2")
rx_teams <- get_teams_from_submissions_dbi(submission_dbi, "x")
r3_teams <- get_teams_from_submissions_dbi(submission_dbi, "3") %>% 
    dplyr::as_tibble() %>% 
    dplyr::inner_join(translation_tbl, by = c("value" = "Insect_name")) %>% 
    dplyr::pull(Bird_name)
    
project_teams <- 
    c(r1_teams, r2_teams, r3_teams, rx_teams) %>% 
    unique()

survey_teams <- get_survey_teams(project_teams, survey_dbi)

if(!all(project_teams %in% project_df$team)) {
    stop("Not all teams have existing projects")
}

team_lists <- list(
    r1_teams,
    r2_teams,
    r3_teams,
    survey_teams,
    project_teams
)

types <- c("round1", "round2", "round3", "survey", "root")

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
        dplyr::inner_join(
            create_wiki_param_df, 
            by = c(
                "project_id" = "owner",
                "parent_id" = "parentWikiId", 
                "title" = "wikiName"
            )
        ) %>% 
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
    dplyr::filter(wikiName != "Report") 

nested_create_wiki_param_df3 <- nested_create_wiki_param_df %>% 
    dplyr::filter(round == "3") %>% 
    dplyr::inner_join(translation_tbl, by = c("team" = "Bird_name")) %>% 
    dplyr::select(-team) %>% 
    dplyr::rename(team = Insect_name)

nested_create_wiki_param_df <- nested_create_wiki_param_df %>% 
    dplyr::filter(round != "3") %>% 
    dplyr::bind_rows(nested_create_wiki_param_df3) %>% 
    tidyr::nest(tbl = -c(team, round, source)) 

nested_create_root_wiki_param_df <- create_wiki_param_df %>% 
    dplyr::filter(wikiName == "Report") %>% 
    dplyr::select(-parentWikiId) %>% 
    tidyr::nest(tbl = -c(team, round, source)) 


# output

if(args$replace_behavior == "replace" && nrow(delete_wiki_param_df) > 0){
    purrr::pmap(delete_wiki_param_df, delete_project_wiki)
}
purrr::pmap(nested_create_wiki_param_df, knit_markdown_by_group)
purrr::pmap(nested_create_root_wiki_param_df, knit_markdown_by_group)





    
