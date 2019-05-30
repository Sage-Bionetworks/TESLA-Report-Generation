
# these functions return a list of teams ----

get_teams_from_submissions_dbi <- function(submission_dbi, teams = c("1", "2")){
    submission_dbi %>%
        dplyr::filter(ROUND %in% teams) %>%
        dplyr::as_tibble() %>% 
        magrittr::use_series(TEAM)
}

get_survey_teams <- function(project_teams, survey_dbi){
    survey_dbi %>%
        dplyr::filter(TEAM %in% project_teams) %>%
        dplyr::as_tibble() %>% 
        magrittr::use_series(TEAM) 
}


# ----

join_markdown_to_project_df <- function(project_df, markdown_df, teams, typ){
    project_df %>%
        filter(team %in% teams) %>%
        tidyr::crossing(markdown_df) %>%
        dplyr::filter(type == typ)
}



knit_markdown_by_group <- function(team, source, round, df, root_wiki, verbose = T){
    synapse_project_id <- df$owner[[1]]
    if(verbose) log_markdown_by_group(team, source, round, synapse_project_id)
    create_config_yaml(team, source, round)
    df %>%
        dplyr::select(file, owner, parentWikiId, wikiName) %>% 
        purrr::pmap(knit2synapse::knitfile2synapse)
    
    clean_up_directory()
    
}


get_or_create_root_wiki_id <- function(project_id){
    wiki_id <- try(synGetWiki(project_id)$id)
    if(class(wiki_id) == "try-error"){
        wiki <- synapser::Wiki(title = "Report", owner = project_id, markdown = "")
        wiki <- synapser::synStore(wiki)
        wiki_id <- wiki$id
    }
    return(wiki_id)
}

get_wiki_df_by_project <- function(project_id){
    wikis <- synapser::synGetWikiHeaders(project_id) 
    df <- 
        dplyr::tibble(
            "id" = as.character(map(wikis, "id")),
            "parent_id" = as.character(map(wikis, "parentId")),
            "title" = as.character(map(wikis, "title")),
            "project_id" = project_id
        ) %>% 
        mutate(parent_id = na_if(parent_id, "NULL"))
}

get_or_create_root_wiki_id <- function(project_id){
    wiki_id <- try(synGetWiki(project_id)$id)
    if(class(wiki_id) == "try-error"){
        wiki <- synapser::Wiki(title = "Report", owner = project_id, markdown = "")
        wiki <- synapser::synStore(wiki)
        wiki_id <- wiki$id
    }
    return(wiki_id)
}

log_markdown_by_group <- function(team, source, round, synapse_project_id){
    print("######")
    print(stringr::str_c("Project id: ", synapse_project_id))
    print(stringr::str_c("Team: ", team))
    print(stringr::str_c("Source: ", source))
    print(stringr::str_c("Round: ", round))
    print("######")
}

get_wiki_df_by_project <- function(project_id){
    wikis <- synapser::synGetWikiHeaders(project_id) 
    df <- 
        dplyr::tibble(
            "id" = as.character(map(wikis, "id")),
            "parent_id" = as.character(map(wikis, "parentId")),
            "title" = as.character(map(wikis, "title")),
            "project_id" = project_id
        ) %>% 
        mutate(parent_id = na_if(parent_id, "NULL"))
}

create_config_yaml <- function(team, source, round){
    lst <-  list(
        "team" = team, 
        "source" = source,
        "round" = round
    ) 
    yaml::write_yaml(lst, "config.yaml")
}

clean_up_directory <- function(){
    file.remove("config.yaml")
    dir(pattern = "*_cache$") %>% 
        str_c("rm -rf ", .) %>% 
        purrr::map(system)
}

delete_project_wiki <- function(project_id, wiki_id){
    res <-
        stringr::str_c("/entity/", project_id, "/wiki/", wiki_id) %>% 
        purrr::map(synapser::synRestDELETE)
}



