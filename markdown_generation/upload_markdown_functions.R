knit_markdown_by_group <- function(team, src, round, df, verbose = T){
    synapse_project_id <- df$owner[[1]]
    if(verbose) log_markdown_by_group(team, source, round, synapse_project_id)
    wiki_id <- get_or_create_root_wiki_id(synapse_project_id)
    wiki_df <- get_wiki_df_by_project(synapse_project_id)
    create_config_yaml(team, src, round)

    df %>%
        dplyr::mutate(parentWikiId = wiki_id) %>%
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

