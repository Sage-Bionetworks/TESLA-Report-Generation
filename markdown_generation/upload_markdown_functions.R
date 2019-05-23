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

knit_markdown_by_team <- function(df, verbose = T){
    synapse_project_id <- df$owner[[1]]
    team_name <- df$team[[1]]
    if(verbose){
        print(stringr::str_c("Project id: ", synapse_project_id))
        print(stringr::str_c("Team: ", team_name))
    }

    wiki_id <- get_or_create_root_wiki_id(synapse_project_id)
    create_config_yaml(team_name)
    
    df %>% 
        dplyr::select(-team) %>% 
        dplyr::mutate(parentWikiId = wiki_id) %>% 
        purrr::pmap(knit2synapse::knitfile2synapse)
    
    clean_up_directory()
}