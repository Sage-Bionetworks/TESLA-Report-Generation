library(magrittr)

devtools::source_url("https://raw.githubusercontent.com/Sage-Bionetworks/synapse_tidy_utils/master/utils.R")
source("../markdown_generation/query_functions.R")

synapser::synLogin()

rounds <- c("1", "2", "X")
admin_ids <- c("1968150", "3346610", "3348278")


project_tbl <- "syn11612493" %>% 
    synapse_file_to_tbl(delim = ",") %>% 
    dplyr::arrange(desc(Report_alias))

insect_tbl <- query_synapse_table("select * from syn20782002")
bird_tbl   <- query_synapse_table("select * from syn8220615") 
translation_tbl <-
    dplyr::inner_join(insect_tbl, bird_tbl, by = "realTeam") %>% 
    dplyr::select(Insect_name = alias.x, Bird_name = alias.y)

bird_teams_with_submissions <- BQ_DBI %>% 
    dplyr::tbl("Submissions") %>% 
    dplyr::filter(ROUND %in% rounds) %>% 
    dplyr::select(TEAM) %>% 
    dplyr::distinct() %>% 
    dplyr::as_tibble()

round3_teams_with_submissions <- BQ_DBI %>% 
    dplyr::tbl("Submissions") %>% 
    dplyr::filter(ROUND %in% "3") %>% 
    dplyr::select(TEAM) %>% 
    dplyr::distinct() %>% 
    dplyr::as_tibble() %>% 
    dplyr::left_join(translation_tbl, by = c("TEAM" =  "Insect_name")) %>% 
    dplyr::select(TEAM = Bird_name)

teams_with_submissions <- 
    list(bird_teams_with_submissions, round3_teams_with_submissions) %>% 
    dplyr::bind_rows() %>% 
    dplyr::distinct()

teams_without_reports <- teams_with_submissions %>% 
    dplyr::filter(!TEAM %in% project_tbl$Bird_alias) %>% 
    dplyr::pull(TEAM)
  
new_project_tbl <- bird_tbl %>% 
    dplyr::rename(Team = realTeam, Bird_alias = alias) %>% 
    dplyr::filter(Bird_alias %in% teams_without_reports) %>% 
    dplyr::mutate(
        Report_number = 1:dplyr::n() + nrow(project_tbl),
        Report_alias = stringr::str_c("team", Report_number),
        Project_name = stringr::str_c("Team ", Report_number, " report"),
        Test_project_name = stringr::str_c("test team ", Report_number, " report")
    ) %>% 
    tidyr::pivot_longer(
        Project_name:Test_project_name, 
        names_to = "type", 
        values_to = "project_name"
    ) %>% 
    dplyr::mutate(id = purrr::map_chr(
        project_name, 
        ~synapser::synStore(synapser::Project(.x))$get("id")
    )) %>% 
    dplyr::select(-c(Report_number, project_name)) %>% 
    tidyr::pivot_wider(names_from = type, values_from = id) %>% 
    dplyr::rename(
        Report_project_id = Project_name, 
        Test_report_project_id = Test_project_name
    )

new_project_tbl %>% 
    tidyr::pivot_longer(Report_project_id:Test_report_project_id) %>% 
    dplyr::select(value) %>% 
    merge(admin_ids) %>% 
    dplyr::as_tibble() %>% 
    dplyr::rename(entity = value, principalId = y) %>% 
    dplyr::mutate(principalId = as.integer(as.character(principalId))) %>% 
    purrr::pmap(synapser::synSetPermissions)

new_project_tbl %>% 
    dplyr::select(principalId = Team, entity = Report_project_id) %>% 
    dplyr::mutate(principalId = purrr::map_int(
        principalId, 
        ~as.integer(synapser::synGetTeam(.x)$get("id"))
    )) %>% 
    purrr::pmap(synapser::synSetPermissions)

x <- 
    dplyr::bind_rows(project_tbl, new_project_tbl) %>% 
    readr::write_csv("team_project_map.csv")

ent <- 
    synapser::File("team_project_map.csv", "syn8123644") %>% 
    synapser::synStore()
