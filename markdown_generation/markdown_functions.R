relevel_df <- function(df){
    new_lvs <- df %>% 
        filter(TEAM == "Measured") %>% 
        arrange(LOG_BINDING) %>% 
        use_series(ALT_EPI_SEQ)
    df$ALT_EPI_SEQ = factor(df$ALT_EPI_SEQ, levels = new_lvs, ordered = TRUE)
    return(df)
}

code_df_by_team <- function(plot_data_df, team){
    mutate(plot_data_df, 
           team_status = ifelse(
               is.na(TEAM), 
               "No team", 
               ifelse(
                   TEAM == team, 
                   "Your team", 
                   "Other teams")))
}



make_dotplot <- function(allele){
    df <- allele_dfs[[allele]]
    if(!is.null(df)){
        make_validation_dotplot(df)
        text <- ""
    } else {
        text <- "None of your teams predicted neoepitopes for this allele were validated."
    }
    return(text)
}



make_scatterplot <- function(patient){
    df <- patient_dfs[[patient]]
    if(!is.null(df)){
        make_validation_scatterplot(df)
        text <- ""
    } else {
        text <- "None of your teams predicted neoepitopes for this patient were validated."
    }
    return(text)
}

create_p_common_matrix <- function(patient, df, max_rank = 20){
    patient_df <- df %>% 
        filter(PATIENT_ID == patient) %>% 
        dplyr::select(PMHC, TEAM)
    teams <- sort(unique(patient_df$TEAM))
    matrix <- matrix(
        0, 
        nrow = length(teams),
        ncol = length(teams), 
        dimnames = list(teams, teams))
    for(team1 in teams){
        for (team2 in teams){
            matrix[team1, team2] <- 
                find_percent_unique(patient_df, team1, team2, max_rank)
        }
    }
    return(matrix)
}

find_percent_unique <- function(df, team1, team2, max_rank){
    e1 <- df %>% 
        filter(TEAM == team1) %>% 
        use_series(PMHC)
    e2 <- df %>% 
        filter(TEAM == team2) %>% 
        use_series(PMHC)
    en <- min(max_rank, length(e1), length(e2))
    if(length(e1) > 0) e1 <- e1[1:en]
    if(length(e2) > 0) e2 <- e2[1:en]
    diff_n  <- length(setdiff(e1, e2))
    if(length(e1) == 0) return(0.0)
    result <- 1 - diff_n / length(e1)
    return(result)
}

create_median_table <- function(matrix){
    matrix %>%
        apply(2, median) %>%
        data.frame %>%
        rownames_to_column("team") %>%
        set_names(c("team", "median")) 
}



