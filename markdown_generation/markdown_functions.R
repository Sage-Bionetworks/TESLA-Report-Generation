code_df_by_team <- function(df, team){
    code_df_by_column(
        df, 
        "TEAM", 
        team, 
        "team_status", 
        "Your team", 
        "No Team", 
        "Other teams")
}

code_df_by_column <- function(
    df, column, matching_values,
    new_column  = "New_column",
    new_value   = "matches", 
    na_value    = NA, 
    other_value = "doesn't match"){
    
    df %>% 
        dplyr::rename(TEMP_COL = column) %>% 
        dplyr::mutate(
            NEW_COL = ifelse(
                is.na(TEMP_COL),
                na_value,
                ifelse(
                    TEMP_COL %in% matching_values,
                    new_value,
                    other_value
                )
            )
        ) %>% 
        dplyr::rename(!!column := TEMP_COL) %>% 
        dplyr::rename(!!new_column := NEW_COL)
}

relevel_df <- function(df){
    new_lvs <- df %>% 
        dplyr::filter(TEAM == "Measured") %>% 
        dplyr::arrange(LOG_BINDING) %>% 
        magrittr::use_series(ALT_EPI_SEQ) %>% 
        unique()
    
    new_df <- dplyr::mutate(
        df, 
        ALT_EPI_SEQ = forcats::fct_relevel(
            ALT_EPI_SEQ, new_lvs))
    
    return(new_df)
}

relevel_df <- function(df, column){
    new_lvs <- df %>% 
        dplyr::filter(TEAM == "Measured") %>% 
        dplyr::arrange(LOG_BINDING) %>% 
        magrittr::use_series(ALT_EPI_SEQ) %>% 
        unique()
    
    new_df <- dplyr::mutate(
        df, 
        ALT_EPI_SEQ = forcats::fct_relevel(
            ALT_EPI_SEQ, new_lvs))
    
    return(new_df)
}

make_submission_plot_dfs <- function(round, source, team){
    prediction_dbi <- make_submission_plot_prediction_dbi(round, source)
    
    log_peptides_df <- make_log_peptides_df(prediction_dbi, team)
    peptide_length_df <- make_peptide_length_df(prediction_dbi, team)
    agretopicity_df <- make_agretopicity_df(prediction_dbi, team)
    overlap_df <- make_overlap_df(prediction_dbi)
    
    patients <- overlap_df %>%
        magrittr::use_series(PATIENT_ID) %>%
        unique %>% unique %>%
        sort
    
    epitope_overlap_df <-
        purrr::map(
            patients,
            create_p_common_matrix,
            overlap_df) %>%
        purrr::map(create_median_table) %>%
        purrr::map2(patients, ~inset(.x, "patient", value = .y)) %>%
        dplyr::bind_rows() %>%
        magrittr::set_names(c("TEAM", "SCORE", "PATIENT_ID")) %>%
        dplyr::mutate(PATIENT_ID = as.factor(PATIENT_ID)) %>%
        code_df_by_team(team)
    
    lst <- list(
        "log_peptides_df" = log_peptides_df,
        "peptide_length_df" = peptide_length_df,
        "agretopicity_df" = agretopicity_df,
        "epitope_overlap_df" = epitope_overlap_df
    )
    return(lst)
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


# functions written br Kristen

make_overlap_list <- function(synapse_id){
    data <- synapse_id %>% 
        synapser::synGet() %>% 
        magrittr::use_series(path) %>% 
        read.delim(header = FALSE)
    patMatrix = makePatientMatrix(data)
    nanVals = apply(patMatrix,MARGIN = 2,function(x){if (sum(x %in% NaN) == nrow(patMatrix)) return(1) else (return(0))})
    invalidVals = apply(patMatrix,MARGIN = 2,function(x){if (sum(x,na.rm = TRUE) == (nrow(patMatrix)-1)) return(1) else (return(0))})
    toRemove = c(which(nanVals == 1), which(invalidVals == 1))
    tempMatrix = patMatrix[-toRemove,-toRemove]
    patMatrix = tempMatrix
    rownames(patMatrix) = as.vector(sapply(rownames(patMatrix),function(x){unlist(strsplit(x = x,split = "_"))[1]}))
    colnames(patMatrix) = as.vector(sapply(colnames(patMatrix),function(x){unlist(strsplit(x = x,split = "_"))[1]}))
    aves = apply(patMatrix,MARGIN = 2,function(x){median(x,na.rm = TRUE)})
}

# Returns a matrix of overlap values for every two team combination.
# Columns are the values wrt team A.
makePatientMatrix=function(patientDataset){
    allbirds = as.list(unique(patientDataset[,1]))
    x = sapply(allbirds,function(x){ birdUnique(bird=x,patientDataset=patientDataset)})
    z = birdUnique(bird=allbirds[[1]],patientDataset = patientDataset)
    colnames(x) = names(z)
    return(x)
}

# Returns a vector of values for all pairs where A = bird.
# Value returned is setdiff(A,B) / size(A)
birdUnique=function(bird,patientDataset,type="records:"){
    smallSet = (patientDataset[patientDataset[,1] == bird & patientDataset[,3] == type,])
    ssDM = data.matrix(smallSet[,4:6])
    x = apply(ssDM,MARGIN = 1,function(x){ x[1] / (x[1] + x[3]) })
    names(x) = smallSet[,2]
    return(x)
}

getbird=function(x){as.vector(strsplit(x,split = "_")[[1]][1])}

# these are functions that use synapse data, but shoudl be turned into BQ queries

make_df_from_synapse_table_query <- function(query){
    result <- synapser::synTableQuery(query, includeRowIdAndRowVersion = F)
    result$asDataFrame() %>% 
        dplyr::as_tibble()
}

make_variant_counts_df <- function(round, team){
    synapser::synLogin()
    if(round == "1"){
        patient_ids <- c("1", "2", "3", "4")
    } else if (round == "2"){
        patient_ids <- c("10", "11", "12", "14", "15", "16")
    } else {
        stop("Data for round not available")
    }
    df <- "select team, patient, SNPs, MNPs, indels, others from syn11465705" %>% 
        make_df_from_synapse_table_query %>% 
        dplyr::rename(PATIENT_ID = patient) %>% 
        dplyr::rename(TEAM = team) %>% 
        dplyr::filter(PATIENT_ID %in% patient_ids) %>% 
        tidyr::drop_na() %>% 
        tidyr::gather(key = "VARIANT", value = "N", -c(TEAM, PATIENT_ID)) %>% 
        dplyr::mutate(LOG_N = log10(N + 1)) %>% 
        code_df_by_team(team) %>% 
        dplyr::select(-c(TEAM, N)) %>% 
        dplyr::mutate(PATIENT_ID = as.character(PATIENT_ID))
}

make_variant_overlap_df <- function(round, team){
    synapser::synLogin()
    if(round == "1"){
        synapse_ids <- c("syn11498546", "syn11498547", "syn11498548", "syn11498549") 
        patient_ids <- c("1", "2", "3", "4")
    } else if (round == "2"){
        synapse_ids <- c("syn11465124", "syn11465125", "syn11465126", "syn11465127", "syn11465128", "syn11465130")
        patient_ids <- c("10", "11", "12", "14", "15", "16")
    } else {
        stop("Data for round not available")
    }
    df <- synapse_ids %>% 
        purrr::map(make_overlap_list) %>% 
        purrr::map(tibble::enframe) %>% 
        purrr::map2(patient_ids, ~dplyr::mutate(.x, patient = .y)) %>% 
        dplyr::bind_rows()  %>%
        dplyr::rename(TEAM = name)  %>% 
        code_df_by_team(team)
    if(!team %in% df$TEAM){
        df <- dplyr::filter(df, TEAM == team)
    }
    return(df)
}

