make_dbis_for_submission_plots <- function(round, src = "fastq"){
    prediction_dbi <- make_prediction_dbi_for_submission_plots(round, src)
    log_peptides_dbi <- make_log_peptides_dbi(prediction_dbi)
    peptide_length_dbi <- make_peptide_length_dbi(prediction_dbi)
    agretopicity_dbi <- make_agretopicity_dbi(prediction_dbi)
    overlap_dbi <- make_overlap_dbi(prediction_dbi)
    lst <- list(
        "log_peptides_dbi" = log_peptides_dbi,
        "peptide_length_dbi" = peptide_length_dbi,
        "agretopicity_dbi" = agretopicity_dbi,
        "overlap_dbi" = overlap_dbi
        )
    return(lst)
}

make_prediction_dbi_for_submission_plots <- function(round, src){
    submission_dbi <- 
        DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
        dplyr::tbl("Submissions") %>% 
        dplyr::filter(ROUND == round) %>% 
        dplyr::select(SUBMISSION_ID, PATIENT_ID, TEAM)
    
    prediction_dbi <-
        DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
        dplyr::tbl("Predictions") %>% 
        dplyr::filter(!is.na(RANK)) %>% 
        dplyr::filter(!is.na(HLA_ALT_BINDING)) %>% 
        dplyr::filter(SOURCE == src) %>% 
        dplyr::inner_join(submission_dbi) %>% 
        dplyr::select(
            PATIENT_ID, 
            HLA_ALLELE, 
            ALT_EPI_SEQ,
            TEAM, 
            HLA_ALT_BINDING, 
            HLA_REF_BINDING, 
            PEP_LEN,
            RANK) 
}

make_log_peptides_dbi <- function(prediction_dbi){
    prediction_dbi %>% 
        dplyr::group_by(TEAM, PATIENT_ID) %>% 
        dplyr::summarise(COUNT = n()) %>% 
        dplyr::mutate(LOG_COUNT = log10(COUNT)) %>% 
        dplyr::select(TEAM, PATIENT_ID, LOG_COUNT)
}

make_peptide_length_dbi <- function(prediction_dbi, max_rank = 20){
    prediction_dbi %>% 
        dplyr::filter(RANK <= max_rank) %>% 
        dplyr::select(TEAM, PEP_LEN)
}

make_agretopicity_dbi <- function(prediction_dbi, max_rank = 20){
    prediction_dbi %>% 
        dplyr::filter(RANK <= max_rank) %>% 
        dplyr::filter(!is.na(HLA_ALT_BINDING)) %>% 
        dplyr::filter(!is.na(HLA_REF_BINDING)) %>% 
        dplyr::mutate(LOG_AGRETOPICITY = log10(HLA_ALT_BINDING /HLA_REF_BINDING)) %>% 
        dplyr::select(TEAM, LOG_AGRETOPICITY)
}

make_overlap_dbi <- function(prediction_dbi, max_rank = 20){
    prediction_dbi %>% 
        dplyr::filter(RANK <= max_rank) %>% 
        dplyr::arrange(RANK) %>% 
        dplyr::select(RANK, TEAM, PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ)
}



make_binding_prediction_dbi <- function(round, src){
    submission_dbi <- 
        DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
        dplyr::tbl("Submissions") %>% 
        dplyr::filter(ROUND == round) %>% 
        dplyr::select(SUBMISSION_ID, PATIENT_ID, TEAM)
    
    prediction_dbi <-
        DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
        dplyr::tbl("Predictions") %>% 
        dplyr::filter(!is.na(RANK)) %>% 
        dplyr::filter(!is.na(HLA_ALT_BINDING)) %>% 
        dplyr::filter(SOURCE == src) %>% 
        dplyr::inner_join(submission_dbi) %>% 
        dplyr::mutate(LOG_PREDICTED_BINDING = log10(HLA_ALT_BINDING + 1)) %>% 
        dplyr::select(
            PATIENT_ID, 
            HLA_ALLELE, 
            ALT_EPI_SEQ,
            TEAM,
            LOG_PREDICTED_BINDING)
}

make_variant_counts_df <- function(patients, team){
    df <- "select team, patient, SNPs, MNPs, indels, others from syn11465705" %>% 
        synapser::synTableQuery() %>% 
        as.data.frame %>% 
        dplyr::as_tibble() %>%
        dplyr::rename(PATIENT_ID = patient) %>% 
        dplyr::rename(TEAM = team) %>% 
        dplyr::filter(PATIENT_ID %in% patients) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(-c(ROW_ID, ROW_VERSION)) %>% 
        tidyr::gather(key = "VARIANT", value = "N", -c(TEAM, PATIENT_ID)) %>% 
        dplyr::mutate(LOG_N = log10(N + 1)) %>% 
        code_df_by_team(team) %>% 
        dplyr::select(-c(TEAM, N)) %>% 
        dplyr::mutate(PATIENT_ID = as.character(PATIENT_ID))
}

make_variant_overlap_df <- function(round, team){
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

make_survey_df <- function(team){
    question_dbi <- 
        DBI::dbConnect(bigrquery::bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
        dplyr::tbl("Survey_Questions") %>% 
        dplyr::select(-QUESTION_PART_1)
    
    survey_dbi <-
        DBI::dbConnect(bigrquery::bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
        dplyr::tbl("Survey_Answers") %>% 
        dplyr::inner_join(question_dbi)
    
    survey_summary_df <- survey_dbi %>% 
        dplyr::as_tibble() %>% 
        dplyr::group_by(QUESTION_NAME) %>% 
        dplyr::summarise(fraction = sum(ANSWER) / n())
    
    plot_df <- survey_dbi %>% 
        dplyr::filter(TEAM == team) %>% 
        dplyr::select(-TEAM) %>% 
        dplyr::as_tibble() %>% 
        dplyr::left_join(survey_summary_df)
}

