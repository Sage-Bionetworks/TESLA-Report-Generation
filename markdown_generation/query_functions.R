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

make_prediction_dbi_for_submission_plots <- function(round, src = "fastq"){
    dbi <- 
        make_prediction_dbi(round, src) %>% 
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



make_prediction_dbi_for_validation_plots <- function(round, src = "fastq"){
    dbi <- 
        make_prediction_dbi(round, src) %>% 
        dplyr::select(
            PATIENT_ID, 
            HLA_ALLELE, 
            ALT_EPI_SEQ,
            TEAM, 
            HLA_ALT_BINDING) %>% 
        dplyr::mutate(LOG_PREDICTED_BINDING = log10(HLA_ALT_BINDING + 1)) %>% 
        dplyr::select(-HLA_ALT_BINDING)
}


make_prediction_dbi <- function(round, src = "fastq"){
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
        dplyr::inner_join(submission_dbi)
}



