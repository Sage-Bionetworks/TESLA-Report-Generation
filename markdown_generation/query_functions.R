## query functions 

## These functions query Bigquery tables and either return a DBI::dbConnect connection, or a dplyr::tibble df

BQ_PROJECT <- "neoepitopes"
BQ_DATASET <- "Version_3"

BQ_DBI     <- DBI::dbConnect(
    bigrquery::bigquery(), 
    project = BQ_PROJECT,
    dataset = BQ_DATASET)

# general ----


make_submission_dbi <- function(rounds = c("1", "2", "x")){
    submission_dbi <- BQ_DBI %>% 
        dplyr::tbl("Submissions") %>% 
        dplyr::filter(ROUND %in% rounds) %>% 
        dplyr::select(SUBMISSION_ID, PATIENT_ID, TEAM)
}

make_prediction_dbi <- function(sources = c("fastq", "vcf"), ranked = T){
    prediction_dbi <- BQ_DBI %>% 
        dplyr::tbl("Predictions") %>% 
        dplyr::filter(SOURCE %in% sources) %>% 
        dplyr::arrange(RANK)
    if(ranked) prediction_dbi <- dplyr::filter(prediction_dbi, !is.na(RANK))
    return(prediction_dbi)
}

# submission ----

make_submission_plot_prediction_dbi <- function(
    rounds = c("1", "2", "x"),
    sources = c("fastq", "vcf"),
    ranked = T){

    submission_dbi <- make_submission_dbi(rounds)
    prediction_dbi <- make_prediction_dbi(sources, ranked) 
    
    dbi <- 
        dplyr::inner_join(prediction_dbi, submission_dbi) %>% 
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
    overlap_dbi <- prediction_dbi %>% 
        dplyr::filter(RANK <= 20) %>% 
        dplyr::mutate(PMHC = CONCAT(HLA_ALLELE, "_", ALT_EPI_SEQ)) %>% 
        dplyr::select(PATIENT_ID, TEAM,  PMHC) 
}

# binding validation ----

make_combined_binding_prediction_dbi <- function(round, src){
    submission_dbi <- make_submission_dbi(round)
    prediction_dbi <- make_prediction_dbi(src) %>% 
        dplyr::filter(!is.na(HLA_ALT_BINDING))
        
    combined_dbi <-
        dplyr::inner_join(prediction_dbi, submission_dbi) %>% 
        dplyr::mutate(LOG_BINDING = log10(HLA_ALT_BINDING + 1)) %>% 
        dplyr::select(
            PATIENT_ID, 
            HLA_ALLELE, 
            ALT_EPI_SEQ,
            TEAM,
            LOG_BINDING)
}

make_binding_validation_dbi <- function(){
    validation_dbi <- BQ_DBI %>%
        dplyr::tbl("Validated_Bindings") %>%
        dplyr::filter(!is.na(LJI_BINDING)) %>%
        dplyr::mutate(LOG_BINDING = log10(LJI_BINDING + 1))  %>%
        dplyr::select(ALT_EPI_SEQ, HLA_ALLELE, PATIENT_ID, LOG_BINDING) 
}

# validation 1 ----

make_binding_dotplot_df <- function(round, src, team){
    prediction_dbi <- make_combined_binding_prediction_dbi(round, src) 
    validation_dbi <- make_binding_validation_dbi() 
    pmhc_dbi       <- make_pmhc_dbi(prediction_dbi, validation_dbi, team)
    
    prediction_df <- prediction_dbi %>% 
        dplyr::inner_join(pmhc_dbi) %>% 
        dplyr::as_tibble() %>% 
        code_df_by_team(team) %>% 
        dplyr::select(-TEAM) %>% 
        dplyr::rename(TEAM = team_status)
    
    validation_df <- validation_dbi %>% 
        dplyr::inner_join(pmhc_dbi) %>% 
        dplyr::mutate(TEAM = "Measured") %>% 
        dplyr::as_tibble()
    
    plot_df <- 
        dplyr::bind_rows(prediction_df, validation_df) %>% 
        relevel_df
}

make_pmhc_dbi <- function(prediction_dbi, validation_dbi, team){
    pmhc_dbi <-  prediction_dbi %>% 
        dplyr::filter(TEAM == team) %>% 
        dplyr::inner_join(
            validation_dbi,
            by = c("PATIENT_ID", "HLA_ALLELE", "ALT_EPI_SEQ")) %>% 
        dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ) %>% 
        dplyr::distinct()
}


# validation 2

make_binding_scatterplot_df <- function(round, src, team){
    prediction_dbi <- make_combined_binding_prediction_dbi(round, src) %>% 
        dplyr::filter(TEAM == team) %>% 
        dplyr::select(-TEAM) %>% 
        dplyr::rename(LOG_PREDICTED_BINDING = LOG_BINDING)
    
    validation_dbi <- make_binding_validation_dbi() %>% 
        dplyr::rename(LOG_MEASURED_BINDING = LOG_BINDING)
    
    scatterplot_df <- 
        dplyr::inner_join(
            prediction_dbi,
            validation_dbi,
            by = c("HLA_ALLELE", "ALT_EPI_SEQ", "PATIENT_ID")) %>% 
        tibble::as_tibble()
}

# validation 3 ----

make_validation_df <- function(round, src, team){
    submission_dbi <- make_submission_dbi(round)
    
    combined_df <-
        make_prediction_dbi(src) %>% 
        dplyr::filter(RANK <= 20) %>%
        dplyr::inner_join(submission_dbi) %>%
        dplyr::select(PATIENT_ID, TEAM, HLA_ALLELE, ALT_EPI_SEQ, RANK) %>% 
        as_tibble()
    
    validated_bindings_df <- 
        make_validated_bindings_df() %>% 
        dplyr::inner_join(combined_df, by = c("PATIENT_ID", "HLA_ALLELE", "ALT_EPI_SEQ"))
        
    validated_epitopes_df <- 
        make_validated_epitopes_df() %>% 
        dplyr::inner_join(combined_df, by = c("PATIENT_ID", "ALT_EPI_SEQ"))
    
    df <- 
        dplyr::bind_rows(validated_epitopes_df, validated_bindings_df) %>% 
        dplyr::select(TEAM, PATIENT_ID, ASSAY, RESULT, RANK) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(ASSAY_NUM = ifelse(RESULT == "+", 1, 0)) %>%
        dplyr::group_by(TEAM, PATIENT_ID, ASSAY) %>%
        dplyr::summarise(COUNT = n(), MEAN_RANK = mean(RANK), RATE = mean(ASSAY_NUM)) %>%
        code_df_by_team(team)
}

make_validated_bindings_df <- function(){
    validated_bindings_df <- BQ_DBI %>% 
        dplyr::tbl("Validated_Bindings") %>% 
        dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ, TCR_NANOPARTICLE, TCR_FLOW_I, TCR_FLOW_II) %>% 
        dplyr::as_tibble() %>% 
        tidyr::gather(key = "ASSAY", value = "RESULT", TCR_NANOPARTICLE, TCR_FLOW_I, TCR_FLOW_II) %>% 
        tidyr::drop_na() 
}

make_validated_epitopes_df <- function(){
    validated_epitopes_df <- BQ_DBI %>% 
        dplyr::tbl("Validated_Epitopes") %>% 
        dplyr::select(PATIENT_ID, ALT_EPI_SEQ, TCELL_REACTIVITY) %>% 
        dplyr::as_tibble() %>% 
        tidyr::gather(key = "ASSAY", value = "RESULT", TCELL_REACTIVITY) %>% 
        tidyr::drop_na()
}



# others -----

make_survey_df <- function(team){
    question_dbi <- BQ_DBI %>% 
        dplyr::tbl("Survey_Questions") %>% 
        dplyr::select(-QUESTION_PART_1)
    
    survey_dbi <- BQ_DBI %>% 
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

