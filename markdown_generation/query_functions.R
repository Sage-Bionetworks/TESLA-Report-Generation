require(magrittr)
require(bigrquery)
require(dplyr)
require(purrr)

## query functions 

## These functions query Bigquery tables and either return a DBI::dbConnect connection, or a dplyr::tibble df

BQ_PROJECT <- "neoepitopes"
BQ_DATASET <- "Version_3"

BQ_DBI     <- DBI::dbConnect(
    bigrquery::bigquery(), 
    project = BQ_PROJECT,
    dataset = BQ_DATASET)

# general ----

make_combined_prediction_dbi <- function(
    round = c("1", "2", "x"),
    source = c("fastq", "vcf"),
    ranked = T,
    max_rank = F,
    hla_binding = F
){
    
    submission_dbi <- make_submission_dbi(round)
    prediction_dbi <- make_prediction_dbi(source, ranked, max_rank) 
    
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
            RANK,
            SOURCE
        ) 
    if(hla_binding){
        dbi <- dbi %>% 
            dplyr::filter(!is.na(HLA_ALT_BINDING)) %>% 
            dplyr::mutate(LOG_BINDING = log10(HLA_ALT_BINDING + 1))
    } 
    return(dbi)
    
}

make_submission_dbi <- function(round = c("1", "2", "x")){
    submission_dbi <- BQ_DBI %>% 
        dplyr::tbl("Submissions") %>% 
        dplyr::filter(ROUND %in% round) %>% 
        dplyr::select(SUBMISSION_ID, PATIENT_ID, TEAM)
}

make_prediction_dbi <- function(
    source = c("fastq", "vcf"), 
    ranked = T,
    max_rank = F
){
    dbi <- BQ_DBI %>% 
        dplyr::tbl("Predictions") %>% 
        dplyr::filter(SOURCE %in% source) %>% 
        dplyr::arrange(RANK)
    if(max_rank) dbi <- dplyr::filter(dbi, RANK <= max_rank)
    if(ranked && !max_rank) dbi <- dplyr::filter(dbi, !is.na(RANK))
    return(dbi)
}

make_lji_binding_dbi <- function(){
    BQ_DBI %>%
        dplyr::tbl("Validated_Bindings") %>%
        dplyr::filter(!is.na(LJI_BINDING)) %>%
        dplyr::mutate(LOG_BINDING = log10(LJI_BINDING + 1))  %>%
        dplyr::select(ALT_EPI_SEQ, HLA_ALLELE, LJI_BINDING, LOG_BINDING) %>% 
        dplyr::distinct()
}


make_binding_assay_dbi <- function(){
    dbi <- 
        c("TCR_NANOPARTICLE", "TCR_FLOW_I", "TCR_FLOW_II") %>% 
        purrr::map(make_binding_assay_dbi_by_col) %>% 
        purrr::reduce(dplyr::union_all) %>% 
        dplyr::filter(!is.na(RESULT)) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(RESULT = dplyr::if_else(RESULT == "+", T, F)) %>% 
        dplyr::group_by_at(vars(-"RESULT")) %>% 
        dplyr::summarise(RESULT = any(RESULT)) %>% 
        dplyr::ungroup()
}

make_binding_assay_dbi_by_col <- function(col){
    BQ_DBI %>% 
        dplyr::tbl("Validated_Bindings") %>% 
        dplyr::select(
            "PATIENT_ID", 
            "HLA_ALLELE", 
            "ALT_EPI_SEQ",
            "RESULT" = col
        ) %>% 
        dplyr::mutate(ASSAY = col) 
}

make_epitope_assay_dbi <- function(){
    BQ_DBI %>% 
        dplyr::tbl("Validated_Epitopes") %>% 
        dplyr::select(
            "PATIENT_ID", 
            "ALT_EPI_SEQ",
            "RESULT" = "TCELL_REACTIVITY"
        ) %>% 
        dplyr::mutate(ASSAY = "TCELL_REACTIVITY") %>%
        dplyr::filter(!is.na(RESULT)) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(RESULT = dplyr::if_else(RESULT == "+", T, F)) %>% 
        dplyr::group_by_at(vars(-"RESULT")) %>% 
        dplyr::summarise(RESULT = any(RESULT)) %>% 
        dplyr::ungroup()
}




# submission ----


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



# validation 1 ----

make_binding_dotplot_df <- function(round, source, team){
    prediction_dbi <- make_combined_prediction_dbi(round, source, hla_binding = T)
    validation_dbi <- make_lji_binding_dbi()
    pmhc_dbi       <- make_pmhc_dbi(prediction_dbi, validation_dbi, team)
    
    prediction_df <- prediction_dbi %>% 
        dplyr::select(
            PATIENT_ID,
            HLA_ALLELE,
            ALT_EPI_SEQ,
            TEAM,
            LOG_BINDING
        ) %>% 
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
        relevel_df()
}

make_pmhc_dbi <- function(prediction_dbi, validation_dbi, team){
    pmhc_dbi <-  prediction_dbi %>% 
        dplyr::filter(TEAM == team) %>% 
        dplyr::inner_join(
            validation_dbi,
            by = c("HLA_ALLELE", "ALT_EPI_SEQ")) %>% 
        dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ) %>% 
        dplyr::distinct()
}


# validation 2

make_binding_scatterplot_df <- function(round, source, team){
    prediction_dbi <- make_combined_prediction_dbi(round, source, hla_binding = T) %>% 
        dplyr::filter(TEAM == team) %>% 
        dplyr::select(-TEAM) %>% 
        dplyr::rename(LOG_PREDICTED_BINDING = LOG_BINDING) 
    
    validation_dbi <- 
        make_lji_binding_dbi () %>% 
        dplyr::rename(LOG_MEASURED_BINDING = LOG_BINDING) 
    
    scatterplot_df <- 
        dplyr::inner_join(
            prediction_dbi,
            validation_dbi,
            by = c("HLA_ALLELE", "ALT_EPI_SEQ")) %>% 
        tibble::as_tibble()
}

# validation 3 ----

make_validation_df <- function(round, source, team){
    combined_dbi <- 
        make_combined_prediction_dbi(round, source, max_rank = 20) %>%
        dplyr::select(PATIENT_ID, TEAM, HLA_ALLELE, ALT_EPI_SEQ, RANK) 
    
    validated_bindings_dbi <- dplyr::inner_join(
        make_binding_assay_dbi(), 
        combined_dbi, 
        by = c("PATIENT_ID", "HLA_ALLELE", "ALT_EPI_SEQ"))
    
    combined_dbi2 <- combined_dbi %>% 
        dplyr::select(-HLA_ALLELE) %>% 
        dplyr::distinct()
        
    validated_epitopes_dbi <- dplyr::inner_join(
        make_epitope_assay_dbi(),
        combined_dbi2, 
        by = c("PATIENT_ID", "ALT_EPI_SEQ"))
    
    df <- 
        dplyr::union_all(validated_epitopes_dbi, validated_bindings_dbi) %>% 
        dplyr::select(TEAM, PATIENT_ID, ASSAY, RESULT, RANK) %>% 
        dplyr::distinct() %>% 
        as_tibble() %>% 
        dplyr::group_by(TEAM, PATIENT_ID, ASSAY) %>%
        dplyr::summarise(
            COUNT = n(), 
            MEAN_RANK = mean(RANK),
            RATE = mean(as.numeric(RESULT))
        ) %>%
        code_df_by_team(team)
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

