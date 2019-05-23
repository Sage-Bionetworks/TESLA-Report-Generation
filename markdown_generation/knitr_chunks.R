## @knitr submissions

submission_dbis <- make_dbis_for_submission_plots(round, src)

log_peptides_df <- submission_dbis[["log_peptides_dbi"]] %>% 
    dplyr::as_tibble() %>% 
    code_df_by_team(CURRENT_TEAM)

peptide_length_df <- submission_dbis[["peptide_length_dbi"]] %>% 
    dplyr::as_tibble() %>% 
    code_df_by_team(CURRENT_TEAM)

agretopicity_df <- submission_dbis[["agretopicity_dbi"]] %>% 
    dplyr::as_tibble() %>% 
    code_df_by_team(CURRENT_TEAM)

overlap_df <- submission_dbis[["overlap_dbi"]] %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(PMHC = stringr::str_c(HLA_ALLELE, "_", ALT_EPI_SEQ)) %>% 
    dplyr::select(TEAM, PATIENT_ID, PMHC, RANK)

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
    code_df_by_team(CURRENT_TEAM)

## @knitr validation1

prediction_dbi <- 
    make_binding_prediction_dbi(round, src) 

team_prediction_dbi <- prediction_dbi %>% 
    dplyr::filter(TEAM == CURRENT_TEAM) %>% 
    dplyr::select(
        PATIENT_ID, 
        HLA_ALLELE, 
        ALT_EPI_SEQ,
        TEAM, 
        LOG_PREDICTED_BINDING) %>% 
    dplyr::mutate(TEAM = "Your team") %>% 
    dplyr::rename(LOG_BINDING = LOG_PREDICTED_BINDING)

validation_dbi <- 
    DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>%
    dplyr::tbl("Validated_Bindings") %>%
    dplyr::filter(!is.na(LJI_BINDING)) %>%
    dplyr::mutate(LOG_MEASURED_BINDING = log10(LJI_BINDING + 1))  %>%
    dplyr::select(ALT_EPI_SEQ, HLA_ALLELE, PATIENT_ID, LOG_MEASURED_BINDING) %>% 
    dplyr::mutate(TEAM = "Measured") %>% 
    dplyr::rename(LOG_BINDING = LOG_MEASURED_BINDING)

validated_epitope_dbi <- 
    dplyr::inner_join(
        team_prediction_dbi, 
        validation_dbi, 
        by = c("PATIENT_ID", "HLA_ALLELE", "ALT_EPI_SEQ")) %>% 
    dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ) %>% 
    dplyr::distinct()

nonteam_prediction_dbi <- prediction_dbi %>% 
    dplyr::filter(TEAM != CURRENT_TEAM) %>% 
    dplyr::select(
        PATIENT_ID, 
        HLA_ALLELE, 
        ALT_EPI_SEQ,
        TEAM, 
        LOG_PREDICTED_BINDING) %>% 
    dplyr::mutate(TEAM = "Other teams") %>% 
    dplyr::rename(LOG_BINDING = LOG_PREDICTED_BINDING)

dotplot_df <-
    list(team_prediction_dbi, nonteam_prediction_dbi, validation_dbi) %>% 
    purrr::map(dplyr::inner_join, validated_epitope_dbi) %>% 
    purrr::map(dplyr::as_tibble) %>% 
    dplyr::bind_rows() %>% 
    relevel_df


## @knitr validation2

prediction_dbi <- 
    make_binding_prediction_dbi(round, src)


validation_dbi <- 
    DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>%
    dplyr::tbl("Validated_Bindings") %>%
    dplyr::filter(!is.na(LJI_BINDING)) %>%
    dplyr::mutate(LOG_MEASURED_BINDING = log10(LJI_BINDING + 1))  %>%
    dplyr::select(ALT_EPI_SEQ, HLA_ALLELE, PATIENT_ID, LOG_MEASURED_BINDING) 


team_prediction_dbi <- prediction_dbi %>% 
    dplyr::filter(TEAM == CURRENT_TEAM) %>% 
    dplyr::select(
        PATIENT_ID, 
        HLA_ALLELE, 
        ALT_EPI_SEQ,
        LOG_PREDICTED_BINDING)

scatterplot_df <- 
    dplyr::inner_join(
        team_prediction_dbi,
        validation_dbi,
        by = c("HLA_ALLELE", "ALT_EPI_SEQ", "PATIENT_ID")) %>% 
    tibble::as_tibble()

## @knitr validation3

submission_dbi <- 
    DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Submissions") %>% 
    dplyr::filter(ROUND == RND) %>% 
    dplyr::select(SUBMISSION_ID, PATIENT_ID, TEAM)

prediction_dbi <-
    DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Predictions") %>% 
    dplyr::filter(RANK <= 20) %>% 
    dplyr::filter(SOURCE == SRC) %>% 
    dplyr::inner_join(submission_dbi) %>% 
    dplyr::select(PATIENT_ID, TEAM, HLA_ALLELE, ALT_EPI_SEQ, RANK)

prediction_dbi2 <- prediction_dbi %>% 
    dplyr::select(-HLA_ALLELE)

validation_bindings_dbi <- 
    DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Validated_Bindings") %>% 
    dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ, TCR_NANOPARTICLE, TCR_FLOW_I, TCR_FLOW_II) 

TCR_NANOPARTICLE_df <- validation_bindings_dbi %>% 
    dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ, TCR_NANOPARTICLE) %>% 
    dplyr::filter(!is.na(TCR_NANOPARTICLE)) %>% 
    dplyr::inner_join(prediction_dbi) %>% 
    dplyr::mutate(ASSAY_NUM = ifelse(TCR_NANOPARTICLE == "+", 1, 0)) %>% 
    dplyr::group_by(TEAM, PATIENT_ID) %>% 
    dplyr::summarise(COUNT = n(), MEAN_RANK = mean(RANK), RATE = mean(ASSAY_NUM)) %>% 
    tibble::as_tibble() %>% 
    code_df_by_team(CURRENT_TEAM)

TCR_FLOW_I_df <- validation_bindings_dbi %>% 
    dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ, TCR_FLOW_I) %>% 
    dplyr::filter(!is.na(TCR_FLOW_I)) %>% 
    dplyr::inner_join(prediction_dbi) %>% 
    dplyr::mutate(ASSAY_NUM = ifelse(TCR_FLOW_I == "+", 1, 0)) %>% 
    dplyr::group_by(TEAM, PATIENT_ID) %>% 
    dplyr::summarise(COUNT = n(), MEAN_RANK = mean(RANK), RATE = mean(ASSAY_NUM)) %>% 
    tibble::as_tibble() %>% 
    code_df_by_team(CURRENT_TEAM)

TCR_FLOW_II_df <- validation_bindings_dbi %>% 
    dplyr::select(PATIENT_ID, HLA_ALLELE, ALT_EPI_SEQ, TCR_FLOW_II) %>% 
    dplyr::filter(!is.na(TCR_FLOW_II)) %>% 
    dplyr::inner_join(prediction_dbi) %>% 
    dplyr::mutate(ASSAY_NUM = ifelse(TCR_FLOW_II == "+", 1, 0)) %>% 
    dplyr::group_by(TEAM, PATIENT_ID) %>% 
    dplyr::summarise(COUNT = n(), MEAN_RANK = mean(RANK), RATE = mean(ASSAY_NUM)) %>% 
    tibble::as_tibble() %>% 
    code_df_by_team(CURRENT_TEAM)

TCELL_REACTIVITY_df <- 
    DBI::dbConnect(bigquery(), project = "neoepitopes", dataset = "Version_3") %>% 
    dplyr::tbl("Validated_Epitopes") %>% 
    dplyr::select(PATIENT_ID, ALT_EPI_SEQ, TCELL_REACTIVITY) %>% 
    dplyr::filter(!is.na( TCELL_REACTIVITY)) %>% 
    dplyr::inner_join(prediction_dbi2) %>%
    dplyr::mutate(ASSAY_NUM = ifelse(TCELL_REACTIVITY == "+", 1, 0)) %>% 
    dplyr::group_by(TEAM, PATIENT_ID) %>% 
    dplyr::summarise(COUNT = n(), MEAN_RANK = mean(RANK), RATE = mean(ASSAY_NUM)) %>% 
    tibble::as_tibble() %>% 
    code_df_by_team(CURRENT_TEAM)

## @knitr variant counts

variant_count_df <- make_variant_counts_df(patients, CURRENT_TEAM)

## @knitr variant overlap

variant_overlap_df <- make_variant_overlap_df(round, CURRENT_TEAM)

## @knitr survey results

survey_df <- make_survey_df(CURRENT_TEAM)
    




    
