## @knitr setup
library(yaml)
library(tidyverse)
library(bigrquery)
library(magrittr)

config <- read_yaml("config.yaml")
ROUND  <- config$round
SOURCE <- config$source
TEAM   <- config$team

source("markdown_functions.R")
source("plot_functions.R")
source("query_functions.R")


## @knitr submissions

submission_dfs <- make_submission_plot_dfs(ROUND, SOURCE)

log_peptides_df <-
    submission_dfs[["log_peptides_df"]] %>%
    code_df_by_team(TEAM)

peptide_length_df <-
    submission_dfs[["peptide_length_df"]] %>%
    code_df_by_team(TEAM)

agretopicity_df <- 
    submission_dfs[["agretopicity_df"]] %>%
    code_df_by_team(TEAM)

overlap_df <- 
    submission_dfs[["overlap_df"]] %>%
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

dotplot_df <- make_binding_dotplot_df(ROUND, SOURCE, TEAM)


## @knitr validation2

scatterplot_df <- make_binding_scatterplot_df(ROUND, SOURCE, TEAM)

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

variant_count_df <- make_variant_counts_df(patients, TEAM)

## @knitr variant overlap

variant_overlap_df <- make_variant_overlap_df(ROUND, TEAM)

## @knitr survey results

survey_df <- make_survey_df(TEAM)
    




    
