require(dplyr)
require(tidyr)
require(magrittr)
require(vctrs)
require(assertthat)


# utils ----

code_df_by_column <- function(
    df, matching_values, column, 
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

code_df_by_team <- purrr::partial(
    code_df_by_column,
    column      = "TEAM",
    new_column  = "team_status",
    new_value   = "Your team", 
    na_value    = "No Team", 
    other_value = "Other teams")

# submission plot functions ----




create_median_overlap_df <- function(df, team){
    df <- dplyr::as_tibble(df)
    assertthat::assert_that(assertthat::has_name(df, c("TEAM"))) 
    assertthat::assert_that(assertthat::has_name(df, c("PATIENT_ID"))) 
    assertthat::assert_that(assertthat::has_name(df, c("HLA_ALLELE"))) 
    assertthat::assert_that(assertthat::has_name(df, c("ALT_EPI_SEQ"))) 
    df %>% 
        dplyr::mutate(PMHC = stringr::str_c(HLA_ALLELE, "_", ALT_EPI_SEQ)) %>%
        dplyr::select(TEAM, PATIENT_ID, PMHC) %>% 
        dplyr::group_by(PATIENT_ID, TEAM) %>% 
        dplyr::summarise(PMHC = list(PMHC)) %>% 
        dplyr::group_by(PATIENT_ID) %>% 
        dplyr::do(create_pmhc_combinations_df(.)) %>% 
        dplyr::mutate(OVERLAP = purrr::map2_dbl(
            PMHC, 
            PMHC2, 
            calc_overlap_perc)) %>% 
        dplyr::select(PATIENT_ID, TEAM, OVERLAP) %>% 
        dplyr::group_by(PATIENT_ID, TEAM) %>% 
        dplyr::summarise(MED_OVERLAP = median(OVERLAP)) %>% 
        dplyr::ungroup() 
}


create_pmhc_combinations_df <- function(df){
    assertthat::assert_that(assertthat::has_name(df, c("TEAM"))) 
    assertthat::assert_that(assertthat::has_name(df, c("PATIENT_ID"))) 
    assertthat::assert_that(assertthat::has_name(df, c("PMHC"))) 
    df %>% 
        tidyr::crossing(., .) %>% 
        dplyr::filter(TEAM != TEAM1) %>% 
        dplyr::select(PATIENT_ID, TEAM, PMHC, PMHC2 = PMHC1)
}

calc_overlap_perc <- function(v1, v2){
    vctrs::vec_assert(v1)
    vctrs::vec_assert(v2)
    vctrs::vec_assert(v1, ptype = character())
    vctrs::vec_assert(v2, ptype = character())
    l1 <- length(v1)
    l2 <- length(v2)
    min_length <- min(l1, l2)
    if(min_length == 0) return(0.0)
    v1 <- v1[1:min_length]
    v2 <- v2[1:min_length]
    n_diff <- length(setdiff(v1, v2))
    result   <-  1 - (n_diff / l1)
}
