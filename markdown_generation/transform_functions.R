require(dplyr)
require(tidyr)
require(forcats)
require(rlang)
require(stringr)
require(magrittr)
require(vctrs)
require(assertthat)


# utils ----

code_df_by_column <- function(
    df, matching_values, column, 
    new_column      = New_column,
    match_value     = "matches", 
    non_match_value = "doesn't match",
    na_value        = NA_character_
){
    dplyr::mutate(df, {{new_column}} := dplyr::case_when(
        is.na({{column}})                ~ na_value,
        {{column}} %in% matching_values  ~ match_value,
        TRUE                             ~ non_match_value
    ))
}

code_df_by_team <- purrr::partial(
    code_df_by_column,
    column      = TEAM,
    new_column  = team_status,
    match_value   = "Your team", 
    non_match_value = "Other teams",
    na_value    = "No Team"
)



relevel_df_column <- function(df, col, lvl_function, ...){
    new_lvls <- lvl_function(df, ...)
    dplyr::mutate(df, {{col}} := forcats::fct_relevel({{col}}, new_lvls))
}

relevel_epitope_column <- purrr::partial(
    relevel_df_column,
    col = ALT_EPI_SEQ,
    lvl_function = get_epitope_levels
)

get_ordered_values <- function(df, arrange_col, pull_col, filter_expr = T){
    df %>% 
        dplyr::filter(!!!rlang::enquos(filter_expr)) %>%
        dplyr::arrange({{arrange_col}}) %>% 
        dplyr::pull({{pull_col}}) %>% 
        unique
}

get_epitope_levels <- purrr::partial(
    get_ordered_values,
    arrange_col = LOG_BINDING,
    pull_col = ALT_EPI_SEQ,
    filter_expr = TEAM == "Measured"
)

# submission plot functions ----


create_median_overlap_df <- function(dbi, team){
    df <- dplyr::as_tibble(dbi)
    assertthat::assert_that(assertthat::has_name(df, c("TEAM"))) 
    assertthat::assert_that(assertthat::has_name(df, c("PATIENT_ID"))) 
    assertthat::assert_that(assertthat::has_name(df, c("PMHC"))) 
    df %>% 
        dplyr::group_by(PATIENT_ID, TEAM) %>%
        dplyr::summarise(PMHCS = list(PMHC)) %>%
        dplyr::ungroup() %>% 
        dplyr::full_join(., ., by = c("PATIENT_ID"), suffix = c("", "2")) %>% 
        dplyr::filter(TEAM != TEAM2) %>%
        dplyr::select(-TEAM2) %>% 
        dplyr::mutate(OVERLAP = purrr::map2_dbl(
            PMHCS, 
            PMHCS2, 
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
