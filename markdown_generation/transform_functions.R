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


# relevel_df_column <- function(
#     df, 
#     new_col, 
#     order_col, 
#     group_col, 
#     filter_col, 
#     filter_vals
# ){
#     new_levels <- get_ordered_column_values(
#         df, group_col, order_col, filter_col, filter_vals)
#     
#     relevel_column(df, new_levels, new_col)
# }
# 
# 
# 
# relevel_epitopes_by_binding <- purrr::partial(
#     relevel_df_column,
#     new_col = "ALT_EPI_SEQ", 
#     order_col = "LOG_BINDING", 
#     group_col = "ALT_EPI_SEQ", 
#     filter_col = "TEAM", 
#     filter_vals = "Measured"
# )
# 
# get_ordered_column_values <- function(
#     df, group_col, order_col, filter_col, filter_vals){
#     
#     order_col_symbol <- rlang::sym(order_col)
#     filter_string <- create_filter_string(filter_col, filter_vals)
#     
#     df %>% 
#         dplyr_df_by_string(filter_string, dplyr::filter) %>% 
#         dplyr::arrange(!!order_col_symbol) %>% 
#         magrittr::extract2(group_col) %>% 
#         unique()
# }
# 
# relevel_column <- function(df, new_levels, column){
#     string <- create_relevel_mutate_string(column, new_levels)
#     dplyr_df_by_string(df, string, dplyr::mutate)
# }
# 
# dplyr_mutate_by_string <- function(
#     df, 
#     string, 
#     dplyr_function = c(dplyr::filter, dplyr::mutate)
# ){
#     expr <- rlang::parse_expr(string)
#     result_df  <- dplyr_function(df, !!expr)
# }
# 
# create_relevel_mutate_string <- function(column, new_levels){
#     mutate_string <- stringr::str_c(
#         column,
#         " = forcats::fct_relevel(", 
#         column, 
#         ", new_levels)") 
# }
# 
# create_filter_string <- function(column, values, exclude = F){
#     filter_string <- values %>% 
#         stringr::str_c("'", ., "'") %>% 
#         stringr::str_c(collapse = ", ") %>% 
#         stringr::str_c(column, " %in% c(", ., ")") 
#     if(exclude) filter_string <- stringr::str_c("!", filter_string)
#     filter_string
# }

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
