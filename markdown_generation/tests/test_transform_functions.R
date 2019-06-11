library(tibble)
library(testthat)
source("../transform_functions.R")
context("Tranform functions")


test_that("create_filter_string", {
    expect_equal(
        create_filter_string("Column", "Value"),
        "Column %in% c('Value')")
    expect_equal(
        create_filter_string("Column", "Value", exclude = T),
        "!Column %in% c('Value')")
    expect_equal(
        create_filter_string("Column", c("Value1", "Value2")),
        "Column %in% c('Value1', 'Value2')")
})


test_that("relevel_df_column", {
    df <- tibble::tribble(
        ~COL1, ~COL2, ~COL3,
        "A", 2, "GROUP1",
        "A", 3, "GROUP2",
        "B", 2, "GROUP1",
        "B", 1, "GROUP2",
        "C", 1, "GROUP1",
        "C", 2, "GROUP2"
    )
    expect_equal(code_df_by_column(input_df, "T1", "COL1")$New_column,
                 c("matches", "doesn't match", "matches", "doesn't match", NA))
    expect_equal(code_df_by_column(input_df, "T2", "COL1")$New_column,
                 c("doesn't match", "matches", "doesn't match", "doesn't match", NA))
    expect_equal(code_df_by_column(input_df, "A", "COL2")$New_column,
                 c("matches", "doesn't match", "doesn't match", NA, "doesn't match"))

})


test_that("code_df_by_column", {
    input_df <- tibble::tribble(
        ~COL1, ~COL2,
        "T1", "A",
        "T2", "B",
        "T1", "C",
        "T3", NA,
        NA,   "D"
    )
    expect_equal(code_df_by_column(input_df, "T1", "COL1")$New_column,
                 c("matches", "doesn't match", "matches", "doesn't match", NA))
    expect_equal(code_df_by_column(input_df, "T2", "COL1")$New_column,
                 c("doesn't match", "matches", "doesn't match", "doesn't match", NA))
    expect_equal(code_df_by_column(input_df, "A", "COL2")$New_column,
                 c("matches", "doesn't match", "doesn't match", NA, "doesn't match"))
    expect_equal(
        code_df_by_column(
            input_df, 
            "A", 
            "COL2",
            new_column  = "RESULT",
            new_value   = T, 
            na_value    = "nothing", 
            other_value = F
        )$RESULT,
        c(T, F, F, "nothing", F))
    
})



test_that("create_pmhc_combinations_df",{
    input_df <- tibble::tribble(
        ~PATIENT_ID, ~TEAM, ~PMHC,
        "PAT1", "T1", c("P1", "P2", "P3"),
        "PAT1", "T2", c("P1", "P2", "P3"),
        "PAT1", "T3", c("P1", "P2", "P4")
    )
    result_df <- create_pmhc_combinations_df(input_df)
    expected_df <- tibble::tibble(
        PATIENT_ID = rep("PAT1", 6),
        TEAM = c("T1", "T1", "T2", "T2", "T3", "T3"),
        PMHC = list(
            c("P1", "P2", "P3"),
            c("P1", "P2", "P3"),
            c("P1", "P2", "P3"),
            c("P1", "P2", "P3"),
            c("P1", "P2", "P4"),
            c("P1", "P2", "P4")
        ),
        PMHC2 = list(
            c("P1", "P2", "P3"),
            c("P1", "P2", "P4"),
            c("P1", "P2", "P3"),
            c("P1", "P2", "P4"),
            c("P1", "P2", "P3"),
            c("P1", "P2", "P3")
        )
    )
    expect_true(identical(result_df, expected_df))
})
test_that("calc_overlap_perc", {
    v1 <- c("A", "B", "C")
    v2 <- c("B", "C", "A")
    v3 <- c("A", "B", "C", "D")
    v4 <- c("D", "B", "C", "A")
    v5 <- as.character(c())
    v6 <- c("X", "Y", "Z")
    v7 <- c(1, 2, 3)
    lst <- list("A", "B", "C")
    expect_equal(calc_overlap_perc(v1, v1), 1.0)
    expect_equal(calc_overlap_perc(v1, v2), 1.0)
    expect_equal(calc_overlap_perc(v1, v3), 1.0)
    expect_equal(calc_overlap_perc(v1, v4), 2/3)
    expect_equal(calc_overlap_perc(v1, v5), 0.0)
    expect_equal(calc_overlap_perc(v1, v6), 0.0)
    expect_equal(calc_overlap_perc(v1, v6), 0.0)
    expect_error(calc_overlap_perc(v1, v7), 
                 "`v2` must be <chr>, not <dbl>.")
    expect_error(calc_overlap_perc(v1, lst), 
                 "`v2` must be <chr>, not <list>.")
})



# df <- tibble::tribble(
#     ~TEAM, ~PMHC,
#     "T1", "P1",
#     "T1", "P2",
#     "T1", "P3",
#     "T2", "P1",
#     "T2", "P2",
#     "T2", "P4"
# )



# df <- tibble::tribble(
#     ~PATIENT_ID, ~TEAM, ~HLA_ALLELE, ~ALT_EPI_SEQ,
#     "PAT1", "T1", "H1", "A1",
#     "PAT1", "T1", "H2", "A2",
#     "PAT1", "T2", "H1", "A1",
#     "PAT1", "T2", "H3", "A3",
#     "PAT1", "T3", "H1", "A1",
#     "PAT1", "T3", "H2", "A2",
#     "PAT2", "T1", "H1", "A1",
#     "PAT2", "T1", "H2", "A2",
#     "PAT2", "T3", "H1", "A1",
#     "PAT2", "T3", "H2", "A2",
#     "PAT2", "T3", "H3", "A3"
# )


# test_that("calc_epitope_overlap_by_teams", {
#     df1 <- tibble::tribble(
#         ~TEAM, ~PMHC,
#         "T1", "P1",
#         "T1", "P2",
#         "T1", "P3",
#         "T2", "P1",
#         "T2", "P2",
#         "T2", "P4"
#     )
#     df2 <- tibble::tribble(
#         ~X, ~PMHC,
#         "T1", "P1"
#     )
#     df3 <- tibble::tribble(
#         ~TEAM, ~Y,
#         "T1", "P1"
#     )
#     expect_equal(calc_epitope_overlap_by_teams(df1, "T1", "T1"), 1.0)
#     expect_equal(calc_epitope_overlap_by_teams(df1, "T1", "T2"), 2/3)
#     expect_equal(calc_epitope_overlap_by_teams(df1, "T2", "T1"), 2/3)
#     expect_equal(calc_epitope_overlap_by_teams(df1, "T1", "T3"), 0.0)
#     expect_error(calc_epitope_overlap_by_teams(df2, "T1", "T2"))
#     expect_error(calc_epitope_overlap_by_teams(df3, "T1", "T2"))
# })

