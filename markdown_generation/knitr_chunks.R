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
source("transform_functions.R")


## @knitr submissions

submission_dfs <- make_submission_plot_dfs(ROUND, SOURCE, TEAM)

log_peptides_df <- submission_dfs[["log_peptides_df"]] 
    
peptide_length_df <- submission_dfs[["peptide_length_df"]]

agretopicity_df <- submission_dfs[["agretopicity_df"]]

epitope_overlap_df <- submission_dfs[["epitope_overlap_df"]] 
    
## @knitr validation1

dotplot_df <- make_binding_dotplot_df(ROUND, SOURCE, TEAM)

## @knitr validation2

scatterplot_df <- make_binding_scatterplot_df(ROUND, SOURCE, TEAM)

## @knitr validation3

validation_df <- make_validation_df(ROUND, SOURCE, TEAM)

## @knitr variant counts

variant_count_df <- make_variant_counts_df(ROUND, TEAM)

## @knitr variant overlap

variant_overlap_df <- make_variant_overlap_df(ROUND, TEAM)

## @knitr survey results

survey_df <- make_survey_df(TEAM)
    




    
