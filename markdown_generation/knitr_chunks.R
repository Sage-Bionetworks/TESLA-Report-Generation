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

validation_dfs <- make_validation_dfs(ROUND, SOURCE, TEAM)

TCR_NANOPARTICLE_df <- validation_dfs[["TCR_NANOPARTICLE_df"]]

TCR_FLOW_I_df <- validation_dfs[["TCR_FLOW_I_df"]]

TCR_FLOW_II_df <- validation_dfs[["TCR_FLOW_II_df"]]

TCELL_REACTIVITY_df <- validation_dfs[["TCELL_REACTIVITY_df"]]


## @knitr variant counts

variant_count_df <- make_variant_counts_df(patients, TEAM)

## @knitr variant overlap

variant_overlap_df <- make_variant_overlap_df(ROUND, TEAM)

## @knitr survey results

survey_df <- make_survey_df(TEAM)
    




    
