relevel_df <- function(df){
    new_lvs <- df %>% 
        filter(TEAM == "Measured") %>% 
        arrange(LOG_BINDING) %>% 
        use_series(ALT_EPI_SEQ) %>% 
        unique()
    df$ALT_EPI_SEQ = factor(df$ALT_EPI_SEQ, levels = new_lvs, ordered = TRUE)
    return(df)
}

code_df_by_team <- function(plot_data_df, team){
    mutate(plot_data_df, 
           team_status = ifelse(
               is.na(TEAM), 
               "No team", 
               ifelse(
                   TEAM == team, 
                   "Your team", 
                   "Other teams")))
}



make_dotplot <- function(allele, df){
    plot_df <- dplyr::filter(df, HLA_ALLELE == allele) 
    if(nrow(plot_df) > 0){
        make_validation_dotplot(plot_df)
        text <- ""
    } else {
        text <- "None of your teams predicted neoepitopes with binding scores for this allele were validated."
    }
    return(text)
}



make_scatterplot <- function(pat, df){
    plot_df <- dplyr::filter(df, PATIENT_ID == pat) 
    if(nrow(plot_df) > 0){
        make_validation_scatterplot(plot_df)
        text <- ""
    } else {
        text <- "None of your teams predicted neoepitopes with binding scores for this patient were validated."
    }
    return(text)
}

create_p_common_matrix <- function(patient, df, max_rank = 20){
    patient_df <- df %>% 
        filter(PATIENT_ID == patient) %>% 
        dplyr::select(PMHC, TEAM)
    teams <- sort(unique(patient_df$TEAM))
    matrix <- matrix(
        0, 
        nrow = length(teams),
        ncol = length(teams), 
        dimnames = list(teams, teams))
    for(team1 in teams){
        for (team2 in teams){
            matrix[team1, team2] <- 
                find_percent_unique(patient_df, team1, team2, max_rank)
        }
    }
    return(matrix)
}

find_percent_unique <- function(df, team1, team2, max_rank){
    e1 <- df %>% 
        filter(TEAM == team1) %>% 
        use_series(PMHC)
    e2 <- df %>% 
        filter(TEAM == team2) %>% 
        use_series(PMHC)
    en <- min(max_rank, length(e1), length(e2))
    if(length(e1) > 0) e1 <- e1[1:en]
    if(length(e2) > 0) e2 <- e2[1:en]
    diff_n  <- length(setdiff(e1, e2))
    if(length(e1) == 0) return(0.0)
    result <- 1 - diff_n / length(e1)
    return(result)
}

create_median_table <- function(matrix){
    matrix %>%
        apply(2, median) %>%
        data.frame %>%
        rownames_to_column("team") %>%
        set_names(c("team", "median")) 
}

make_variant_boxplot <- function(pat, df){
    plot_df <- dplyr::filter(df, PATIENT_ID == pat) 
    if(!is.null(plot_df)){
        make_variant_boxplot_obj(plot_df)
        text <- ""
    } else {
        text <- "Your team had none of the plotted variants for this patient."
    }
    return(text)
}

make_variant_histogram <- function(pat, df){
    plot_df <- dplyr::filter(df, patient == pat) 
    if("Your team" %in% plot_df$team_status){
        make_variant_histogram_obj(plot_df)
        text <- ""
    } else {
        text <- "Your team has no variant data for this patient."
    }
    return(text)
}


make_survey_barchart <- function(cat, df){
    plot_df <- dplyr::filter(df, CATEGORY == cat) 
    if(nrow(plot_df) > 0){
        make_survey_barchart_obj(plot_df)
        text <- ""
    } else {
        text <- "Your team did not answer questions in this category."
    }
    return(text)
}


# functions written br Kristen

make_overlap_list <- function(synapse_id){
    data <- synapse_id %>% 
        synapser::synGet() %>% 
        magrittr::use_series(path) %>% 
        read.delim(header = FALSE)
    patMatrix = makePatientMatrix(data)
    nanVals = apply(patMatrix,MARGIN = 2,function(x){if (sum(x %in% NaN) == nrow(patMatrix)) return(1) else (return(0))})
    invalidVals = apply(patMatrix,MARGIN = 2,function(x){if (sum(x,na.rm = TRUE) == (nrow(patMatrix)-1)) return(1) else (return(0))})
    toRemove = c(which(nanVals == 1), which(invalidVals == 1))
    tempMatrix = patMatrix[-toRemove,-toRemove]
    patMatrix = tempMatrix
    rownames(patMatrix) = as.vector(sapply(rownames(patMatrix),function(x){unlist(strsplit(x = x,split = "_"))[1]}))
    colnames(patMatrix) = as.vector(sapply(colnames(patMatrix),function(x){unlist(strsplit(x = x,split = "_"))[1]}))
    aves = apply(patMatrix,MARGIN = 2,function(x){median(x,na.rm = TRUE)})
}

# Returns a matrix of overlap values for every two team combination.
# Columns are the values wrt team A.
makePatientMatrix=function(patientDataset){
    allbirds = as.list(unique(patientDataset[,1]))
    x = sapply(allbirds,function(x){ birdUnique(bird=x,patientDataset=patientDataset)})
    z = birdUnique(bird=allbirds[[1]],patientDataset = patientDataset)
    colnames(x) = names(z)
    return(x)
}

# Returns a vector of values for all pairs where A = bird.
# Value returned is setdiff(A,B) / size(A)
birdUnique=function(bird,patientDataset,type="records:"){
    smallSet = (patientDataset[patientDataset[,1] == bird & patientDataset[,3] == type,])
    ssDM = data.matrix(smallSet[,4:6])
    x = apply(ssDM,MARGIN = 1,function(x){ x[1] / (x[1] + x[3]) })
    names(x) = smallSet[,2]
    return(x)
}

getbird=function(x){as.vector(strsplit(x,split = "_")[[1]][1])}




