# The functions take a df and print a plot, or return a text message 
## submission plots ----

make_submissions_boxplot <- function(df){
    if(nrow(df) > 0){
        print(make_submissions_boxplot_obj(df))
        text <- ""
    } else {
        text <- "Empty input df"
    }
    return(text)
}

make_agretopicity_boxplot <- function(df){
    if(nrow(df) > 0){
        print(make_agretopicity_boxplot_obj(df))
        text <- ""
    } else {
        text <- "Empty input df"
    }
    return(text)
}

make_pep_length_barchart <- function(df){
    if(nrow(df) > 0){
        print(make_pep_length_barchart_obj(df))
        text <- ""
    } else {
        text <- "Empty input df"
    }
    return(text)
}

make_epitope_overlap_boxplot <- function(df){
    if(nrow(df) > 0){
        print(make_epitope_overlap_boxplot_obj(df))
        text <- ""
    } else {
        text <- "Empty input df"
    }
    return(text)
}

## validation plots ----
make_validation_dotplot <- function(allele, df){
    plot_df <- dplyr::filter(df, HLA_ALLELE == allele) 
    if(nrow(plot_df) > 0){
        print(make_validation_dotplot_obj(plot_df))
        text <- ""
    } else {
        text <- "None of your teams predicted neoepitopes with binding scores for this allele were validated."
    }
    return(text)
}

make_validation_scatterplot <- function(pat, df){
    plot_df <- dplyr::filter(df, PATIENT_ID == pat) 
    if(nrow(plot_df) > 0){
        print(make_validation_scatterplot_obj(plot_df))
        text <- ""
    } else {
        text <- "None of your teams predicted neoepitopes with binding scores for this patient were validated."
    }
    return(text)
}

make_validated_submissions_boxplot <- function(df, assay, column, ylab){
    plot_df <- dplyr::filter(df, ASSAY == assay) 
    if(nrow(plot_df) > 0){
        print(make_validated_submissions_boxplot_obj(plot_df, column, ylab))
        text <- ""
    } else {
        text <- "Empty input df"
    }
    return(text)
}

# variant plots ----

make_variant_boxplot <- function(pat, df){
    plot_df <- dplyr::filter(df, PATIENT_ID == pat) 
    if(!is.null(plot_df)){
        print(make_variant_boxplot_obj(plot_df))
        text <- ""
    } else {
        text <- "Your team had none of the plotted variants for this patient."
    }
    return(text)
}

make_variant_histogram <- function(pat, df){
    plot_df <- dplyr::filter(df, patient == pat) 
    if("Your team" %in% plot_df$team_status){
        print(make_variant_histogram_obj(plot_df))
        text <- ""
    } else {
        text <- "Your team has no variant data for this patient."
    }
    return(text)
}

# survey

make_survey_barchart <- function(cat, df){
    plot_df <- dplyr::filter(df, CATEGORY == cat) 
    if(nrow(plot_df) > 0){
        print(make_survey_barchart_obj(plot_df))
        text <- ""
    } else {
        text <- "Your team did not answer questions in this category."
    }
    return(text)
}




# These functions take a df and return a plot_object 

## submission plots ----


make_submissions_boxplot_obj <- function(df){
    df %>% 
        ggplot(aes(x = PATIENT_ID, y = LOG_COUNT)) +
        geom_boxplot(color = "black", fill = "white", outlier.shape = NA) +
        geom_jitter(aes(color = team_status, shape = team_status), size = 4) + 
        labs(shape = "") +
        labs(color = "") +
        scale_color_manual(values = c("black", "red")) +
        ylab("Log10 count submitted peptides") +
        xlab("Patient") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
}

make_agretopicity_boxplot_obj <- function(df){
    df %>% 
        ggplot(aes(x = TEAM,  y = LOG_AGRETOPICITY)) +
        geom_boxplot(aes(fill = team_status)) +
        scale_fill_manual(values=c("white", "red")) +
        labs(fill = "") +
        xlab("Team") +
        ylab("Log10 Agretopicity Index") +
        theme_bw() +
        theme(axis.text.x = element_blank())
}

make_pep_length_barchart_obj <- function(df){
    df %>% 
        ggplot(aes(x = TEAM)) +
        geom_bar(aes(fill = as.factor(PEP_LEN), color = team_status), size = 2) +
        labs(color = 'Bar outline(Team)') +
        labs(fill  = 'Bar color(Peptide Length)') +
        scale_color_manual(values = c("black", "red")) +
        xlab("Team") +
        ylab("Peptide Length Counts") +
        theme_bw() +
        theme(axis.text.x = element_blank()) 
}

make_epitope_overlap_boxplot_obj <- function(df){
    df %>% 
        ggplot(aes(x = PATIENT_ID, y = SCORE)) +
        geom_boxplot(color = "black", fill = "white", outlier.shape = NA) +
        geom_jitter(aes(color = team_status, shape = team_status), size = 4) + 
        labs(shape = "") +
        labs(color = "") +
        scale_color_manual(values = c("black", "red")) +
        ylab("Median overlap of top 20 epitopes with other teams") +
        xlab("Patient") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
}

## validation 1 ----

make_validation_dotplot_obj <- function(df){
    lines_vector <- make_dotplot_lines_vector(df)
    df %>% 
        ggplot(aes(
            x = interaction(ALT_EPI_SEQ, PATIENT_ID, sep = ";Pat:"),
            y = LOG_BINDING)) +
        geom_jitter(aes(color = TEAM, size = TEAM, alpha = TEAM), width = .4) +
        scale_alpha_manual(values = c(0.7, 1, 0.7)) +
        scale_size_manual(values = c(4, 2, 4)) +
        scale_color_manual(values = c("green", "black", "red")) +
        labs(size = "", color = "", alpha = "") +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1), 
            text = element_text(size=20)) +
        ylab("Log 10 (Binding + 1)") +
        xlab("Tested Epitope; Patient") +
        geom_vline(xintercept = lines_vector)
}

make_dotplot_lines_vector <- function(df){
    plot_epitopes <- df %>%
        use_series(ALT_EPI_SEQ) %>%
        unique %>%
        sort
    lines <- 1.5:length(plot_epitopes)
}

## validation 2 ----

make_validation_scatterplot_obj <- function(df){
    title <- make_validation_scatterplot_title(df)
    plot <- df %>% 
        ggplot(aes(x = LOG_PREDICTED_BINDING,  y = LOG_MEASURED_BINDING)) +
        geom_point(
            size = 4, 
            aes(x = LOG_PREDICTED_BINDING, 
                y = LOG_MEASURED_BINDING, 
                color = HLA_ALLELE, 
                shape = HLA_ALLELE)) +
        labs(color = "") +
        labs(shape = "") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
        ylab("LOG 10 (Measured Binding + 1)") +
        xlab("LOG 10 (Predicted Binding + 1)") +
        ggtitle(title) +
        geom_smooth(method = 'lm')
}

make_validation_scatterplot_title <- function(df){
    if(nrow(df) == 1){
        title <- ""
    } else{
        spearman_cor <- get_spearman_cor(df)
        spearman_p <- get_spearman_p_value(df)
        
        title <- str_c(
            "Spearman correlation = ",
            spearman_cor, 
            ", pvalue = ",
            spearman_p)
    }
    return(title)
}

get_spearman_cor <- function(df){
    spearman_cor <- 
        cor(
            df$LOG_PREDICTED_BINDING, 
            df$LOG_MEASURED_BINDING, 
            method = "spearman") %>% 
        round(4) %>% 
        as.character
}

get_spearman_p_value <- function(df){
    spearman_p <- 
        cor.test(
            df$LOG_PREDICTED_BINDING, 
            df$LOG_MEASURED_BINDING, 
            method = "spearman") %>% 
        use_series(p.value) %>% 
        round(4) %>% 
        as.character
}

## validation 3 ----


make_validated_submissions_boxplot_obj <- function(df, column, ylab){
    df %>% 
        ggplot(aes_string(x = "PATIENT_ID", y = column)) +
        geom_boxplot(color = "black", fill = "white", outlier.shape = NA) +
        geom_jitter(
            aes(color = team_status, shape = team_status), 
            size = 4,
            height = 0) + 
        labs(shape = "", color = "") +
        scale_color_manual(values = c("black", "red")) +
        ylab(ylab) +
        xlab("Patient") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
}

# variant counts ----

make_variant_boxplot_obj <- function(df){
    df %>% 
        ggplot(aes(x = VARIANT, y = LOG_N)) +
        geom_boxplot(color = "black", fill = "white", outlier.shape = NA) +
        geom_jitter(
            aes(color = team_status, shape = team_status), 
            size = 4,
            height = 0) + 
        labs(shape = "", color = "") +
        scale_color_manual(values = c("black", "red")) +
        ylab("Log 10(Number of variants + 1)") +
        xlab("Variant type") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
}

# variant overlap ----

make_variant_histogram_obj <- function(df){
    point_x <- get_variant_overlap_point_value(df)
    df %>% 
        ggplot(aes(x = value)) +
        geom_histogram(binwidth = .05, fill = "white", color = "black") +
        geom_point(x = point_x, y = 0, color = "red", size = 4) +
        xlim(0.5, 1.05) +
        ylab("Number of teams") +
        xlab("Median per-team overlap") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
}

get_variant_overlap_point_value <- function(df){
    df %>% 
        dplyr::filter(team_status == "Your team") %>% 
        magrittr::extract2("value")
}

# survey ----

make_survey_barchart_obj <- function(df){
    df %>% 
        ggplot() +
        geom_bar(
            stat = "identity", 
            aes(x = QUESTION_PART_2, y = fraction, color = ANSWER),
            fill = "white") +
        scale_color_manual(values = c("black", "red"), name = "Your teams answer") +
        ylab("Fraction of teams answering yes") +
        xlab("Survey question") +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size = 20)) +
        coord_flip()
}


