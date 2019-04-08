make_submissions_boxplot <- function(df){
    plot <- df %>% 
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
    print(plot)
}

make_agretopicity_boxplot <- function(df){
    plot <- df %>% 
        ggplot(aes(x = TEAM,  y = LOG_AGRETOPICITY)) +
        geom_boxplot(aes(fill = team_status)) +
        scale_fill_manual(values=c("white", "red")) +
        labs(fill = "") +
        xlab("Team") +
        ylab("Log10 Agretopicity Index") +
        theme_bw() +
        theme(axis.text.x = element_blank())
    print(plot)
}

make_pep_length_barchart <- function(df){
    plot <- df %>% 
        ggplot(aes(x = TEAM)) +
        geom_bar(aes(fill = as.factor(PEP_LEN), color = team_status), size = 2) +
        labs(color = 'Bar outline(Team)') +
        labs(fill  = 'Bar color(Peptide Length)') +
        scale_color_manual(values = c("black", "red")) +
        xlab("Team") +
        ylab("Peptide Length Counts") +
        theme_bw() +
        theme(axis.text.x = element_blank()) 
    print(plot)
}

make_epitope_overlap_boxplot <- function(df){
    plot <- df %>% 
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
    print(plot)
}



make_validation_dotplot <- function(df){
    lines_vector <- make_dotplot_lines_vector(df)
    plot_object <- make_validation_dotplot_obj(df, lines_vector)
    print(plot_object)
}

make_dotplot_lines_vector <- function(df){
    plot_epitopes <- df %>%
        use_series(ALT_EPI_SEQ) %>%
        unique %>%
        sort
    lines <- 1.5:length(plot_epitopes)
}


make_validation_dotplot_obj <- function(df, lines_vector){
    plot <- df %>% 
        ggplot(aes(x = ALT_EPI_SEQ,  y = LOG_BINDING)) +
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
        xlab("Tested Epitope") +
        geom_vline(xintercept = lines_vector)
}

make_validation_scatterplot <- function(df){
    title <- make_validation_dotplot_title(df)
    plot_object <- make_validation_scatterplot_obj(df, title)
    print(plot_object)
}

make_validation_dotplot_title <- function(df){
    spearman_cor <- 
        cor(
            df$LOG_PREDICTED_BINDING, 
            df$LOG_MEASURED_BINDING, 
            method = "spearman") %>% 
        round(4) %>% 
        as.character
    
    spearman_p <- 
        cor.test(
            df$LOG_PREDICTED_BINDING, 
            df$LOG_MEASURED_BINDING, 
            method = "spearman") %>% 
        use_series(p.value) %>% 
        round(4) %>% 
        as.character
    
    title <- str_c(
        "Spearman correlation = ",
        spearman_cor, 
        ", pvalue = ",
        spearman_p)
}

make_validation_scatterplot_obj <- function(df, title){
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


make_validated_submissions_boxplot <- function(df, ylab){
    plot <- df %>% 
        ggplot(aes_string(x = "PATIENT_ID", y = "STAT")) +
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
    print(plot)
}