# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(lme4)
library(ggplot2)
library(gghighlight)
library(ggthemes)


# Data I/O ---------------------------------------------------------------------
# ABCD data path
# abcd_path <- "/path/to/abcd/folder"
abcd_folder_path <- "/Users/tywong/GoogleDrive/opendata/abcd/abcd-data-release-5.1/core"

master_preprocessed_df <- here("data", "processed", "master_preprocessed_df.rds") %>% 
    read_rds()

health_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Health")

mh_vars <- health_variables_df %>% 
    filter(category == "Mental Health") %>% 
    pull(variable_name)
mh_vars <- mh_vars[1:11]

mh_vars_label <- c(
    "01. Anxious/Depressed", "02. Withdrawn/Depressed", "03. Somatic", 
    "04. Social", "05. Thought",  "06. Attention", "07.Rule-breaking", 
    "08. Aggressive", "09. Internalizing", "10. Externalizing", 
    "11. Total Problems"
)

mh_y3_df <- here(abcd_folder_path, "mental-health", "mh_p_cbcl.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(src_subject_id, eventname, any_of(mh_vars)) %>% 
    filter(eventname == "3_year_follow_up_y_arm_1") %>% 
    rename_at(vars(-(1:2)), ~ paste0(., '_3yfu')) %>% 
    select(-eventname)

here(abcd_folder_path, "mental-health", "mh_p_cbcl.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(src_subject_id, eventname, any_of(mh_vars)) %>% 
    filter(eventname == "3_year_follow_up_y_arm_1") %>% 
    rename_at(vars(-(1:2)), ~ paste0(., '_3yfu')) %>% 
    select(-eventname)


mh_2ymean_df <- here(abcd_folder_path, "mental-health", "mh_p_cbcl.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(src_subject_id, eventname, any_of(mh_vars)) %>% 
    filter(eventname %in% c(
        "baseline_year_1_arm_1",
        "1_year_follow_up_y_arm_1",
        "2_year_follow_up_y_arm_1"
    )) %>% 
    group_by(src_subject_id) %>% 
    summarize_at(
        vars(cbcl_scr_syn_anxdep_t:cbcl_scr_syn_totprob_t),
        mean, 
        na.rm = TRUE
    ) %>% 
    rename_at(vars(-(1)), ~ paste0(., '_2ymean')) %>% 
    ungroup()


master_preprocessed_3yfu <- master_preprocessed_df %>% 
    left_join(mh_y3_df, by = "src_subject_id") %>% 
    left_join(mh_2ymean_df, by = "src_subject_id") %>%


res_list <- list()
for (ith in 1:(length(mh_vars))) {
    mh_var <- mh_vars[ith]
    
    f1 <- glue::glue(
        "{mh_var}_3yfu ~ loneliness_2y + cbcl_scr_syn_totprob_t_2ymean + demo_sex_v2 + \
        interview_age + race_ethnicity + (1 | site_id_l)"
    ) %>% 
        as.formula()
    f2 <- glue::glue(
        "{mh_var}_3yfu ~ loneliness_2y + cbcl_scr_syn_totprob_t_2ymean + demo_sex_v2 + \
        interview_age + race_ethnicity + demo_comb_income_v2 + (1 | site_id_l)"
    ) %>% 
        as.formula()
    f3 <- glue::glue(
        "{mh_var}_3yfu ~ loneliness_2y + cbcl_scr_syn_totprob_t_2ymean + demo_sex_v2 + \
        interview_age + race_ethnicity + demo_comb_income_v2 + \
        asr_scr_totprob_t + (1 | site_id_l)"
    ) %>% 
        as.formula()
    f4 <- glue::glue(
        "{mh_var}_3yfu ~ loneliness_2y + cbcl_scr_syn_totprob_t_2ymean + demo_sex_v2 + \
        interview_age + race_ethnicity + demo_comb_income_v2 + ksads_ptsd_sum +\
        asr_scr_totprob_t + (1 | site_id_l)"
    ) %>% 
        as.formula()
    f5 <- glue::glue(
        "{mh_var}_3yfu ~ loneliness_2y + cbcl_scr_syn_totprob_t_2ymean + demo_sex_v2 + \
        interview_age + race_ethnicity + demo_comb_income_v2 + ksads_ptsd_sum +\
        asr_scr_totprob_t + fes_p_ss_fc + (1 | site_id_l)"
    ) %>% 
        as.formula()
    f6 <- glue::glue(
        "{mh_var}_3yfu ~ loneliness_2y + cbcl_scr_syn_totprob_t_2ymean + demo_sex_v2 + \
        interview_age + race_ethnicity + demo_comb_income_v2 + ksads_ptsd_sum +\
        asr_scr_totprob_t + fes_p_ss_fc + srpf_y_ss_iiss + srpf_y_ss_ses + 
        (1 | site_id_l)"
    ) %>% 
        as.formula()
    
    if (ith %in% 1:8) {
        
    }
    
    
    f1_res <- master_preprocessed_3yfu %>% 
        lmer(formula = f1) %>% 
        parameters::model_parameters(
            method = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = mh_var,
            variable_label = mh_vars_label[ith]
        ) %>% 
        mutate(
            model = "model 1",
            cohend = effectsize::t_to_d(t, df_error), .before = "p"
        )
    f2_res <- master_preprocessed_3yfu %>% 
        lmer(formula = f2) %>% 
        parameters::model_parameters(
            method = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = mh_var,
            variable_label = mh_vars_label[ith]
        ) %>% 
        mutate(
            model = "model 2: model 1 + income",
            cohend = effectsize::t_to_d(t, df_error), .before = "p"
        ) 
    f3_res <- master_preprocessed_3yfu %>% 
        lmer(formula = f3) %>% 
        parameters::model_parameters(
            method = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = mh_var,
            variable_label = mh_vars_label[ith]
        ) %>% 
        mutate(
            model = "model 3: model 2 + parental psychopathology",
            cohend = effectsize::t_to_d(t, df_error), .before = "p"
        ) 
    f4_res <- master_preprocessed_3yfu %>% 
        lmer(formula = f4) %>% 
        parameters::model_parameters(
            method = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = mh_var,
            variable_label = mh_vars_label[ith]
        ) %>% 
        mutate(
            model = "model 4: model 3 + total traumatic events",
            cohend = effectsize::t_to_d(t, df_error), 
            .before = "p"
        ) 
    f5_res <- master_preprocessed_3yfu %>% 
        lmer(formula = f5) %>% 
        parameters::model_parameters(
            method = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = mh_var,
            variable_label = mh_vars_label[ith]
        ) %>% 
        mutate(
            model = "model 5: model 4 + family conflict",
            cohend = effectsize::t_to_d(t, df_error), 
            .before = "p"
        )
    f6_res <- master_preprocessed_3yfu %>% 
        lmer(formula = f6) %>% 
        parameters::model_parameters(
            method = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = mh_var,
            variable_label = mh_vars_label[ith]
        ) %>% 
        mutate(
            model = "model 6: model 5 + school environment",
            cohend = effectsize::t_to_d(t, df_error), 
            .before = "p"
        ) 
    
    res <- bind_rows(f1_res, f2_res, f3_res, f4_res, f5_res, f6_res) %>% 
        select(model, variable_name, -Parameter, everything())

    res_list[[ith]] <- res
}

mh_long_res <- reduce(res_list, bind_rows)

fu_mh_fig <- mh_long_res %>% 
    group_by(model) %>% 
    mutate(p.adj = p.adjust(p, method = "holm")) %>% 
    ungroup() %>% 
    mutate(
        p.adj = p.adjust(p, method = "fdr"),
        sig = if_else(p.adj < 0.05, "Yes", "No"),
        name = factor(variable_label)
    ) %>% 
    ggplot(aes(
        x = cohend$d, y = name, color = model
    )) +
    geom_point(size = 4) +
    gghighlight(
        p.adj < 0.05,
        unhighlighted_params = list(size = 3.5, color = "grey70"), 
        use_direct_label=FALSE,
        keep_scales = TRUE
    ) +
    geom_vline(xintercept = 0) +
    scale_y_discrete(limits = rev) +
    labs(x = "Effect Size (Cohen's d)", y = "", color = "") +
    theme_pander() +
    theme(axis.text.y = element_text(hjust = 0))

fu_mh_fig

# Save -------------------------------------------------------------------------
ggpubr::ggexport(
    fu_mh_fig, 
    filename = here("outputs", "figs", "fu_mh_fig.pdf"), 
    width = 8, 
    height = 4
)