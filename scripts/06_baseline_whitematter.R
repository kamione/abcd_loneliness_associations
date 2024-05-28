# Environment ------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(ggseg)
library(patchwork)
library(ggsegJHU)

# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) %>% 
    read_rds() %>% 
    filter(
        # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1,
        imgincl_dmri_include == 1
    ) # 6,612

wmt_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "White Matter Tracts") %>% 
    filter(!variable_name %in% c("dmdtifp1_38", "dmdtifp1_80", "dmdtifp1_132", "dmdtifp1_164"))
wmt_roinames <- pull(wmt_variables_df, variable_name)
wmt_label <- pull(wmt_variables_df, region)
wmt_hemi <- pull(wmt_variables_df, hemi)

res_list <- list()
for (ith in 1:length(wmt_roinames)) {
    tract <- wmt_roinames[ith]
    f1 <- glue::glue(
        "{tract} ~ loneliness_bl + demo_sex_v2 + interview_age + \
        race_ethnicity_3 + dmri_meanmotion + \
        (1 | mri_info_deviceserialnumber)"
    ) %>% 
        as.formula()
    
    res <- master_preprocessed_df %>% 
        lmer(formula = f1) %>% 
        parameters::model_parameters(
            standardize = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = wmt_roinames[ith],
            hemi = wmt_hemi[ith],
            region = wmt_label[ith],
            cohend = effectsize::t_to_d(t, df_error)$d,
            cohend_low = effectsize::t_to_d(t, df_error)$CI_low,
            cohend_high = effectsize::t_to_d(t, df_error)$CI_high,
            .before = Parameter
        )
    res_list[[ith]] <- res
}

fa_wmt_res <- reduce(res_list, bind_rows) %>% 
    mutate(
        p.adj = p.adjust(p, method = "fdr"),
        modality = c(rep("fa", 20), rep("md", 20), rep("ad", 20), rep("rd", 20))
    ) %>% 
    filter(p < 0.05) %>% 
    select(modality, hemi, region, cohend, p, p.adj)

# save table
reduce(res_list, bind_rows) %>% 
    mutate(
        p.adj = p.adjust(p, method = "fdr"),
        modality = c(rep("fa", 20), rep("md", 20), rep("ad", 20), rep("rd", 20))
    ) %>% 
    select(variable_name, hemi, region, Coefficient, t, df_error, cohend, 
           cohend_low, cohend_high, p, p.adj) %>% 
    rename(
        "beta" = "Coefficient",
        "cohen's d" = "cohend", 
        "95% CI low cohen's d" = "cohend_low",
        "95% CI high cohen's d" = "cohend_high"
    ) %>% 
    write_csv(here("outputs", "tables", "lmm_results_baseline_wmt.csv"))
