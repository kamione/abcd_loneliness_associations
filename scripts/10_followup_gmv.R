# Environment ------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(ggseg)
library(patchwork)

source(here("src", "R", "plot_ggseg_brain.R"))

# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) %>% 
    read_rds() %>% 
    filter(
        # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1
    )

gmv_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Gray Matter Volume") %>% 
    filter(!variable_name %in% c(
        "smri_vol_scs_intracranialv", 
        "smri_vol_cdk_total", 
        "smri_vol_scs_subcorticalgv"
    ))
gmv_roinames <- pull(gmv_variables_df, variable_name)
gmv_ggseg_roinames <- pull(gmv_variables_df, ggseg_label)


res_list <- list()
for (ith in 1:length(gmv_roinames)) {
    roi <- gmv_roinames[ith]
    f1 <- glue::glue(
        "loneliness_followup_bin ~ {roi} + loneliness_bl + demo_sex_v2 + \
         interview_age + race_ethnicity_3 + smri_vol_scs_intracranialv + \
         (1 | mri_info_deviceserialnumber)"
    ) %>% 
        as.formula()
    
    res <- master_preprocessed_df %>% 
        glmer(
            formula = f1, 
            family = binomial(link = "logit"),
            control = glmerControl(tolPwrss = 1e-3)
        ) %>% 
        parameters::model_parameters(
            standardize = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            variable_name = roi,
            label = gmv_ggseg_roinames[ith],
            cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
            cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
            cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE),
            .before = Parameter
        ) 
    res_list[[ith]] <- res
}

cgmv_res <- reduce(res_list, bind_rows) %>% 
    slice(1:68) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    select(label, cohend, p, p.adj) %>% 
    filter(p < 0.05) %>% 
    right_join(as_tibble(dk), by = "label")

cgmv_fig <- plot_ggseg_brain(
    dat = cgmv_res,
    atlas = "dk", 
    fill = "cohend", 
    min = -0.08, 
    max = 0.08, 
    break_int = 0.04
)
cgmv_fig

# no significant results
scgmv_res <- reduce(res_list, bind_rows) %>% 
    slice(69:82) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    select(label, cohend, cohend_low, cohend_high, p, p.adj) %>% 
    drop_na(label) %>% 
    filter(p < 0.05)

# Save -------------------------------------------------------------------------
reduce(res_list, bind_rows) %>% 
    slice(1:68) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    select(variable_name, label, Coefficient, z, cohend, cohend_low, 
           cohend_high, p, p.adj) %>% 
    rename(
        "beta" = "Coefficient",
        "cohen's d" = "cohend", 
        "95% CI low cohen's d" = "cohend_low",
        "95% CI high cohen's d" = "cohend_high"
    ) %>% 
    write_csv(here("outputs", "tables", "lmm_results_followup_gmv.csv"))

ggpubr::ggexport(
    cgmv_fig, 
    filename = here("outputs", "figs", "cgmv_followup.pdf"), 
    width = 8, 
    height = 4
)

