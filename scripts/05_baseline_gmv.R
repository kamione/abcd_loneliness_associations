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
        "{roi} ~ loneliness_bl + demo_sex_v2 + interview_age + \
         race_ethnicity_3 + smri_vol_scs_intracranialv + \
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
            variable_name = roi,
            label = gmv_ggseg_roinames[ith],
            cohend = effectsize::t_to_d(t, df_error)$d,
            cohend_low = effectsize::t_to_d(t, df_error)$CI_low,
            cohend_high = effectsize::t_to_d(t, df_error)$CI_high,
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
    min = -0.06, 
    max = 0.06, 
    break_int = 0.03
)
cgmv_fig

scgmv_res <- reduce(res_list, bind_rows) %>% 
    slice(69:82) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    select(label, cohend, cohend_low, cohend_high, p, p.adj) %>% 
    drop_na(label) %>% 
    filter(p < 0.05) %>% 
    right_join(as_tibble(aseg), by = "label")

scgmv_fig <- plot_ggseg_brain(
    dat = scgmv_res,
    atlas = "aseg", 
    fill = "cohend", 
    min = -0.06, 
    max = 0.06, 
    break_int = 0.03
)
scgmv_fig
gmv_fig <- cgmv_fig + scgmv_fig + plot_layout(widths = c(5, 2), guides = "collect")
gmv_fig


# Save -------------------------------------------------------------------------
# result tables
reduce(res_list, bind_rows) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    select(variable_name, label, Coefficient, t, df_error, cohend, cohend_low, cohend_high, p, 
           p.adj) %>% 
    rename(
        "beta" = "Coefficient",
        "cohen's d" = "cohend", 
        "95% CI low cohen's d" = "cohend_low",
        "95% CI high cohen's d" = "cohend_high"
    ) %>% 
    write_csv(here("outputs", "tables", "lmm_results_baseline_gmv.csv"))

# result figure
ggpubr::ggexport(
    gmv_fig, 
    filename = here("outputs", "figs", "gmv.pdf"), 
    width = 6, 
    height = 4
)

