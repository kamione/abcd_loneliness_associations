# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(lme4)
library(ggseg)
library(patchwork)
library(ggpubr)



# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here("data", "processed", "master_preprocessed_df.rds") %>% 
    read_rds()

gmv_master_preprocessed_df <- master_preprocessed_df %>% 
    filter(
        # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1
    )
wmt_master_preprocessed_df <- master_preprocessed_df %>% 
    filter(
        # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1,
        imgincl_dmri_include == 1
    )

# total cortical volume
cgmv_model <- as.formula(
    "loneliness_followup_bin ~ smri_vol_cdk_total + loneliness_bl + \
    demo_sex_v2 + interview_age + race_ethnicity_3 + \
    (1 | mri_info_deviceserialnumber)"
)
cgmv_res <- master_preprocessed_df %>% 
    glmer(
        formula = cgmv_model, 
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
        cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
        cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
        cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE)
    ) %>% 
    mutate(name = "1. Cortical GMV", .before = "Parameter")

# total subcortical volume
scgmv_model <- as.formula(
    "loneliness_followup_bin ~ smri_vol_scs_subcorticalgv + loneliness_bl + \
    demo_sex_v2 + interview_age + race_ethnicity_3 + \
    (1 | mri_info_deviceserialnumber)"
)
scgmv_res <- master_preprocessed_df %>% 
    glmer(
        formula = scgmv_model, 
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
        cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
        cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
        cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE)
    ) %>% 
    mutate(name = "2. Subcortical GMV", .before = "Parameter")

# DTI fa
wmt_fa_model <- as.formula(
    "loneliness_followup_bin ~ dmdtifp1_38 + loneliness_bl + demo_sex_v2 + \
    interview_age + race_ethnicity_3 + dmri_meanmotion + \
    (1 | mri_info_deviceserialnumber)"
)
wmt_fa_res <- wmt_master_preprocessed_df %>% 
    glmer(
        formula = wmt_fa_model, 
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
        cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
        cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
        cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE)
    ) %>% 
    mutate(name = "3. All WM tract FA", .before = "Parameter")

wmt_md_model <- as.formula(
    "loneliness_followup_bin ~ dmdtifp1_80 + loneliness_bl + demo_sex_v2 + \
    interview_age + race_ethnicity_3 + dmri_meanmotion + \
    (1 | mri_info_deviceserialnumber)"
)
wmt_md_res <- wmt_master_preprocessed_df %>% 
    glmer(
        formula = wmt_md_model, 
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
        cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
        cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
        cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE)
    ) %>% 
    mutate(name = "4. All WM tract MD", .before = "Parameter")

wmt_ad_model <- as.formula(
    "loneliness_followup_bin ~ dmdtifp1_122 + loneliness_bl + demo_sex_v2 + \
    interview_age + race_ethnicity_3 + dmri_meanmotion + \
    (1 | mri_info_deviceserialnumber)"
)
wmt_ad_res <- wmt_master_preprocessed_df %>% 
    glmer(
        formula = wmt_ad_model, 
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
        cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
        cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
        cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE)
    ) %>% 
    mutate(name = "5. All WM tract AD", .before = "Parameter")

wmt_rd_model <- as.formula(
    "loneliness_followup_bin ~ dmdtifp1_164 + loneliness_bl + demo_sex_v2 + \
    interview_age + race_ethnicity_3 + dmri_meanmotion + \
    (1 | mri_info_deviceserialnumber)"
)
wmt_rd_res <- wmt_master_preprocessed_df %>% 
    glmer(
        formula = wmt_rd_model, 
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
        cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
        cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
        cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE)
    ) %>% 
    mutate(name = "6. All WM tract RD", .before = "Parameter")

brain_effectsize <- bind_rows(
    cgmv_res, scgmv_res, wmt_fa_res, wmt_md_res,
    wmt_ad_res, wmt_rd_res
) %>% 
    mutate(
        p.adj = p.adjust(p, method = "fdr"),
        sig = if_else(p.adj < 0.05, "Yes", "No"),
        name = factor(name)
    ) %>% 
    ggplot(aes(x = cohend, y = name, color = sig)) +
        geom_point(size = 4) +
        geom_errorbar(
            aes(xmin = cohend_low, xmax = cohend_high, color = sig),
            linewidth = 1,
            width = 0.3
        ) +
        scale_color_manual(values = c("grey60", "tomato3")) +
        geom_vline(xintercept = 0) +
        scale_y_discrete(limits = rev) +
        labs(x = "Effect Size (Cohen's d)", y = "") +
        ggthemes::theme_pander() +
        guides(color = "none")

brain_effectsize

ggexport(
    brain_effectsize,
    filename = here("outputs", "figs", "brain_effectsize_followup.pdf"),
    width = 5,
    height = 5
)
