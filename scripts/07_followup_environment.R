# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(lme4)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(gghighlight)
library(ggrepel)

source(here("src", "R", "is_binary_column.R"))


# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) %>% 
    read_rds() %>% 
    mutate(demo_comb_income_v2 = as.numeric(demo_comb_income_v2)) %>% 
    drop_na(loneliness_followup_bin)

exposome_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Exposome")
exposome_variables <- master_preprocessed_df %>% 
    select(any_of(pull(exposome_variables_df, variable_name))) %>% 
    colnames()

res_list <- list()
for (ith in 1:(length(exposome_variables))) {
    
    env <- exposome_variables[ith]
    
    if (is_binary_column(master_preprocessed_df[[env]]) == FALSE) {
        master_preprocessed_df[[env]] <- master_preprocessed_df[[env]] %>% 
            scale() %>% 
            as.numeric()
    }
    
    f1 <- glue::glue(
        "loneliness_followup_bin ~ {env} + loneliness_bl + demo_sex_v2 + \
        interview_age + race_ethnicity_3 + (1 | site_id_l)"
    ) %>% 
        as.formula()
    
    res <- master_preprocessed_df %>% 
        glmer(
            formula = f1, 
            family = binomial(link = "logit"),
            control = glmerControl(tolPwrss=1e-3)
        ) %>% 
        parameters::model_parameters(
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        dplyr::slice(2) %>% 
        mutate(
            cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
            cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
            cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE),
            variable_name = env,
            .before = "p"
        )
    
    res_list[[ith]] <- res

}


exposome_res <- reduce(res_list, bind_rows) %>% 
    mutate(p.adj = p.adjust(p, method = "holm")) %>% 
    left_join(exposome_variables_df, by = "variable_name") %>% 
    arrange(category) %>% 
    mutate(x = factor(1:n()))

exposome_res %>% 
    mutate(p_log = -log10(p)) %>% 
    ggplot(aes(x = x, y = p_log, color = category)) +
    geom_point(size = 4, alpha = 0.75) +
    labs(x = "", y = "-log10") +
    gghighlight(
        p.adj < 0.05,
        unhighlighted_params = list(size = 3, colour = NULL, alpha = 0.1), 
        use_direct_label = FALSE,
        keep_scales = TRUE
    ) +
    theme_clean() +
    scale_y_continuous(limits = c(-5, 50)) +
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(colour = NA, fill = NA),
        legend.position = c(0.5, 0.05),
        legend.direction = "horizontal",
        legend.key.size = unit(1, "mm"),
        plot.background = element_rect(fill = NA, color = NA)
    )

exposome_fig <- exposome_res %>% 
    ggplot(aes(x = x, y = cohend, color = category)) +
        geom_point(size = 4, alpha = 0.75) +
        labs(x = "", y = "Effect Size (Cohen's d)") +
        scale_x_discrete(expand = c(0.05, 0)) +
        scale_y_continuous(limits = c(-0.38, 0.38),
                           breaks = c(-0.2, -0.1, 0, 0.1, 0.2)) +
        theme_clean() +
        theme(
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.line.y = element_blank(), 
            axis.line.x = element_blank(),
            legend.title = element_blank(),
            legend.background = element_rect(colour = NA, fill = NA),
            legend.position = c(0.5, 0.05),
            legend.direction = "horizontal",
            legend.key.size = unit(1, "mm"),
            plot.background = element_rect(fill = NA, color = NA)
        ) + 
        gghighlight(
            p.adj < 0.05,
            unhighlighted_params = list(size = 2, color = "grey70"), 
            use_direct_label=FALSE,
            keep_scales = TRUE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Developmental History") %>% 
                arrange(desc(abs(cohend))) %>% 
                slice(1:3) %>% 
                mutate(label = c(
                    "having a twin", 
                    "planned pregnancy"
                )),
            nudge_x = c(-10, -20),
            nudge_y = 0,
            hjust = 1,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Parental Psychopathogy") %>% 
                arrange(desc(abs(cohend))) %>% 
                slice(1:7) %>% 
                mutate(label = c(
                    "total",
                    "internalizing",
                    "externalizing", 
                    "attention", 
                    "anxious/depressed", 
                    "thought",
                    "withdrawn"
                )),
            nudge_x = c(-20, 35, -30, -30, -30, 40, -35),
            nudge_y = c(0.05, -0.05, 0.1, 0.05, 0, 0, -0.05),
            hjust = 1,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Life Adversity") %>% 
                arrange(desc(abs(cohend))) %>%
                mutate(label = "total traumatic events"),
            nudge_y = 0.07,
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "School Environment") %>% 
                arrange(desc(abs(cohend))) %>%
                slice(1:3) %>% 
                mutate(label = c(
                    "school environment", 
                    "school involvement",
                    "school disengagement"
            )),
            nudge_x = -10,
            nudge_y = c(-0.17, -0.1, -0.25),
            hjust = 1,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Family and Social Community") %>% 
                arrange(desc(abs(cohend))) %>%
                slice(1:3) %>% 
                mutate(label = c(
                    "family conflict", 
                    "family support", 
                    "family as referent"
                )),
            nudge_x = c(10, 10, 45),
            nudge_y = c(0.1, -0.08, -0.02),
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Neighbourhood and Built Environment") %>% 
                arrange(desc(abs(cohend))) %>%
                slice(1:2) %>% 
                mutate(label = c(
                    "safe from crime", 
                    "violence not a problem"
                )),
            nudge_x = -70,
            nudge_y = c(-0.1, -0.2),
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        )
exposome_fig


# Save -------------------------------------------------------------------------
ggpubr::ggexport(
    exposome_fig, 
    filename = here("outputs", "figs", "exposome_longitudinal.pdf"), 
    width = 10, 
    height = 5
)

