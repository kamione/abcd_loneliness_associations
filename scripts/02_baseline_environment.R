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
    mutate(demo_comb_income_v2 = as.numeric(demo_comb_income_v2))

exposome_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Exposome")
exposome_variables <- master_preprocessed_df %>% 
    select(any_of(pull(exposome_variables_df, variable_name))) %>% 
    colnames()

res_list <- list()
for (ith in 1:(length(exposome_variables))) {
    env <- exposome_variables[ith]
    
    f1 <- glue::glue(
        "{env} ~ loneliness_bl + demo_sex_v2 + interview_age + \
        race_ethnicity_3 + (1 | site_id_l)"
    ) %>% 
        as.formula()
    
    if (is_binary_column(master_preprocessed_df[[env]])) {
        res <- master_preprocessed_df %>% 
            glmer(
                formula = f1, 
                family = binomial(link = "logit"),
                control = glmerControl(tolPwrss=1e-3)
            ) %>% 
            parameters::model_parameters(
                method = "refit",
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
    } else {
        res <- master_preprocessed_df %>% 
            lmer(formula = f1) %>% 
            parameters::model_parameters(
                method = "refit", 
                effect = "fixed", 
                verbose = FALSE
            ) %>%
            as_tibble() %>% 
            slice(2) %>% 
            mutate(
                cohend = effectsize::t_to_d(t, df_error)$d, 
                cohend_low = effectsize::t_to_d(t, df_error)$CI_low,
                cohend_high = effectsize::t_to_d(t, df_error)$CI_high,
                variable_name = env,
                .before = "p"
            )
    }
    
    res_list[[ith]] <- res
}

exposome_res <- reduce(res_list, bind_rows) %>% 
    mutate(p.adj = p.adjust(p, method = "holm")) %>% 
    left_join(exposome_variables_df, by = "variable_name") %>% 
    arrange(category) %>% 
    mutate(x = factor(1:n()))

# save results to csv
exposome_res %>% 
    select(category, table_name, variable_name, variable_label, Coefficient,
           t, df_error, cohend, cohend_low, cohend_high, p, p.adj) %>% 
    rename(
        "beta" = "Coefficient",
        "cohen's d" = "cohend", 
        "95% CI low cohen's d" = "cohend_low",
        "95% CI high cohen's d" = "cohend_high"
    ) %>% 
    write_csv(here("outputs", "tables", "lmm_results_baseline_environment.csv"))

# Visualization ----------------------------------------------------------------
exposome_fig <- exposome_res %>% 
    ggplot(aes(x = x, y = cohend, color = category)) +
        geom_point(size = 4, alpha = 0.75) +
        labs(x = "", y = "Effect Size (Cohen's d)") +
        scale_x_discrete(expand = c(0.05, 0.05)) +
        scale_y_continuous(limits = c(-0.7, 0.7),
                           breaks = c(-0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5)) +
        theme_clean() +
        theme(
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.line.y = element_blank(), 
            axis.line.x = element_blank(),
            legend.title = element_blank(),
            legend.background = element_rect(colour = NA, fill = NA),
            legend.position = c(0.5, 0.1),
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
                    "bedwetting", 
                    "planned pregnancy", 
                    "tobacco per days"
                )),
            nudge_x = c(-10, -20, -25),
            nudge_y = c(0.2, -0.1, 0.15),
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
                    "total problems",
                    "internalizing",
                    "anxious/depressed", 
                    "withdrawn", 
                    "externalizing",
                    "attention", 
                    "thought"
                )),
            nudge_x = c(-20, 35, -30, -30, -30, 50, -35),
            nudge_y = c(0.06, -0.05, 0.1, 0.05, 0, 0, -0.05),
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
                subset(p.adj < 0.05 & category == "Laws and Bias") %>% 
                arrange(desc(abs(cohend))) %>%
                mutate(label = "CBD/THC law became effective"),
            nudge_y = -0.06,
            nudge_x = 30,
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Socioeconomic Status") %>% 
                arrange(desc(abs(cohend))) %>%
                mutate(label = c("income", "has partner", "parental education")),
            nudge_y = c(-0.1, 0.15, -0.3),
            hjust = 0,
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
                    "educational attainment"
            )),
            nudge_x = c(-20, -20, -10),
            nudge_y = c(-0.06, -0.11, -0.2),
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
                    "parental monitoring",
                    "caregiver acceptance"
                )),
            nudge_x = c(10, -10, 30),
            nudge_y = c(0.2, -0.16, -0.12),
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Community Health Care") %>% 
                arrange(desc(abs(cohend))) %>%
                slice(1:3) %>% 
                mutate(label = c("mental health", "smoking", "colorectal screening")),
            nudge_x = -20,
            nudge_y = c(0.08, 0, -0.1),
            hjust = 1,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Neighbourhood Socioeconomic Status") %>% 
                arrange(desc(abs(cohend))) %>%
                slice(1:4) %>% 
                mutate(label = c(
                    "neighborhood affluence",
                    "income greater than 75K proportion", 
                    "income greater than 100K proportion", 
                    "college proportion"
            )),
            nudge_x = 20,
            nudge_y = c(-0.2, -0.28, -0.32, -0.24),
            hjust = 1,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = exposome_res %>% 
                subset(p.adj < 0.05 & category == "Neighbourhood and Built Environment") %>% 
                arrange(desc(abs(cohend))) %>%
                slice(1:3) %>% 
                mutate(label = c(
                    "safe walking", 
                    "safe from crime", 
                    "violence not a problem"
                )),
            nudge_x = c(-40, 30, -5),
            nudge_y = c(0.42, 0.4, 0.45),
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        )
exposome_fig

ggpubr::ggexport(
    exposome_fig, 
    filename = here("outputs", "figs", "exposome.pdf"), 
    width = 10, 
    height = 5
)