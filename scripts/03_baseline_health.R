# Environment ------------------------------------------------------------------
library(here)
library(tidyverse)
library(lme4)
library(ggplot2)
library(ggthemes)
library(gghighlight)
library(ggrepel)

# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) %>% 
    read_rds()

health_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Health")
health_variables <- master_preprocessed_df %>% 
    select(any_of(pull(health_variables_df, variable_name))) %>% 
    colnames()

source(here("src", "R", "is_binary_column.R"))


res_list <- list()
for (ith in 1:(length(health_variables))) {
    health_var <- health_variables[ith]
    
    f1 <- glue::glue(
        "{health_var} ~ loneliness_bl + demo_sex_v2 + interview_age + \
         race_ethnicity_3 + (1 | site_id_l)"
    ) %>% 
        as.formula()
    
    is_col_bin <- is_binary_column(master_preprocessed_df[[health_var]])
    
    if (is_col_bin) {
        res <- master_preprocessed_df %>% 
            glmer(formula = f1, family = binomial(link = "logit")) %>% 
            parameters::model_parameters(
                method = "refit", 
                effect = "fixed", 
                verbose = FALSE
            ) %>%
            as_tibble() %>% 
            slice(2) %>% 
            mutate(
                cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE), 
                variable_name = health_var,
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
                variable_name = health_var,
                .before = "p"
            )
    }
    
    res_list[[ith]] <- res
}

res_df <- reduce(res_list, bind_rows) %>% 
    mutate(p.adj = p.adjust(p, method = "holm")) %>% 
    left_join(health_variables_df, by = "variable_name") %>% 
    arrange(category) %>% 
    mutate(x = factor(1:n()))

# save results to csv
res_df %>% 
    select(category, table_name, variable_name, variable_label, Coefficient,
           t, df_error, cohend, cohend_low, cohend_high, p, p.adj) %>% 
    rename(
        "beta" = "Coefficient",
        "cohen's d" = "cohend", 
        "95% CI low cohen's d" = "cohend_low",
        "95% CI high cohen's d" = "cohend_high"
    ) %>% 
    write_csv(here("outputs", "tables", "lmm_results_baseline_health.csv"))


# Visualization ----------------------------------------------------------------
health_fig <- res_df %>% 
    ggplot(aes(x = x, y = cohend, color = category)) +
        geom_point(size = 4, alpha = 0.75) +
        labs(x = "", y = "Effect Size (Cohen's d)") +
        scale_x_discrete(expand = c(0.05, 0.05)) +
        scale_y_continuous(limits = c(-0.85, 1.3),
                           breaks = c(-0.8, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 0.8)) +
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
            data = res_df %>% 
                subset(p.adj < 0.05 & category == "Mental Health") %>% 
                arrange(desc(cohend)) %>% 
                slice(1:8) %>% 
                mutate(label = c(
                    "social", "total", "anxious/depressed", "internalizing", 
                    "aggressive", "withdrawn/depressed", "externalizing",
                    "thought"
                )),
            nudge_x = c(5, 4, 4, 8, 6, -2, 5, 0),
            nudge_y = c(-0.05, 0.1, 0.2, 0, 0.13, 0.1, -0.17, -0.27),
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            mapping = aes(label = label),
            data = res_df %>% 
                subset(p.adj < 0.05 & category == "Neurocognition") %>% 
                arrange(desc(abs(cohend))) %>% 
                slice(1:5) %>% 
                mutate(label = c(
                    "oral reading recognition", "total cognition", "fluid", 
                    "long-delayed recall", "short-delayed recall"
                )),
            nudge_x = c(0, -4, -1, -1, -3),
            nudge_y = c(0.2, -0.23, -0.15, -0.3, 0.2),
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        geom_text_repel(
            data = res_df %>% 
                subset(p.adj < 0.05 & category == "Physical Health") %>% 
                arrange(desc(cohend)) %>% 
                slice(1, 7, 9, 10, 11) %>% 
                mutate(label = c(
                    "sleep diificulties", "total medical problems", 
                    "puberty", "waist", "BMI"
                )),
            aes(label = label),
            nudge_x = c(-2, 2, 2, -4, -4),
            nudge_y = c(0.15, -0.32, 0.2, -0.18, -0.05),
            hjust = 0.5,
            segment.size = 0.25,
            show.legend = FALSE
        ) +
        ggforce::geom_mark_ellipse(
            aes(group = category, filter = str_detect(variable_name, "sds")),
            expand = unit(2, "mm"),
            show.legend = FALSE
        )
health_fig

ggpubr::ggexport(
    health_fig, 
    filename = here("outputs", "figs", "health.pdf"), 
    width = 7, 
    height = 5
)
