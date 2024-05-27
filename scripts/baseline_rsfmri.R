# Environment ------------------------------------------------------------------
library(here)
library(tidyverse)
library(igraph)
library(ggraph)
library(ggsegGordon)
library(lme4)
library(ggpubr)


# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here("data", "processed", "master_preprocessed_df.rds") %>% 
    read_rds() %>% 
    filter(
        # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1,
        rsfmri_ntpoints > 375 & imgincl_rsfmri_include == 1
    ) # 6,120


rsfmri_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Resting-state FC")
    
fc_connnames <- rsfmri_variables_df %>% 
    pull(variable_name) # 306 FC

res_list <- list()
for (ith in 1:length(fc_connnames)) {
    fc <- fc_connnames[ith]
    f1 <- glue::glue(
        "{fc} ~ loneliness_bl + demo_sex_v2 + interview_age + \
        race_ethnicity_3 + rsfmri_meanmotion + \
        (1 | mri_info_deviceserialnumber)"
    ) %>% 
        as.formula()
    
    res_fit <- lmer(formula = f1, data = master_preprocessed_df)
    res <- res_fit %>% 
        parameters::model_parameters(
            standardize = "refit", 
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        slice(2) %>% 
        mutate(
            label = fc,
            from = rsfmri_variables_df %>% filter(variable_name == fc) %>% pull(from),
            to = rsfmri_variables_df %>% filter(variable_name == fc) %>% pull(to),
            cohend = effectsize::t_to_d(t, df_error)$d,
            cohend_low = effectsize::t_to_d(t, df_error)$CI_low,
            cohend_high = effectsize::t_to_d(t, df_error)$CI_high,
            .before = Parameter
        )
    res_list[[ith]] <- res
}

gordon_fmri_res <- reduce(res_list, bind_rows) %>% 
    slice(1:169) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    filter(p < 0.05)

gordon_subcortical_fmri_res <- reduce(res_list, bind_rows) %>% 
    slice(170:306) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    filter(p < 0.05)


# cortical-to-cortical
networkname <- rsfmri_variables_df %>% 
    filter(from != "none") %>% 
    pull(from) %>% 
    unique()

within_network_connection <- gordon_fmri_res %>% 
    filter(from == to) %>% 
    select(from, cohend) %>% 
    rename(name = from, size = cohend)

hierarchy <- data.frame(from = "origin", to = networkname)

gordon_network_color <- c(
    "#A963F4", "#6F30B6", "#F8D196", "#EA463B",
    "#68D04D", "#FEFF5E", "#F3BCFB", "#2E2E2E",
    "#8EF8FE", "#F2A448", "#56A7A9", "#2F2DB3"
)

vertices <- data.frame(
    name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))
) %>% 
    left_join(within_network_connection, by = "name")

mygraph <- graph_from_data_frame(hierarchy, vertices = vertices)

connect <- gordon_fmri_res %>% 
    filter(from != "none") %>% 
    filter(to != "none") %>% 
    select(from, to, cohend, p.adj) %>% 
    drop_na() %>% 
    mutate(linetype = if_else(p.adj < 0.05, 1, 3))

from <- match(connect$from, vertices$name)
to <- match(connect$to, vertices$name)

gordon_fc_map <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    geom_node_point(
        aes(filter = leaf, x = x*1.1, y = y*1.1, size = 5, color = name)
    ) +
    scale_color_manual(values = gordon_network_color) +
    geom_conn_bundle(
        data = get_con(from = from, to = to, col = connect$cohend, type = connect$linetype), 
        aes(color = col),
        alpha = 0.6,
        width = 3,
        tension = 1
    ) +
    scale_edge_linetype_manual(values = 3) +
    scale_edge_color_gradientn(
        colours = c("navy", "grey90", "tomato1"),
        limits = c(-0.1, 0.1),
        name = "Cohen's d"
    ) +
    guides(color = "none", size = "none") +
    theme_void() + 
    theme(legend.position = "right")
gordon_fc_map

gordon_data <- gordon %>% 
    as_tibble() %>% 
    mutate(network = word(region, 1, sep = "_")) %>%
    filter(network != "NA", network != "None") %>% 
    mutate(network = recode(network,
        FrontoParietal = "Frontoparietal",
        CinguloOperc = "Cingulo-Opercular",
        DorsalAttn = "Dorsal Attention",
        VentralAttn = "Ventral Attention",
        SMhand = "Sensorimotor Hand",
        SMmouth = "Sensorimotor Mouth",
        ParietoOccip = "Retrosplenial Temporal",
        MedialParietal = "Cinguloparietal"
    ))

gordon_brain <- gordon_data %>% 
    ggplot() +
        geom_brain(
            atlas = gordon,
            position = position_brain(side ~ hemi),
            mapping = aes(fill = network),
            color = "grey80"
        )  + 
        scale_fill_manual(
            values = gordon_network_color,
            na.translate = FALSE
        ) +
        theme_void() +
        labs(fill = "") +
        theme(legend.position = "bottom")

gordon_brain
gordon_fc_map

ggexport(
    gordon_brain, 
    filename = here("outputs", "figs", "gordon_brain.pdf"),
    width = 6,
    height = 5
)
ggexport(
    gordon_fc_map, 
    filename = here("outputs", "figs", "gordon_fc_map.pdf"),
    width = 4,
    height = 3.5
)

# cortical-to-subcortical

gs_networkname <- rsfmri_variables_df %>% 
    filter(to != "none") %>% 
    pull(to) %>% 
    unique()

vertices <- data.frame(
    name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))
)
gs_hierarchy <- data.frame(from = "origin", to = gs_networkname)

gs_graph <- graph_from_data_frame(gs_hierarchy, vertices = vertices)

gs_connect <- gordon_subcortical_fmri_res %>% 
    filter(from != "none") %>% 
    filter(to != "none") %>% 
    select(from, to, cohend, p.adj) %>% 
    drop_na() %>% 
    mutate(linetype = if_else(p.adj < 0.05, 1, 3))

gs_from <- c(from, match(gs_connect$from, vertices$name))
gs_to <- c(to, match(gs_connect$to, vertices$name))


gs_fc_map <- ggraph(gs_graph, layout = 'dendrogram', circular = TRUE) + 
    geom_node_point(
        aes(filter = leaf, x = x*1.1, y = y*1.1, size = 5, color = name)
    ) +
    #scale_color_manual(values = gordon_network_color) +
    geom_conn_bundle(
        data = get_con(from = gs_from, to = gs_to, 
                       col = gs_connect$cohend, type = gs_connect$linetype), 
        aes(color = col),
        alpha = 0.6,
        width = 3,
        tension = 0
    ) +
    scale_edge_linetype_manual(values = 3) +
    scale_edge_color_gradientn(
        colours = c("navy", "grey90", "tomato1"),
        limits = c(-0.1, 0.1),
        name = "Cohen's d"
    ) +
    guides(color = "none", size = "none") +
    theme_void() + 
    theme(legend.position = "right")
gs_fc_map
