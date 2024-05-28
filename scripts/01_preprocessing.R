# Environment ------------------------------------------------------------------
library(here)
library(tidyverse)
library(gtsummary)
library(glue)
library(flextable)
library(officer)

source(here("src", "R", "is_binary_column.R"))

# Data I/O ---------------------------------------------------------------------
# ABCD data path
# abcd_path <- "/path/to/abcd/folder"
abcd_folder_path <- "/Users/tywong/GoogleDrive/opendata/abcd/abcd-data-release-5.1/core"

# general
abcd_y_lt <- here(abcd_folder_path, "abcd-general", "abcd_y_lt.csv") %>% 
    read_csv(show_col_types = FALSE)
abcd_p_demo <- here(abcd_folder_path, "abcd-general", "abcd_p_demo.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_adm_info <- here(abcd_folder_path, "imaging", "mri_y_adm_info.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_qc_clfind <- here(abcd_folder_path, "imaging", "mri_y_qc_clfind.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_qc_motion <- here(abcd_folder_path, "imaging", "mri_y_qc_motion.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_qc_incl <- here(abcd_folder_path, "imaging", "mri_y_qc_incl.csv") %>% 
    read_csv(show_col_types = FALSE)

# loneliness data
loneliness_baseline <- here(abcd_folder_path, "mental-health", "mh_p_cbcl.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(src_subject_id, eventname, cbcl_q12_p) %>% 
    drop_na(cbcl_q12_p) %>% 
    filter(eventname == "baseline_year_1_arm_1") %>% 
    mutate(cbcl_q12_p = if_else(cbcl_q12_p == 0, 0, 1)) %>% 
    rename(loneliness_bl = cbcl_q12_p)
loneliness_followup <- here(abcd_folder_path, "mental-health", "mh_p_cbcl.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(src_subject_id, eventname, cbcl_q12_p) %>% 
    drop_na(cbcl_q12_p) %>% 
    filter(
        eventname == "1_year_follow_up_y_arm_1" | 
        eventname == "2_year_follow_up_y_arm_1" | 
        eventname == "3_year_follow_up_y_arm_1"
    ) %>% 
    mutate(loneliness_bin = if_else(cbcl_q12_p > 0, 1, 0)) %>% 
    group_by(src_subject_id) %>% 
    mutate(n = n()) %>% 
    filter(n == 3) %>% 
    mutate(
        loneliness_followup = sum(loneliness_bin),
        loneliness_followup_bin = if_else(loneliness_followup == 0, 0, 1)
    ) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(src_subject_id, eventname, loneliness_followup, loneliness_followup_bin) %>% 
    mutate(eventname = "baseline_year_1_arm_1")

# read in selected variables
health_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Health")
exposome_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Exposome")
gmv_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "Gray Matter Volume")
wmt_variables_df <- here("data", "raw", "included_variables.xlsx") %>% 
    readxl::read_excel(sheet = "White Matter Tracts")

selected_filenames <- c(
    unique(pull(exposome_variables_df, table_name)), 
    unique(pull(health_variables_df, table_name)), 
    unique(pull(gmv_variables_df, table_name)),
    unique(pull(wmt_variables_df, table_name))
) %>% 
    discard(is.na) 
selected_variables <- c(
    pull(exposome_variables_df, variable_name),
    pull(health_variables_df, variable_name),
    pull(gmv_variables_df, variable_name),
    pull(wmt_variables_df, variable_name)
)

selected_files <- list()
for (ith in 1:length(selected_filenames)) {
    selected_files[[ith]] <- list.files(
        path = abcd_folder_path, 
        full.names = TRUE,
        pattern = glue("{selected_filenames[ith]}.csv"),
        recursive = TRUE
    ) %>% 
        read_csv(show_col_types = FALSE) %>% 
        select(src_subject_id, eventname, any_of(selected_variables))
}

# combine
master_df <- selected_files %>% 
    reduce(left_join, by = c("src_subject_id", "eventname")) %>% 
    right_join(mri_y_qc_motion, by = c("src_subject_id", "eventname")) %>% 
    right_join(mri_y_qc_clfind, by = c("src_subject_id", "eventname")) %>% 
    right_join(mri_y_qc_incl, by = c("src_subject_id", "eventname")) %>% 
    right_join(mri_y_adm_info, by = c("src_subject_id", "eventname")) %>% 
    right_join(loneliness_baseline, by = c("src_subject_id", "eventname")) %>% 
    right_join(loneliness_followup, by = c("src_subject_id", "eventname")) %>% 
    right_join(abcd_p_demo, by = c("src_subject_id", "eventname"), suffix = c("", ".y")) %>% 
    right_join(abcd_y_lt, by = c("src_subject_id", "eventname"))

sample_with_seed <- function(seed_id, x) {
    set.seed(seed_id)
    res <- sample(x)
    return(res)
}

# Preprocessing ----------------------------------------------------------------
master_preprocessed_df <- master_df %>% 
    filter(eventname == "baseline_year_1_arm_1") %>% # 11,868 subjects
    group_by(rel_family_id) %>% 
    mutate(
        random_order_id = sample_with_seed(rel_family_id, n())
    ) %>% 
    filter(random_order_id == 1) %>% # 9,807 subjects
    ungroup() %>% 
    filter(
        demo_sex_v2 %in% c(1, 2), # 9,804 subjects
    ) %>% 
    mutate(
        across(where(is.numeric), ~na_if(., 777)),
        across(where(is.numeric), ~na_if(., 999))
    ) %>%
    drop_na(
        src_subject_id, site_id_l, interview_age, demo_sex_v2, race_ethnicity,
        loneliness_bl
    ) %>% # 7,968
    mutate(
        demo_sex_v2 = factor(demo_sex_v2, labels = c("Male", "Female")),
        race_white = if_else(demo_race_a_p___10 == 1, 1, 0),
        race_black = if_else(demo_race_a_p___11 == 1, 1, 0),
        race_asian = if_else(
            (demo_race_a_p___18 + demo_race_a_p___19 + demo_race_a_p___20 +
                 demo_race_a_p___21 + demo_race_a_p___22 + demo_race_a_p___23 +
                 demo_race_a_p___24) > 0, 1, 0
        ),
        race_aian = if_else((demo_race_a_p___12 + demo_race_a_p___13) > 0, 1, 0),
        race_nhpi = if_else(
            (demo_race_a_p___14 + demo_race_a_p___15 + demo_race_a_p___17 +
                 demo_race_a_p___17) > 0, 1, 0
        ),
        race_other = if_else(demo_race_a_p___25 == 1, 1, 0),
        race_mixed = if_else(
            (race_white + race_black + race_asian + race_aian + race_aian + 
                 race_other) > 1, 1, 0
        ),
        race_ethnicity = case_when(
            race_white == 1 & race_mixed == 0 & demo_ethn_v2 == 2 ~ 1,
            race_black == 1 & race_mixed == 0 ~ 2,
            race_asian == 1 & race_mixed == 0 ~ 3,
            demo_ethn_v2 == 1 & race_mixed == 0 ~ 4,
            race_mixed == 1 ~ 5,
            demo_race_a_p___77 | demo_race_a_p___99 ~ NA,
            TRUE ~ 6
        ),
        race_ethnicity = factor(
            race_ethnicity, 
            levels = 1:6,
            label = c("Non-Hispanic White", "Black", "Asian", "Hispanic", "Mixed", "Other")
        ),
        race_ethnicity_3 = case_when(
            race_white == 1 & race_mixed == 0 & demo_ethn_v2 == 2 ~ 1,
            race_black == 1 & race_mixed == 0 ~ 2,
            demo_race_a_p___77 | demo_race_a_p___99 ~ NA,
            TRUE ~ 3
        ),
        race_ethnicity_3 = factor(
            race_ethnicity_3, 
            levels = 1:3,
            label = c("White", "Black", "Other")
        ),
        demo_comb_income_v2 = case_when(
            demo_comb_income_v2 <= 4 ~ 1,
            demo_comb_income_v2 > 4 & demo_comb_income_v2 <= 6 ~ 2,
            demo_comb_income_v2 == 7 ~ 3,
            demo_comb_income_v2 == 8 ~ 4,
            demo_comb_income_v2 == 9 ~ 5,
            demo_comb_income_v2 == 10 ~ 6,
            TRUE ~ NA
        ),
        demo_comb_income_v2 = ordered(
            demo_comb_income_v2,
            levels = 1:6,
            label = c("1. <$25K", "2. $25K-$49K", "3. $50K-$74K", 
                      "4. $75K-$99K", "5. $100K-$199K", "6. $200K+")
        ),
        medhx_total = medhx_2a + medhx_2b + medhx_2d + medhx_2e + medhx_2f +
            medhx_2g + medhx_2h + medhx_2i + medhx_2j + medhx_2k + medhx_2l +
            medhx_2m + medhx_2n + medhx_2o + medhx_2p + medhx_2q + medhx_2r +
            medhx_2s,
        birth_weight_kg = birth_weight_lbs * 0.45359237 + birth_weight_oz * 0.0283495231,
        famhx_ss_fath_mh_history = famhx_ss_fath_prob_alc_p + 
            famhx_ss_fath_prob_dg_p + famhx_ss_fath_prob_dprs_p + 
            famhx_ss_fath_prob_ma_p + famhx_ss_fath_prob_vs_p + 
            famhx_ss_fath_prob_trb_p + famhx_ss_fath_prob_nrv_p + 
            famhx_ss_fath_prob_prf_p + famhx_ss_fath_prob_hspd_p + 
            famhx_ss_fath_prob_scd_p,
        famhx_ss_moth_mh_history = famhx_ss_moth_prob_alc_p + 
            famhx_ss_moth_prob_dg_p + famhx_ss_moth_prob_dprs_p + 
            famhx_ss_moth_prob_ma_p + famhx_ss_moth_prob_vs_p + 
            famhx_ss_moth_prob_trb_p + famhx_ss_moth_prob_nrv_p + 
            famhx_ss_moth_prob_prf_p + famhx_ss_moth_prob_hspd_p + 
            famhx_ss_moth_prob_scd_p,
        devhx_10sum3 = devhx_10a3_p + devhx_10b3_p + devhx_10c3_p +
            devhx_10d3_p + devhx_10e3_p + devhx_10f3_p + devhx_10g3_p +
            devhx_10h3_p + devhx_10i3_p + devhx_10j3_p + devhx_10k3_p +
            devhx_10l3_p + devhx_10m3_p,
        devhx_14sum3 = devhx_14a3_p + devhx_14b3_p + devhx_14c3_p + 
            devhx_14d3_p + devhx_14e3_p + devhx_14f3_p + devhx_14g3_p +
            devhx_14h3_p,
        ksads_ptsd_sum = ksads_ptsd_raw_754_p + ksads_ptsd_raw_755_p + 
            ksads_ptsd_raw_756_p + ksads_ptsd_raw_757_p + ksads_ptsd_raw_758_p +
            ksads_ptsd_raw_759_p + ksads_ptsd_raw_760_p + ksads_ptsd_raw_761_p +
            ksads_ptsd_raw_762_p + ksads_ptsd_raw_763_p + ksads_ptsd_raw_764_p +
            ksads_ptsd_raw_765_p + ksads_ptsd_raw_766_p + ksads_ptsd_raw_767_p +
            ksads_ptsd_raw_768_p + ksads_ptsd_raw_769_p + ksads_ptsd_raw_770_p,
        devhx_15 = as.numeric(devhx_15),
        devhx_16_p = as.numeric(devhx_16_p),
        devhx_16_p = if_else(devhx_16_p == 9990, NA, devhx_16_p),
        devhx_17_p = as.numeric(devhx_17_p),
        reshist_addr1_nanda_parks_tc10 = if_else(
            reshist_addr1_nanda_parks_tc10 == "10 or more", "10", 
            reshist_addr1_nanda_parks_tc10
        ),
        reshist_addr1_nanda_parks_tc10 = as.numeric(reshist_addr1_nanda_parks_tc10),
        smri_vol_scs_intracranialv = as.numeric(
            scale(smri_vol_scs_intracranialv, center = TRUE, scale = TRUE)
        ),
        bmi = anthroweightcalc / (anthroheightcalc^2) * 703,
        bmi = if_else(bmi > 100 | bmi < 5, NA, bmi),
        pds_ss_male_category = (pds_p_ss_male_category + pds_y_ss_male_category) / 2,
        pds_ss_female_category = (pds_p_ss_female_category + pds_y_ss_female_category) / 2,
        pds_ss_category = if_else(
            demo_sex_v2 == "Male", pds_ss_male_category, pds_ss_female_category
        )
    ) %>% 
    select(
        -c(birth_weight_lbs, birth_weight_oz, 
           anthroweightcalc, anthroheightcalc,
           pds_p_ss_male_category, pds_y_ss_male_category,
           pds_p_ss_female_category, pds_y_ss_female_category),
        -matches("medhx_2"),
        -matches("famhx_ss_fath_prob_"),
        -matches("famhx_ss_moth_prob_"),
        -matches("devhx_10[abcdefghijklm]3"),
        -matches("devhx_14[abcdefgh]3"),
        -matches("ksads_ptsd_raw_7")
    ) %>% 
    select(
        src_subject_id, site_id_l, mri_info_deviceserialnumber, mrif_score,
        imgincl_t1w_include, imgincl_rsfmri_include, rsfmri_ntpoints, 
        rsfmri_meanmotion, imgincl_dmri_include, dmri_meanmotion, interview_age, 
        demo_sex_v2, race_ethnicity, race_ethnicity_3, bmi, pds_ss_category, loneliness_bl,
        loneliness_followup, loneliness_followup_bin, demo_comb_income_v2, 
        birth_weight_kg, famhx_ss_fath_mh_history, famhx_ss_moth_mh_history, 
        medhx_total, devhx_10sum3, devhx_14sum3, ksads_ptsd_sum,
        smri_vol_cdk_total, smri_vol_scs_subcorticalgv,
        any_of(selected_variables)
    ) %>% 
    mutate_if(is_binary_column, factor)

# pull binary variables with less than 1:9
remove_binary_variables <- master_preprocessed_df %>% 
    select_if(is_binary_column) %>% 
    map_dfr(table, .id = "variable_name") %>% 
    mutate(sum = `0` + `1`) %>% 
    mutate(ratio = `0` / sum) %>% 
    mutate(low = if_else(ratio > 0.9 | ratio < 0.1, 1, 0)) %>%
    filter(low == 1) %>% 
    slice(-1) %>% # qc variable
    pull(variable_name)

# pull participants with more than 10% missing values
remove_id <- master_preprocessed_df %>% 
    is.na() %>% 
    rowSums() %>% 
    as_tibble() %>% 
    mutate(
        src_subject_id = master_preprocessed_df$src_subject_id,
        prop_missing = value / dim(.)[1]
    ) %>% 
    filter(prop_missing > 0.1) %>% 
    pull(src_subject_id)


master_preprocessed_df <- master_preprocessed_df %>% 
    filter(!(src_subject_id %in% remove_id)) %>% # 7,73
    select(-all_of(remove_binary_variables)) 


# visually inspect the data
master_preprocessed_df %>% 
    select(any_of(exposome_variables_df$variable_name)) %>% 
    skimr::skim()
    


# Table 1 ----------------------------------------------------------------------
table1 <- master_preprocessed_df %>% 
    select(loneliness_bl, loneliness_followup_bin, interview_age,
           demo_sex_v2, race_ethnicity_3) %>%
    mutate(
        loneliness_bl = factor(
            loneliness_bl, levels = c(0, 1), labels = c("No", "Yes")),
        loneliness_followup_bin = factor(
            loneliness_followup_bin, levels = c(0, 1), labels = c("No", "Yes"))
    ) %>%
    tbl_summary(
        by = loneliness_bl,
        label = list(
            loneliness_followup_bin ~ "Loneliness in the followups",
            interview_age ~ "Age (months)",
            demo_sex_v2 ~ "Sex at birth",
            race_ethnicity_3 ~ "Race"
        ),
        missing = "no",
        statistic = list(
            all_continuous() ~ "{mean} ± {sd}"
        ),
        digits = list(
            all_continuous() ~ 2
        )
    ) %>% 
    add_p() %>%
    add_q(method = "bonferroni") %>% 
    bold_p(q = TRUE) %>% 
    modify_header(label = "**Characteristics**") %>% 
    modify_spanning_header(
        all_stat_cols() ~ "**Loneliness at baseline**"
    )

# docx page setup
sect_properties <- prop_section(
    page_size = page_size(),
    type = "continuous",
    page_margins = page_mar(
        bottom = 0.5, top = 0.5, right = 0.5, left = 0.5, gutter = 0
    )
)

table1 %>% 
    as_flex_table() %>% 
    fontsize(size = 8.5, part = "all") %>% 
    padding(padding.top = 1, padding.bottom = 1, part = "all") %>% 
    bold(part = "header") %>% 
    set_table_properties(width = 1, layout = "autofit") %>% 
    save_as_docx(
        path = here("outputs", "tables", "table1.docx"),
        pr_section = sect_properties
    )

# proportion
(master_preprocessed_df$loneliness_followup %>% table()) / dim(master_preprocessed_df)[1] * 100


# Save -------------------------------------------------------------------------
master_preprocessed_df %>% 
    write_csv(file = here("data", "processed", "master_preprocessed_df.csv"))
master_preprocessed_df %>% 
    write_rds(file = here("data", "processed", "master_preprocessed_df.rds"))

