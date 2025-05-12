rm(list = ls())

library(haven)
library(tidyverse)
library(lavaan)
library(lme4)
library(here)

# years of education, secondary schooling, tertiary education
df_master <- read_dta(file = here("Daten/SUF_9-0-0/ZA6701_master_v9-0-0.dta"),
                      col_select = c(pid, zyg0112, cgr))

df_f2f1 <- read_dta(file = here("Daten/SUF_9-0-0/ZA6701_person_wid1_v9-0-0.dta"),
                    col_select = c(pid, fid, ptyp, age0100, igf0182, igf0282, igf0382, igf0482, eca0100)) %>% 
    rename_with(.cols = c(age0100, eca0100),
                .fn = \(x) paste0(x, "_1"))
df_f2f2 <- read_dta(file = here("Daten/SUF_9-0-0/ZA6701_person_wid3_v9-0-0.dta"),
                    col_select = c(pid, fid, ptyp, age0100, eca0100, starts_with("bmi"))) %>% 
    rename_with(.cols = c(age0100, eca0100),
                .fn = \(x) paste0(x, "_2"))
df_f2f3 <- read_dta(file = here("Daten/SUF_9-0-0/ZA6701_person_wid5_v9-0-0.dta"),
                    col_select = c(pid, fid, ptyp, age0100, starts_with(c("sdo","rwa","bmi")))) %>% 
    rename_with(.cols = age0100,
                .fn = \(x) paste0(x, "_3"))

df_full <- df_f2f1 %>%
    full_join(df_f2f2) %>% 
    full_join(df_f2f3) %>% 
    inner_join(df_master)

df_full <- df_full %>% 
    mutate(across(.cols = everything(),
                  .fns = \(x) case_when(x < 0 ~ NA, .default = x))) %>% 
    filter(ptyp %in% c(1,2) & cgr > 1 & !is.na(zyg0112)) %>% 
    mutate(rwa0101 = 6 - rwa0101,
           rwa0102 = 6 - rwa0102,
           sdo0102 = 6 - sdo0102)

df_rwa <- df_full %>% 
    select(starts_with(c("sdo","rwa")))
cor(df_rwa, use = "pairwise.complete.obs")

create_IQ_scores <- function(df, c, id, IQ = TRUE) {
    
    df_cfa <- df %>% 
        filter(ptyp %in% c(1,2) & cgr == c)
    df_cfa <- zap_labels(df_cfa)
    
    if (IQ == TRUE) {
        if (c == 1) {
            df_cfa <- df_cfa %>% 
                filter(age0100 > 4)
            cfa_model <- paste0('IQ =~ NA*igf0582 + igf0682 + igf0782 
                IQ ~~ 1*IQ')
        } else {
            cfa_model <- paste0('IQ =~ NA*igf0182 + igf0282 + igf0382 + igf0482 
                IQ ~~ 1*IQ')
        }
    } 
    cfa_result <- cfa(cfa_model, data = df_cfa, missing = "fiml")
    cfa_predict <- lavPredict(cfa_result)
    df_cfa <- df_cfa %>% 
        select(!!sym(id)) %>% 
        bind_cols(cfa_predict) 
    
    return(df_cfa)
}

residualize_var <- function(df, model) {
    var_res_label <- model %>% 
        str_extract(".*~") %>% 
        str_replace("~", "") %>% 
        str_trim() %>% 
        paste0("_res")
    
    df <- df %>% 
        mutate(row_id = row_number())
    
    var_res <- lm(formula = model, data = df)[["residuals"]]
    
    df_res <- bind_cols(!!sym(var_res_label) := as.numeric(var_res),
                        row_id = as.numeric(names(var_res)))
    
    df <- df %>% 
        full_join(df_res) %>% 
        select(-row_id)
    
    return(df)
}

create_latent_scores <- function(df, cfa_model, filter_cond, id) {
    df_cfa <- df %>% 
        filter(eval(parse(text = filter_cond)))
    df_cfa <- zap_labels(df_cfa)
    factor <- cfa_model %>% str_extract(".*=~") %>% str_replace("=~","") %>% str_trim()
    cfa_result <- cfa(cfa_model, data = df_cfa, missing = "fiml")
    cfa_predict <- lavPredict(cfa_result)
    
    df_cfa <- df_cfa %>% 
        select(!!sym(id)) %>% 
        bind_cols(cfa_predict) 
    
    return(df_cfa)
}

rwa_model <- "RWA =~ NA*rwa0100+rwa0101+rwa0102+rwa0103\nRWA~~1*RWA"
sdo_model <- "SDO =~ NA*sdo0100+sdo0101+sdo0102+sdo0103\nSDO~~1*SDO"
conditions <- c("cgr == 2", "cgr == 3", "cgr == 4")

# Create latent scores
df_full <- c(2:4) %>% 
    map(create_IQ_scores, df = df_full, id = "pid") %>% 
    bind_rows() %>% 
    full_join(df_full)

df_full <- residualize_var(df = df_full, model = "IQ ~ age0100_1")

df_full <- conditions %>% 
    map(create_latent_scores, cfa_model = rwa_model, df = df_full, id = "pid") %>% 
    bind_rows() %>% 
    full_join(df_full)

df_full <- conditions %>% 
    map(create_latent_scores, cfa_model = sdo_model, df = df_full, id = "pid") %>% 
    bind_rows() %>% 
    full_join(df_full)

summary(df_full)

cor(df_full[,c("RWA","SDO","IQ_res")], use = "pairwise.complete.obs")


df_full <- df_full %>% 
    group_by(fid) %>% 
    mutate(across(.cols = c("RWA","SDO","IQ_res"),
                  .fns = \(x) mean(x, na.rm = TRUE),
                  .names = "{col}_mean"),
           across(.cols = c("RWA","SDO","IQ_res"),
                  .fns = \(x) (x - mean(x, na.rm = TRUE)),
                  .names = "{col}_dev")) %>% 
    ungroup() %>% 
    mutate(zygosity = case_when(zyg0112 == 1 ~ "MZ",
                                zyg0112 == 2 ~ "DZ", 
                                TRUE ~ NA),
           cohort = case_when(cgr == 2 ~ "C2",
                              cgr == 3 ~ "C3",
                              cgr == 4 ~ "C4",
                              .default = NA))

# check case numbers
df_full <- df_full %>% 
    mutate(validIQ = case_when(!is.na(IQ) ~ 1, .default = 0),
           validRWA = case_when(!is.na(RWA) ~ 1, .default = 0),
           validSDO = case_when(!is.na(SDO) ~ 1, .default = 0)) %>% 
    group_by(fid) %>% 
    mutate(N = n(),
           validIQ = sum(validIQ),
           validRWA = sum(validRWA),
           validSDO = sum(validSDO)) %>% 
    ungroup() %>% 
    mutate(full_pair_RWA = case_when(validIQ == 2 & validRWA == 2 ~ 1, .default = 0),
           full_pair_SDO = case_when(validIQ == 2 & validSDO == 2 ~ 1, .default = 0))

View(a[,c("fid","pid", "IQ","RWA","SDO","validIQ","validRWA","validSDO")])
df_rwa %>% 
    filter(!is.na(IQ) & !is.na(SDO)) %>% 
    group_by(fid) %>% 
    mutate(N = n()) %>% 
    ungroup() %>% 
    mutate(N1 = case_when(N == 1 ~ 1, .default = 0),
           N2 = case_when(N == 2 ~ 1, .default = 0)) %>% 
    group_by(zygosity) %>% 
    summarise(
        N1_sum = sum(N1),
        N2_sum = sum(N2)/2
    )

# m0 <- lmer("RWA ~ IQ_res + (1|fid)", data = df_full)
# summary(m0)
#
# m1 <- lmer("RWA ~ IQ_res_dev + IQ_res_mean + (1|fid)", data = df_full)
# summary(m1)
#
m3 <- lmer("RWA ~ IQ_res_dev*zygosity + IQ_res_mean + (1|fid)", data = df_full)
summary(m3)
#
# m0 <- lmer("SDO ~ IQ_res + (1|fid)", data = df_full)
# summary(m0)
#
# m1 <- lmer("SDO ~ IQ_res_dev + IQ_res_mean + (1|fid)", data = df_full)
# summary(m1)
#
m2 <- lmer("SDO ~ IQ_res_dev*zygosity + IQ_res_mean + (1|fid)", data = df_full %>% filter(full_pair_SDO == 1))
summary(m2)
