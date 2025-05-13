rm(list = ls())

library(MASS)
library(lme4)
library(broom.mixed)
library(tidyverse)
library(marginaleffects)
library(lmerTest)
library(doParallel)
library(foreach)
library(paletteer)

# Input-Parameter
a11 <- c(sqrt(0.5), sqrt(0.3), sqrt(0.6), sqrt(0.7))
c11 <- c(sqrt(0.2), sqrt(0.2), sqrt(0.2), sqrt(0.3))
e11 <- c(sqrt(0.3), sqrt(0.5), sqrt(0.2), sqrt(0.2))
ace11 <- tibble(a11 = a11, c11 = c11, e11 = e11)
n_ace <- nrow(ace11)

params <- expand_grid(
    b = c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
    a21 = c(0, 0.2, 0.4, 0.5, 0.6),
    c21 = c(0, 0.2, 0.4, 0.5, 0.6),
    e21 = c(0, 0.2, 0.4, 0.5, 0.6),
    a22_weights = c(0.4),
    c22_weights = c(0.2),
    e22_weights = c(0.4),
    N_mz = c(572, 573),
    N_dz = 602,
    alpha = 0.05
)
n_params <- nrow(params)

ace_expanded <- ace11[rep(seq_len(n_ace), times = n_params), ]
params_expanded <- params[rep(seq_len(n_params), each = n_ace), ]
par_grid <- bind_cols(ace_expanded, params_expanded)
par_grid <- data.frame(par_grid)

get_info_lme <- function(model, variable, info) {
    tidy(model) %>% 
        as_tibble() %>% 
        dplyr::filter(term == variable) %>% 
        dplyr::select(all_of(info)) %>% 
        pull()
}

twinfe_sim <- function(round = 1, alpha, N_mz, N_dz, a11, c11, e11, b, a21, c21, e21, a22_weights, c22_weights, e22_weights) {
    
    set.seed(round)
    
    # Residual Variance Y
    res_variance <- 1-(b^2 + a21^2 + c21^2 + e21^2)
    if (res_variance > 0) {
        a22 <- sqrt(a22_weights*res_variance)
        c22 <- sqrt(c22_weights*res_variance)
        e22 <- sqrt(e22_weights*res_variance)
        
        SigmaA_mz <- matrix(c(1,1,0,0,
                              1,1,0,0,
                              0,0,1,1,
                              0,0,1,1),
                            byrow = TRUE,
                            ncol = 4)
        
        SigmaA_dz <- matrix(c(1,0.5,0,0,
                              0.5,1,0,0,
                              0,0,1,0.5,
                              0,0,0.5,1),
                            byrow = TRUE,
                            ncol = 4)
        
        SigmaC <- matrix(c(1,0,
                           0,1),
                         byrow = TRUE,
                         ncol = 2)
        SigmaE <- matrix(c(1,0,0,0,
                           0,1,0,0,
                           0,0,1,0,
                           0,0,0,1),
                         byrow = TRUE,
                         ncol = 4)
        
        A_mz <- mvrnorm(n = N_mz,
                        mu = rep(0, 4),
                        empirical = TRUE,
                        Sigma = SigmaA_mz) %>% 
            as_tibble() %>% 
            set_names(c("Ax1","Ax2","Ay1","Ay2"))
        
        A_dz <- mvrnorm(n = N_dz,
                        mu = rep(0, 4),
                        empirical = TRUE,
                        Sigma = SigmaA_dz) %>% 
            as_tibble() %>% 
            set_names(c("Ax1","Ax2","Ay1","Ay2"))
        
        A <- bind_rows(A_mz, A_dz)
        
        C <- mvrnorm(n = N_mz + N_dz,
                     mu = rep(0, 2),
                     empirical = TRUE,
                     Sigma = SigmaC) %>% 
            as_tibble() %>% 
            set_names(c("Cx","Cy"))
        
        E <- mvrnorm(n = N_mz + N_dz,
                     mu = rep(0, 4),
                     empirical = TRUE,
                     Sigma = SigmaE) %>% 
            as_tibble() %>% 
            set_names(c("Ex1","Ex2","Ey1","Ey2"))
        
        df <- bind_cols(A, C, E) %>% 
            mutate(zyg = case_when(row_number() <= N_mz ~ 0,
                                   row_number() > N_mz ~ 1),
                   fid = 1:(N_mz+N_dz),
                   X1 = a11*Ax1 + c11*Cx + e11*Ex1,
                   X2 = a11*Ax2 + c11*Cx + e11*Ex2,
                   Y1 = X1*b + a21*Ax1 + c21*Cx + a22*Ay1 + c22*Cy + e22*Ey1,
                   Y2 = X2*b + a21*Ax2 + c21*Cx + a22*Ay2 + c22*Cy + e22*Ey2)
        
        df_long <- df %>% 
            pivot_longer(cols = ends_with(c("1","2")),
                         names_to = c(".value","sib"),
                         names_pattern = "(.*)(\\d)") %>% 
            group_by(fid) %>% 
            mutate(X_mean = mean(X),
                   X_dev = X - X_mean)
        
        m1 <- lmer(Y ~ X + (1|fid), data = df_long)
        m2 <- lmer(Y ~ X_dev + X_mean + (1|fid), data = df_long)
        m3 <- lmer(Y ~ X_dev*zyg + X_mean + (1|fid), data = df_long)
        
        m3_by_zyg <- avg_slopes(m3,
                                variables = "X_dev",
                                by = "zyg")
        
        
        
        m1_est <- get_info_lme(model = m1, variable = "X", info = "estimate")
        m1_p <- get_info_lme(model = m1, variable = "X", info = "p.value")
        
        m2_est <- get_info_lme(model = m2, variable = "X_dev", info = "estimate")
        m2_p <- get_info_lme(model = m2, variable = "X_dev", info = "p.value")
        
        m3_est_mz <- m3_by_zyg %>% 
            filter(zyg == 0) %>% 
            dplyr::select(estimate) %>% 
            pull()
        m3_est_dz <- m3_by_zyg %>% 
            filter(zyg == 1) %>% 
            dplyr::select(estimate) %>% 
            pull()
        m3_p_mz <- m3_by_zyg %>% 
            filter(zyg == 0) %>% 
            dplyr::select(p.value) %>% 
            pull()
        m3_p_dz <- m3_by_zyg %>% 
            filter(zyg == 1) %>% 
            dplyr::select(p.value) %>% 
            pull()
        
        return_df <- tibble(Run = round,
                            Model = c("M1","M2","M3","M3"),
                            Group = c("POLS",
                                      "Within-Pooled",
                                      "Within-DZ",
                                      "Within-MZ"),
                            Est = c(m1_est, m2_est, m3_est_dz, m3_est_mz),
                            P = c(m1_p, m2_p, m3_p_dz, m3_p_mz),
                            Sig = P <= alpha,
                            Note = NA)
    } else {
        return_df <- tibble(Run = round,
                            Model = c("M1","M2","M3","M3"),
                            Group = c("POLS",
                                      "Within-Pooled",
                                      "Within-DZ",
                                      "Within-MZ"),
                            Est = rep(NA, 4),
                            P = rep(NA, 4),
                            Sig = rep(NA, 4),
                            Note = "Res Var < 0")
    }
    
    
    return(return_df)
}

twinfe_power <- function(n_sim = 500, alpha = alpha, N_mz, N_dz, a11, c11, e11, b, a21, c21, e21, a22_weights, c22_weights, e22_weights) {
    power_result <- list()
    for (i in 1:n_sim) {
        cat(".")
        power_result[[i]] <- twinfe_sim(round = i, alpha = alpha, N_mz = N_mz, N_dz = N_dz, a11 = a11, c11 = c11, e11 = e11, 
                                        b = b, a21 = a21, c21 = c21, e21 = e21, 
                                        a22_weights = a22_weights, c22_weights = c22_weights, e22_weights = e22_weights)
    }
    power_df <- power_result %>% 
        bind_rows() %>% 
        mutate(b = b,
               alpha = alpha, N_mz = N_mz, N_dz = N_dz, 
               a11 = a11, c11 = c11, e11 = e11,
               a21 = a21, c21 = c21, e21 = e21,
               a22_weights = a22_weights, c22_weights = c22_weights, e22_weights = e22_weights) %>% 
        group_by(Group) %>% 
        mutate(Power = mean(Sig)) %>% 
        ungroup() %>%
        dplyr::select(-c(Run, Est, P, Sig)) %>% 
        distinct()
    return(power_df)
}




# Parallel Processing Config
parallel::detectCores()
n_cores <- 12
my_cluster <- parallel::makeCluster(
    n_cores, 
    type = "PSOCK"
)
doParallel::registerDoParallel(cl = my_cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

# Do the power analysis
power_result <- foreach (k = 1:nrow(par_grid),
                         .packages = c("tidyverse","MASS", "marginaleffects", "broom.mixed", "lme4", "lmerTest")) %dopar% {
                             twinfe_power(n_sim = 500,
                                          alpha = par_grid[k,"alpha"], N_mz = par_grid[k,"N_mz"], N_dz = par_grid[k,"N_dz"], 
                                          a11 = par_grid[k,"a11"], c11 = par_grid[k,"c11"], e11 = par_grid[k,"e11"], 
                                          b = par_grid[k,"b"], a21 = par_grid[k,"a21"], c21 = par_grid[k,"c21"], e21 = par_grid[k,"e21"], 
                                          a22_weights = par_grid[k,"a22_weights"], c22_weights = par_grid[k,"c22_weights"], e22_weights = par_grid[k,"e22_weights"])
                         }
stopCluster(my_cluster)

saveRDS(object = power_result, file = "power_result.rds")

power_result <- read_rds("Authoritarianism/power_result.rds") %>% 
    dplyr::select(-b) %>% 
    distinct() %>% 
    bind_rows() %>% 
    mutate(Outcome = case_when(N_mz == 572 ~ "RWA",
                               N_mz == 573 ~ "SDO"))

# one plot for every combination of a11, c11 and e11
for (i in 1:nrow(ace11)) {
    mz_power <- power_result %>% 
        filter(is.na(Note) & Group == "Within-MZ" & a11 == as.numeric(ace11[i, "a11"]) & c11 == as.numeric(ace11[i, "c11"]) & e11 == as.numeric(ace11[i, "e11"]) & e21 == 0) 
    
    power_plot <- ggplot(mz_power, 
                         aes(x = b, y = Power, color = factor(Outcome), group = Outcome)) +
        geom_line(size = 1.5, alpha = 0.6) + 
        geom_hline(yintercept = 0.8, linetype = 'dashed') + 
        geom_vline(xintercept = 0.3, linetype = 'dashed', color = paletteer_d("nationalparkcolors::Badlands")[3]) + 
        theme_minimal() + 
        scale_y_continuous(labels = scales::percent) + 
        scale_color_paletteer_d("nationalparkcolors::Badlands") +
        theme(plot.caption = element_text(hjust = 0),
              plot.title = element_text(hjust = 0.5)) +
        labs(x = "Linear Effect Size",
             y = 'Power', 
             #title = "Results of Power Analysis",
             caption = paste0("Note: Power calculated based on Monte Carlo Simulations with M=500 rounds. \nEffect size reported by meta analysis by Onraet et al. (2015) marked by dashed red line. \nFixed parameters: alpha=",
                              as.numeric(mz_power[1,"alpha"]),", a11 = sqrt(",
                              as.numeric(mz_power[1,"a11"])**2,"), c11 = sqrt(",
                              as.numeric(mz_power[1,"c11"])**2,"), e11 = sqrt(",
                              as.numeric(mz_power[1,"e11"])**2,"), \na22 = sqrt(",
                              as.numeric(mz_power[1,"a22_weights"]),"*e(Y)), c22 = sqrt(",
                              as.numeric(mz_power[1,"c22_weights"]),"*e(Y)), e22 = sqrt(",
                              as.numeric(mz_power[1,"e22_weights"]),"*e(Y)), \ne(Y) = 1-(b^2 + a21^2 + c21^2 + e21^2)"
             ),
             color = 'Outcome') +
        facet_grid(
            a21 ~ c21,
            labeller = labeller(
                a21 = c("0" = "a21 = 0", "0.2" = "a21 = 0.2", "0.4" = "a21 = 0.4", "0.5" = "a21 = 0.5"),
                c21 = c("0" = "c21 = 0", "0.2" = "c21 = 0.2", "0.4" = "c21 = 0.4", "0.5" = "c21 = 0.5")
            )
        )
    power_plot
    ggsave(filename = paste0("Authoritarianism/power_plot",i,".png"), bg = "white",
           width = 8, height = 8,
           plot = power_plot, dpi = 800)
}







