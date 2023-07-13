rm(list = ls())
library(tidyverse)
library(cmdstanr)
library(lubridate)
library(bayesplot)
library(gt)
library(flextable)
library(posterior)
library(glue)
library(MetBrewer)
library(ggtext)
library(bayestestR)
library(officer)
library(fastDummies)
set_flextable_defaults(table.layout = "autofit")

#Turn true to create the figure of RD-STDC vs Igg vs cases in Mexico City
#you need a MARIADB installation because the City's cases are too much. 
#See the packages vignettes for more info
generate_plot <- T
if (generate_plot){
  library(covidmx) #devtools::install_github("RodrigoZepeda/covidmx@0732467ff19fb3bae883c8c0e995552e6197e17e")
#  datos_covid <- descarga_datos_abiertos(language = "Español") #datos_covid <- read_datos_abiertos()
}

if (!dir.exists("plots_do_not_share")){
  dir.create("plots_do_not_share")
}

if (!dir.exists("results")){
  dir.create("results")
}

if (!dir.exists("paper_tables")){
  dir.create("paper_tables")
}

complete_dataset       <- read_rds("final_datasets/annonymized_dataset_paper_4_cycles.rds") %>%
  ungroup()

#Check repeated per cycle
complete_dataset |>
  filter(Seleccion == "ALEATORIO") |> 
  group_by(ID) |>
  mutate(idnum = 1:n()) |>
  ungroup() |>
  filter(idnum > 1) |>
  group_by(Ciclo) |>
  tally()

complete_itt           <- read_rds("final_datasets/annonimized_twd_dataset.rds") %>%
  ungroup()

for (filtertype in c("RANDOM", "CONVENIENCE", "EVERYONE")){
    
  if (exists("filtertype") && filtertype == "RANDOM"){ 
    mydata <- complete_dataset %>%
      filter(Seleccion == "ALEATORIO")
    
    ittdata <- complete_itt %>%
      filter(Seleccion == "ALEATORIO")
  } else if (exists("filtertype") && filtertype == "CONVENIENCE"){
    mydata <- complete_dataset %>%
      filter(Seleccion != "ALEATORIO")

    ittdata <- complete_itt %>%
      filter(Seleccion != "ALEATORIO")
  } else {
    mydata     <- complete_dataset
    ittdata    <- complete_itt 
    filtertype <- "EVERYONE"
  }
  
  #Get business name which is not kept in this file for annonimity
  empresa      <- "[_EMPRESA_]"
  if (empresa == ""){
    empresa <- "CENSORED-NAME"
  }
  
  #STAN Options----
  options(mc.cores = max(parallel::detectCores() - 2,1))
  
  #In my machine I have a special compiler
  user <- Sys.info()
  if (user["user"] == "rod"){
    cpath       <- "/usr/bin/gcc"
  } else {
    cpath       <- NULL
  }
  
  cxx_flags   <- "-g -O3 -Wall -pedantic -mtune=native -pipe"
  if (is.null(cpath)){
    cpp_options <- list(cxx_flags = cxx_flags, stan_threads = T)
  } else {
    cpp_options <- list(cxx_flags = cxx_flags, stan_threads = T,
                        cpath = cpath)
  }
  
  #Weights---------------------------------------------
  weights <- mydata %>%
    group_by(Ciclo, strata, pondef, N_Total, N_Strata) %>%
    tally() %>%
    dplyr::select(-n)
  
  stan_weights <- weights %>% 
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = pondef,
                values_fill = 0) %>% 
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #LOGISTIC (MODEL)-------------
  if (filtertype != "CONVENIENCE"){
    
    Ncovariates <- 7
    
    fdata <- mydata %>%
      filter(Ciclo < 4) %>% #Last cycle messes everything up
      filter(!is.na(IGG) & !is.na(edad_cat) & !is.na(sexo) & !is.na(covid_sintoma) &
               covid_sintoma != "No sabe" & covid_sintoma != "No deseo responder") %>% 
      mutate(IGG = if_else(IGG == "Positivo", 1, 0))  %>% 
      mutate(sexo = if_else(sexo == "Hombre", 1, 0)) %>% 
      mutate(covid_sintoma = if_else(covid_sintoma == "Si", 1, 0))  %>% 
      mutate(edad_35 = if_else(edad_cat == "Menores a 35", 1, 0)) %>% 
      mutate(edad_44 = if_else(edad_cat == "35 a 44", 1, 0)) %>%
      mutate(edad_59 = if_else(edad_cat == "45 a 59", 1, 0))  %>%
      mutate(edad_60 = if_else(edad_cat == "60 y más", 1, 0)) %>%
      select(Ciclo, strata, IGG, sexo, covid_sintoma, edad_35, edad_44, edad_59, edad_60, `RD-STDC`)
    
    cts <- fdata %>% group_by(Ciclo, strata) %>% tally()
    
    #Get sample size
    nciclos       <- length(unique(fdata$Ciclo))
    nstrata       <- length(unique(fdata$strata))
    max_obs_cycle <- max(cts$n)
    
    #IGG
    stan_groups <- array(0, dim = c(nstrata, nciclos, max_obs_cycle))
    stan_groups <- provideDimnames(stan_groups , sep = "_", 
                                   base = list("Ciclo","strata","individuo"))
    
    for (ciclo in 1:nciclos){
      for (strata in 1:nstrata){
        estrato <- unique(fdata$strata)[strata]
        nstrat  <- fdata %>% filter(Ciclo == ciclo & strata == !!estrato) 
        
        if (nrow(nstrat) > 0){
          
          stan_groups[strata, ciclo, 1:nrow(nstrat)] <- unlist(nstrat$IGG)
          
        }
      }
    }
    
    #Covariates
    stan_cov <- array(0, dim = c(nstrata, nciclos, Ncovariates, max_obs_cycle))
    stan_cov <- provideDimnames(stan_cov , sep = "_", 
                                   base = list("Ciclo","strata","cov","individuo"))
    
    for (ciclo in 1:nciclos){
      for (strata in 1:nstrata){
        
        estrato <- unique(fdata$strata)[strata]
        nstrat  <- fdata %>% filter(Ciclo == ciclo & strata == !!estrato) 
        
        if (nrow(nstrat) > 0){
          stan_cov[strata, ciclo, 1, 1:nrow(nstrat)] <- nstrat$sexo
          stan_cov[strata, ciclo, 2, 1:nrow(nstrat)] <- nstrat$edad_35
          stan_cov[strata, ciclo, 3, 1:nrow(nstrat)] <- nstrat$edad_44
          stan_cov[strata, ciclo, 4, 1:nrow(nstrat)] <- nstrat$edad_59
          stan_cov[strata, ciclo, 5, 1:nrow(nstrat)] <- nstrat$edad_60
          stan_cov[strata, ciclo, 6, 1:nrow(nstrat)] <- nstrat$covid_sintoma
          stan_cov[strata, ciclo, 7, 1:nrow(nstrat)] <- nstrat$`RD-STDC`
        }
      }
    }
    
    cp_model    <- cmdstan_model("code/additional_functions/logit_sens_spec.stan",
                                 include_paths = ".", stanc_options = list("O1"), 
                                 force_recompile = F)
    
    datos  <- list(
      Ncovariates   = Ncovariates,
      Ncycles       = nciclos,
      Ncenters      = nstrata,
      igg           = stan_groups,
      cov           = stan_cov,
      pondef        = stan_weights[,1:nciclos],
      N_obs_cycle   = cts %>% 
        pivot_wider(id_cols = strata, values_from = n, names_from = Ciclo) %>%
        select(-strata) %>% as.matrix,
      max_obs_cycle = dim(stan_cov)[4],
      sens        = 0.98,
      spec        = 0.907,
      sigma_prior = 0.01
    )
    
    #Valores iniciales 
    initf2 <- function(chain_id = 1) {
      list(
        beta_baseline = as.matrix(rnorm(datos$Ncenters*datos$Ncycles), ncol = datos$Ncycles),
        beta_sexo = as.matrix(rnorm(datos$Ncenters*datos$Ncycles), ncol = datos$Ncycles),
        beta_edad_35_44 = as.matrix(rnorm(datos$Ncenters*datos$Ncycles), ncol = datos$Ncycles),
        beta_edad_45_60 = as.matrix(rnorm(datos$Ncenters*datos$Ncycles), ncol = datos$Ncycles),
        beta_edad_60_plus = as.matrix(rnorm(datos$Ncenters*datos$Ncycles), ncol = datos$Ncycles),
        beta_sintoma = as.matrix(rnorm(datos$Ncenters*datos$Ncycles), ncol = datos$Ncycles),
        beta_itt = as.matrix(rnorm(datos$Ncenters*datos$Ncycles), ncol = datos$Ncycles)
      )
    }
    
    init_ll      <- lapply(1:4, function(id) initf2(chain_id = id))
    
    model_sample <- cp_model$sample(data              = datos,
                                    chains            = 4,
                                    seed              = 454565456,
                                    refresh           = 100,
                                    iter_warmup       = 500,
                                    iter_sampling     = 500,
                                    init              = init_ll,
                                    adapt_delta       = 0.95)
    
    model_draws  <- model_sample$draws() %>% 
      as_draws_df() %>%
      dplyr::select(matches("average_OR")) 
    
    ORs      <- model_draws %>% 
      summarise_draws(~quantile(., c(0.025, 0.5, 0.975))) %>%
      mutate(Ciclo = str_remove_all(variable,"average_OR.*\\[|\\]")) %>%
      mutate(variable = str_remove_all(variable,"average_OR_|\\[[0-9]\\]")) %>%
      rename(`Estimacion_puntual`=`50%` ) %>%
      rename(`CI_95%_low` = `2.5%`) %>%
      rename(`CI_95%_high` = `97.5%`) %>%
      mutate(variable = str_replace_all(variable, "sexo","male")) %>%
      mutate(variable = str_replace_all(variable, "60",">=60")) %>%
      mutate(variable = str_replace_all(variable, "45","45 to 60")) %>%
      mutate(variable = str_replace_all(variable, "35","35 to 45")) %>%
      mutate(variable = str_replace_all(variable, "baseline","intercept")) %>%
      mutate(variable = str_replace_all(variable, "sintoma","covid symptom (yes)")) %>%
      mutate(variable = str_replace_all(variable, "itt","RDSTDC (yes)")) %>%
      mutate(centro = "Todos")
    
    # model_draws  <- model_sample$draws() %>% 
    #   as_draws_df() %>%
    #   dplyr::select(matches("\\bbeta_baseline\\b|\\bbeta_sexo\\b|\\bbeta_edad_35_44\\b|\\bbeta_edad_45_60\\b|\\bbeta_edad_60_plus\\b|\\bbeta_sintoma\\b|\\bbeta_itt\\b")) 
    # 
    # betas_2      <- model_draws %>% 
    #   summarise_draws(~quantile(exp(.), c(0.025, 0.5, 0.975))) %>%
    #   mutate(Ciclo  = str_remove_all(variable,"beta.*\\[[0-9]+,|\\]")) %>%
    #   mutate(Centro = as.numeric(str_remove_all(variable,"beta.*\\[|,[0-9]+\\]"))) %>%
    #   mutate(variable = str_remove_all(variable,"beta_|\\[.*\\]")) %>%
    #   rename(`Estimacion_puntual`=`50%` ) %>%
    #   rename(`CI_95%_low` = `2.5%`) %>%
    #   rename(`CI_95%_high` = `97.5%`) %>%
    #   mutate(variable = str_replace_all(variable, "sexo","hombre")) %>%
    #   left_join(
    #     cts %>% 
    #       pivot_wider(id_cols = strata, values_from = n, names_from = Ciclo) %>% 
    #       select(strata) %>%
    #       mutate(Centro = 1:n())
    #   ) %>%
    #   select(-Centro) %>%
    #   rename(centro = strata)
      
      beta <- ORs #%>% bind_rows(betas_2)
    
      beta %>% 
        mutate(CI = 
                    paste0(scales::comma(as.numeric(`CI_95%_low`),0.01),"-",
                           scales::comma(as.numeric(`CI_95%_high`),0.01))) %>%
        mutate(Estimacion_puntual = scales::comma(as.numeric(Estimacion_puntual), 0.01)) %>%
        select(variable, Estimacion_puntual, CI, Ciclo) %>%
        filter(variable != "intercept") %>%
        mutate(variable = case_when(
          variable == "35 to 45" ~ "35 - 44",
          variable == "45 to 60" ~ "45 - 59",
          variable == ">=60" ~ "60 over",
          variable == "covid symptom (yes)" ~ "covid symptom",
          variable == "RDSTDC (yes)" ~ "RD-STDC (in the last 6 months)",
          TRUE ~ variable
        )) %>%
        pivot_wider(id_cols = variable, names_from = Ciclo,
                    values_from = c(Estimacion_puntual, CI)) %>%
        select(variable, Estimacion_puntual_1, CI_1,
               Estimacion_puntual_2, CI_2,
               Estimacion_puntual_3, CI_3) %>%
        flextable() %>%
        set_header_labels(values = list(
          variable = "Variable",
          Estimacion_puntual_1 = "OR",
          Estimacion_puntual_2 = "OR",
          Estimacion_puntual_3 = "OR",
          CI_1 = "95% CI",
          CI_2 = "95% CI",
          CI_3 = "95% CI"
        )) %>%
        add_header_row(colwidths = c(1,2,2,2),
                       values = c("","Cycle 1", "Cycle 2", "Cycle 3")) %>%
        save_as_docx(path = glue::glue("paper_tables/{filtertype}_OR_3_ciclos_logitstica_positivo_igg_por_ciclo_{today()}.docx"))
      
      
  }
  
  
  
  #AGE (CONNTINUOUS)---------------------------------------------------------------------------------------
  
  #Get counts
  age <- mydata %>%
    filter(!is.na(edad)) %>%
    select(edad, strata, Ciclo) %>%
    mutate(strata = as.numeric(as.factor(strata))) %>%
    as.matrix
  
  cp_model    <- cmdstan_model("code/additional_functions/age.stan",
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    N           = nrow(age),
    Ncycles     = length(unique(mydata$Ciclo)),
    Ncenters    = length(unique(mydata$strata)),
    ages        = age,
    prior_mu    = 29,
    prior_sigma = 1,
    sigma_hiper = 1,
    pondef      = stan_weights
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 454565456,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.95)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    dplyr::select(matches("mu_all_centers")) 
  
  #Medias
  mean      <- model_draws %>% 
    summarise_draws(~quantile(., c(0.025, 0.5, 0.975))) %>%
    mutate(Ciclo = str_remove_all(variable,"mu_all_centers\\[|\\]")) %>%
    mutate(variable = "Media de edad") %>%
    rename(`Estimacion_puntual`=`50%` ) %>%
    rename(`CI_95%_low` = `2.5%`) %>%
    rename(`CI_95%_high` = `97.5%`) %>%
    write_excel_csv(glue("results/{filtertype}_medias_edad_ciclo_{today()}.csv"))
  
  #Medias
  mean      <- model_sample$draws() %>% 
    as_draws_df() %>%
    dplyr::select(matches("mu_centro_cycle"))  %>% 
    summarise_draws(~quantile(., c(0.025, 0.5, 0.975))) %>%
    mutate(Centro = as.numeric(str_remove_all(variable,"mu_centro_cycle\\[|,[0-9]\\]"))) %>%
    mutate(Ciclo  = as.numeric(str_remove_all(variable,"mu_centro_cycle\\[[0-9],|\\]"))) %>%
    left_join(
      mydata %>% 
        distinct(strata) %>%
        mutate(Centro      = as.numeric(as.factor(strata)))
    ) %>%
    select(-Centro) %>%
    mutate(variable = "Media de edad") %>%
    rename(`Estimacion_puntual`=`50%` ) %>%
    rename(`CI_95%_low` = `2.5%`) %>%
    rename(`CI_95%_high` = `97.5%`) %>%
    write_excel_csv(glue("results/{filtertype}_medias_edad_ciclo_centro{today()}.csv"))
  
  #PCR---------------------------------------------------------------------------------------
  
  #Get counts
  pcr <- mydata %>%
    filter(!is.na(PCR) & PCR != "Inconcluso" & PCR != "SIN PRUEBA") %>%
    group_by(Ciclo, strata, pondef, PCR) %>%
    tally()
  
  #Complete with 0's
  pcr_values <- tibble(PCR = unique(pcr$PCR))
  pcr_values <- expand_grid(pcr_values, weights)
  
  #Expand grid
  pcr_values <- pcr_values %>%
    left_join(pcr, by = c("PCR", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, PCR)
  
  #Get sample size
  stan_ns <- pcr_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n,
                values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get sample size
  stan_positive <- pcr_values %>% 
    filter(PCR == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_sens_spec.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(pcr_values$Ciclo)),
    Ncenters    = length(unique(pcr_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 1.0,
    sigma_hiper = 1.0 / 1000,
    pondef      = stan_weights,
    ns          = stan_ns,
    cases       = stan_positive
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 454565456,
                                  refresh           = 100,
                                  iter_warmup       = 1000,
                                  iter_sampling     = 1000,
                                  adapt_delta       = 0.999,
                                  max_treedepth     = 12)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    dplyr::select(matches("proportion_positive_cycle|prevalence_positive_cycle")) 
  
  #Medias
  mean      <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95) %>% 
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"proportion_positive_cycle\\[|prevalence_positive_cycle\\[|\\]")) %>%
    mutate(variable = if_else(
      str_detect(variable,"proportion_positive_cycle"), "Unadjusted prevalence","Adjusted prevalence"))
  
  ggplot(mean) +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = variable),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = variable),
               position=position_dodge(width=0.5),
               size = 2) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Incidence",
      title = glue("Incidencia de SARS-CoV-2 estimada en {empresa}"),
      subtitle = "Ajuste considerando sensibilidad del 98% y especificidad del 100%",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = c("firebrick","deepskyblue4")) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_PCR_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_PCR_{today()}.csv"))
  
  
  #IGG-----
  
  #Get counts
  igg <- mydata %>%
    filter(!is.na(IGG)) %>%
    group_by(Ciclo, strata, pondef, IGG) %>%
    tally()
  
  #Complete with 0's
  igg_values <- tibble(IGG = unique(igg$IGG))
  igg_values <- expand_grid(igg_values, weights)
  
  #Expand grid
  igg_values <- igg_values %>%
    left_join(igg, by = c("IGG", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, IGG)
  
  #Get sample size
  stan_ns <- igg_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get sample size
  stan_positive <- igg_values %>% 
    filter(IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_sens_spec.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(igg_values$Ciclo)),
    Ncenters    = length(unique(igg_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    sigma_hiper = 1.0 / 1000.0,
    pondef      = stan_weights,
    ns          = stan_ns,
    cases       = stan_positive
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 1000,
                                  iter_sampling     = 1000,
                                  adapt_delta       = 0.99,
                                  max_treedepth     = 14,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    dplyr::select(matches(glue("proportion_positive_cycle|prevalence_positive_cycle",
                        "|prevalence_positive|proportion_positive")))
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95) %>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable, glue("proportion_positive_cycle\\[|",
                                                 "proportion_positive\\[[1-9],|",
                                                 "prevalence_positive\\[[1-9],|",
                                                 "prevalence_positive_cycle\\[|",
                                                 "\\]"))) %>%
    mutate(Ciclo = as.numeric(Ciclo)) %>%
    mutate(workcenter = str_remove_all(variable, glue("proportion_positive\\[|",
                                                  "prevalence_positive\\[|",
                                                 ",[1-9]\\]"))) %>%
    mutate(workcenter = as.numeric(workcenter)) %>%
    left_join(
      complete_dataset %>% ungroup() %>% dplyr::select(Ciclo, strata) %>% 
        distinct(strata) %>% arrange() %>% mutate(workcenter = 1:n()),
      by = c("workcenter")
    ) %>%
    mutate(workcenter = if_else(is.na(workcenter),"Overall Population",strata)) %>%
    dplyr::select(-strata) %>%
    mutate(variable = if_else(
      str_detect(variable,"proportion_positive"), "Unadjusted prevalence","Adjusted prevalence")
    )
  
  ggplot(mean) +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = workcenter, linetype = variable),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = workcenter, shape = variable),
               position=position_dodge(width=0.5),
               size = 2) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Anticuerpos (IgG)",
      title = glue("Prevalencia de anticuerpos para SARS-CoV-2 estimada en {empresa} "),
      subtitle = "Ajuste considerando sensibilidad del 98% y especificidad del 90.7%",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = met.brewer("Egypt", n = 4)) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_IgG_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_IgG_{today()}.csv"))
  
  #SEX-----
  
  #Get counts
  sexo <- mydata %>%
    filter(!is.na(sexo)) %>%
    group_by(Ciclo, strata, pondef, sexo) %>%
    tally()
  
  #Complete with 0's
  sexo_values <- tibble(sexo = unique(sexo$sexo))
  sexo_values <- expand_grid(sexo_values, weights)
  
  #Expand grid
  sexo_values <- sexo_values %>%
    left_join(sexo, by = c("sexo", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, sexo)
  
  #Get sample size
  stan_ns <- sexo_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get sample size
  stan_positive <- sexo_values %>% 
    filter(sexo == "Hombre") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_sens_spec.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(sexo$Ciclo)),
    Ncenters    = length(unique(sexo_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 1.0,
    spec        = 1.0,
    sigma_hiper = 1.0 / 1000,
    pondef      = stan_weights,
    ns          = stan_ns,
    cases       = stan_positive
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 1000,
                                  iter_sampling     = 1000,
                                  adapt_delta       = 0.99,
                                  max_treedepth     = 12,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_positive_cycle|prevalence_positive_cycle")) 
  
  #Medias
  mean      <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"proportion_positive_cycle\\[|prevalence_positive_cycle\\[|\\]")) %>%
    mutate(variable = if_else(
      str_detect(variable,"proportion_positive_cycle"), "Unadjusted prevalence","Adjusted prevalence"))
  
  mean %>%
    filter(str_detect(variable,"Unadjusted")) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high),
                  color = "firebrick",
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median),
               color = "firebrick",
               size = 2) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Porcentaje de hombres",
      title = glue("Proporción de hombres estimada en {empresa} "),
      subtitle = "Estimación a partir de la muestra",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_Hombres_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% filter(str_detect(variable,"Unadjusted")) %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>%
    write_excel_csv(glue("results/{filtertype}_Hombres_{today()}.csv"))
  
  #RD-STDC----
  
  #Get counts
  itt <- mydata %>%
    filter(!is.na(`RD-STDC`)) %>%
    group_by(Ciclo, strata, pondef, `RD-STDC`) %>%
    tally()
  
  #Complete with 0's
  itt_values <- tibble(`RD-STDC` = unique(itt$`RD-STDC`))
  itt_values <- expand_grid(itt_values, weights)
  
  #Expand grid
  itt_values <- itt_values %>%
    left_join(itt, by = c("RD-STDC", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, `RD-STDC`)
  
  #Get sample size
  stan_ns <- itt_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get sample size
  stan_positive <- itt_values %>% 
    filter(`RD-STDC` == 1) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(itt$Ciclo)),
    Ncenters    = length(unique(itt_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 1.0,
    spec        = 1.0,
    pondef      = stan_weights,
    sigma_hiper = 1.0 / 1000,
    ns          = stan_ns,
    cases       = stan_positive
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 54556,
                                  refresh           = 100,
                                  iter_warmup       = 1000,
                                  iter_sampling     = 1000,
                                  adapt_delta       = 0.99,
                                  max_treedepth     = 12,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_positive_cycle|prevalence_positive_cycle")) 
  
  #Medias
  mean      <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"proportion_positive_cycle\\[|prevalence_positive_cycle\\[|\\]")) %>%
    mutate(variable = if_else(
      str_detect(variable,"proportion_positive_cycle"), "Unadjusted prevalence","Adjusted prevalence"))
  
  mean %>%
    filter(str_detect(variable,"Unadjusted")) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high),
                  color = "firebrick",
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median),
               color = "firebrick",
               size = 2) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Porcentaje de incapacidades",
      title = glue("Proporción de incapacidades temporales en el trabajo estimadas en {empresa} "),
      subtitle = "Estimación a partir de la muestra",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_ITT_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% filter(str_detect(variable,"Unadjusted")) %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>%
    write_excel_csv(glue("results/{filtertype}_ITT_{today()}.csv"))
  
  #VAC1-----------------
  
  #Get counts
  vacuna_1 <- mydata %>%
    mutate(Vacuna1 = !str_detect(`Vacuna 1ra`,"^NO\\b")) %>%
    filter(!is.na(Vacuna1)) %>%
    group_by(Ciclo, strata, pondef, Vacuna1) %>%
    tally()
  
  #Complete with 0's
  vacuna_1_values <- tibble(`Vacuna1` = unique(vacuna_1$`Vacuna1`))
  vacuna_1_values <- expand_grid(vacuna_1_values, weights)
  
  #Expand grid
  vacuna_1_values <- vacuna_1_values %>%
    left_join(vacuna_1, by = c("Vacuna1", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, `Vacuna1`)
  
  #Get sample size
  stan_ns <- vacuna_1_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  stan_ns_1_3 <- mydata %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = n()) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  stan_ns[,1:(length(unique(complete_dataset$Ciclo)) - 1)] <- 
    stan_ns_1_3[,1:(length(unique(complete_dataset$Ciclo)) - 1)]
  
  #Get sample size
  stan_positive <- vacuna_1_values %>% 
    filter(Vacuna1) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(vacuna_1_values$Ciclo)),
    Ncenters    = length(unique(vacuna_1_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 1.0,
    spec        = 1.0,
    sigma_hiper = 1.0 / 100,
    pondef      = stan_weights,
    ns          = stan_ns,
    cases       = stan_positive
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 2345,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.9999,
                                  max_treedepth     = 16,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_positive_cycle|prevalence_positive_cycle")) 
  
  #Medias
  mean      <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"proportion_positive_cycle\\[|prevalence_positive_cycle\\[|\\]")) %>%
    mutate(variable = if_else(
      str_detect(variable,"proportion_positive_cycle"), "Unadjusted prevalence","Adjusted prevalence"))
  
  mean %>%
    filter(str_detect(variable,"Unadjusted")) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high),
                  color = "firebrick",
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median),
               color = "firebrick",
               size = 2) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Porcentaje de vacunados",
      title = glue("Proporción de vacunados (al menos 1 dosis) estimadas en {empresa} "),
      subtitle = "Estimación a partir de la muestra",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_Vacunas_AUTORREPPRTE_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    filter(str_detect(variable,"Unadjusted")) %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ floor(. * 100)/100)) %>%
    write_excel_csv(glue("results/{filtertype}_Vacunas_AUTORREPPRTE_{today()}.csv"))
  
  #VAC2----
  
  #Get counts
  vacuna_2 <- mydata %>%
    mutate(Vacuna2 = !str_detect(`Vacuna 2da`,"^NO\\b")) %>%
    filter(!is.na(Vacuna2)) %>%
    group_by(Ciclo, strata, pondef, Vacuna2) %>%
    tally()
  
  #Complete with 0's
  vacuna_2_values <- tibble(Vacuna2 = unique(vacuna_2$Vacuna2))
  vacuna_2_values <- expand_grid(vacuna_2_values, weights)
  
  #Expand grid
  vacuna_2_values <- vacuna_2_values %>%
    left_join(vacuna_2, by = c("Vacuna2", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, `Vacuna2`)
  
  #Get sample size
  stan_ns <- vacuna_2_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  stan_ns_1_3 <- mydata %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = n()) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  stan_ns[,1:(length(unique(complete_dataset$Ciclo)) - 1)] <- 
    stan_ns_1_3[,1:(length(unique(complete_dataset$Ciclo)) - 1)]
  
  
  #Get sample size
  stan_positive <- vacuna_2_values %>% 
    filter(Vacuna2) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(vacuna_1_values$Ciclo)),
    Ncenters    = length(unique(vacuna_1_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 1.0,
    spec        = 1.0,
    sigma_hiper = 1.0 / 100,
    pondef      = stan_weights,
    ns          = stan_ns,
    cases       = stan_positive
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 1000,
                                  iter_sampling     = 1000,
                                  adapt_delta       = 0.99999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_positive_cycle|prevalence_positive_cycle")) 
  
  #Medias
  mean      <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"proportion_positive_cycle\\[|prevalence_positive_cycle\\[|\\]")) %>%
    mutate(variable = if_else(
      str_detect(variable,"proportion_positive_cycle"), "Unadjusted prevalence","Adjusted prevalence"))
  
  
  mean %>%
    filter(str_detect(variable,"Unadjusted")) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high),
                  color = "firebrick",
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median),
               color = "firebrick",
               size = 2) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Porcentaje de vacunados",
      title = glue("Proporción de vacunados (2a dosis) estimadas en {empresa} "),
      subtitle = "Estimación a partir de la muestra",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_Vacunas_2dosis_AUTORREPPRTE_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    filter(str_detect(variable,"Unadjusted")) %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ floor(. * 100)/100)) %>%
    write_excel_csv(glue("results/{filtertype}_Vacunas_2dosis_AUTORREPPRTE_{today()}.csv"))
  
  #IGG | RD-STDC----------------------
  #Get counts
  igg_itt <- mydata %>%
    filter(!is.na(IGG) & !is.na(`RD-STDC`)) %>%
    group_by(Ciclo, strata, pondef, IGG, `RD-STDC`) %>%
    tally()
  
  #Complete with 0's
  igg_values <- tibble(IGG = unique(igg_itt$IGG))
  itt_values <- tibble(`RD-STDC` = unique(igg_itt$`RD-STDC`))
  igg_values <- expand_grid(igg_values, itt_values, weights)
  
  #Expand grid
  igg_values <- igg_values %>%
    left_join(igg_itt, by = c("IGG","RD-STDC","Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, IGG, `RD-STDC`)
  
  #Get sample size for primary
  stan_ns <- igg_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary
  stan_positive_primary <- igg_values %>% 
    filter(`RD-STDC` == 1) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary and secondary
  stan_positive_primary_positive_secondary <- igg_values %>% 
    filter(`RD-STDC` == 1 & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get negative for primary and positive for secondary
  stan_negative_primary_positive_secondary <- igg_values %>% 
    filter(`RD-STDC` != 1 & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_conditional_sens_spec.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(igg_values$Ciclo)),
    Ncenters    = length(unique(igg_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    sigma_hiper = 1.0 / 1000,
    pondef      = stan_weights,
    ns          = stan_ns,
    primary_cases       = stan_positive_primary,
    secondary_positive_cases = stan_positive_primary_positive_secondary,
    secondary_negative_cases = stan_negative_primary_positive_secondary
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 1000,
                                  iter_sampling     = 1000,
                                  adapt_delta       = 0.999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_secondary_cycle_positive|proportion_secondary_cycle_negative"),
           matches("prevalence_secondary_cycle_positive|prevalence_secondary_cycle_negative"),
           matches("proportion_primary_cycle|prevalence_secondary_cycle|proportion_secondary_cycle"),
           matches("conditional_proportion|conditional_prevalence"))
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"^(.*?)\\[|\\]")) %>%
    mutate(Adjusted = if_else(str_detect(variable,"proportion"), "Unadjusted","Adjusted")) %>%
    mutate(variable = case_when(
      str_detect(variable,"proportion_secondary_cycle_positive") ~ "P(IgG|ITT)",
      str_detect(variable,"proportion_secondary_cycle_negative") ~ "P(IgG|¬ ITT)",
      str_detect(variable,"prevalence_secondary_cycle_positive") ~ "P(IgG|ITT)",
      str_detect(variable,"prevalence_secondary_cycle_negative") ~ "P(IgG|¬ ITT)",
      str_detect(variable,"prevalence_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_primary_cycle") ~ "P(ITT)",
      str_detect(variable,"conditional_proportion_primary_secondary_positive") ~ "P(ITT|IgG)",
      str_detect(variable,"conditional_proportion_primary_secondary_negative") ~ "P(ITT|¬ IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_positive") ~ "P(ITT|IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_negative") ~ "P(ITT|¬ IgG)",
    )) %>%
    mutate(conditional = case_when(
      str_detect(variable, "IgG(.*?)ITT") ~ "IgG",
      str_detect(variable, "ITT(.*?)IgG") ~ "ITT",
      str_detect(variable, "P\\(IgG\\)") ~ "IgG",
      str_detect(variable, "P\\(ITT\\)") ~ "ITT"
    )) %>%
    mutate(negation = case_when(
      str_detect(variable, "¬ ITT") ~ "dado que no tuvo ITT",
      str_detect(variable, "¬ IgG") ~ "dado que no tuvo IgG",
      str_detect(variable, "P\\(IgG\\)") ~ "general",
      str_detect(variable, "P\\(ITT\\)") ~ "general",
      str_detect(variable, "ITT\\)") ~ "dado que tuvo ITT",
      str_detect(variable, "IgG\\)") ~ "dado que tuvo IgG"
    )) 
  
  mean %>%
    filter(Adjusted == "Adjusted" | variable == "P(ITT)") %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = negation),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = negation, shape = negation),
               position=position_dodge(width=0.5),
               size = 2) +
    facet_wrap(~conditional, ncol = 1) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Prevalencia",
      title = glue("Prevalencia de anticuerpos e ITT estimada en {empresa} "),
      subtitle = "Ajuste considerando sensibilidad del 98% y especificidad del 90.7%",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = met.brewer("Degas", n = 5)) +
    scale_y_continuous(labels = scales::percent) +
    scale_shape_manual("Prevalencia", values = c(15, 15, 17, 17, 16))
  ggsave(glue("plots_do_not_share/{filtertype}_IgG_ITT_{today()}.pdf"), width = 7, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_IgG_conditional_itt_{today()}.csv"))
  
  
  #IGG | SEX-------------------------
  
  #Get counts
  igg_sexo <- mydata %>%
    filter(!is.na(IGG) & !is.na(sexo)) %>%
    group_by(Ciclo, strata, pondef, IGG, sexo) %>%
    tally()
  
  #Complete with 0's
  igg_values <- tibble(IGG = unique(igg_sexo$IGG))
  sex_values <- tibble(sexo = unique(igg_sexo$sexo))
  igg_values <- expand_grid(igg_values, sex_values, weights)
  
  #Expand grid
  igg_values <- igg_values %>%
    left_join(igg_sexo, by = c("IGG","sexo","Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, IGG, `sexo`)
  
  #Get sample size for primary
  stan_ns <- igg_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary
  stan_positive_primary <- igg_values %>% 
    filter(sexo == "Hombre") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary and secondary
  stan_positive_primary_positive_secondary <- igg_values %>% 
    filter(sexo == "Hombre" & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get negative for primary and positive for secondary
  stan_negative_primary_positive_secondary <- igg_values %>% 
    filter(sexo != "Hombre" & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_conditional_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(igg_values$Ciclo)),
    Ncenters    = length(unique(igg_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    sigma_hiper = 1.0 / 100,
    pondef      = stan_weights,
    ns          = stan_ns,
    primary_cases       = stan_positive_primary,
    secondary_positive_cases = stan_positive_primary_positive_secondary,
    secondary_negative_cases = stan_negative_primary_positive_secondary
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.99999,
                                  max_treedepth     = 16,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_secondary_cycle_positive|proportion_secondary_cycle_negative"),
           matches("prevalence_secondary_cycle_positive|prevalence_secondary_cycle_negative"),
           matches("proportion_primary_cycle|prevalence_secondary_cycle|proportion_secondary_cycle"),
           matches("conditional_proportion|conditional_prevalence")) 
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"^(.*?)\\[|\\]")) %>%
    mutate(Adjusted = if_else(str_detect(variable,"proportion"), "Unadjusted","Adjusted")) %>%
    mutate(variable = case_when(
      str_detect(variable,"proportion_secondary_cycle_positive") ~ "P(IgG|Hombre)",
      str_detect(variable,"proportion_secondary_cycle_negative") ~ "P(IgG|Mujer)",
      str_detect(variable,"prevalence_secondary_cycle_positive") ~ "P(IgG|Hombre)",
      str_detect(variable,"prevalence_secondary_cycle_negative") ~ "P(IgG|Mujer)",
      str_detect(variable,"prevalence_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_primary_cycle") ~ "P(Hombre)",
      str_detect(variable,"conditional_proportion_primary_secondary_positive") ~ "P(Hombre|IgG)",
      str_detect(variable,"conditional_proportion_primary_secondary_negative") ~ "P(Hombre|¬ IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_positive") ~ "P(Hombre|IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_negative") ~ "P(Hombre|¬ IgG)",
    )) %>%
    mutate(conditional = case_when(
      str_detect(variable, "P\\(IgG") ~ "IgG",
      str_detect(variable, "Hombre") ~ "Hombre",
    )) %>%
    mutate(negation = case_when(
      str_detect(variable, "Mujer") ~ "dado que es mujer",
      str_detect(variable, "¬ IgG") ~ "dado que no tuvo IgG",
      str_detect(variable, "P\\(IgG\\)") ~ "general",
      str_detect(variable, "P\\(Hombre\\)") ~ "general",
      str_detect(variable, "Hombre\\)") ~ "dado que es hombre",
      str_detect(variable, "IgG\\)") ~ "dado que tuvo IgG"
    )) 
  
  mujer <- mean %>%
    filter(str_detect(variable, "P\\(Hombre")) %>%
    mutate(across(c(median,CI_low,CI_high), ~ 1 - .)) %>%
    mutate(across(c(variable,conditional), ~ str_replace(.,"Hombre","Mujer")))
  
  mean <- mean %>%
    bind_rows(mujer)
  
  mean %>%
    filter(Adjusted == "Adjusted" | variable == "P(Hombre)" | variable == "P(Mujer)") %>%
    mutate(conditional = factor(conditional, levels = c("Hombre","Mujer","IgG"))) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = negation),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = negation, shape = negation),
               position=position_dodge(width=0.5),
               size = 2) +
    facet_wrap(~conditional, ncol = 1) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Prevalencia",
      title = glue("Prevalencia de anticuerpos por sexo estimada en {empresa} "),
      subtitle = "Ajuste considerando sensibilidad del 98% y especificidad del 90.7%",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = met.brewer("Degas", n = 5)) +
    scale_y_continuous(labels = scales::percent) +
    scale_shape_manual("Prevalencia", values = c(15, 15, 17, 17, 16))
  ggsave(glue("plots_do_not_share/{filtertype}_IgG_sexo_{today()}.pdf"), width = 8, height = 6)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_IgG_conditional_sexo_{today()}.csv"))
  
  #IGG | COVID SYMPTOM-------------------------
    
  #Get counts
  igg_covid_sintoma <- mydata %>%
    filter(!is.na(IGG) & !is.na(covid_sintoma) & 
             covid_sintoma != "No sabe" & covid_sintoma != "No deseo responder") %>%
    group_by(Ciclo, strata, pondef, IGG, covid_sintoma) %>%
    tally()
  
  #Complete with 0's
  igg_values           <- tibble(IGG = unique(igg_covid_sintoma$IGG))
  covid_sintoma_values <- tibble(covid_sintoma = unique(igg_covid_sintoma$covid_sintoma))
  igg_values           <- expand_grid(igg_values, covid_sintoma_values, weights)
  
  #Expand grid
  igg_values <- igg_values %>%
    left_join(igg_covid_sintoma,
              by = c("IGG","covid_sintoma","Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, IGG, `covid_sintoma`)
  
  #Get sample size for primary
  stan_ns <- igg_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary
  stan_positive_primary <- igg_values %>% 
    filter(covid_sintoma == "Si") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary and secondary
  stan_positive_primary_positive_secondary <- igg_values %>% 
    filter(covid_sintoma == "Si" & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get negative for primary and positive for secondary
  stan_negative_primary_positive_secondary <- igg_values %>% 
    filter(covid_sintoma != "Si" & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_conditional_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(igg_values$Ciclo)),
    Ncenters    = length(unique(igg_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    sigma_hiper = 1.0 / 100,
    pondef      = stan_weights,
    ns          = stan_ns,
    primary_cases       = stan_positive_primary,
    secondary_positive_cases = stan_positive_primary_positive_secondary,
    secondary_negative_cases = stan_negative_primary_positive_secondary
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.99999,
                                  max_treedepth     = 16,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_secondary_cycle_positive|proportion_secondary_cycle_negative"),
           matches("prevalence_secondary_cycle_positive|prevalence_secondary_cycle_negative"),
           matches("proportion_primary_cycle|prevalence_secondary_cycle|proportion_secondary_cycle"),
           matches("conditional_proportion|conditional_prevalence")) 
  
  #Medias
  mean         <- model_draws     %>%
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"^(.*?)\\[|\\]")) %>%
    mutate(Adjusted = if_else(str_detect(variable,"proportion"), "Unadjusted","Adjusted")) %>%
    mutate(variable = case_when(
      str_detect(variable,"proportion_secondary_cycle_positive") ~ "P(IgG|Tuvo Enfermedad Respiratoria)",
      str_detect(variable,"proportion_secondary_cycle_negative") ~ "P(IgG|No tuvo Enfermedad Respiratoria)",
      str_detect(variable,"prevalence_secondary_cycle_positive") ~ "P(IgG|Tuvo Enfermedad Respiratoria)",
      str_detect(variable,"prevalence_secondary_cycle_negative") ~ "P(IgG|No tuvo Enfermedad Respiratoria)",
      str_detect(variable,"prevalence_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_primary_cycle") ~ "P(Tuvo Enfermedad Respiratoria)",
      str_detect(variable,"conditional_proportion_primary_secondary_positive") ~ "P(Tuvo Enfermedad Respiratoria|IgG)",
      str_detect(variable,"conditional_proportion_primary_secondary_negative") ~ "P(Tuvo Enfermedad Respiratoria|¬ IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_positive") ~ "P(Tuvo Enfermedad Respiratoria|IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_negative") ~ "P(Tuvo Enfermedad Respiratoria|¬ IgG)",
    )) %>%
    mutate(conditional = case_when(
      str_detect(variable, "P\\(IgG") ~ "IgG",
      str_detect(variable, "Tuvo Enfermedad Respiratoria") ~ "Tuvo Enfermedad Respiratoria",
    )) %>%
    mutate(negation = case_when(
      str_detect(variable, "No tuvo Enfermedad Respiratoria") ~ "dado que no tuvo Enfermedad Respiratoria",
      str_detect(variable, "¬ IgG") ~ "dado que no tuvo IgG",
      str_detect(variable, "P\\(IgG\\)") ~ "general",
      str_detect(variable, "P\\(Tuvo Enfermedad Respiratoria\\)") ~ "general",
      str_detect(variable, "Tuvo Enfermedad Respiratoria\\)") ~ "dado que tuvo Enfermedad Respiratoria",
      str_detect(variable, "IgG\\)") ~ "dado que tuvo IgG"
    )) 
  
  mujer <- mean %>%
    filter(str_detect(variable, "P\\(Tuvo Enfermedad Respiratoria")) %>%
    mutate(across(c(median,CI_low,CI_high), ~ 1 - .)) %>%
    mutate(across(c(variable,conditional), ~ str_replace(.,"Tuvo Enfermedad Respiratoria","No tuvo Enfermedad Respiratoria")))
  
  mean <- mean %>%
    bind_rows(mujer)
  
  mean %>%
    filter(Adjusted == "Adjusted" | variable == "P(Tuvo Enfermedad Respiratoria)" | variable == "P(No tuvo Enfermedad Respiratoria)") %>%
    mutate(conditional = factor(conditional, levels = c("Tuvo Enfermedad Respiratoria","No tuvo Enfermedad Respiratoria","IgG"))) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = negation),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = negation, shape = negation),
               position=position_dodge(width=0.5),
               size = 2) +
    facet_wrap(~conditional, ncol = 1) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Prevalencia",
      title = glue("Prevalencia de anticuerpos por autorreporte de enfermedad respiratoria estimada en {empresa} "),
      subtitle = "Ajuste considerando sensibilidad del 98% y especificidad del 90.7%",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = met.brewer("Degas", n = 5)) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_shape_manual("Prevalencia", values = c(15, 15, 17, 17, 16))
  ggsave(glue("plots_do_not_share/{filtertype}_IgG_autorreporte_Enfermedad Respiratoria_{today()}.pdf"), 
         width = 10, height = 6)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_IgG_autorreporte_covid_{today()}.csv"))
  
  #ITT | COVID SYMPTOM-------------------------
  
  #Get counts
  itt_covid_sintoma <- mydata %>%
    filter(!is.na(`RD-STDC`) & !is.na(covid_sintoma) & 
             covid_sintoma != "No sabe" & covid_sintoma != "No deseo responder") %>%
    group_by(Ciclo, strata, pondef, `RD-STDC`, covid_sintoma) %>%
    tally()
  
  #Complete with 0's
  itt_values           <- tibble(`RD-STDC` = unique(itt_covid_sintoma$`RD-STDC`))
  covid_sintoma_values <- tibble(covid_sintoma = unique(itt_covid_sintoma$covid_sintoma))
  itt_values           <- expand_grid(itt_values, covid_sintoma_values, weights)
  
  #Expand grid
  itt_values <- itt_values %>%
    left_join(itt_covid_sintoma,
              by = c("RD-STDC","covid_sintoma","Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, `RD-STDC`, `covid_sintoma`)
  
  #Get sample size for primary
  stan_ns <- itt_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary
  stan_positive_primary <- itt_values %>% 
    filter(covid_sintoma == "Si") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary and secondary
  stan_positive_primary_positive_secondary <- itt_values %>% 
    filter(covid_sintoma == "Si" & `RD-STDC` == 1) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get negative for primary and positive for secondary
  stan_negative_primary_positive_secondary <- itt_values %>% 
    filter(covid_sintoma != "Si" & `RD-STDC` == 1) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_conditional_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(igg_values$Ciclo)),
    Ncenters    = length(unique(igg_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 1.0,
    spec        = 1.0,
    sigma_hiper = 1.0 / 100,
    pondef      = stan_weights,
    ns          = stan_ns,
    primary_cases       = stan_positive_primary,
    secondary_positive_cases = stan_positive_primary_positive_secondary,
    secondary_negative_cases = stan_negative_primary_positive_secondary
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.99999,
                                  max_treedepth     = 16,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_secondary_cycle_positive|proportion_secondary_cycle_negative"),
           matches("prevalence_secondary_cycle_positive|prevalence_secondary_cycle_negative"),
           matches("proportion_primary_cycle|prevalence_secondary_cycle|proportion_secondary_cycle"),
           matches("conditional_proportion|conditional_prevalence"))
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"^(.*?)\\[|\\]")) %>%
    mutate(Adjusted = if_else(str_detect(variable,"proportion"), "Unadjusted","Adjusted")) %>%
    mutate(variable = case_when(
      str_detect(variable,"proportion_secondary_cycle_positive") ~ "P(ITT|Tuvo Enfermedad Respiratoria)",
      str_detect(variable,"proportion_secondary_cycle_negative") ~ "P(ITT|No tuvo Enfermedad Respiratoria)",
      str_detect(variable,"prevalence_secondary_cycle_positive") ~ "P(ITT|Tuvo Enfermedad Respiratoria)",
      str_detect(variable,"prevalence_secondary_cycle_negative") ~ "P(ITT|No tuvo Enfermedad Respiratoria)",
      str_detect(variable,"prevalence_secondary_cycle") ~ "P(ITT)",
      str_detect(variable,"proportion_secondary_cycle") ~ "P(ITT)",
      str_detect(variable,"proportion_primary_cycle") ~ "P(Tuvo Enfermedad Respiratoria)",
      str_detect(variable,"conditional_proportion_primary_secondary_positive") ~ "P(Tuvo Enfermedad Respiratoria|ITT)",
      str_detect(variable,"conditional_proportion_primary_secondary_negative") ~ "P(Tuvo Enfermedad Respiratoria|¬ ITT)",
      str_detect(variable,"conditional_prevalence_primary_secondary_positive") ~ "P(Tuvo Enfermedad Respiratoria|ITT)",
      str_detect(variable,"conditional_prevalence_primary_secondary_negative") ~ "P(Tuvo Enfermedad Respiratoria|¬ ITT)",
    )) %>%
    mutate(conditional = case_when(
      str_detect(variable, "P\\(ITT") ~ "ITT",
      str_detect(variable, "Tuvo Enfermedad Respiratoria") ~ "Tuvo Enfermedad Respiratoria",
    )) %>%
    mutate(negation = case_when(
      str_detect(variable, "No tuvo Enfermedad Respiratoria") ~ "dado que no tuvo Enfermedad Respiratoria",
      str_detect(variable, "¬ ITT") ~ "dado que no tuvo ITT",
      str_detect(variable, "P\\(ITT\\)") ~ "general",
      str_detect(variable, "P\\(Tuvo Enfermedad Respiratoria\\)") ~ "general",
      str_detect(variable, "Tuvo Enfermedad Respiratoria\\)") ~ "dado que tuvo Enfermedad Respiratoria",
      str_detect(variable, "ITT\\)") ~ "dado que tuvo ITT"
    )) 
  
  mujer <- mean %>%
    filter(str_detect(variable, "P\\(Tuvo Enfermedad Respiratoria")) %>%
    mutate(across(c(median,CI_low,CI_high), ~ 1 - .)) %>%
    mutate(across(c(variable,conditional), ~ str_replace(.,"Tuvo Enfermedad Respiratoria","No tuvo Enfermedad Respiratoria")))
  
  mean <- mean %>%
    bind_rows(mujer)
  
  mean %>%
    filter(Adjusted == "Adjusted" | variable == "P(Tuvo Enfermedad Respiratoria)" | variable == "P(No tuvo Enfermedad Respiratoria)") %>%
    mutate(conditional = factor(conditional, levels = c("Tuvo Enfermedad Respiratoria","No tuvo Enfermedad Respiratoria","ITT"))) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = negation),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = negation, shape = negation),
               position=position_dodge(width=0.5),
               size = 2) +
    facet_wrap(~conditional, ncol = 1) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Prevalencia",
      title = glue("Prevalencia de incapacidades por autorreporte de enfermedad respiratoria estimada en {empresa} "),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = met.brewer("Degas", n = 5)) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_shape_manual("Prevalencia", values = c(15, 15, 17, 17, 16))
  ggsave(glue("plots_do_not_share/{filtertype}_ITT_autorreporte_Enfermedad Respiratoria_{today()}.pdf"), 
         width = 10, height = 6)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_ITT_autorreporte_covid_{today()}.csv"))
  
  #IGG | COVID SYMPTOM & ITT-------------------------
  
  #Get counts
  igg_itt_covid_sintoma <- mydata %>%
    filter(!is.na(`RD-STDC`) & !is.na(IGG) & !is.na(covid_sintoma) & 
             covid_sintoma != "No sabe" & covid_sintoma != "No deseo responder") %>%
    mutate(symptomatic_itt = `RD-STDC` == 1 & covid_sintoma == "Si") %>%
    group_by(Ciclo, strata, pondef, symptomatic_itt, IGG) %>%
    tally()
  
  #Complete with 0's
  igg_values           <- tibble(IGG = unique(igg_itt_covid_sintoma$IGG))
  sympyitt_values      <- tibble(symptomatic_itt = unique(igg_itt_covid_sintoma$symptomatic_itt))
  igg_values           <- expand_grid(igg_values, sympyitt_values, weights)
  
  #Expand grid
  igg_values <- igg_values %>%
    left_join(igg_itt_covid_sintoma,
              by = c("symptomatic_itt","IGG","Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, `symptomatic_itt`, `IGG`)
  
  #Get sample size for primary
  stan_ns <- igg_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary
  stan_positive_primary <- igg_values %>% 
    filter(symptomatic_itt) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary and secondary
  stan_positive_primary_positive_secondary <- igg_values %>% 
    filter(symptomatic_itt & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get negative for primary and positive for secondary
  stan_negative_primary_positive_secondary <- igg_values %>% 
    filter(!symptomatic_itt & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_conditional_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(igg_values$Ciclo)),
    Ncenters    = length(unique(igg_values$strata)),
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    sigma_hiper = 1.0 / 100,
    pondef      = stan_weights,
    ns          = stan_ns,
    primary_cases       = stan_positive_primary,
    secondary_positive_cases = stan_positive_primary_positive_secondary,
    secondary_negative_cases = stan_negative_primary_positive_secondary
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 23576,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.99999,
                                  max_treedepth     = 16,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_secondary_cycle_positive|proportion_secondary_cycle_negative"),
           matches("prevalence_secondary_cycle_positive|prevalence_secondary_cycle_negative"),
           matches("proportion_primary_cycle|prevalence_secondary_cycle|proportion_secondary_cycle"),
           matches("conditional_proportion|conditional_prevalence")) 
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"^(.*?)\\[|\\]")) %>%
    mutate(Adjusted = if_else(str_detect(variable,"proportion"), "Unadjusted","Adjusted")) %>%
    mutate(variable = case_when(
      str_detect(variable,"proportion_secondary_cycle_positive") ~ "P(IgG|Tuvo Enfermedad Respiratoria e ITT)",
      str_detect(variable,"proportion_secondary_cycle_negative") ~ "P(IgG|No Tuvo Enfermedad Respiratoria ni ITT)",
      str_detect(variable,"prevalence_secondary_cycle_positive") ~ "P(IgG|Tuvo Enfermedad Respiratoria e ITT)",
      str_detect(variable,"prevalence_secondary_cycle_negative") ~ "P(IgG|No Tuvo Enfermedad Respiratoria ni ITT)",
      str_detect(variable,"prevalence_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_primary_cycle") ~ "P(Tuvo Enfermedad Respiratoria e ITT)",
      str_detect(variable,"conditional_proportion_primary_secondary_positive") ~ "P(Tuvo Enfermedad Respiratoria e ITT|IgG)",
      str_detect(variable,"conditional_proportion_primary_secondary_negative") ~ "P(Tuvo Enfermedad Respiratoria e ITT|¬ IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_positive") ~ "P(Tuvo Enfermedad Respiratoria e ITT|IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_negative") ~ "P(Tuvo Enfermedad Respiratoria e ITT|¬ IgG)",
    )) %>%
    mutate(conditional = case_when(
      str_detect(variable, "P\\(IgG") ~ "IgG",
      str_detect(variable, "Tuvo Enfermedad Respiratoria e ITT") ~ "Tuvo Enfermedad Respiratoria e ITT",
    )) %>%
    mutate(negation = case_when(
      str_detect(variable, "No Tuvo Enfermedad Respiratoria ni ITT") ~ "dado que No Tuvo Enfermedad Respiratoria ni ITT",
      str_detect(variable, "¬ IgG") ~ "dado que no tuvo IgG",
      str_detect(variable, "P\\(IgG\\)") ~ "general",
      str_detect(variable, "P\\(Tuvo Enfermedad Respiratoria e ITT\\)") ~ "general",
      str_detect(variable, "Tuvo Enfermedad Respiratoria e ITT\\)") ~ "dado que Tuvo Enfermedad Respiratoria e ITT",
      str_detect(variable, "IgG\\)") ~ "dado que tuvo IgG"
    )) 
  
  mujer <- mean %>%
    filter(str_detect(variable, "P\\(Tuvo Enfermedad Respiratoria e ITT")) %>%
    mutate(across(c(median,CI_low,CI_high), ~ 1 - .)) %>%
    mutate(across(c(variable,conditional), ~ str_replace(.,"Tuvo Enfermedad Respiratoria e ITT","No Tuvo Enfermedad Respiratoria ni ITT")))
  
  mean <- mean %>%
    bind_rows(mujer)
  
  mean %>%
    filter(Adjusted == "Adjusted" | variable == "P(Tuvo Enfermedad Respiratoria e ITT)" | variable == "P(No Tuvo Enfermedad Respiratoria ni ITT)") %>%
    mutate(conditional = factor(conditional, levels = c("Tuvo Enfermedad Respiratoria e ITT","No Tuvo Enfermedad Respiratoria ni ITT","IgG"))) %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = negation),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = negation, shape = negation),
               position=position_dodge(width=0.5),
               size = 2) +
    facet_wrap(~conditional, ncol = 1) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Prevalencia",
      title = glue("Prevalencia de incapacidades por autorreporte de enfermedad respiratoria estimada en {empresa} "),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = met.brewer("Degas", n = 5)) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_shape_manual("Prevalencia", values = c(15, 15, 17, 17, 16))
  ggsave(glue("plots_do_not_share/{filtertype}_ITT_autorreporte_Enfermedad Respiratoria_itt_{today()}.pdf"), 
         width = 10, height = 6)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_ITT_autorreporte_covid_itt_{today()}.csv"))
  
  
  #IGG | VACCINE 1---------------------------
  
  
  #Get counts
  igg_vacuna <- mydata %>%
    filter(Ciclo == !!max(mydata$Ciclo)) %>%
    mutate(Vacuna1 = !str_detect(`Vacuna 1ra`,"^NO\\b")) %>%
    filter(!is.na(IGG) & !is.na(Vacuna1)) %>%
    group_by(Ciclo, strata, pondef, IGG, `Vacuna1`) %>%
    tally()
  
  #Complete with 0's
  igg_values <- tibble(IGG = unique(igg_vacuna$IGG))
  vac_values <- tibble(`Vacuna1` = unique(igg_vacuna$Vacuna1))
  igg_values <- expand_grid(igg_values, vac_values, weights %>% filter(Ciclo == !!max(mydata$Ciclo)))
  
  #Expand grid
  igg_values <- igg_values %>%
    left_join(igg_vacuna, by = c("IGG","Vacuna1","Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, IGG, `Vacuna1`)
  
  #Get sample size for primary
  stan_ns <- igg_values %>% 
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary
  stan_positive_primary <- igg_values %>% 
    filter(`Vacuna1`) %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get positive for primary and secondary
  stan_positive_primary_positive_secondary <- igg_values %>% 
    filter(`Vacuna1` == 1 & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Get negative for primary and positive for secondary
  stan_negative_primary_positive_secondary <- igg_values %>% 
    filter(`Vacuna1` != 1 & IGG == "Positivo") %>%
    group_by(Ciclo, strata) %>%
    summarise(n = sum(n)) %>%
    pivot_wider(id_cols = strata, names_from = Ciclo, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    dplyr::select(-`strata`) %>%
    as.matrix
  
  #Model
  cp_model    <- cmdstan_model("code/additional_functions/betabernoulli_conditional_sens_spec.stan", 
                                 include_paths = ".", stanc_options = list("O1"), force_recompile = F)
  
  datos  <- list(
    Ncycles     = length(unique(igg_values$Ciclo)),
    Ncenters    = length(unique(igg_values$strata)),
    prior_alpha = 10,
    prior_beta  = 0.5,
    sens        = 0.98,
    spec        = 0.907,
    sigma_hiper = 1.0,
    pondef      = stan_weights[,ncol(stan_weights)] %>% as.matrix,
    ns          = stan_ns,
    primary_cases       = stan_positive_primary,
    secondary_positive_cases = stan_positive_primary_positive_secondary,
    secondary_negative_cases = stan_negative_primary_positive_secondary
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 6757,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.99999,
                                  max_treedepth     = 16,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion_secondary_cycle_positive|proportion_secondary_cycle_negative"),
           matches("prevalence_secondary_cycle_positive|prevalence_secondary_cycle_negative"),
           matches("proportion_primary_cycle|prevalence_secondary_cycle|proportion_secondary_cycle"),
           matches("conditional_proportion|conditional_prevalence")) 
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = "4") %>%
    mutate(Adjusted = if_else(str_detect(variable,"proportion"), "Unadjusted","Adjusted")) %>%
    mutate(variable = case_when(
      str_detect(variable,"proportion_secondary_cycle_positive") ~ "P(IgG|Vac)",
      str_detect(variable,"proportion_secondary_cycle_negative") ~ "P(IgG|¬ Vac)",
      str_detect(variable,"prevalence_secondary_cycle_positive") ~ "P(IgG|Vac)",
      str_detect(variable,"prevalence_secondary_cycle_negative") ~ "P(IgG|¬ Vac)",
      str_detect(variable,"prevalence_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_secondary_cycle") ~ "P(IgG)",
      str_detect(variable,"proportion_primary_cycle") ~ "P(Vac)",
      str_detect(variable,"conditional_proportion_primary_secondary_positive") ~ "P(Vac|IgG)",
      str_detect(variable,"conditional_proportion_primary_secondary_negative") ~ "P(Vac|¬ IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_positive") ~ "P(Vac|IgG)",
      str_detect(variable,"conditional_prevalence_primary_secondary_negative") ~ "P(Vac|¬ IgG)",
    )) %>%
    mutate(conditional = case_when(
      str_detect(variable, "IgG(.*?)Vac") ~ "IgG",
      str_detect(variable, "Vac(.*?)IgG") ~ "Vacunado (al menos 1 dosis)",
      str_detect(variable, "P\\(IgG\\)") ~ "IgG",
      str_detect(variable, "P\\(Vac\\)") ~ "Vacunado (al menos 1 dosis)"
    )) %>%
    mutate(negation = case_when(
      str_detect(variable, "¬ Vac") ~ "dado que no está vacunado (al menos 1 dosis)",
      str_detect(variable, "¬ IgG") ~ "dado que no tuvo IgG",
      str_detect(variable, "P\\(IgG\\)") ~ "general",
      str_detect(variable, "P\\(Vac\\)") ~ "general",
      str_detect(variable, "Vac\\)") ~ "dado que está vacunado (al menos 1 dosis)",
      str_detect(variable, "IgG\\)") ~ "dado que tuvo IgG"
    )) 
  
  mean %>%
    filter(Adjusted == "Adjusted" | variable == "P(Vac)") %>%
    ggplot() +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = negation),
                  position=position_dodge(width=0.5),
                  width = 0.3) +
    geom_point(aes(x = Ciclo, y = median, color = negation, shape = negation),
               position=position_dodge(width=0.5),
               size = 2) +
    facet_wrap(~conditional, ncol = 1) +
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Prevalencia",
      title = glue("Prevalencia de anticuerpos y vacunación (al menos 1 dosis) estimada en {empresa} "),
      subtitle = "Ajuste considerando sensibilidad del 98% y especificidad del 90.7%",
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Prevalencia", values = met.brewer("Degas", n = 5)) +
    scale_y_continuous(labels = scales::percent) +
    scale_shape_manual("Prevalencia", values = c(15, 15, 17, 17, 16))
  ggsave(glue("plots_do_not_share/{filtertype}_IgG_Vac_{today()}.pdf"), width = 8, height = 5)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_IgG_conditional_vaccine_{today()}.csv"))
  
  #AGE---------------------
  
  #Get counts
  edad <- mydata %>%
    filter(!is.na(edad_cat)) %>%
    group_by(Ciclo, strata, pondef, edad_cat) %>%
    tally()
  
  #Complete with 0's
  edad_values <- tibble(edad_cat = unique(complete_dataset$edad_cat))
  edad_values <- expand_grid(edad_values, weights) %>%
    filter(!is.na(edad_cat))
  
  #Expand grid
  edad_values <- edad_values %>%
    left_join(edad, by = c("edad_cat", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, edad_cat)
  
  #Get sample size
  nciclos     <- length(unique(edad_values$Ciclo))
  ngrupos     <- sum(!is.na(unique(complete_dataset$edad_cat)))
  nstrata     <- length(unique(edad_values$strata))
  
  stan_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_groups <- provideDimnames(stan_groups , sep = "_", 
                                 base = list("Ciclo","strata","edadcat"))
  
  for (ciclo in 1:max(edad_values$Ciclo)){
    stan_groups[ciclo,,] <- edad_values %>% 
      filter(Ciclo == ciclo) %>%
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = edad_cat, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(`Menores a 35`,`35 a 44`,`45 a 59`,`60 y más`) %>%
      as.matrix
  }
  
  cp_model    <- cmdstan_model("code/additional_functions/multinomial.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = nciclos,
    Ncenters    = nstrata,
    Groups      = ngrupos,
    pondef      = stan_weights,
    cases       = stan_groups,
    kappa_prior = 1.0,
    k_sigma = 2.5,
    k_mu = 0.0,
    sigma_1 = 0.0,
    sigma_2 = 2.5
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 43534,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion")) 
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"proportions_per_cycle\\[|,(.*)\\]")) %>%
    select(-variable) %>%
    arrange(Ciclo) %>%
    mutate(Edad_cat = rep(c("Menores a 35","35 a 44","45 a 59","60 y más"), nciclos)) %>%
    mutate(Edad_cat = factor(Edad_cat, levels = c("Menores a 35","35 a 44","45 a 59","60 y más")))
  
  ggplot(mean) +
    geom_errorbar(aes(x = Edad_cat, ymin = CI_low, ymax = CI_high, color = Edad_cat),
                  width = 0.3) +
    geom_point(aes(x = Edad_cat, y = median, color = Edad_cat),
               size = 2) +
    facet_wrap(~Ciclo) + 
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Porcentaje en\ngrupo de edad",
      title = glue("Composición etaria en {empresa} "),
      subtitle = glue("Ciclos 1 a {nciclos}"),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Grupo etario", values = met.brewer("Degas", n = 4)) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_edad_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_edad_{today()}.csv"))
  
  #VACCINES---------------------
  
  #Get counts
  vacunas_data <- mydata %>%
    filter(Ciclo == 4) %>%
    mutate(Vaccine = if_else(!is.na(vacuna), vacuna, `Vacuna 1ra`)) %>%
    mutate(Vaccine = str_replace_all(toupper(Vaccine)," ", "")) %>%
    filter(!is.na(Vaccine))
  
  vacunas <- vacunas_data %>%
    group_by(strata, pondef, Vaccine) %>%
    tally()
  
  #Complete with 0's
  vacuna_values <- tibble(Vaccine = unique(toupper(c(unique(complete_dataset$vacuna), 
                                    unique(complete_dataset$`Vacuna 1ra`))))) %>%
                    mutate(Vaccine = str_replace_all(toupper(Vaccine), " ", "")) %>%
    distinct()
  vacuna_values <- expand_grid(vacuna_values, weights) %>%
    filter(!is.na(Vaccine)) %>%
    filter(Ciclo == 4)
  
  #Expand grid
  vacuna_values <- vacuna_values %>%
    left_join(vacunas, by = c("Vaccine", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(strata, Vaccine)
  
  #Get sample size
  nciclos     <- 1
  ngrupos     <- sum(!is.na(unique(vacuna_values$Vaccine)))
  nstrata     <- length(unique(vacuna_values$strata))
  
  stan_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_groups <- provideDimnames(stan_groups , sep = "_", 
                                 base = list("Ciclo","strata","vacuna_values"))
  
  
  stan_groups[1,,] <- vacuna_values %>% 
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = Vaccine, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(sort(unique(vacuna_values$Vaccine))) %>%
      as.matrix
  
  
  cp_model    <- cmdstan_model("code/additional_functions/multinomial.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = nciclos,
    Ncenters    = nstrata,
    Groups      = ngrupos,
    pondef      = as.matrix(stan_weights[,4]),
    cases       = stan_groups,
    kappa_prior = 1.0,
    k_sigma = 2.5,
    k_mu = 0.0,
    sigma_1 = 0.0,
    sigma_2 = 2.5
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 43534,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion")) 
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = 4) %>%
    select(-variable) %>%
    arrange(Ciclo) %>%
    mutate(Vaccines = sort(unique(vacuna_values$Vaccine))) %>%
    mutate(Vaccines = if_else(str_detect(Vaccines,"VACUNADEMARCANOESPECIFICADA"),
                              "Unspecified", Vaccines)) %>%
    mutate(Vaccines = if_else(str_detect(Vaccines,"SPUTNIKV"),
                              "SPUTNIK V", Vaccines)) %>%
    mutate(Vaccines = str_to_sentence(Vaccines)) %>%
    mutate(Vaccines = factor(Vaccines, 
                             levels = c("Astrazeneca","Cansino","Janssen","Moderna",
                                        "Pfizer","Sinovac","Sputnik v","Unspecified", "No")))
  
  ggplot(mean) +
    geom_errorbar(aes(x = Vaccines, ymin = CI_low, ymax = CI_high, color = Vaccines),
                  width = 0.3) +
    geom_point(aes(x = Vaccines, y = median, color = Vaccines),
               size = 2) +
    facet_wrap(~Ciclo) + 
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Porcentaje por grupo de vacuna",
      title = glue("Composición de vacunas en {empresa} "),
      subtitle = glue("Ciclo 4"),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Grupo etario", values = met.brewer("Degas", n = 9)) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_vacunatipo_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_vacunatipo_{today()}.csv"))
  
  #ICD-10 CODES---------------------
  mydata <- mydata %>%
    mutate(ICD10 = case_when(
      ULT_COD_CIE %in% c("U070", "U071", "U072", "U07E", "U07S","U07D","U099", "B342", "B972",
                         "U129","U099") ~ "COVID-19",
      str_detect(ULT_COD_CIE,"J10|J11") ~ "Influenza",
      str_detect(ULT_COD_CIE,"J01|J04|J05|J06|J20|J21") ~ "Acute Respiratory Illness",
      str_detect(ULT_COD_CIE,"J12|J13|J14|J15|J16|J17|J18") ~ "Pneumonia",
      ULT_COD_CIE %in% c("J029", "J00X", "J02X", "J039", "J22X") ~ "Other related respiratory illnesses",
      is.na(ULT_COD_CIE) ~ NA_character_,
      TRUE ~ ULT_COD_CIE
    ))
  
  #Get counts
  icd10 <- mydata %>%
    filter(!is.na(ICD10)) %>%
    group_by(Ciclo, strata, pondef, ICD10) %>%
    tally()
  
  #Complete with 0's
  icd10_values <- tibble(ICD10 = c("COVID-19","Influenza","Acute Respiratory Illness",
                                   "Pneumonia","Other related respiratory illnesses"))
  icd10_values <- expand_grid(icd10_values, weights) %>%
    filter(!is.na(ICD10))
  
  #Expand grid
  icd10_values <- icd10_values %>%
    left_join(icd10, by = c("ICD10", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, ICD10)
  
  #Get sample size
  nciclos     <- length(unique(icd10_values$Ciclo))
  ngrupos     <- length(unique(icd10_values$ICD10))
  nstrata     <- length(unique(icd10_values$strata))
  
  stan_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_groups <- provideDimnames(stan_groups , sep = "_", 
                                 base = list("Ciclo","strata","ICD10"))
  
  for (ciclo in 1:max(icd10_values$Ciclo)){
    stan_groups[ciclo,,] <- icd10_values %>% 
      filter(Ciclo == ciclo) %>%
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = ICD10, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(`Acute Respiratory Illness`, `COVID-19`, Influenza, 
             `Other related respiratory illnesses`, Pneumonia) %>%
      as.matrix
  }
  
  cp_model    <- cmdstan_model("code/additional_functions/multinomial.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = nciclos,
    Ncenters    = nstrata,
    Groups      = ngrupos,
    pondef      = stan_weights,
    cases       = stan_groups,
    kappa_prior = 1.0,
    k_sigma = 2.5,
    k_mu = 0.0,
    sigma_1 = 0.0,
    sigma_2 = 2.5
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 43534,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("proportion")) 
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = str_remove_all(variable,"proportions_per_cycle\\[|,(.*)\\]")) %>%
    select(-variable) %>%
    arrange(Ciclo) %>%
    mutate(ICD_10 = rep(sort(c("COVID-19","Influenza","Acute Respiratory Illness",
                            "Pneumonia","Other related respiratory illnesses")), nciclos)) %>%
    mutate(ICD_10 = factor(ICD_10, levels = c("COVID-19","Influenza","Acute Respiratory Illness",
                                              "Pneumonia","Other related respiratory illnesses")))
  
  ggplot(mean) +
    geom_errorbar(aes(x = ICD_10, ymin = CI_low, ymax = CI_high, color = ICD_10),
                  width = 0.3) +
    geom_point(aes(x = ICD_10, y = median, color = ICD_10),
               size = 2) +
    facet_wrap(~Ciclo) + 
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Porcentaje de itt\npor diagnóstico",
      title = glue("Composición etaria en {empresa} "),
      subtitle = glue("Ciclos 1 a {nciclos}"),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    theme(axis.text.x = element_blank()) +
    scale_color_manual("Diagnóstico", values = met.brewer("Degas", n = 5)) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_itt_diagnostico_icd10_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_itt_diagnostico_icd10_{today()}.csv"))
  
  #IGG | AGE---------------------
  edad <- mydata %>%
    filter(!is.na(edad_cat)) %>%
    group_by(Ciclo, strata, pondef, edad_cat) %>%
    tally()
  
  igg_edad <- mydata %>%
    filter(!is.na(edad_cat) & IGG == "Positivo") %>%
    group_by(Ciclo, strata, pondef, edad_cat) %>%
    tally()
  
  #Complete with 0's
  edad_values_grid <- tibble(edad_cat = unique(complete_dataset$edad_cat))
  edad_values_grid <- expand_grid(edad_values_grid, weights) %>%
    filter(!is.na(edad_cat))
  
  #Expand grid
  edad_values <- edad_values_grid %>%
    left_join(edad, by = c("edad_cat", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, edad_cat)
  
  igg_edad_values <- edad_values_grid %>%
    left_join(igg_edad, by = c("edad_cat", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, edad_cat)
  
  #Get sample size
  nciclos     <- length(unique(edad_values$Ciclo))
  ngrupos     <- sum(!is.na(unique(complete_dataset$edad_cat)))
  nstrata     <- length(unique(edad_values$strata))
  
  stan_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_groups <- provideDimnames(stan_groups , sep = "_", 
                                 base = list("Ciclo","strata","edadcat"))
  
  for (ciclo in 1:max(edad_values$Ciclo)){
    stan_groups[ciclo,,] <- edad_values %>% 
      filter(Ciclo == ciclo) %>%
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = edad_cat, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(`Menores a 35`,`35 a 44`,`45 a 59`, `60 y más`) %>%
      as.matrix
  }
  
  stan_positive_in_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_positive_in_groups <- provideDimnames(stan_groups , sep = "_", 
                                             base = list("Ciclo","strata","edadcat"))
  
  for (ciclo in 1:max(edad_values$Ciclo)){
    stan_positive_in_groups[ciclo,,] <- ciclo
    stan_positive_in_groups[ciclo,,] <- igg_edad_values %>% 
      filter(Ciclo == ciclo) %>%
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = edad_cat, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(`Menores a 35`,`35 a 44`,`45 a 59`, `60 y más`) %>%
      as.matrix
  }
  
  cp_model    <- cmdstan_model("code/additional_functions/multi_2.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = nciclos,
    Ncenters    = nstrata,
    Groups      = ngrupos,
    pondef      = stan_weights,
    cases       = stan_groups,
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    k_sigma = 2.5, #1/100
    k_mu = 0.0, #150
    sigma_1 = 0.0,
    sigma_hiper = 1.0 / 1000.0,
    sigma_2 = 2.5,
    sigma_3 = 1.0 ,
    sigma_4 = 1.0 ,
    cases_in_group = stan_positive_in_groups,
    kappa_prior = 1.0
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 2,
                                  seed              = 43534,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.999999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("p_adjusted_cycle"))
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = as.numeric(str_remove_all(variable,"p_adjusted_cycle\\[|,(.*?)\\]"))) %>%
    select(-variable) %>%
    arrange(Ciclo) %>%
    mutate(Edad_cat = rep(c("Menores a 35","35 a 44","45 a 59","60 y más"), nciclos)) %>%
    mutate(Edad_cat = factor(Edad_cat, 
                             levels = c("Menores a 35","35 a 44","45 a 59","60 y más")))
  
  ggplot(mean) +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = Edad_cat),
                  width = 0.2, alpha = 0.5) +
    geom_line(aes(x = Ciclo, y = median, color = Edad_cat),
              size = 0.75) +
    geom_point(aes(x = Ciclo, y = median, color = Edad_cat),
               size = 2) +
    geom_point(aes(x = Ciclo, y = median),
               size = 1, color = "white") +
    facet_wrap(~Edad_cat) + 
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Anticuerpos (IgG)",
      title = glue("Proporción por grupo de edad en {empresa} "),
      subtitle = glue("Ciclos 1 a {nciclos}"),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Grupo etario", values = met.brewer("Degas", n = 4)) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_edad_igg_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_edad_igg_{today()}.csv"))
  
  #IGG | VACCINE TYPE---------------------
  
  #Get counts
  vacunas_data <- mydata %>%
    filter(Ciclo == 4) %>%
    filter(!is.na(`Vacuna 1ra`)) %>%
    mutate(Vaccine = if_else(!is.na(vacuna), vacuna, `Vacuna 1ra`)) %>%
    mutate(Vaccine = str_replace_all(toupper(Vaccine)," ", "")) %>%
    filter(!is.na(Vaccine))
  
  vacunas <- vacunas_data %>%
    group_by(strata, pondef, Vaccine) %>%
    tally()
  
  #Complete with 0's
  vacuna_values_grid <- tibble(Vaccine = unique(toupper(c(unique(complete_dataset$vacuna), 
                                                     unique(complete_dataset$`Vacuna 1ra`))))) %>%
    mutate(Vaccine = str_replace_all(toupper(Vaccine), " ", "")) %>%
    distinct()
  vacuna_values_grid <- expand_grid(vacuna_values_grid, weights) %>%
    filter(!is.na(Vaccine)) %>%
    filter(Ciclo == 4)
  
  #Expand grid
  vacuna_values <- vacuna_values_grid %>%
    left_join(vacunas, by = c("Vaccine", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(strata, Vaccine)
  
  igg_vacuna <- mydata %>%
    filter(Ciclo == 4) %>%
    filter(!is.na(`Vacuna 1ra`)) %>%
    mutate(Vaccine = if_else(!is.na(vacuna), vacuna, `Vacuna 1ra`)) %>%
    mutate(Vaccine = str_replace_all(toupper(Vaccine)," ", "")) %>%
    filter(!is.na(Vaccine)) %>% 
    filter(IGG == "Positivo") %>%
    group_by(strata, pondef, Vaccine) %>%
    tally()
  
  igg_vaccine_values <- vacuna_values_grid %>%
    left_join(igg_vacuna, by = c("Vaccine", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(strata, Vaccine)
  
  #Get sample size
  nciclos     <- 1
  ngrupos     <- sum(!is.na(unique(vacuna_values_grid$Vaccine)))
  nstrata     <- length(unique(mydata$strata))
  
  stan_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_groups <- provideDimnames(stan_groups , sep = "_", 
                                 base = list("Ciclo","strata","Vaccine"))
  
  stan_groups[1,,] <- vacuna_values %>% 
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = Vaccine, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(sort(unique(vacuna_values$Vaccine))) %>%
      as.matrix
  
  stan_positive_in_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_positive_in_groups <- provideDimnames(stan_groups , sep = "_", 
                                             base = list("Ciclo","strata","Vaccine"))
  
  
  stan_positive_in_groups[1,,] <- igg_vaccine_values %>% 
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = Vaccine, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(sort(unique(vacuna_values$Vaccine))) %>%
      as.matrix
  
  
  cp_model    <- cmdstan_model("code/additional_functions/multi_2.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = nciclos,
    Ncenters    = nstrata,
    Groups      = ngrupos,
    pondef      = as.matrix(stan_weights[,4]),
    cases       = stan_groups,
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    k_sigma = 2.5, #1/100
    k_mu = 0.0, #150
    sigma_1 = 0.0,
    sigma_hiper = 1.0 / 1000.0,
    sigma_2 = 2.5,
    sigma_3 = 1.0 ,
    sigma_4 = 1.0 ,
    cases_in_group = stan_positive_in_groups,
    kappa_prior = 1.0
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 2,
                                  seed              = 43534,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.999999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("p_adjusted_cycle"))
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = 4) %>%
    select(-variable) %>%
    arrange(Ciclo) %>%
    mutate(Vaccines = sort(unique(vacuna_values$Vaccine))) %>%
    mutate(Vaccines = if_else(str_detect(Vaccines,"VACUNADEMARCANOESPECIFICADA"),
                              "Unspecified", Vaccines)) %>%
    mutate(Vaccines = if_else(str_detect(Vaccines,"SPUTNIKV"),
                              "SPUTNIK V", Vaccines)) %>%
    mutate(Vaccines = str_to_sentence(Vaccines)) %>%
    mutate(Vaccines = factor(Vaccines, 
                             levels = c("Astrazeneca","Cansino","Janssen","Moderna",
                                        "Pfizer","Sinovac","Sputnik v","Unspecified", "No")))
  ggplot(mean) +
    geom_errorbar(aes(x = Vaccines, ymin = CI_low, ymax = CI_high, color = Vaccines),
                  width = 0.2, alpha = 0.5) +
    geom_line(aes(x = Vaccines, y = median, color = Vaccines),
              size = 0.75) +
    geom_point(aes(x = Vaccines, y = median, color = Vaccines),
               size = 2) +
    geom_point(aes(x = Vaccines, y = median),
               size = 1, color = "white") +
    theme_minimal() +
    labs(
      x = "Vacuna",
      y = "Anticuerpos (IgG)",
      title = glue("Proporción por vacuna en {empresa} "),
      subtitle = glue("Ciclo 4"),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Vacuna", values = met.brewer("Degas", n = 9)) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_vac_igg_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_vac_igg_{today()}.csv"))
  
  #IGG | AGE SEX---------------------
  
  #Get counts
  sexedad <- mydata %>%
    filter(!is.na(edad_cat) & !is.na(sexo)) %>%
    mutate(sexdad = paste0(sexo, "_", edad_cat)) %>%
    group_by(Ciclo, strata, pondef, sexdad) %>%
    tally()
  
  igg_sexedad <- mydata %>%
    filter(!is.na(edad_cat)  & !is.na(sexo) & IGG == "Positivo") %>%
    mutate(sexdad = paste0(sexo, "_", edad_cat)) %>%
    group_by(Ciclo, strata, pondef, sexdad) %>%
    tally()
  
  #Complete with 0's
  age_unique <- unique(complete_dataset$edad_cat[!is.na(complete_dataset$edad_cat)])
  sex_unique <- unique(complete_dataset$sexo[!is.na(complete_dataset$sexo)])
  edad_values_grid <- tibble(sexdad = do.call(paste0, 
                                              expand.grid(sex_unique, paste0("_",age_unique))))
  
  edad_values_grid <- expand_grid(edad_values_grid, weights)
  
  #Expand grid
  sexedad_values <- edad_values_grid %>%
    left_join(sexedad, by = c("sexdad", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, sexdad)
  
  igg_sexedad_values <- edad_values_grid %>%
    left_join(igg_sexedad, by = c("sexdad", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, sexdad)
  
  #Get sample size
  nciclos     <- length(unique(sexedad_values$Ciclo))
  ngrupos     <- length(age_unique)*length(sex_unique)
  nstrata     <- length(unique(sexedad_values$strata))
  
  stan_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_groups <- provideDimnames(stan_groups , sep = "_", 
                                 base = list("Ciclo","strata","sexedad"))
  
  for (ciclo in 1:max(sexedad_values$Ciclo)){
    stan_groups[ciclo,,] <- sexedad_values %>% 
      filter(Ciclo == ciclo) %>%
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = sexdad, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(`Hombre_Menores a 35`,`Hombre_35 a 44`,`Hombre_45 a 59`, `Hombre_60 y más`,
             `Mujer_Menores a 35`,`Mujer_35 a 44`,`Mujer_45 a 59`, `Mujer_60 y más`) %>%
      as.matrix
  }
  
  stan_positive_in_groups <- array(NA, dim = c(nciclos, nstrata, ngrupos))
  stan_positive_in_groups <- provideDimnames(stan_groups , sep = "_", 
                                             base = list("Ciclo","strata","sexedad"))
  
  for (ciclo in 1:max(sexedad_values$Ciclo)){
    stan_positive_in_groups[ciclo,,] <- ciclo
    stan_positive_in_groups[ciclo,,] <- igg_sexedad_values %>% 
      filter(Ciclo == ciclo) %>%
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      pivot_wider(id_cols = strata, names_from = sexdad, values_from = n, values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-`strata`) %>%
      select(`Hombre_Menores a 35`,`Hombre_35 a 44`,`Hombre_45 a 59`, `Hombre_60 y más`,
             `Mujer_Menores a 35`,`Mujer_35 a 44`,`Mujer_45 a 59`, `Mujer_60 y más`) %>%
      as.matrix
  }
  
  cp_model    <- cmdstan_model("code/additional_functions/multi_2.stan", 
                               include_paths = ".", stanc_options = list("O1"), 
                               force_recompile = F)
  
  datos  <- list(
    Ncycles     = nciclos,
    Ncenters    = nstrata,
    Groups      = ngrupos,
    pondef      = stan_weights,
    cases       = stan_groups,
    prior_alpha = 1.0,
    prior_beta  = 1.0,
    sens        = 0.98,
    spec        = 0.907,
    k_sigma = 2.5, #1/100
    k_mu = 0.0, #150
    sigma_1 = 0.0,
    sigma_hiper = 1.0 / 1000.0,
    sigma_2 = 2.5,
    sigma_3 = 1.0 ,
    sigma_4 = 1.0 ,
    cases_in_group = stan_positive_in_groups,
    kappa_prior = 1.0
  )
  
  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 2,
                                  seed              = 43534,
                                  refresh           = 100,
                                  iter_warmup       = 500,
                                  iter_sampling     = 500,
                                  adapt_delta       = 0.999999,
                                  max_treedepth     = 13,
                                  threads_per_chain = 4)
  
  model_draws  <- model_sample$draws() %>% 
    as_draws_df() %>%
    select(matches("p_adjusted_cycle")) 
  
  #Medias
  mean         <- model_draws %>% 
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(Ciclo = as.numeric(str_remove_all(variable,"p_adjusted_cycle\\[|,(.*?)\\]"))) %>%
    select(-variable) %>%
    arrange(Ciclo) %>%
    mutate(Edad_cat = rep(c("Hombre_Menores a 35","Hombre_35 a 44","Hombre_45 a 59", 
                            "Hombre_60 y más", "Mujer_Menores a 35","Mujer_35 a 44",
                            "Mujer_45 a 59", "Mujer_60 y más"), nciclos)) %>%
    mutate(Edad_cat = factor(Edad_cat, 
                             levels = c("Hombre_Menores a 35","Hombre_35 a 44","Hombre_45 a 59", 
                                        "Hombre_60 y más", "Mujer_Menores a 35","Mujer_35 a 44",
                                        "Mujer_45 a 59", "Mujer_60 y más")))
  
  ggplot(mean) +
    geom_errorbar(aes(x = Ciclo, ymin = CI_low, ymax = CI_high, color = Edad_cat),
                  width = 0.2, alpha = 0.5) +
    geom_line(aes(x = Ciclo, y = median, color = Edad_cat),
              size = 0.75) +
    geom_point(aes(x = Ciclo, y = median, color = Edad_cat),
               size = 2) +
    geom_point(aes(x = Ciclo, y = median),
               size = 1, color = "white") +
    facet_wrap(~Edad_cat) + 
    theme_minimal() +
    labs(
      x = "Cycle",
      y = "Anticuerpos (IgG)",
      title = glue("Proporción por grupo de edad y sexo en {empresa} "),
      subtitle = glue("Ciclos 1 a {nciclos}"),
      caption = paste0("Elaborada el ", lubridate::today())
    ) + 
    scale_color_manual("Grupo sexo_edad", values = met.brewer("Degas", n = 8)) +
    scale_y_continuous(labels = scales::percent)
  ggsave(glue("plots_do_not_share/{filtertype}_edad_igg_sex_{today()}.pdf"), width = 8, height = 4)
  
  mean %>% 
    mutate(across(c(`median`,CI_low, CI_high), ~ round(., 4))) %>% 
    write_excel_csv(glue("results/{filtertype}_edad_igg_sex_{today()}.csv"))
  
  #RD-STDC ----
  itt <- ittdata %>%
    mutate(epiweek = epiweek(COLAPS_FEC_INICIO)) %>%
    mutate(epiyear = epiyear(COLAPS_FEC_INICIO)) %>%
    group_by(Ciclo, strata, pondef, epiweek, epiyear) %>%
    summarise(n = n())

  #Complete with 0's
  itt_values <- expand_grid(
      epiyear = unique(epiyear(complete_itt$COLAPS_FEC_INICIO)),
      epiweek = unique(epiweek(complete_itt$COLAPS_FEC_INICIO)),
    ) %>%
    filter_all(~!is.na(.)) %>%
    arrange(epiweek, epiyear) %>%
    filter(!(epiweek == 53 & epiyear != 2020)) %>%
    filter(epiyear < 2022 | epiweek < 14) %>%
    mutate(weeknum = epiweek + 52*(epiyear - 2020) + as.numeric(epiyear != 2020))

  itt_values <- expand_grid(itt_values, weights)

  #Expand grid
  itt_values <- itt_values %>%
    left_join(itt, by = c("epiweek","epiyear", "Ciclo", "strata", "pondef")) %>%
    mutate(n = replace_na(n, 0)) %>%
    arrange(Ciclo, strata, epiyear, epiweek, weeknum)

  #Get sample size
  nciclos     <- length(unique(itt_values$Ciclo))
  nweeks      <- max(itt_values$weeknum)
  nstrata     <- length(unique(itt_values$strata))

  stan_ns     <- ittdata %>%
    select(Ciclo, strata, n_Strata) %>%
    distinct() %>%
    pivot_wider(names_from = Ciclo, values_from = n_Strata, values_fill = 0) %>%
    arrange(strata) %>%
    select(-strata) %>%
    as.matrix

  #Number of people with RD-STDC
  stan_cases <- array(NA, dim = c(nweeks, nstrata, nciclos))
  stan_cases <- provideDimnames(stan_cases , sep = "_",
                             base = list("week","strata","cycle"))

  for (ciclo in 1:nciclos){
    stan_cases[,,ciclo] <- itt_values %>%
      filter(Ciclo == ciclo) %>%
      full_join(complete_dataset %>% select(strata) %>% distinct(), by = "strata") %>%
      select(weeknum, strata, n) %>%
      pivot_wider(id_cols = weeknum, names_from = strata, values_from = n, values_fill = 0) %>%
      filter(!is.na(weeknum)) %>%
      arrange(weeknum) %>%
      dplyr::select(-`weeknum`) %>%
      select(`Centro de Distribución`, `Centro de ventas`, `Planta de producción`) %>%
      as.matrix
  }

  #Model
  cp_model    <- cmdstan_model(
    "code/additional_functions/betabernoulli_conditional_sens_spec_time.stan",
                               include_paths = ".", stanc_options = list("O1"), 
    force_recompile = F)

  datos  <- list(
    Ncycles     = nciclos,
    Ncenters    = nstrata,
    Nweeks      = nweeks,
    prior_alpha = 1,
    prior_beta  = 1,
    pondef      = stan_weights,
    sigma_hiper = 1/100,
    ns          = stan_ns,
    cases       = stan_cases
  )

  model_sample <- cp_model$sample(data              = datos,
                                  chains            = 4,
                                  seed              = 54556,
                                  refresh           = 100,
                                  iter_warmup       = 250,
                                  iter_sampling     = 250,
                                  adapt_delta       = 0.8,
                                  max_treedepth     = 12,
                                  threads_per_chain = 4)

  model_draws  <- model_sample$draws() %>%
    as_draws_df() %>%
    select(starts_with("proportion_week[")) 

  #Medias
  mean         <- model_draws %>%
    bayestestR::hdi(., ci = 0.95)%>%
    as_tibble() %>% 
    rename(variable = Parameter) %>%
    left_join( model_draws %>% summarise_draws("median")) %>%
    mutate(semana = as.numeric(str_remove_all(variable,"proportion_week\\[|\\]")))

  meanN <- ittdata %>%
    distinct(N_Total) %>%
    summarise(mean = mean(N_Total)) %>%
    unlist

  igg <- read_csv(glue("results/{filtertype}_IgG_{today()}.csv"))
  igg <- igg %>%
    filter(str_detect(variable,"Adjusted")) %>%
    filter(workcenter == "Overall Population")

  igg_c1 <- igg %>% filter(Ciclo == 1) %>% select(median) %>% unlist %>% scales::percent(accuracy = 0.01)
  igg_c2 <- igg %>% filter(Ciclo == 2) %>% select(median) %>% unlist %>% scales::percent(accuracy = 0.01)
  igg_c3 <- igg %>% filter(Ciclo == 3) %>% select(median) %>% unlist %>% scales::percent(accuracy = 0.01)
  igg_c4 <- igg %>% filter(Ciclo == 4) %>% select(median) %>% unlist %>% scales::percent(accuracy = 0.01)

  #> Plot for RD-STDC ------
  if (generate_plot){
    
    fechas    <- tibble(
      fecha = seq(ymd("2020/01/01"), today(), by = "1 day")
    ) %>%
      mutate(semana = epiweek(fecha)) %>%
      mutate(año    = epiyear(fecha)) %>%
      distinct(semana, año, .keep_all = T)
    
    casos_umf <- read_csv("final_datasets/casos_per_capita_en_unidades.csv") %>%
      mutate(semana = epiweek(COLAPS_FEC_INICIO)) %>%
      mutate(año    = epiyear(COLAPS_FEC_INICIO)) %>%
      group_by(semana, año) %>%
      summarise(n_per_capita = sum(n_per_capita)) %>%
      left_join(fechas) %>%
      filter(año >= 2020) %>%
      filter(fecha <= ymd("2021/11/22"))
    
    casos_empresa <- read_csv("datasets/empresa_by_itself.csv") %>%
      mutate(COLAPS_FEC_INICIO = dmy(COLAPS_FEC_INICIO)) %>%
      mutate(semana = epiweek(COLAPS_FEC_INICIO)) %>%
      mutate(año    = epiyear(COLAPS_FEC_INICIO)) %>%
      group_by(semana, año) %>%
      summarise(n_per_capita = sum(n)) %>%
      mutate(n_per_capita = n_per_capita/13397) %>%
      left_join(fechas) %>%
      filter(año >= 2020) %>%
      filter(fecha <= ymd("2021/11/22"))
    
    mean %>%
      filter(ymd("2020/01/01") + weeks(semana) <= ymd("2021/11/24")) %>%
      ggplot() +
      geom_rect(aes(xmin = ymd("2020/12/20"), xmax = ymd("2021/02/14"),
                    ymin = 0, ymax = 0.025), alpha = 0.75, fill = "gray95") +
      geom_rect(aes(xmin = ymd("2021/02/14"), xmax = ymd("2021/05/01"),
                    ymin = 0, ymax = 0.025), alpha = 0.75, fill = "floralwhite") +
      geom_rect(aes(xmin = ymd("2021/05/01"), xmax = ymd("2021/06/01"),
                    ymin = 0, ymax = 0.025), alpha = 0.75, fill = "gray95") +
      geom_rect(aes(xmin = ymd("2021/06/01"), xmax = ymd("2021/07/01"),
                    ymin = 0, ymax = 0.025), alpha = 0.75, fill = "floralwhite") +
      geom_rect(aes(xmin = ymd("2021/07/01"), xmax = ymd("2021/08/01"),
                    ymin = 0, ymax = 0.025), alpha = 0.75, fill = "gray95") +
      geom_rect(aes(xmin = ymd("2021/08/01"), xmax = ymd("2021/11/22"),
                    ymin = 0, ymax = 0.025), alpha = 0.75, fill = "floralwhite") +
      geom_ribbon(aes(x = ymd("2020/01/01") + weeks(semana),
                      ymin = CI_low, ymax = CI_high), alpha = 0.15,
                  fill = "#006778") +
      geom_line(aes(x = ymd("2020/01/01") + weeks(semana), y = median, 
                    linetype = "RD-STDC in the work centers under study\n(weighted sample with 95% CI)",
                    color = "RD-STDC in the work centers under study\n(weighted sample with 95% CI)"), size = 1) +
      geom_line(aes(x = fecha, y = n_per_capita,
                    linetype = "RD-STDC in same healthcare facilities as\nstudy participants (Source: IMSS)", 
                    color = "RD-STDC in same healthcare facilities as\nstudy participants (Source: IMSS)"),
                size = 1, data = casos_umf) +
      geom_line(aes(x = fecha, y = n_per_capita,
                    linetype = "RD-STDC in the company across all of\ntheir work centers (Source: IMSS)", 
                    color = "RD-STDC in the company across all of\ntheir work centers (Source: IMSS)"),
                size = 1, data = casos_empresa) +
      geom_vline(aes(xintercept =  ymd("2020/09/22")), linetype = "dotted",
                 color = "#00AFC1") +
      geom_vline(aes(xintercept =  ymd("2020/11/09")), linetype = "dotted",
                 color = "#00AFC1") +
      geom_vline(aes(xintercept =  ymd("2021/01/04")), linetype = "dotted",
                 color = "#00AFC1") +
      geom_vline(aes(xintercept = ymd("2021/11/22")), linetype = "dotted",
                 color = "#00AFC1") +
      geom_errorbar(aes(x = ymd("2020/09/22"),
                        ymin = CI_low*0.025, ymax = CI_high*0.025),
                    width = 6,
                    color = "#EE5007",
                    data = igg %>% filter(Ciclo == 1)) +
      geom_point(aes(x = ymd("2020/09/22"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 1), size = 3,
                 color = "#EE5007") +
      annotate("label", x = ymd("2020/09/22"), color = "#EE5007", 
               fill = "white", alpha = 0.95, label.size = NA,
               y = (igg %>% filter(Ciclo == 1) %>% select(median) %>% unlist)*0.025,
               label = glue("IgG+: \n{igg_c1} "), hjust = 1.1, size = 3) +
      annotate("label", x = ymd("2020/11/09"), color = "#EE5007",
               fill = "white", alpha = 0.95, label.size = NA,
               y = (igg %>% filter(Ciclo == 2) %>% select(median) %>% unlist)*0.025,
               label = glue("IgG+: \n{igg_c2} "), hjust = 1.1, size = 3) +
      annotate("label", x = ymd("2021/01/04"), color = "#EE5007", 
               fill = "white", alpha = 0.95, label.size = NA,
               y = (igg %>% filter(Ciclo == 3) %>% select(median) %>% unlist)*0.025,
               label = glue("IgG+: \n{igg_c3} "), hjust = 1.1, size = 3) +
      annotate("label", x = ymd("2021/11/22"), color = "#EE5007",
               fill = "white", alpha = 0.95, label.size = NA,
               y = (igg %>% filter(Ciclo == 4) %>% select(median) %>% unlist)*0.025,
               label = glue("IgG+: {igg_c4} "), hjust = 1.1, size = 3) +
      geom_point(aes(x = ymd("2020/09/22"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 1), size = 1, color = "white") +
      geom_errorbar(aes(x = ymd("2020/11/09"),
                        ymin = CI_low*0.025, ymax = CI_high*0.025),
                    width = 6,
                    color = "#EE5007",
                    data = igg %>% filter(Ciclo == 2)) +
      geom_point(aes(x = ymd("2020/11/09"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 2), size = 3,
                 color = "#EE5007") +
      geom_point(aes(x = ymd("2020/11/09"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 2), size = 1, color = "white") +
      geom_errorbar(aes(x = ymd("2021/11/22"),
                        ymin = CI_low*0.025, ymax = CI_high*0.025),
                    width = 6,
                    color = "#EE5007",
                    data = igg %>% filter(Ciclo == 4)) +
      geom_point(aes(x = ymd("2021/11/22"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 4), size = 3,
                 color = "#EE5007") +
      geom_point(aes(x = ymd("2021/11/22"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 4), size = 1, color = "white") +
      geom_errorbar(aes(x = ymd("2021/01/04"),
                        ymin = CI_low*0.025, ymax = CI_high*0.025),
                    width = 6,
                    color = "#EE5007",
                    data = igg %>% filter(Ciclo == 3)) +
      geom_point(aes(x = ymd("2021/01/04"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 3), size = 3,
                 color = "#EE5007") +
      geom_point(aes(x = ymd("2021/01/04"), y = median*0.025),
                 data = igg %>% filter(Ciclo == 3), size = 1, color = "white") +
      scale_x_date(breaks = "1 month", date_labels = "%b-%y", expand = c(0,3)) +
      theme_classic() +
      labs(
        x = "",
        y = "_Per capita_ incidence"#,
        #caption = glue("[*] RD-STDC from the same healthcare facilities considered the total RD-STDC ",
        #               "registered at IMSS in the same healthcare facilities that the company's workers\n",
        #               "went to. Per capita estimation was done using the total number of ",
        #               "affiliated workers in the units by June 2021."),
        #title = glue("Weekly incidence of <span style='color:#006778'>",
        #             "**Respiratory Disease Short-Term Disability Claims (RD-STDC)**</span> in the company",
        #             " submitted to<br>IMSS (as percent of workers) and RD-STDC",
        #             " registered at IMSS in the same healthcare facilities*."),
      ) +
      theme(axis.text.x        = element_text(angle = 90, hjust = 1),
            plot.title         = element_markdown(size = 16, lineheight = 1.2),
            axis.line.y        = element_line(color = "black"),
            axis.line.y.right  = element_line(color = "#EE5007"),
            axis.text.y        = element_text(color = "black"),
            axis.text.y.right  = element_text(color = "#EE5007"),
            axis.ticks.y       = element_line(color = "black"),
            axis.ticks.y.right = element_line(color = "#EE5007"),
            axis.title.y       = element_markdown(),
            axis.title.y.right = element_markdown()) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 0.025), expand = c(0, 0.001),
                         minor_breaks = seq(0,1, by = 0.1),
                         sec.axis = sec_axis( trans=~.*40,
                                              name=glue("<span style='color:#EE5007'>",
                                              "Percent of the company's workers with ",
                                              "SARS-CoV-2 antibodies</span>"),
                                              labels = scales::percent)) +
      annotate("text", x = ymd("2020/12/20") + 1.1*(ymd("2021/02/14") - ymd("2020/12/20"))/2,
               y = 0.95*0.025,  label = "Vaccination of\nFrontline\nHealthcare\nWorkers", angle = 90,
               hjust = 1, size = 3) +
      annotate("text", x = ymd("2021/02/14") + (ymd("2021/05/01") - ymd("2021/02/14"))/2,
               y = 0.95*0.025, label = "Vaccination of Senior Citizens 60+", angle = 90, hjust = 1,
               size = 3) +
      annotate("text", x = ymd("2021/05/01") + (ymd("2021/06/01") - ymd("2021/05/01"))/2,
               y = 0.95*0.025, label = "Vaccination of Citizens 50 - 59", angle = 90, hjust = 1,
               size = 3) +
      annotate("text", x = ymd("2021/06/01") + (ymd("2021/07/01") - ymd("2021/06/01"))/2,
               y = 0.95*0.025, label = "Vaccination of Citizens 40 - 49", angle = 90, hjust = 1,
               size = 3) +
      annotate("text", x = ymd("2021/07/01") + (ymd("2021/08/01") - ymd("2021/07/01"))/2,
               y = 0.95*0.025, label = "Vaccination of Citizens 30 - 39", angle = 90, hjust = 1,
               size = 3) +
      annotate("text", x = ymd("2021/08/01") + days(15),
             y = 0.95*0.025, label = "Vaccination of\nRemaining Citizens 18+", angle = 90, hjust = 1,
             size = 3) +
      annotate("label", label = "Cycle 2", x =  ymd("2020/11/09"), y = 0.7*0.025, hjust = 1.1,
               color = "black") +
      annotate("label", label = "Cycle 4", x =  ymd("2021/11/22"), y = 0.7*0.025, hjust = 1.1,
               color = "black") +
      annotate("label", label = "Cycle 3", x =  ymd("2021/01/04"), y = 0.7*0.025, hjust = 1.1,
               color = "black") +
      annotate("label", label = "Cycle 1", x =  ymd("2020/09/22"), y = 0.7*0.025, hjust = 1.1,
               color = "black") +
      scale_linetype_manual("Weekly Incidence", values = c("RD-STDC in the work centers under study\n(weighted sample with 95% CI)" = "solid", 
                                                           "RD-STDC in the company across all of\ntheir work centers (Source: IMSS)" = "solid",
                                                           "RD-STDC in same healthcare facilities as\nstudy participants (Source: IMSS)" = "dotted")) +
      scale_color_manual("Weekly Incidence", values = c("RD-STDC in the work centers under study\n(weighted sample with 95% CI)" = "#006778", 
                                                        "RD-STDC in the company across all of\ntheir work centers (Source: IMSS)" = "deepskyblue3",
                                                        "RD-STDC in same healthcare facilities as\nstudy participants (Source: IMSS)" = "black")) +
      theme(legend.position = "bottom") +
      geom_segment(aes(x = ymd("2020/03/23"), xend = ymd("2020/05/31"),
                       y = 0.02, yend = 0.02)) +
      geom_segment(aes(x = ymd("2020/03/23"), xend = ymd("2020/03/23"),
                       y = 0.0201, yend = 0.0199)) +
      geom_segment(aes(x = ymd("2020/05/31"), xend = ymd("2020/05/31"),
                       y = 0.0201, yend = 0.0199)) +
      geom_richtext(aes(x = ymd("2020/03/23") + (ymd("2020/05/31") - ymd("2020/03/23"))/2, 
               y = 0.022), hjust = 0.5, fill = NA, label.color = NA,
               label = "Jornada Nacional<br>de Sana Distancia<br><span style='color:gray25'>_(Voluntary lockdown for<br>non-essential workers)_</span>", size = 3)
    
  
    ggsave(glue("plots_share/{filtertype}_ITT_inc_{today()}.pdf"), width = 13, height = 6.5)
    ggsave(glue("plots_share/{filtertype}_ITT_inc_{today()}.png"), width = 13, height = 6.5, dpi = 750)
  }

  mean %>%
    mutate(across(c(`median`,CI_low, CI_high), ~ round(. * meanN, 4))) %>%
    write_excel_csv(glue("results/{filtertype}_ITT_inc_{today()}.csv"))

  
  #TABLE 1----
  
  # |> Characteristics of sample ------
  total <- mydata %>% 
    group_by(Ciclo) %>% 
    tally() %>% 
    full_join(tibble(Ciclo = unique(complete_dataset$Ciclo)), by = "Ciclo") %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Ciclo", .," \n Sample (as is)"), .cols = everything()) %>%
    mutate("Variable" = "Total") %>%
    mutate(across(everything(), ~ as.character(.)))
  
  strata <- mydata %>% 
    group_by(Ciclo, strata) %>%
    tally() %>% 
    full_join(expand_grid(Ciclo  = unique(complete_dataset$Ciclo),
                          strata = unique(complete_dataset$strata)), by = c("Ciclo","strata")) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Ciclo", .," \n Sample (as is)"), .cols = everything()) %>%
    rename("Variable" = 1) %>%
    mutate(across(everything(), ~ as.character(.)))
  
  age <- mydata %>% 
    group_by(Ciclo, edad_cat) %>% tally() %>% 
    full_join(expand_grid(Ciclo  = unique(complete_dataset$Ciclo),
                          edad_cat = unique(complete_dataset$edad_cat)), by = c("Ciclo","edad_cat")) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Ciclo", .," \n Sample (as is)"), .cols = everything())  %>%
    rename("Variable" = 1) %>%
    filter(!is.na(Variable))  %>%
    mutate(across(everything(), ~ as.character(.)))
  
  sex <- mydata %>% 
    group_by(Ciclo, sexo) %>% tally() %>% 
    full_join(expand_grid(Ciclo  = unique(complete_dataset$Ciclo),
                          sexo = unique(complete_dataset$sexo)), by = c("Ciclo","sexo")) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Ciclo", .," \n Sample (as is)"), .cols = everything())  %>%
    rename("Variable" = 1) %>%
    filter(!is.na(Variable))  %>%
    mutate(across(everything(), ~ as.character(.))) %>%
    mutate(Variable = if_else(str_detect(Variable,"Hombre"), "Men","Women"))
  
  pcr <- mydata %>% 
    group_by(Ciclo, PCR) %>% tally() %>% 
    full_join(expand_grid(Ciclo  = unique(complete_dataset$Ciclo),
                          PCR = unique(complete_dataset$PCR)), by = c("Ciclo","PCR")) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Ciclo", .," \n Sample (as is)"), .cols = everything())  %>%
    rename("Variable" = 1)  %>%
    filter(Variable == "Positivo")  %>%
    mutate("Variable" = "PCR positive") %>%
    mutate(across(everything(), ~ as.character(.)))
  
  igg <- mydata %>% 
    group_by(Ciclo, IGG) %>% tally() %>% 
    full_join(expand_grid(Ciclo  = unique(complete_dataset$Ciclo),
                          IGG = unique(complete_dataset$IGG)), by = c("Ciclo","IGG")) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Ciclo", .," \n Sample (as is)"), .cols = everything())  %>%
    rename("Variable" = 1)  %>%
    filter(Variable != "Negativo")  %>%
    mutate("Variable" = "IgG positive") %>%
    mutate(across(everything(), ~ as.character(.)))
  
  rdstdc <- mydata %>% 
    group_by(Ciclo, `RD-STDC`) %>% tally() %>% 
    full_join(expand_grid(Ciclo  = unique(complete_dataset$Ciclo),
                          `RD-STDC` = unique(complete_dataset$`RD-STDC`)), by = c("Ciclo","RD-STDC")) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Ciclo", .," \n Sample (as is)"), .cols = everything())  %>%
    rename("Variable" = 1)  %>%
    filter(Variable != 0) %>%
    mutate(Variable = "Yes")  %>%
    mutate("Variable" = "RD-STDC in last 6 months") %>%
    mutate(across(everything(), ~ as.character(.)))
  
  # |> Expansion-----
  total_exp <- mydata %>% 
    group_by(Ciclo) %>% distinct(N_Total) %>% 
    rename(n = N_Total) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Cycle", .,"Sample (weighted)"), .cols = everything()) %>%
    mutate("Variable" = "Total") %>%
    mutate(across(everything(), ~ as.character(.)))
  
  strata_exp <- mydata %>% 
    group_by(Ciclo, strata) %>% distinct(N_Strata) %>% 
    rename(n = N_Strata) %>%
    pivot_wider(names_from = Ciclo, values_from = n, values_fill = 0) %>%
    rename_with( ~ paste("Cycle", .,"Sample (weighted)"), .cols = everything()) %>%
    rename("Variable" = 1) %>%
    mutate(across(everything(), ~ as.character(.)))
  
  sexo_exp <- read_csv(glue("results/{filtertype}_Hombres_{today()}.csv")) %>%
    mutate(Variable = "Men") %>%
    bind_rows(read_csv(glue("results/{filtertype}_Hombres_{today()}.csv")) %>%
                mutate(median = 1 - median) %>%
                mutate(`2.5%_a` = 1 - CI_high) %>%
                mutate(CI_high = 1 - CI_low) %>%
                mutate(CI_low = `2.5%_a`)  %>%
                mutate(Variable = "Women")) %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo, Variable) %>%
    pivot_wider(names_from = Ciclo, values_from = c(estim, CI), values_fill = "")  %>%
    rename_with( ~ paste("Ciclo", .,"Sample (weighted)"), .cols = everything()) %>%
    rename("Variable" = 1)  
  
  edad_exp <- read_csv(glue("results/{filtertype}_edad_{today()}.csv")) %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo, Edad_cat) %>%
    pivot_wider(names_from = Ciclo, values_from = c(estim, CI), values_fill = "")  %>%
    rename_with( ~ paste("Ciclo", .,"Sample (weighted)"), .cols = everything()) %>%
    rename("Variable" = 1)
  
  igg_exp <- read_csv(glue("results/{filtertype}_IgG_{today()}.csv")) %>%
    filter(str_detect(variable,"Adjusted")) %>%
    filter(str_detect(workcenter,"Overall")) %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo) %>%
    pivot_wider(names_from = Ciclo, values_from = c(estim, CI), values_fill = "")  %>%
    rename_with( ~ paste("Ciclo", .,"Sample (weighted)"), .cols = everything()) %>%
    mutate("Variable" = "IgG positive")
  
  itt_exp <- read_csv(glue("results/{filtertype}_ITT_{today()}.csv")) %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo) %>%
    pivot_wider(names_from = Ciclo, values_from = c(estim, CI), values_fill = "")  %>%
    rename_with( ~paste("Ciclo", .,"Sample (weighted)"), .cols = everything()) %>%
    mutate("Variable" = "RD-STDC in last 6 months")
  
  pcr_exp <- read_csv(glue("results/{filtertype}_PCR_{today()}.csv")) %>%
    filter(str_detect(variable,"Adjusted")) %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo) %>%
    pivot_wider(names_from = Ciclo, values_from = c(estim, CI), values_fill = "")  %>%
    rename_with( ~ paste("Ciclo", .,"Sample (weighted)"), .cols = everything()) %>%
    mutate("Variable" = "PCR positive")
  
  sample <- (total %>% left_join(total_exp)) %>%
    bind_rows(strata %>% left_join(strata_exp)) %>%
    bind_rows(age %>% left_join(edad_exp)) %>%
    bind_rows(pcr %>% left_join(pcr_exp)) %>%
    bind_rows(igg %>% left_join(igg_exp)) %>%
    bind_rows(rdstdc %>% left_join(itt_exp)) %>%
    bind_rows(sex %>% left_join(sexo_exp)) %>%
    select(Variable, everything()) %>%
    mutate(Category = c("Total","Strata","Strata","Strata",
                        rep("Age", 4), "PCR","IgG","RD-STDC (last 6 months)", rep("Sex",2))) %>%
    select(Category, everything()) %>%
    mutate(across(c(3:6), ~ scales::comma(as.numeric(.)))) %>%
    mutate(across(c(7:10), ~ if_else(str_detect(Category,"Strata|Total"),
                                     scales::comma(as.numeric(.)), .))) %>%
    mutate(across(everything(), ~ replace_na(., "0")))
  
  sample %>% 
    filter(Category != "Strata") %>%
    select(-c(7:10)) %>%
    select(-Category) %>%
    mutate(Variable = case_when(
      Variable == "35 a 44" ~ "35 - 44",
      Variable == "45 a 59" ~ "45 - 59",
      Variable == "60 y más" ~ "60 over",
      Variable == "Menores a 35" ~ "Under 35",
      Variable == "covid symptom (yes)" ~ "covid symptom",
      TRUE ~ Variable
    )) %>%
    mutate(across(6:13, ~if_else(Variable == "Total", "100", .))) %>%
    mutate(Variable = factor(Variable,
                             levels = c("Total","Men","Women",
                                        "Under 35","35 - 44","45 - 59",
                                        "60 over","PCR positive","IgG positive",
                                        "RD-STDC in last 6 months"))) %>%
    arrange(Variable) %>%
    select(Variable,matches("1"),matches("2"),matches("3"),matches("4")) %>%
    flextable() %>%
    set_header_labels(values = list(
      variable = "Variable",
      `Ciclo 1  \n Sample (as is)` = "Crude",
      `Ciclo 2  \n Sample (as is)` = "Crude",
      `Ciclo 3  \n Sample (as is)` = "Crude",
      `Ciclo 4  \n Sample (as is)` = "Crude",
      `Ciclo estim_1 Sample (weighted)` = "Estimate (weighted)",
      `Ciclo estim_2 Sample (weighted)` = "Estimate (weighted)",
      `Ciclo estim_3 Sample (weighted)` = "Estimate (weighted)",
      `Ciclo estim_4 Sample (weighted)` = "Estimate (weighted)",
      `Ciclo CI_1 Sample (weighted)` = "95% CI",
      `Ciclo CI_2 Sample (weighted)` = "95% CI",
      `Ciclo CI_3 Sample (weighted)` = "95% CI",
      `Ciclo CI_4 Sample (weighted)` = "95% CI"
    )) %>%
    add_header_row(colwidths = c(1,3,3,3,3),
                   values = c("","Cycle 1", "Cycle 2", "Cycle 3","Cycle 4")) %>%
    save_as_docx(path = glue::glue("paper_tables/{filtertype}_Table1_{today()}.docx"),
                 pr_section =  prop_section(page_size = page_size(orient = "landscape"),
                                            type = "continuous"))
  
  
  suptable <- sample %>% 
    filter(Category == "Strata" | Category == "Total") %>%
    select(-Category)
  
  colnames(suptable) <- str_replace_all(colnames(suptable), "Ciclo","Cycle")
  colnames(suptable) <- str_replace_all(colnames(suptable), "as is","crude")
  
  suptable %>%
    mutate(Variable = case_when(
      str_detect(Variable, "ventas") ~ "Sales center",
      str_detect(Variable, "Distrib") ~ "Distribution center",
      str_detect(Variable, "Planta") ~ "Production plant",
      TRUE ~ Variable
    )) %>%
    rename(`Work center` = Variable) %>%
    select(1:9) %>%
    write_excel_csv(glue("paper_tables/{filtertype}_Suplement_Table1_{today()}.csv"))
  
  #Plot 1 for paper ----
  cycle_labs <- as_labeller(
    c(`1` = "**Cycle 1**<br>------------<br><i style='color:#7f7f7f'>_September 22-29, 2020_</i>", 
      `2` = "**Cycle 2**<br>------------<br><i style='color:#7f7f7f'>_November 9-13, 2020_</i>",
      `3` = "**Cycle 3**<br>------------<br><i style='color:#7f7f7f'>_January 4-9, 2021_</i>",
      `4` = "**Cycle 4**<br>------------<br><i style='color:#7f7f7f'>_November 22-26, 2021_</i>"))
  
  iggsexage <- read_csv(glue("results/{filtertype}_edad_igg_sex_{today()}.csv"))
  iggsexage <- iggsexage %>%
    mutate(sex = if_else(str_detect(Edad_cat,"Hombre"),"Male","Female")) %>%
    mutate(Age = case_when(
      str_detect(Edad_cat, "Menores a 35") ~ "< 35",
      str_detect(Edad_cat, "35 a 44") ~ "35 - 44",
      str_detect(Edad_cat, "45 a 59") ~ "45 - 59",
      str_detect(Edad_cat, "60 y más") ~ "≥ 60",
    ))  %>%
    mutate(Age =
             factor(Age,
                    c("< 35","35 - 44", "45 - 59", "≥ 60"))) %>%
    filter(Age != "≥ 60")
  
  ggplot(iggsexage) +
    geom_point(aes(x = sex, y = median, color = Age, shape = sex), 
               position = position_dodge(width=0.75), size = 2) +
    geom_errorbar(aes(x = sex, ymin = CI_low,
                      width = 0.25,
                      ymax = CI_high, color = Age), 
                  position = position_dodge(width=0.75)) +
    facet_wrap(~Ciclo, labeller = cycle_labs) +
    theme_minimal() +
    labs(
      x = "",
      y = "Antibody prevalence"#,
      #title = "SARS-CoV-2 antibody prevalence by sex and age category per survey cycle",
      #subtitle = "Analysis adjusted by test's specificity (0.98) and sensitivity (0.907)"
    ) +
    scale_color_manual("Age", values = c("#785EF0","#DC267F","#648FFF")) +
    theme(
      strip.text.x = element_markdown()
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_shape_manual("Sex", values = c("Male" = 16, "Female" = 15))
  ggsave(glue("plots_share/{filtertype}_edad_igg_sex_{today()}.pdf"), width = 8, height = 5)
  
  
  #TABLE 2 (IGG) ----
  igg_age      <- read_csv(glue("results/{filtertype}_edad_igg_{today()}.csv")) %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo, Edad_cat) %>%
    rename(Variable = Edad_cat)
  
  igg_age <- mydata %>%
    filter(!is.na(edad_cat)) %>%
    group_by(Ciclo, edad_cat) %>%
    tally() %>%
    rename(Variable = edad_cat) %>%
    right_join(igg_age, by = c("Ciclo", "Variable")) %>%
    rename(`Sample size` = n) %>%
    mutate(Variable = case_when(
      Variable == "Menores a 35" ~ "< 35",
      Variable == "35 a 44" ~ "35 - 44",
      Variable == "45 a 59" ~ "45 - 59",
      Variable == "60 y más" ~ "60+",
    )) %>%
    arrange(Variable) %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  igg_sex <- read_csv(glue("results/{filtertype}_IgG_conditional_sexo_{today()}.csv")) %>%
    filter(Adjusted != "Unadjusted") %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo, variable) %>%
    rename(Variable = variable) %>%
    filter(str_detect(Variable,"IgG\\|Hombre|IgG\\|Mujer")) %>%
    mutate(Variable = if_else(str_detect(Variable,"Hombre"),"Men","Women"))
  
  igg_sex <- mydata %>%
    filter(!is.na(sexo)) %>%
    group_by(Ciclo, sexo) %>%
    tally() %>%
    rename(Variable = sexo) %>%
    mutate(Variable = if_else(str_detect(Variable,"Hombre"),"Men","Women")) %>%
    right_join(igg_sex, by = c("Ciclo", "Variable")) %>%
    rename(`Sample size` = n) 
  
  igg_rdstdc   <- read_csv(glue("results/{filtertype}_IgG_conditional_itt_{today()}.csv")) %>%
    filter(Adjusted != "Unadjusted") %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo, variable) %>%
    rename(Variable = variable) %>%
    filter(str_detect(Variable,"IgG\\|ITT|IgG\\|¬ ITT")) %>%
    mutate(Variable = if_else(str_detect(Variable,"¬ ITT"),
                              "No RD-STDC in the last 6 months",
                              "RD-STDC in the last 6 months"))
  
  igg_rdstdc <- mydata %>%
    filter(!is.na(`RD-STDC`)) %>%
    group_by(Ciclo, `RD-STDC`) %>%
    tally() %>%
    rename(Variable = `RD-STDC`) %>%
    mutate(Variable = if_else(Variable == 0,
                              "No RD-STDC in the last 6 months",
                              "RD-STDC in the last 6 months")) %>%
    right_join(igg_rdstdc, by = c("Ciclo", "Variable")) %>%
    rename(`Sample size` = n) 
  
  igg_symptoms <- read_csv(glue("results/{filtertype}_IgG_autorreporte_covid_{today()}.csv")) %>%
    filter(Adjusted != "Unadjusted") %>%
    mutate(estim = scales::comma(100*median, accuracy = 0.01)) %>%
    mutate(CI    = paste0(
      scales::comma(100*CI_low, accuracy = 0.01),"-",
      scales::comma(100*CI_high, accuracy = 0.01))) %>%
    select(estim, CI, Ciclo, variable) %>%
    rename(Variable = variable) %>%
    filter(str_detect(Variable,"\\|Tuvo|\\|No")) %>%
    mutate(Variable = if_else(str_detect(Variable,"No"),"No self-reported respiratory illness",
                              "Self-reported respiratory illness"))
  
  igg_symptoms <- mydata %>%
    filter(covid_sintoma != "No sabe" & covid_sintoma != "No deseo responder" & !is.na(covid_sintoma)) %>%
    group_by(Ciclo, covid_sintoma) %>%
    tally() %>%
    rename(Variable = covid_sintoma) %>%
    mutate(Variable = if_else(Variable == "No","No self-reported respiratory illness",
                              "Self-reported respiratory illness")) %>%
    right_join(igg_symptoms, by = c("Ciclo", "Variable")) %>%
    rename(`Sample size` = n) 
  
  igg_sex %>%
    bind_rows(igg_symptoms) %>%
    bind_rows(igg_rdstdc) %>%
    bind_rows(igg_age) %>%
    mutate(`Sample size` = scales::comma(`Sample size`)) %>%
    pivot_wider(id_cols = c("Variable"),
                values_from = c("Sample size", estim, CI),
                names_from = "Ciclo", 
                values_fill = "0",
                names_glue = "Cycle {Ciclo}: {.value}") %>%
    mutate(Variable = case_when(
      Variable == "35 a 44" ~ "35 - 44",
      Variable == "45 a 59" ~ "45 - 59",
      Variable == "60+" ~ "60 over",
      Variable == "< 35" ~ "Under 35",
      Variable == "covid symptom (yes)" ~ "covid symptom",
      TRUE ~ Variable
    )) %>%
    mutate(Variable = factor(Variable,
                             levels = c("Men","Women",
                                        "Under 35","35 - 44","45 - 59",
                                        "60 over","Self-reported respiratory illness",
                                        "No self-reported respiratory illness",
                                        "RD-STDC in last 6 months",
                                        "No RD-STDC in the last 6 months" ))) %>%
    arrange(Variable) %>%
    select(Variable,matches("1"),matches("2"),matches("3"),matches("4")) %>%
    flextable() %>%
    set_header_labels(values = list(
      variable = "Variable",
      `Cycle 1: Sample size` = "Sample size",
      `Cycle 2: Sample size` = "Sample size",
      `Cycle 3: Sample size` = "Sample size",
      `Cycle 4: Sample size` = "Sample size",
      `Cycle 1: estim` = "Estimate (weighted)",
      `Cycle 2: estim` = "Estimate (weighted)",
      `Cycle 3: estim` = "Estimate (weighted)",
      `Cycle 4: estim` = "Estimate (weighted)",
      `Cycle 1: CI` = "95% CI",
      `Cycle 2: CI` = "95% CI",
      `Cycle 3: CI` = "95% CI",
      `Cycle 4: CI` = "95% CI"
    )) %>%
    add_header_row(colwidths = c(1,3,3,3,3),
                   values = c("","Cycle 1", "Cycle 2", "Cycle 3","Cycle 4")) %>%
    fontsize(size = 9, part = "all") %>%
    save_as_docx(path = glue::glue("paper_tables/{filtertype}_Table2_{today()}.docx"),
                 pr_section =  prop_section(page_size = page_size(orient = "landscape"),
                                            type = "continuous"))
  
  #TABLE 3 (Vaccines) ----
  vaccine_data <- read_csv(glue("results/{filtertype}_vacunatipo_{today()}.csv"))
  vaccine_igg  <- read_csv(glue("results/{filtertype}_vac_igg_{today()}.csv"))
  vaccine_igg  <- vaccine_igg %>%
    rename_with(~paste0("IGG_", .), c(starts_with("CI", ignore.case = F), median))
  
  vaccine_data <- vaccine_data %>%
    left_join(vaccine_igg, by = c("Vaccines","Ciclo")) %>%
    select(-CI)
  
  #Get counts
  vacunas_data <- mydata %>%
    filter(Ciclo == 4) %>%
    mutate(Vaccine = if_else(!is.na(vacuna), vacuna, `Vacuna 1ra`)) %>%
    mutate(Vaccine = str_replace_all(toupper(Vaccine)," ", "")) %>%
    filter(!is.na(Vaccine)) %>%
    group_by(Vaccine) %>%
    tally()  %>%
    mutate(Vaccine = if_else(str_detect(Vaccine,"VACUNADEMARCANOESPECIFICADA"),
                             "Unspecified", Vaccine)) %>%
    mutate(Vaccine = if_else(str_detect(Vaccine,"SPUTNIKV"),
                             "SPUTNIK V", Vaccine)) %>%
    mutate(Vaccine = str_to_sentence(Vaccine)) 
  
  vacunas_data_positivo <- mydata %>%
    filter(Ciclo == 4) %>%
    filter(IGG == "Positivo") %>%
    mutate(Vaccine = if_else(!is.na(vacuna), vacuna, `Vacuna 1ra`)) %>%
    mutate(Vaccine = str_replace_all(toupper(Vaccine)," ", "")) %>%
    filter(!is.na(Vaccine)) %>%
    group_by(Vaccine) %>%
    tally()  %>%
    mutate(Vaccine = if_else(str_detect(Vaccine,"VACUNADEMARCANOESPECIFICADA"),
                             "Unspecified", Vaccine)) %>%
    mutate(Vaccine = if_else(str_detect(Vaccine,"SPUTNIKV"),
                             "SPUTNIK V", Vaccine)) %>%
    mutate(Vaccine = str_to_sentence(Vaccine)) 
  
  vaccine_data %>%
    left_join(vacunas_data, by = c("Vaccines" = "Vaccine")) %>%
    rename(`Total workers within group (sample)` = n) %>%
    left_join(vacunas_data_positivo, by = c("Vaccines" = "Vaccine")) %>%
    rename(`Total workers with antibodies within group (sample)` = n) %>%
    mutate(Vaccines = factor(Vaccines, 
                             levels = c("Astrazeneca","Cansino","Janssen","Moderna",
                                        "Pfizer","Sinovac","Sputnik v","Unspecified", "No"))) %>%
    arrange(Vaccines) %>%
    select(-Ciclo, -IGG_CI) %>%
    mutate(CI_vac = paste0(
      scales::comma(100*CI_low, 0.01), "-",
      scales::comma(100*CI_high, 0.01)
    )) %>%
    mutate(`Weighted estimate (%)` = scales::comma(100*median, 0.01)) %>%
    mutate(CI_igg = paste0(
      scales::comma(100*IGG_CI_low, 0.01), "-",
      scales::comma(100*IGG_CI_high, 0.01)
    )) %>%
    mutate(`Within vaccinated weighted estimate (%)` = scales::comma(100*IGG_median, 0.01)) %>%
    select(Vaccines, `Total workers within group (sample)`, 
           `Weighted estimate (%)`,
           CI_vac, `Total workers with antibodies within group (sample)`,
           `Within vaccinated weighted estimate (%)`, CI_igg) %>%
    rename(`Total workers vaccinated` = `Total workers within group (sample)`) %>%
    rename(`Total workers IGG` = `Total workers with antibodies within group (sample)`) %>%
    flextable() %>%
    set_header_labels(values = list(
      Vaccines = "Vaccine",
      `Total workers vaccinated` = "Sample size (crude)",
      `Total workers IGG`        = "Sample size (crude)",
      `Weighted estimate (%)`    = "Weighted estimate (%)",
      `Within vaccinated weighted estimate (%)` = "Weighted estimate (%)",
      `CI_vac` = "95% CI",
      `CI_igg` = "95% CI"
    )) %>%
    add_header_row(colwidths = c(1,3,3),
                   values = c("","Vaccinated workers","Antibodies within vaccine group")) %>%
    footnote(value = as_paragraph(
      paste0("Estimation assuming each percent is distributed Beta with a standard uniform prior. ", 
      "For small sample sizes the prior is strong enough that the posterior estimate does not ",
      "move far enough from 0.5.")),
      ref_symbols = c("*"), i = 2, j = 1) %>%
    save_as_docx(path = glue::glue("paper_tables/{filtertype}_Table4_{today()}.docx"))
    
  
}

#Delete complied files-----
if (user["user"] == "rod" | user["user"] == "rodrigo"){
  message("Deleting compiled files")
  compiled_files <- list.files("code/additional_functions/", "*threads", full.names = T)
  file.remove(compiled_files)
}

