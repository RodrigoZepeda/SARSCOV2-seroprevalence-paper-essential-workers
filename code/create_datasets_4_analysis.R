rm(list = ls())
pacman::p_load(tidyverse, readxl, tidybayes, rlang, scales, purrr, openssl, 
               lubridate, janitor)
setwd("~/Dropbox/EPCOVID/paper_SPM")

#DB CREATION-------------
ciclo1   <- read_delim("datasets/EPCOVID_c1_20210622.txt")
ciclo2   <- read_delim("datasets/EPCOVID_c2_20210622.txt")
ciclo3   <- read_delim("datasets/EPCOVID_c3_20210622_n.txt")
ciclo4   <- read_delim("datasets/EP_Bimbo_220101.txt")

rps <- c(
  unique(ciclo1$registro_patronal_bimbo),
  unique(ciclo2$registro_patronal_bimbo),
  unique(ciclo3$registro_patronal_bimbo),
  unique(ciclo4$cve_reg_patronal)
) |> unique() |>
  write_lines("rps_epcovid.txt")

#Create dataset for TWD
enf    <- read_delim("datasets/BASE_ENF_RESP_16042022(1).zip", delim = "|")


#RPS
rps <- c(unique(ciclo1$registro_patronal_bimbo),
         unique(ciclo2$registro_patronal_bimbo),
         unique(ciclo3$registro_patronal_bimbo),
         unique(ciclo4$cve_reg_patronal))
rps <- rps[!is.na(rps)]
rps <- unique(str_remove_all(rps,"-"))

vaccines <- read_excel("datasets/CURP_Encuesta_Panel_3_mayo_2022-vacuna.xlsx")
vaccines <- vaccines %>%
  select(CURPIMSS, vacuna, EsquemaCompleto) %>%
  mutate(across(.cols = c(vacuna, EsquemaCompleto), 
                ~ if_else(. == "NULL", NA_character_, .)))

#Curps para MHA
#ciclo4 %>%
#  select(curp, cve_curp) %>%
#  mutate(curp = if_else(nchar(curp) == 18, curp, NA_character_)) %>%
#  mutate(cve_curp = if_else(nchar(cve_curp) == 18, cve_curp, NA_character_)) %>%
#  mutate(curp = if_else(is.na(curp), cve_curp, curp)) %>%
#  mutate(cve_curp = if_else(is.na(cve_curp), curp, cve_curp)) %>%
#  pivot_longer(curp:cve_curp) %>%
#  distinct(value) %>%
#  rename(CURP = value) %>%
#  write_excel_csv("CURP_Encuesta_Panel_3_mayo_2022.csv")

ciclo1 <- ciclo1 %>% 
  select(num_afiliacion, nombre_completo, sexo, edad, clasificacion, covid_sintoma,
         resultado_sarscov2, itt_uno, resultado_cualit_sero, folio_laboratorio_sero) %>%
  filter(!is.na(clasificacion)) %>%
  mutate(Ciclo = 1)

ciclo2 <- ciclo2 %>% 
  select(num_afiliacion, nombre_completo, sexo, edad, clasificacion, covid_sintoma,
         resultado_sarscov2, itt_uno, resultado_cualit_sero, folio_laboratorio_sero) %>%
  filter(!is.na(clasificacion))  %>%
  mutate(Ciclo = 2)

ciclo3 <- ciclo3 %>% 
  select(num_afiliacion, nombre_completo, sexo, edad, clasificacion, covid_sintoma,
         resultado_sarscov2, itt_uno, resultado_cualit_sero, folio_laboratorio_sero) %>%
  filter(!is.na(clasificacion))  %>%
  mutate(Ciclo = 3)

ciclo4 <- ciclo4 %>% 
  mutate(Nombre = if_else(is.na(ap_paterno) | is.na(ap_materno) | 
                            str_detect(nombre, ap_paterno), nombre, 
                          paste(ap_paterno, ap_materno, nombre))) %>%
  left_join(vaccines, by = c("curp" = "CURPIMSS")) %>%
  select(num_afiliacion, Nombre, sexo, edad, Clasificación, covid_sintoma,
         resultadopcr, itt_uno, resultadoanticuerpo, foliomuestrainmegen,
         vacuna, EsquemaCompleto) %>%
  #filter(!is.na(Clasificación))  %>%
  mutate(edad = as.numeric(edad)) %>%
  filter(!is.na(Nombre)) %>%
  mutate(Ciclo = 4)

#Bind cycle 4 other workplace names
ciclo4_bis <- read_excel("datasets/BD 4to ciclo_Nov21.xlsx", sheet = 2) %>%
  mutate(Ciclo = 4) %>%
  select(Nombre, Clasificación, Ciclo) %>%
  rename(Class_4 = Clasificación)

ciclo4 <- ciclo4 %>% 
  left_join(ciclo4_bis, by = c("Nombre", "Ciclo")) %>%
  mutate(Clasificación = if_else(is.na(Clasificación) & Ciclo == 4, Class_4, Clasificación)) %>%
  filter(!is.na(Clasificación)) %>%
  select(-Class_4) %>%
  distinct_at(vars(Nombre), .keep_all = T) %>%
  filter(!is.na(resultadoanticuerpo)) %>%
  rename(nombre_completo = Nombre)

ciclo4 <- ciclo4 %>% 
  rename(clasificacion = Clasificación) %>%
  rename(resultado_sarscov2 = resultadopcr) %>%
  rename(resultado_cualit_sero = resultadoanticuerpo) %>%
  rename(folio_laboratorio_sero = foliomuestrainmegen) 

panel_datasets <- ciclo1 %>%
  bind_rows(ciclo2) %>%
  bind_rows(ciclo3) %>%
  bind_rows(ciclo4)

panel_datasets <- panel_datasets %>%
  filter(!is.na(nombre_completo)) %>%
  mutate(identifier = nombre_completo) 

dats <- data.frame(identifier = unique(panel_datasets$identifier)) %>%
  mutate(ID = 1:n())

panel_datasets <- panel_datasets %>% 
  left_join(dats, by = "identifier")

panel_datasets <- panel_datasets %>%
  mutate(folio_laboratorio_sero = str_replace_all(folio_laboratorio_sero, "N", "")) %>%
  mutate(folio_laboratorio_sero = str_replace_all(folio_laboratorio_sero, "S", "")) %>%
  mutate(folio_laboratorio_sero = str_replace_all(folio_laboratorio_sero, "M", "")) %>%
  mutate(nombre_completo = str_replace_all(nombre_completo, "  ", " "))

#RANDOM vs CONVENIENCE------------------------------------------------------------------------------------
ciclo1 <- read_excel("datasets/Informacion cierre ciclo 1.xlsx") %>%
  mutate(nombre = str_replace(nombre, ",", "")) %>%
  mutate(Ciclo = 1) %>%
  mutate(Seleccion1 = if_else(Folio == "Sin folio", "CONVENIENCIA", "ALEATORIA")) %>%
  select(nombre, Ciclo, Seleccion1, Folio)

panel_datasets <- panel_datasets %>%
  left_join(ciclo1, by = c("nombre_completo" = "nombre", "Ciclo" = "Ciclo"))
#write_excel_csv(david, "david.csv")

ciclo2 <- read_excel("datasets/lista participantes bimbo ciclo 2.xlsx") %>%
  mutate(Ciclo = 2) %>%
  mutate(Seleccion = "ALEATORIA") %>%
  mutate(Nombre = str_replace(Nombre, ",", "")) %>%
  mutate(`Num. Folio` = str_replace_all(`Num. Folio`, "N", "")) %>%
  mutate(`Num. Folio` = str_replace_all(`Num. Folio`, "S", "")) %>%
  mutate(`Num. Folio` = str_replace_all(`Num. Folio`, " ", "")) %>%
  select(Nombre, Ciclo, Seleccion, "Num. Folio") %>%
  rename("Folio" = "Num. Folio") 

ciclo22 <- read_excel("datasets/match_si_c2.xlsx") %>%
  select(idn_match, nombre_completo) %>%
  rename("Folio" = "idn_match") %>%
  rename("Nombre" = "nombre_completo") %>%
  mutate(Ciclo = 2) %>%
  mutate(Seleccion = "ALEATORIA") %>%
  mutate(Folio = as.character(Folio))

ciclo2 <- bind_rows(ciclo2, ciclo22) %>%
  distinct() %>%
  rename(Seleccion2 = Seleccion)

panel_datasets <- panel_datasets %>%
  left_join(ciclo2, by = c("folio_laboratorio_sero" = "Folio", "Ciclo"))

ciclo3 <- read_excel("datasets/Bimbo ciclo 3 participantes ENERO 28.xlsx", 
                     sheet = "TOMA DE MUESTRA REAL COMPLETO") %>%
  select(NOMBRE, Selecccion) %>%
  mutate(NOMBRE = str_replace(NOMBRE, ",", "")) %>%
  rename(Seleccion3 = Selecccion) %>%
  mutate(Ciclo = 3)

panel_datasets <- panel_datasets %>%
  left_join(ciclo3, by = c("nombre_completo" = "NOMBRE", "Ciclo"))

ciclo4 <- read_excel("datasets/Concentrado General de Ciclo 4 Bimbo.xlsx",
                     sheet = "Hoja1") %>%
  select(Clasificación, NOMBRE, Selecccion, `Vacuna 1ra`,`Vacuna 2da`) %>%
  rename(category = Clasificación) %>%
  mutate(Ciclo = 4) %>%
  rename(Seleccion4 = Selecccion) %>%
  mutate(Seleccion4 = if_else(
    #FIXME ARREGLAR LOS ALEATORTIOS
    str_detect(category,"Distribuc") & Seleccion4 == "tbd", "ALEATORIO", Seleccion4))

panel_datasets <- panel_datasets %>%
  left_join(ciclo4, by = c("nombre_completo" = "NOMBRE", "Ciclo"))

panel_datasets <- panel_datasets %>%
  mutate(Seleccion = case_when(
    !is.na(Seleccion1) ~ Seleccion1,
    !is.na(Seleccion2) ~ Seleccion2,
    !is.na(Seleccion3) ~ Seleccion3,
    !is.na(Seleccion4) ~ Seleccion4,
    TRUE ~ NA_character_
  )) %>%
  mutate(Seleccion = toupper(Seleccion)) %>%
  mutate(Seleccion = case_when(
    str_detect(Seleccion,"ALEATOR") ~ "ALEATORIO",
    TRUE ~ "CONVENIENCIA"
  ))

panel_datasets <- panel_datasets  %>%
  select(ID, sexo, edad, itt_uno, clasificacion, resultado_sarscov2, resultado_cualit_sero,
         Seleccion, Ciclo, `Vacuna 1ra`, `Vacuna 2da`, covid_sintoma, num_afiliacion,
         vacuna, EsquemaCompleto) %>%
  mutate(edad_cat = case_when(
    edad < 35 ~ "Menores a 35",
    edad >= 35 & edad < 45 ~ "35 a 44",
    edad >= 45 & edad < 60 ~ "45 a 59",
    edad >= 60 ~ "60 y más",
  )) %>% 
  mutate(Seleccion = if_else(is.na(Seleccion),"CONVENIENCIA", Seleccion)) %>%
  mutate(resultado_sarscov2 = if_else(is.na(resultado_sarscov2), "SIN PRUEBA", resultado_sarscov2)) %>%
  rename(PCR = resultado_sarscov2) %>%
  mutate(resultado_cualit_sero = if_else(is.na(resultado_cualit_sero), "SIN PRUEBA", resultado_cualit_sero)) %>%
  rename(IGG = resultado_cualit_sero) %>%
  rename(`RD-STDC` = itt_uno) %>%
  mutate(sexo = case_when(
    str_detect(sexo, "Hombre|Masculino|MASCULINO") ~ "Hombre",
    str_detect(sexo, "Mujer|Femenino|FEMENINO") ~ "Mujer",
    TRUE ~ NA_character_
  )) %>%
  mutate(
    strata = case_when(
      str_detect(clasificacion, "ventas|VENTAS|Ventas") ~ "Centro de ventas",
      str_detect(clasificacion, "producción|Producción|PRODUCCIÓN") ~ "Planta de producción",
      str_detect(clasificacion, "Distribución|Distribucion|DISTRIBUCIÓN") ~ "Centro de Distribución",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    IGG = case_when(
      str_detect(IGG, "Positiv") ~ "Positivo",
      str_detect(IGG, "Negativ") ~ "Negativo",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-clasificacion) 

#RD-STDC in base
panel_datasets <- panel_datasets %>% select(-`RD-STDC`)
rdstd_search   <- panel_datasets %>% 
  distinct(num_afiliacion, Ciclo) %>%
  mutate(fecha = case_when(
    Ciclo == 1 ~ ymd("2020/09/22"),
    Ciclo == 2 ~ ymd("2020/11/09"),
    Ciclo == 3 ~ ymd("2021/01/04"),
    Ciclo == 4 ~ ymd("2021/11/22")
  ))

rdstd_search <- rdstd_search %>% 
  left_join(enf, by = c("num_afiliacion" = "NUM_AFILIACION")) %>%
  filter(dmy(COLAPS_FEC_INICIO) >= fecha - months(6) & dmy(COLAPS_FEC_INICIO) < fecha) %>%
  select(num_afiliacion, Ciclo) %>%
  mutate(`RD-STDC` = 1) %>%
  distinct(num_afiliacion, Ciclo, .keep_all = T)

panel_datasets <- panel_datasets %>%
  left_join(rdstd_search, by = c("num_afiliacion", "Ciclo")) %>%
  mutate(`RD-STDC` = if_else(is.na(`RD-STDC`), 0, 1))

#WEIGHTS---------------------------------------------------------
Ciclo     <- 1:4
strata    <- c("Centro de ventas", "Planta de producción", "Centro de Distribución")
weights   <- expand_grid(Ciclo, strata) 
weights   <- weights %>%
  mutate(N_Strata = c(7017,3171,231,6312,2874,172,7414,3290,262,7760,3525,317))

weights   <- weights %>%
  group_by(Ciclo) %>%
  summarise(N_Total = sum(N_Strata)) %>%
  right_join(weights, by = "Ciclo") %>%
  mutate(pondef = N_Strata/N_Total)

#Create sample size
ssize <- panel_datasets %>%
  group_by(Ciclo, strata) %>%
  summarise(n_Strata = n(), 
            igg_positive_Strata = sum(IGG == "Positivo"),
            .groups = "keep")

ssize <- ssize %>% 
  group_by(Ciclo) %>%
  summarise(n_total = sum(n_Strata)) %>%
  right_join(ssize, by = "Ciclo")

weights <- weights %>%
  left_join(ssize, by = c("Ciclo", "strata"))

panel_datasets <- panel_datasets %>%
  left_join(weights, by = c("strata","Ciclo")) %>%
  mutate(`Vacuna 1ra` = if_else(str_detect(`Vacuna 1ra`,"na"), NA_character_, `Vacuna 1ra`)) %>%
  mutate(`Vacuna 2da` = if_else(str_detect(`Vacuna 2da`,"na"), NA_character_, `Vacuna 2da`)) %>%
  mutate(`Vacuna 1ra` = if_else(str_detect(`Vacuna 1ra`,"\\bSI\\b"), "VACUNA DE MARCA NO ESPECIFICADA", `Vacuna 1ra`)) %>%
  mutate(`Vacuna 2da` = if_else(str_detect(`Vacuna 2da`,"\\bSI\\b"), "VACUNA DE MARCA NO ESPECIFICADA", `Vacuna 2da`)) 

panel_datasets <- panel_datasets %>%
  mutate(ID = md5(paste0(as.character(ID), sexo), key = "IMSSPIRACION"))

itt_matching   <- panel_datasets %>%
  select(`RD-STDC`, ID, Ciclo, Seleccion, strata, num_afiliacion, sexo, pondef,
         N_Total, n_total, N_Strata, n_Strata) %>%
  ungroup()

panel_datasets <- panel_datasets %>%
  select(-num_afiliacion)
  
# panel_datasets %>% 
#   filter(Ciclo %in% c(1,4)) %>%
#   write_rds("final_datasets/annonymized_dataset_paper.rds") %>%
#   write_excel_csv("final_datasets/annonymized_dataset_paper.csv")
# 
# panel_datasets %>% 
#   write_rds("final_datasets/annonymized_dataset_paper_4_cycles.rds") %>%
#   write_excel_csv("final_datasets/annonymized_dataset_paper_4_cycles.csv")


bimbo  <- enf %>%
  filter(CVE_REG_PATRONAL %in% rps) %>%
  group_by(COLAPS_FEC_INICIO) %>%
  tally() %>%
  mutate(n_per_buisness_capita = n/13397) %>%
  write_excel_csv("datasets/empresa_by_itself.csv")

enf    <- enf %>% 
  select(COLAPS_FEC_TERMINO, COLAPS_FEC_INICIO, NUM_AFILIACION, 
         ULT_COD_CIE, CVE_UNIDAD_EXPEDICION) %>%
  mutate(across(matches("_FEC_"), ~ dmy(.)))

itt_matching <- itt_matching %>%
  left_join(enf, by = c("num_afiliacion" = "NUM_AFILIACION"))

umfs   <- unique(unlist(itt_matching$CVE_UNIDAD_EXPEDICION))
umfs   <- umfs[!is.na(umfs)]

enf    <- enf %>% 
  filter(CVE_UNIDAD_EXPEDICION %in% umfs) %>%
  group_by(COLAPS_FEC_INICIO) %>% 
  tally()

adsc   <- read_excel("datasets/2021-12-30_SIAIS_PAMF_Junio_2021_ver1.xlsx",
                     sheet = 2, skip = 11) %>% clean_names()
adsc   <- adsc %>% filter(consultorio_9999_total_unidad != "TOTAL DELEGACIONAL")
adsc   <- adsc %>% select(cve_presupuestal, asegurados_derechohabientes_dir)

#Cargar para leer la clabe priei
unidades <- read_excel("datasets/CUUMSP_ENERO_2022.xlsx",
                       sheet = 2, skip = 9) %>% clean_names()
adsc     <- adsc %>% 
  left_join(unidades, by = c("cve_presupuestal" = "clave_presupuestal")) %>%
  select(unidad_de_informacion_prei, asegurados_derechohabientes_dir)

asegs <- tibble(unidad = umfs) %>%
  mutate(unidad = str_sub(unidad, 1, 6)) %>%
  left_join(adsc, by = c("unidad" = "unidad_de_informacion_prei")) %>%
  summarise(asegurados = sum(asegurados_derechohabientes_dir, na.rm = T)) %>%
  unlist() %>%
  as.numeric()

enf    <- enf %>% 
  mutate(n_per_capita = n/!!asegs)

enf    <- enf %>%
  write_excel_csv("final_datasets/casos_per_capita_en_unidades.csv") %>%
  write_rds("final_datasets/casos_per_capita_en_unidades.rds")

itt_matching <- itt_matching %>%
  mutate(fecha_ciclo = case_when(
    Ciclo == 1 ~ ymd("2020/09/22"),
    Ciclo == 2 ~ ymd("2020/11/09"),
    Ciclo == 3 ~ ymd("2021/01/04"),
    Ciclo == 4 ~ ymd("2021/11/23"),
  )) %>%
  mutate(itt_ciclo = as.numeric(days(fecha_ciclo - COLAPS_FEC_INICIO)) < 180)

itt_matching <- itt_matching %>%
  select(-CVE_UNIDAD_EXPEDICION) %>%
  group_by(ID) %>%
  mutate(idnum = 1:n()) %>%
  ungroup() %>%
  arrange(ID, Ciclo, desc(itt_ciclo))  %>%
  group_by(ID, Ciclo) %>%
  arrange(Ciclo, desc(itt_ciclo)) %>%
  mutate(duplicate = 1:n()) %>%
  filter(duplicate == 1) %>%
  ungroup() 

panel_datasets_2 <- itt_matching %>%
  right_join(panel_datasets)

panel_datasets_2 <- panel_datasets_2 %>%
  mutate(itt_ciclo = if_else(is.na(itt_ciclo), FALSE, itt_ciclo)) %>%
  select(-duplicate, -idnum, -num_afiliacion) %>%
  #select(-`RD-STDC`) %>%
  #rename(`RD-STDC` = itt_ciclo) %>%
  #mutate(`RD-STDC` = as.numeric(`RD-STDC`)) %>% 
   write_rds("final_datasets/annonymized_dataset_paper_4_cycles.rds") %>%
   write_excel_csv("final_datasets/annonymized_dataset_paper_4_cycles.csv")

itt_matching <- itt_matching %>%
  mutate(ID = md5(paste0(as.character(ID), sexo), key = "IMSSPIRACION")) %>%
  select(-num_afiliacion) %>%
  #filter(Ciclo %in% c(1,4)) %>%
  write_excel_csv("final_datasets/annonimized_twd_dataset.csv") %>%
  write_rds("final_datasets/annonimized_twd_dataset.rds")
