install.packages("httpgd")
install.packages("tidyr")
install.packages("gridExtra")
install.packages("forcats")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("minpack.lm")
install.packages("rsq") 
while(TRUE){
  
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(gridExtra)
  library(forcats)
  library(httpgd)
  library(minpack.lm)
  library(rsq) 
  break
}
# ordenando datos 
while(TRUE){


# todos VC
(hplc_vc <- tibble(read.csv("RAW_HPLC_VC.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal") %>%
    pivot_longer(cols = 2:18, names_to = "patient", values_to = "concentration", values_drop_na = FALSE, names_repair = "minimal")  %>%
        mutate(sol = c("")))
  pos0 <- 1
  for (i in hplc_vc$patient){
    if(i == "Paciente 1") hplc_vc$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 2") hplc_vc$sol[pos0]   <- "P" 
    if(i == "Paciente 3") hplc_vc$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 4") hplc_vc$sol[pos0]   <- "P" 
    if(i == "Paciente 5") hplc_vc$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 6") hplc_vc$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 7") hplc_vc$sol[pos0]   <- "P" 
    if(i == "Paciente 8") hplc_vc$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 9") hplc_vc$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 10") hplc_vc$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 11") hplc_vc$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 12") hplc_vc$sol[pos0]  <- "P"
    if(i == "Paciente 13") hplc_vc$sol[pos0]  <- "P"
    if(i == "Paciente 14") hplc_vc$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 15") hplc_vc$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 17") hplc_vc$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 18") hplc_vc$sol[pos0]  <- "TAC 2"
    pos0 <- pos0 + 1
  }
hplc_vc

# promedios VC hplc
(hplc_vc_mean <- tibble(read.csv("RAW_HPLC_VC_Mean.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal") %>%
    pivot_longer(cols = 2:18, names_to = "patient", values_to = "concentration", values_drop_na = FALSE, names_repair = "minimal")  %>%
        mutate(sol = c("")) 
    )

  pos0 <- 1
  for (i in hplc_vc_mean$patient){
    if(i == "P1") hplc_vc_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P2") hplc_vc_mean$sol[pos0]   <- "P" 
    if(i == "P3") hplc_vc_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P4") hplc_vc_mean$sol[pos0]   <- "P" 
    if(i == "P5") hplc_vc_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P6") hplc_vc_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P7") hplc_vc_mean$sol[pos0]   <- "P" 
    if(i == "P8") hplc_vc_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P9") hplc_vc_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P10") hplc_vc_mean$sol[pos0]  <- "TAC 2"
    if(i == "P11") hplc_vc_mean$sol[pos0]  <- "TAC 2"
    if(i == "P12") hplc_vc_mean$sol[pos0]  <- "P"
    if(i == "P13") hplc_vc_mean$sol[pos0]  <- "P"
    if(i == "P14") hplc_vc_mean$sol[pos0]  <- "TAC 2"
    if(i == "P15") hplc_vc_mean$sol[pos0]  <- "TAC 2"
    if(i == "P17") hplc_vc_mean$sol[pos0]  <- "TAC 2"
    if(i == "P18") hplc_vc_mean$sol[pos0]  <- "TAC 2"
    pos0 <- pos0 + 1
  }
hplc_vc_mean

# NAC individual
(hplc_nac <- tibble(read.csv("RAW_HPLC_NAC.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal") %>%
    pivot_longer(cols = 2:13, names_to = "patient", values_to = "concentration", values_drop_na = FALSE, names_repair = "minimal")  %>%
        mutate(sol = c("")))

  pos0 <- 1
  for (i in hplc_nac$patient){
    if(i == "Paciente 1") hplc_nac$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 2") hplc_nac$sol[pos0]   <- "P" 
    if(i == "Paciente 3") hplc_nac$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 4") hplc_nac$sol[pos0]   <- "P" 
    if(i == "Paciente 5") hplc_nac$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 6") hplc_nac$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 7") hplc_nac$sol[pos0]   <- "P" 
    if(i == "Paciente 8") hplc_nac$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 9") hplc_nac$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 10") hplc_nac$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 11") hplc_nac$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 12") hplc_nac$sol[pos0]  <- "P"
    if(i == "Paciente 13") hplc_nac$sol[pos0]  <- "P"
    if(i == "Paciente 14") hplc_nac$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 15") hplc_nac$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 17") hplc_nac$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 18") hplc_nac$sol[pos0]  <- "TAC 2"
    pos0 <- pos0 + 1
  }
hplc_nac

# promedios NAC hplc
(hplc_nac_mean <- tibble(read.csv("RAW_HPLC_NAC_Mean.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal") %>%
    pivot_longer(cols = 2:13, names_to = "patient", values_to = "concentration", values_drop_na = FALSE, names_repair = "minimal")  %>%
        mutate(sol = c("")) 
    )

  pos0 <- 1
  for (i in hplc_nac_mean$patient){
    if(i == "P1") hplc_nac_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P2") hplc_nac_mean$sol[pos0]   <- "P" 
    if(i == "P3") hplc_nac_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P4") hplc_nac_mean$sol[pos0]   <- "P" 
    if(i == "P5") hplc_nac_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P6") hplc_nac_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P7") hplc_nac_mean$sol[pos0]   <- "P" 
    if(i == "P8") hplc_nac_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P9") hplc_nac_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P10") hplc_nac_mean$sol[pos0]  <- "TAC 2"
    if(i == "P11") hplc_nac_mean$sol[pos0]  <- "TAC 2"
    if(i == "P12") hplc_nac_mean$sol[pos0]  <- "P"
    if(i == "P13") hplc_nac_mean$sol[pos0]  <- "P"
    if(i == "P14") hplc_nac_mean$sol[pos0]  <- "TAC 2"
    if(i == "P15") hplc_nac_mean$sol[pos0]  <- "TAC 2"
    if(i == "P17") hplc_nac_mean$sol[pos0]  <- "TAC 2"
    if(i == "P18") hplc_nac_mean$sol[pos0]  <- "TAC 2"
    pos0 <- pos0 + 1
  }
hplc_nac_mean


# DFO

(hplc_dfo <- tibble(read.csv("RAW_HPLC_DFO.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal") %>%
    pivot_longer(cols = 2:13, names_to = "patient", values_to = "concentration", values_drop_na = FALSE, names_repair = "minimal")  %>%
        mutate(sol = c("")))
  pos0 <- 1
  for (i in hplc_dfo$patient){
    if(i == "Paciente 1") hplc_dfo$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 2") hplc_dfo$sol[pos0]   <- "P" 
    if(i == "Paciente 3") hplc_dfo$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 4") hplc_dfo$sol[pos0]   <- "P" 
    if(i == "Paciente 5") hplc_dfo$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 6") hplc_dfo$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 7") hplc_dfo$sol[pos0]   <- "P" 
    if(i == "Paciente 8") hplc_dfo$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 9") hplc_dfo$sol[pos0]   <- "TAC 1" 
    if(i == "Paciente 10") hplc_dfo$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 11") hplc_dfo$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 12") hplc_dfo$sol[pos0]  <- "P"
    if(i == "Paciente 13") hplc_dfo$sol[pos0]  <- "P"
    if(i == "Paciente 14") hplc_dfo$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 15") hplc_dfo$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 17") hplc_dfo$sol[pos0]  <- "TAC 2"
    if(i == "Paciente 18") hplc_dfo$sol[pos0]  <- "TAC 2"
    pos0 <- pos0 + 1
  }
hplc_dfo

# promedios dfo hplc
hplc_dfo_mean <- tibble(read.csv("RAW_HPLC_DFO_Mean.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal")
(hplc_dfo_mean <- tibble(read.csv("RAW_HPLC_DFO_Mean.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal") %>%
    pivot_longer(cols = 2:13, names_to = "patient", values_to = "concentration", values_drop_na = FALSE, names_repair = "minimal")  %>%
        mutate(sol = c("")) 
    )

  pos0 <- 1
  for (i in hplc_dfo_mean$patient){
    if(i == "P1") hplc_dfo_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P2") hplc_dfo_mean$sol[pos0]   <- "P" 
    if(i == "P3") hplc_dfo_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P4") hplc_dfo_mean$sol[pos0]   <- "P" 
    if(i == "P5") hplc_dfo_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P6") hplc_dfo_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P7") hplc_dfo_mean$sol[pos0]   <- "P" 
    if(i == "P8") hplc_dfo_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P9") hplc_dfo_mean$sol[pos0]   <- "TAC 1" 
    if(i == "P10") hplc_dfo_mean$sol[pos0]  <- "TAC 2"
    if(i == "P11") hplc_dfo_mean$sol[pos0]  <- "TAC 2"
    if(i == "P12") hplc_dfo_mean$sol[pos0]  <- "P"
    if(i == "P13") hplc_dfo_mean$sol[pos0]  <- "P"
    if(i == "P14") hplc_dfo_mean$sol[pos0]  <- "TAC 2"
    if(i == "P15") hplc_dfo_mean$sol[pos0]  <- "TAC 2"
    if(i == "P17") hplc_dfo_mean$sol[pos0]  <- "TAC 2"
    if(i == "P18") hplc_dfo_mean$sol[pos0]  <- "TAC 2"
    pos0 <- pos0 + 1
  }
hplc_dfo_mean

break
}
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

summary(filter(hplc_vc_mean, t_min == 30 | t_min == 60, sol == "TAC 1"))
filter(hplc_vc_mean, t_min == 0, sol == "TAC 2") %>%
    summary()
filter(hplc_vc_mean, t_min == 0, sol == "TAC 2") %>%
    pull(concentration) %>%
        sd()
filter(hplc_vc_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    summary()
filter(hplc_vc_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    pull(concentration) %>%
        sd()
densidad <- filter(hplc_vc_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    pull(concentration) %>%
        density()
plot(densidad)

filter(hplc_nac_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    summary()
filter(hplc_nac_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    pull(concentration) %>%
        sd()
densidad <- filter(hplc_nac_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    pull(concentration) %>%
        density()
plot(densidad)

filter(hplc_dfo_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    summary()
filter(hplc_dfo_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    pull(concentration) %>%
        sd()
densidad <- filter(hplc_dfo_mean, t_min == 30 | t_min == 60, sol == "TAC 1") %>%
    pull(concentration) %>%
        density()
plot(densidad)


filter(hplc_vc_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    summary()
filter(hplc_vc_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    pull(concentration) %>%
        sd()
densidad <- filter(hplc_vc_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    pull(concentration) %>%
        density()
plot(densidad)

filter(hplc_nac_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    summary()
filter(hplc_nac_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    pull(concentration) %>%
        sd()
densidad <- filter(hplc_nac_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    pull(concentration) %>%
        density()
plot(densidad)

filter(hplc_dfo_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    summary()
filter(hplc_dfo_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    pull(concentration) %>%
        sd()
densidad <- filter(hplc_dfo_mean, t_min == 30 | t_min == 60, sol == "TAC 2") %>%
    pull(concentration) %>%
        density()
plot(densidad)



####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

# análisis VC
# tabla parámetros por paciente VC

# mg/ml * ml/min * 60min/h / 176.12g/mol  / 1000
# R0 para TAC 1
# 2475/250*6*60/176.12
# 20.2362
# R0 para TAC 2
# 2250/250*3*60/176.12
# 9.198274

(parametros_pk_VC <- tibble(
  paciente = c("P1", "P3", "P5", "P6", "P8", "P9", "P10", "P11", "P14", "P15", "P17", "P18"),
  V = 0, Cl= 0, Cph = 0
  )
)
# modelo PK 

modelo_ajuste_vc_tac1 <- function(V, Cl, Cph) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- c(30, 90)
  
  R <- 20.2362
  conc <- t
  
  conc[t <= t_cambio[1]] <-Cph + (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60)) 
  conc[t <= t_cambio[2] & t > t_cambio[1]] <- Cph +
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60))
  
  conc[t > t_cambio[2]] <- Cph +
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * exp(-k*(t[t > t_cambio[2]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*((t_cambio[2] - t_cambio[1])/60))) * exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60))
  
  return(conc)
}

modelo_simulacion_vc_tac1 <- function(V, Cl, Cph) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- c(30, 90)
  R <- 20.2362
  conc <- t

  conc[t <= t_cambio[1]] <- Cph + (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60))
  conc[t <= t_cambio[2] & t > t_cambio[1]] <- Cph + 
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * 
                            exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60))
  
  conc[t > t_cambio[2]] <- Cph + 
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * 
                            exp(-k*(t[t > t_cambio[2]] - t_cambio[1])/60) + # OJO, antes estaba restando el t_cambio[2] en vez de el t_cambio[1]

                  ((R/6)/Cl) * (1-exp(-k*((t_cambio[2] - t_cambio[1])/60))) * 
                            exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60) +

                  (0/Cl) * (1-exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60))
  
  return(conc)
}





  # TAC 1
while(TRUE){  
  concentracion_p1 <- filter(hplc_vc_mean, patient == "P1") %>% pull(concentration)/1000
  concentracion_p3 <- filter(hplc_vc_mean, patient == "P3") %>% pull(concentration)/1000
  concentracion_p5 <- filter(hplc_vc_mean, patient == "P5") %>% pull(concentration)/1000
  concentracion_p6 <- filter(hplc_vc_mean, patient == "P6") %>% pull(concentration)/1000
  concentracion_p8 <- filter(hplc_vc_mean, patient == "P8") %>% pull(concentration)/1000
  concentracion_p9 <- filter(hplc_vc_mean, patient == "P9") %>% pull(concentration)/1000
  
  fit_p1 <- nlsLM(concentracion_p1 ~ modelo_ajuste_vc_tac1(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p3 <- nlsLM(concentracion_p3 ~ modelo_ajuste_vc_tac1(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p5 <- nlsLM(concentracion_p5 ~ modelo_ajuste_vc_tac1(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p6 <- nlsLM(concentracion_p6 ~ modelo_ajuste_vc_tac1(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p8 <- nlsLM(concentracion_p8 ~ modelo_ajuste_vc_tac1(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p9 <- nlsLM(concentracion_p9 ~ modelo_ajuste_vc_tac1(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
parametros_pk_VC$V [1] <- summary(fit_p1)$coefficients[1,1] # V
parametros_pk_VC$V [2] <- summary(fit_p3)$coefficients[1,1] # V
parametros_pk_VC$V [3] <- summary(fit_p5)$coefficients[1,1] # V
parametros_pk_VC$V [4] <- summary(fit_p6)$coefficients[1,1] # V
parametros_pk_VC$V [5] <- summary(fit_p8)$coefficients[1,1] # V
parametros_pk_VC$V [6] <- summary(fit_p9)$coefficients[1,1] # V



parametros_pk_VC$Cl[1] <- summary(fit_p1)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[2] <- summary(fit_p3)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[3] <- summary(fit_p5)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[4] <- summary(fit_p6)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[5] <- summary(fit_p8)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[6] <- summary(fit_p9)$coefficients[2,1] # Cl

parametros_pk_VC$Cph[1] <- summary(fit_p1)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[2] <- summary(fit_p3)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[3] <- summary(fit_p5)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[4] <- summary(fit_p6)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[5] <- summary(fit_p8)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[6] <- summary(fit_p9)$coefficients[3,1] # Cph

simulacion_tac1_p1 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac1(V = parametros_pk_VC$V[1], Cl = parametros_pk_VC$Cl[1], Cph = parametros_pk_VC$Cph[1])) %>%
    mutate(paciente = "P1")
simulacion_tac1_p3 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac1(V = parametros_pk_VC$V[2], Cl = parametros_pk_VC$Cl[2], Cph = parametros_pk_VC$Cph[2])) %>%
    mutate(paciente = "P3")
simulacion_tac1_p5 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac1(V = parametros_pk_VC$V[3], Cl = parametros_pk_VC$Cl[3], Cph = parametros_pk_VC$Cph[3])) %>%
    mutate(paciente = "P5")
simulacion_tac1_p6 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac1(V = parametros_pk_VC$V[4], Cl = parametros_pk_VC$Cl[4], Cph = parametros_pk_VC$Cph[4])) %>%
    mutate(paciente = "P6")
simulacion_tac1_p8 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac1(V = parametros_pk_VC$V[5], Cl = parametros_pk_VC$Cl[5], Cph = parametros_pk_VC$Cph[5])) %>%
    mutate(paciente = "P8")
simulacion_tac1_p9 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac1(V = parametros_pk_VC$V[6], Cl = parametros_pk_VC$Cl[6], Cph = parametros_pk_VC$Cph[6])) %>%
    mutate(paciente = "P9")

simulacion_tac1 <- bind_rows(
  simulacion_tac1_p1,
  simulacion_tac1_p3,
  simulacion_tac1_p5,
  simulacion_tac1_p6,
  simulacion_tac1_p8,
  simulacion_tac1_p9
  ) %>%
      mutate(sol = "TAC 1")
  break
}
plot(parametros_pk_VC$V, parametros_pk_VC$Cl)
plot(parametros_pk_NAC$V, parametros_pk_NAC$Cl)
plot(parametros_pk_DFO$V, parametros_pk_DFO$Cl)

deviance(fit_p1)
deviance(fit_p3)
deviance(fit_p5)
deviance(fit_p6)
deviance(fit_p8)
deviance(fit_p9)
summary(fit_p1)
summary(fit_p3)
summary(fit_p5)
summary(fit_p6)
summary(fit_p8)
summary(fit_p9)
summary(parametros_pk_VC)


# ggplot()+
#   geom_line( data =  simulacion_tac1, aes(x=tiempo_continuo, y=concentracion*1000, color = paciente),size = 0.7)+
#   geom_point(data = hplc_vc_mean, aes(x=t_min, y=concentration, col = patient),size = 2, alpha = 0.8)+
#   # geom_line(data = gpt_hat, aes(x=tiempo, y=concentracion),size = 0.7, col="#113b19")+
#   theme_bw()+
#   facet_grid(cols = vars(sol))+
#   guides(title=NULL)+
#   theme(legend.position = "bottom", legend.title = element_blank(), 
#         legend.text = element_text(size=8),
#         axis.text = element_text(size = 8))
# 


  #TAC 2

modelo_ajuste_vc_tac2 <- function(V, Cl, Cph) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- 90
  R <- 9.198274
  conc <- t
  conc[t <= t_cambio] <- Cph + (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <- Cph + 
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}
modelo_simulacion_vc_tac2 <- function(V, Cl, Cph) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- 90
  R <- 9.198274
  conc <- t
  conc[t <= t_cambio] <- Cph + (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <- Cph +
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}


while(TRUE){  

  concentracion_p10 <- filter(hplc_vc_mean, patient == "P10") %>% pull(concentration)/1000
  concentracion_p11 <- filter(hplc_vc_mean, patient == "P11") %>% pull(concentration)/1000
  concentracion_p14 <- filter(hplc_vc_mean, patient == "P14") %>% pull(concentration)/1000
  concentracion_p15 <- filter(hplc_vc_mean, patient == "P15") %>% pull(concentration)/1000
  concentracion_p17 <- filter(hplc_vc_mean, patient == "P17") %>% pull(concentration)/1000
  concentracion_p18 <- filter(hplc_vc_mean, patient == "P18") %>% pull(concentration)/1000

  
  fit_p10 <- nlsLM(concentracion_p10 ~ modelo_ajuste_vc_tac2(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p11 <- nlsLM(concentracion_p11 ~ modelo_ajuste_vc_tac2(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p14 <- nlsLM(concentracion_p14 ~ modelo_ajuste_vc_tac2(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p15 <- nlsLM(concentracion_p15 ~ modelo_ajuste_vc_tac2(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p17 <- nlsLM(concentracion_p17 ~ modelo_ajuste_vc_tac2(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))
  fit_p18 <- nlsLM(concentracion_p18 ~ modelo_ajuste_vc_tac2(V, Cl, Cph), start = list(V = 6, Cl = 5, Cph = 0))

deviance(fit_p10)
deviance(fit_p11)
deviance(fit_p14)
deviance(fit_p15)
deviance(fit_p17)
deviance(fit_p18)
summary(fit_p10)
summary(fit_p11)
summary(fit_p14)
summary(fit_p15)
summary(fit_p17)
summary(fit_p18)


parametros_pk_VC$V [7]  <- summary(fit_p10)$coefficients[1,1] # V
parametros_pk_VC$V [8]  <- summary(fit_p11)$coefficients[1,1] # V
parametros_pk_VC$V [9]  <- summary(fit_p14)$coefficients[1,1] # V
parametros_pk_VC$V [10] <- summary(fit_p15)$coefficients[1,1] # V
parametros_pk_VC$V [11] <- summary(fit_p17)$coefficients[1,1] # V
parametros_pk_VC$V [12] <- summary(fit_p18)$coefficients[1,1] # V

parametros_pk_VC$Cl[7]  <- summary(fit_p10)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[8]  <- summary(fit_p11)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[9]  <- summary(fit_p14)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[10] <- summary(fit_p15)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[11] <- summary(fit_p17)$coefficients[2,1] # Cl
parametros_pk_VC$Cl[12] <- summary(fit_p18)$coefficients[2,1] # Cl

parametros_pk_VC$Cph[7]  <- summary(fit_p10)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[8]  <- summary(fit_p11)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[9]  <- summary(fit_p14)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[10] <- summary(fit_p15)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[11] <- summary(fit_p17)$coefficients[3,1] # Cph
parametros_pk_VC$Cph[12] <- summary(fit_p18)$coefficients[3,1] # Cph

simulacion_tac2_p10 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac2(V = parametros_pk_VC$V[7],  Cl = parametros_pk_VC$Cl[7],  Cph = parametros_pk_VC$Cph[7]))   %>% mutate(paciente = "P10")
simulacion_tac2_p11 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac2(V = parametros_pk_VC$V[8],  Cl = parametros_pk_VC$Cl[8],  Cph = parametros_pk_VC$Cph[8]))   %>% mutate(paciente = "P11")
simulacion_tac2_p14 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac2(V = parametros_pk_VC$V[9],  Cl = parametros_pk_VC$Cl[9],  Cph = parametros_pk_VC$Cph[9]))   %>% mutate(paciente = "P14")
simulacion_tac2_p15 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac2(V = parametros_pk_VC$V[10], Cl = parametros_pk_VC$Cl[10], Cph = parametros_pk_VC$Cph[10])) %>% mutate(paciente = "P15")
simulacion_tac2_p17 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac2(V = parametros_pk_VC$V[11], Cl = parametros_pk_VC$Cl[11], Cph = parametros_pk_VC$Cph[11])) %>% mutate(paciente = "P17")
simulacion_tac2_p18 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_vc_tac2(V = parametros_pk_VC$V[12], Cl = parametros_pk_VC$Cl[12], Cph = parametros_pk_VC$Cph[12])) %>% mutate(paciente = "P18")

simulacion_tac2 <- bind_rows(
  simulacion_tac2_p10,
  simulacion_tac2_p11,
  simulacion_tac2_p14,
  simulacion_tac2_p15,
  simulacion_tac2_p17,
  simulacion_tac2_p18
  ) %>%
      mutate(sol = "TAC 2")
  break
}

simulacion_ambas_sol_vc <- bind_rows(simulacion_tac1, simulacion_tac2) 
datos_vc <- hplc_vc_mean %>%
    filter(sol != "P")

grafico_vc <- ggplot()+
  # geom_line(data = datos_vc, aes(x=t_min, y=concentration, col = patient),size = 1, alpha = 0.5, linetype = "solid")+
  geom_boxplot(data = datos_vc, aes(x=t_min, y=concentration, group = t_min),size = 1, alpha = 0.5, col = "#00000090")+
  geom_point(data = datos_vc, aes(x=t_min, y=concentration, col = patient),size = 4, alpha = 0.5)+
  # geom_boxplot(data = simulacion_ambas_sol_vc, aes(x=tiempo_continuo, y=concentracion*1000, group = tiempo_continuo),size = 1, fill = "#0000003e", col = "#772929b8")+
  geom_line( data =  simulacion_ambas_sol_vc, aes(x=tiempo_continuo, y=concentracion*1000, color = paciente),size = 0.7)+
  # geom_line(data = datos_vc, aes(x=t_min, y=concentration, col = patient),size = 1, alpha = 1, linetype = "dashed")+
  # geom_point(data = datos_vc, aes(x=t_min, y=concentration, col = patient),size = 4, alpha = 0.5)+
  # geom_line(data = gpt_hat, aes(x=tiempo, y=concentracion),size = 0.7, col="#113b19")+
  scale_x_continuous(breaks = seq(0, 180, 30))+
  facet_grid(cols = vars(sol))+
  # scale_color_viridis_d()+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grafico_vc

View(parametros_pk_VC)

summary(parametros_pk_VC)
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

# análisis NAC
# tabla parámetros por paciente NAC

(parametros_pk_NAC <- tibble(
  paciente = c("P1", "P3", "P5", "P6", "P8", "P9", "P10", "P11", "P14", "P15", "P17", "P18"),
  V = 0, Cl= 0
  )
)
# modelo PK 

# peso molecular NAC : 163.19 g/mol 

# mg/ml * ml/min * 60min/h /163.191 g/mol  / 1000
# R0 para TAC 1
# 2000/250*6*60/163.191
# 17.64803
# R0 para TAC 2
# 4000/250*3*60/163.191
# 17.64803

modelo_ajuste_nac_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- c(30, 90)
  
  R <- 17.64814
  conc <- t
  
  conc[t <= t_cambio[1]] <- (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60)) 
  conc[t <= t_cambio[2] & t > t_cambio[1]] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60))
  
  conc[t > t_cambio[2]] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * exp(-k*(t[t > t_cambio[2]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*((t_cambio[2] - t_cambio[1])/60))) * exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60))
  
  return(conc)
}

modelo_simulacion_nac_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- c(30, 90)
  R <- 17.64814
  conc <- t

  conc[t <= t_cambio[1]] <-  (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60))
  conc[t <= t_cambio[2] & t > t_cambio[1]] <-  
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * 
                            exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60))
  
  conc[t > t_cambio[2]] <-  
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * 
                            exp(-k*(t[t > t_cambio[2]] - t_cambio[1])/60) + # OJO, antes estaba restando el t_cambio[2] en vez de el t_cambio[1]

                  ((R/6)/Cl) * (1-exp(-k*((t_cambio[2] - t_cambio[1])/60))) * 
                            exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60) +

                  (0/Cl) * (1-exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60))
  
  return(conc)
}





  # TAC 1
while(TRUE){  
  concentracion_p1 <- filter(hplc_nac_mean, patient == "P1") %>% pull(concentration)/1000
  concentracion_p3 <- filter(hplc_nac_mean, patient == "P3") %>% pull(concentration)/1000
  concentracion_p5 <- filter(hplc_nac_mean, patient == "P5") %>% pull(concentration)/1000
  concentracion_p6 <- filter(hplc_nac_mean, patient == "P6") %>% pull(concentration)/1000
  concentracion_p8 <- filter(hplc_nac_mean, patient == "P8") %>% pull(concentration)/1000
  concentracion_p9 <- filter(hplc_nac_mean, patient == "P9") %>% pull(concentration)/1000
  
  fit_p1 <- nlsLM(concentracion_p1 ~ modelo_ajuste_nac_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p3 <- nlsLM(concentracion_p3 ~ modelo_ajuste_nac_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p5 <- nlsLM(concentracion_p5 ~ modelo_ajuste_nac_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p6 <- nlsLM(concentracion_p6 ~ modelo_ajuste_nac_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p8 <- nlsLM(concentracion_p8 ~ modelo_ajuste_nac_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p9 <- nlsLM(concentracion_p9 ~ modelo_ajuste_nac_tac1(V, Cl), start = list(V = 6, Cl = 5))




parametros_pk_NAC$V [1] <- summary(fit_p1)$coefficients[1,1] # V
parametros_pk_NAC$V [2] <- summary(fit_p3)$coefficients[1,1] # V
parametros_pk_NAC$V [3] <- summary(fit_p5)$coefficients[1,1] # V
parametros_pk_NAC$V [4] <- summary(fit_p6)$coefficients[1,1] # V
parametros_pk_NAC$V [5] <- summary(fit_p8)$coefficients[1,1] # V
parametros_pk_NAC$V [6] <- summary(fit_p9)$coefficients[1,1] # V

parametros_pk_NAC$Cl[1] <- summary(fit_p1)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[2] <- summary(fit_p3)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[3] <- summary(fit_p5)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[4] <- summary(fit_p6)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[5] <- summary(fit_p8)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[6] <- summary(fit_p9)$coefficients[2,1] # Cl



simulacion_tac1_p1 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac1(V = parametros_pk_NAC$V[1], Cl = parametros_pk_NAC$Cl[1])) %>%
    mutate(paciente = "P1")
simulacion_tac1_p3 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac1(V = parametros_pk_NAC$V[2], Cl = parametros_pk_NAC$Cl[2])) %>%
    mutate(paciente = "P3")
simulacion_tac1_p5 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac1(V = parametros_pk_NAC$V[3], Cl = parametros_pk_NAC$Cl[3])) %>%
    mutate(paciente = "P5")
simulacion_tac1_p6 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac1(V = parametros_pk_NAC$V[4], Cl = parametros_pk_NAC$Cl[4])) %>%
    mutate(paciente = "P6")
simulacion_tac1_p8 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac1(V = parametros_pk_NAC$V[5], Cl = parametros_pk_NAC$Cl[5])) %>%
    mutate(paciente = "P8")
simulacion_tac1_p9 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac1(V = parametros_pk_NAC$V[6], Cl = parametros_pk_NAC$Cl[6])) %>%
    mutate(paciente = "P9")

simulacion_tac1 <- bind_rows(
  simulacion_tac1_p1,
  simulacion_tac1_p3,
  simulacion_tac1_p5,
  simulacion_tac1_p6,
  simulacion_tac1_p8,
  simulacion_tac1_p9
  ) %>%
      mutate(sol = "TAC 1")
  break
}


  #TAC 2

modelo_ajuste_nac_tac2 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- 90
  R <- 17.64814
  conc <- t
  conc[t <= t_cambio] <-  (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <-  
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}
modelo_simulacion_nac_tac2 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- 90
  R <- 17.64814
  
  conc <- t
  conc[t <= t_cambio] <-  (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}


while(TRUE){  

  concentracion_p10 <- filter(hplc_nac_mean, patient == "P10") %>% pull(concentration)/1000
  concentracion_p11 <- filter(hplc_nac_mean, patient == "P11") %>% pull(concentration)/1000
  concentracion_p14 <- filter(hplc_nac_mean, patient == "P14") %>% pull(concentration)/1000
  concentracion_p15 <- filter(hplc_nac_mean, patient == "P15") %>% pull(concentration)/1000
  concentracion_p17 <- filter(hplc_nac_mean, patient == "P17") %>% pull(concentration)/1000
  concentracion_p18 <- filter(hplc_nac_mean, patient == "P18") %>% pull(concentration)/1000

  
  fit_p10 <- nlsLM(concentracion_p10 ~ modelo_ajuste_nac_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p11 <- nlsLM(concentracion_p11 ~ modelo_ajuste_nac_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p14 <- nlsLM(concentracion_p14 ~ modelo_ajuste_nac_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p15 <- nlsLM(concentracion_p15 ~ modelo_ajuste_nac_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p17 <- nlsLM(concentracion_p17 ~ modelo_ajuste_nac_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p18 <- nlsLM(concentracion_p18 ~ modelo_ajuste_nac_tac2(V, Cl), start = list(V = 6, Cl = 5))

parametros_pk_NAC$V [7]  <- summary(fit_p10)$coefficients[1,1] # V
parametros_pk_NAC$V [8]  <- summary(fit_p11)$coefficients[1,1] # V
parametros_pk_NAC$V [9]  <- summary(fit_p14)$coefficients[1,1] # V
parametros_pk_NAC$V [10] <- summary(fit_p15)$coefficients[1,1] # V
parametros_pk_NAC$V [11] <- summary(fit_p17)$coefficients[1,1] # V
parametros_pk_NAC$V [12] <- summary(fit_p18)$coefficients[1,1] # V

parametros_pk_NAC$Cl[7]  <- summary(fit_p10)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[8]  <- summary(fit_p11)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[9]  <- summary(fit_p14)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[10] <- summary(fit_p15)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[11] <- summary(fit_p17)$coefficients[2,1] # Cl
parametros_pk_NAC$Cl[12] <- summary(fit_p18)$coefficients[2,1] # Cl


simulacion_tac2_p10 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac2(V = parametros_pk_NAC$V[7],  Cl = parametros_pk_NAC$Cl[7]))   %>% mutate(paciente = "P10")
simulacion_tac2_p11 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac2(V = parametros_pk_NAC$V[8],  Cl = parametros_pk_NAC$Cl[8]))   %>% mutate(paciente = "P11")
simulacion_tac2_p14 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac2(V = parametros_pk_NAC$V[9],  Cl = parametros_pk_NAC$Cl[9]))   %>% mutate(paciente = "P14")
simulacion_tac2_p15 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac2(V = parametros_pk_NAC$V[10], Cl = parametros_pk_NAC$Cl[10])) %>% mutate(paciente = "P15")
simulacion_tac2_p17 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac2(V = parametros_pk_NAC$V[11], Cl = parametros_pk_NAC$Cl[11])) %>% mutate(paciente = "P17")
simulacion_tac2_p18 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_nac_tac2(V = parametros_pk_NAC$V[12], Cl = parametros_pk_NAC$Cl[12])) %>% mutate(paciente = "P18")

simulacion_tac2 <- bind_rows(
  simulacion_tac2_p10,
  simulacion_tac2_p11,
  simulacion_tac2_p14,
  simulacion_tac2_p15,
  simulacion_tac2_p17,
  simulacion_tac2_p18
  ) %>%
      mutate(sol = "TAC 2")
  break
}


deviance(fit_p1)
deviance(fit_p3)
deviance(fit_p5)
deviance(fit_p6)
deviance(fit_p8)
deviance(fit_p9)
deviance(fit_p10) # ojo
deviance(fit_p11)
deviance(fit_p14) # ojo
deviance(fit_p15)
deviance(fit_p17)
deviance(fit_p18)
summary(fit_p1)
summary(fit_p3)
summary(fit_p5)
summary(fit_p6)
summary(fit_p8)
summary(fit_p9)
summary(fit_p10)
summary(fit_p11)
summary(fit_p14)
summary(fit_p15)
summary(fit_p17)
summary(fit_p18)

simulacion_ambas_sol_nac <- bind_rows(simulacion_tac1, simulacion_tac2) 
datos_nac <- hplc_nac_mean %>%
    filter(sol != "P")

grafico_nac <- ggplot()+
  # geom_line(data = datos_nac, aes(x=t_min, y=concentration, col = patient),size = 1, alpha = 0.5, linetype = "solid")+
  geom_boxplot(data = datos_nac, aes(x=t_min, y=concentration, group = t_min),size = 1, alpha = 0.5, col = "#00000090")+
  # geom_boxplot(data = simulacion_ambas_sol_nac, aes(x=tiempo_continuo, y=concentracion*1000, group = tiempo_continuo),size = 1, fill = "#0000003e", col = "#772929b8")+
  geom_line( data =  simulacion_ambas_sol_nac, aes(x=tiempo_continuo, y=concentracion*1000, color = paciente),size = 0.7)+
  # geom_line(data = datos_nac, aes(x=t_min, y=concentration, col = patient),size = 1, linetype = "dashed")+
  geom_point(data = datos_nac, aes(x=t_min, y=concentration, col = patient),size = 4, alpha = 0.5)+
  # geom_line(data = gpt_hat, aes(x=tiempo, y=concentracion),size = 0.7, col="#113b19")+
  scale_x_continuous(breaks = seq(0, 180, 30))+
  facet_grid(cols = vars(sol))+
  # scale_color_viridis_d()+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grafico_nac
View(parametros_pk_NAC)

summary(parametros_pk_NAC)
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

# análisis DFO
# tabla parámetros por paciente NAC

(parametros_pk_DFO <- tibble(
  paciente = c("P1", "P3", "P5", "P6", "P8", "P9", "P10", "P11", "P14", "P15", "P17", "P18"),
  V = 0, Cl= 0
  )
)
# modelo PK 

# peso molecular DFO : 656.8 g/mol 

# mg/ml * ml/min * 60min/h /656.8 g/mol  / 1000
# R0 para TAC 1
#1000/250*6*60/656.8
#2.192448
# R0 para TAC 2
#1600/250*3*60/656.8
#1.753959


modelo_ajuste_dfo_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- c(30, 90)
  
  R <- 2.192448
  conc <- t
  
  conc[t <= t_cambio[1]] <- (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60)) 
  conc[t <= t_cambio[2] & t > t_cambio[1]] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60))
  
  conc[t > t_cambio[2]] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * exp(-k*(t[t > t_cambio[2]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*((t_cambio[2] - t_cambio[1])/60))) * exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60))
  
  return(conc)
}

modelo_simulacion_dfo_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- c(30, 90)
  R <- 2.192448
  conc <- t

  conc[t <= t_cambio[1]] <-  (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60))
  conc[t <= t_cambio[2] & t > t_cambio[1]] <-  
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * 
                            exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60) +
                  ((R/6)/Cl) * (1-exp(-k*(t[t <= t_cambio[2] & t > t_cambio[1]] - t_cambio[1])/60))
  
  conc[t > t_cambio[2]] <-  
                  (R/Cl) * (1-exp(-k*(t_cambio[1]/60))) * 
                            exp(-k*(t[t > t_cambio[2]] - t_cambio[1])/60) + # OJO, antes estaba restando el t_cambio[2] en vez de el t_cambio[1]

                  ((R/6)/Cl) * (1-exp(-k*((t_cambio[2] - t_cambio[1])/60))) * 
                            exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60) +

                  (0/Cl) * (1-exp(-k*(t[t > t_cambio[2]] - t_cambio[2])/60))
  
  return(conc)
}





  # TAC 1
while(TRUE){  
  concentracion_p1 <- filter(hplc_dfo_mean, patient == "P1") %>% pull(concentration)/1000
  concentracion_p3 <- filter(hplc_dfo_mean, patient == "P3") %>% pull(concentration)/1000
  concentracion_p5 <- filter(hplc_dfo_mean, patient == "P5") %>% pull(concentration)/1000
  concentracion_p6 <- filter(hplc_dfo_mean, patient == "P6") %>% pull(concentration)/1000
  concentracion_p8 <- filter(hplc_dfo_mean, patient == "P8") %>% pull(concentration)/1000
  concentracion_p9 <- filter(hplc_dfo_mean, patient == "P9") %>% pull(concentration)/1000
  
  fit_p1 <- nlsLM(concentracion_p1 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p3 <- nlsLM(concentracion_p3 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p5 <- nlsLM(concentracion_p5 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p6 <- nlsLM(concentracion_p6 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p8 <- nlsLM(concentracion_p8 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
  fit_p9 <- nlsLM(concentracion_p9 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))

parametros_pk_DFO$V [1] <- summary(fit_p1)$coefficients[1,1] # V
parametros_pk_DFO$V [2] <- summary(fit_p3)$coefficients[1,1] # V
parametros_pk_DFO$V [3] <- summary(fit_p5)$coefficients[1,1] # V
parametros_pk_DFO$V [4] <- summary(fit_p6)$coefficients[1,1] # V
parametros_pk_DFO$V [5] <- summary(fit_p8)$coefficients[1,1] # V
parametros_pk_DFO$V [6] <- summary(fit_p9)$coefficients[1,1] # V

parametros_pk_DFO$Cl[1] <- summary(fit_p1)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[2] <- summary(fit_p3)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[3] <- summary(fit_p5)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[4] <- summary(fit_p6)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[5] <- summary(fit_p8)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[6] <- summary(fit_p9)$coefficients[2,1] # Cl



simulacion_tac1_p1 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac1(V = parametros_pk_DFO$V[1], Cl = parametros_pk_DFO$Cl[1])) %>%
    mutate(paciente = "P1")
simulacion_tac1_p3 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac1(V = parametros_pk_DFO$V[2], Cl = parametros_pk_DFO$Cl[2])) %>%
    mutate(paciente = "P3")
simulacion_tac1_p5 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac1(V = parametros_pk_DFO$V[3], Cl = parametros_pk_DFO$Cl[3])) %>%
    mutate(paciente = "P5")
simulacion_tac1_p6 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac1(V = parametros_pk_DFO$V[4], Cl = parametros_pk_DFO$Cl[4])) %>%
    mutate(paciente = "P6")
simulacion_tac1_p8 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac1(V = parametros_pk_DFO$V[5], Cl = parametros_pk_DFO$Cl[5])) %>%
    mutate(paciente = "P8")
simulacion_tac1_p9 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac1(V = parametros_pk_DFO$V[6], Cl = parametros_pk_DFO$Cl[6])) %>%
    mutate(paciente = "P9")

simulacion_tac1 <- bind_rows(
  simulacion_tac1_p1,
  simulacion_tac1_p3,
  simulacion_tac1_p5,
  simulacion_tac1_p6,
  simulacion_tac1_p8,
  simulacion_tac1_p9
  ) %>%
      mutate(sol = "TAC 1")
  break
}

nlsLM(concentracion_p1 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p1 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 100, Cl = 200))
nlsLM(concentracion_p1 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 200, Cl = 5))
nlsLM(concentracion_p1 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 4, Cl = 500))

nlsLM(concentracion_p3 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p3 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 100, Cl = 200))
nlsLM(concentracion_p3 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 200, Cl = 80))
nlsLM(concentracion_p3 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 1, Cl = 500))

nlsLM(concentracion_p5 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p5 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 100, Cl = 200))
nlsLM(concentracion_p5 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 200, Cl = 5))
nlsLM(concentracion_p5 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 5.2, Cl = 500))

nlsLM(concentracion_p6 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p6 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 100, Cl = 200))
nlsLM(concentracion_p6 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 500, Cl = 5))
nlsLM(concentracion_p6 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 4, Cl = 500))

nlsLM(concentracion_p8 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p8 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 500, Cl = 500))
nlsLM(concentracion_p8 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 500, Cl = 5))
nlsLM(concentracion_p8 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 4, Cl = 500))

nlsLM(concentracion_p9 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p9 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 500, Cl =500))
nlsLM(concentracion_p9 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 500, Cl = 5))
nlsLM(concentracion_p9 ~ modelo_ajuste_dfo_tac1(V, Cl), start = list(V = 4, Cl = 500))

  #TAC 2

modelo_ajuste_dfo_tac2 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- 90
  R <- 1.753959
  conc <- t
  conc[t <= t_cambio] <-  (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <-  
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}
modelo_simulacion_dfo_tac2 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- 90
  R <- 1.753959
  
  conc <- t
  conc[t <= t_cambio] <-  (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}


while(TRUE){  

  concentracion_p10 <- filter(hplc_dfo_mean, patient == "P10") %>% pull(concentration)/1000
  concentracion_p11 <- filter(hplc_dfo_mean, patient == "P11") %>% pull(concentration)/1000
  concentracion_p14 <- filter(hplc_dfo_mean, patient == "P14") %>% pull(concentration)/1000
  concentracion_p15 <- filter(hplc_dfo_mean, patient == "P15") %>% pull(concentration)/1000
  concentracion_p17 <- filter(hplc_dfo_mean, patient == "P17") %>% pull(concentration)/1000
  concentracion_p18 <- filter(hplc_dfo_mean, patient == "P18") %>% pull(concentration)/1000

  
  fit_p10 <- nlsLM(concentracion_p10 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p11 <- nlsLM(concentracion_p11 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p14 <- nlsLM(concentracion_p14 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p15 <- nlsLM(concentracion_p15 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p17 <- nlsLM(concentracion_p17 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
  fit_p18 <- nlsLM(concentracion_p18 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))

parametros_pk_DFO$V [7]  <- summary(fit_p10)$coefficients[1,1] # V
parametros_pk_DFO$V [8]  <- summary(fit_p11)$coefficients[1,1] # V
parametros_pk_DFO$V [9]  <- summary(fit_p14)$coefficients[1,1] # V
parametros_pk_DFO$V [10] <- summary(fit_p15)$coefficients[1,1] # V
parametros_pk_DFO$V [11] <- summary(fit_p17)$coefficients[1,1] # V
parametros_pk_DFO$V [12] <- summary(fit_p18)$coefficients[1,1] # V

parametros_pk_DFO$Cl[7]  <- summary(fit_p10)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[8]  <- summary(fit_p11)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[9]  <- summary(fit_p14)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[10] <- summary(fit_p15)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[11] <- summary(fit_p17)$coefficients[2,1] # Cl
parametros_pk_DFO$Cl[12] <- summary(fit_p18)$coefficients[2,1] # Cl


simulacion_tac2_p10 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac2(V = parametros_pk_DFO$V[7],  Cl = parametros_pk_DFO$Cl[7]))   %>% mutate(paciente = "P10")
simulacion_tac2_p11 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac2(V = parametros_pk_DFO$V[8],  Cl = parametros_pk_DFO$Cl[8]))   %>% mutate(paciente = "P11")
simulacion_tac2_p14 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac2(V = parametros_pk_DFO$V[9],  Cl = parametros_pk_DFO$Cl[9]))   %>% mutate(paciente = "P14")
simulacion_tac2_p15 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac2(V = parametros_pk_DFO$V[10], Cl = parametros_pk_DFO$Cl[10])) %>% mutate(paciente = "P15")
simulacion_tac2_p17 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac2(V = parametros_pk_DFO$V[11], Cl = parametros_pk_DFO$Cl[11])) %>% mutate(paciente = "P17")
simulacion_tac2_p18 <- tibble(tiempo_continuo = c(0:180),concentracion = modelo_simulacion_dfo_tac2(V = parametros_pk_DFO$V[12], Cl = parametros_pk_DFO$Cl[12])) %>% mutate(paciente = "P18")

simulacion_tac2 <- bind_rows(
  simulacion_tac2_p10,
  simulacion_tac2_p11,
  simulacion_tac2_p14,
  simulacion_tac2_p15,
  simulacion_tac2_p17,
  simulacion_tac2_p18
  ) %>%
      mutate(sol = "TAC 2")
  break
}


nlsLM(concentracion_p15 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p15 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 500, Cl = 500))
nlsLM(concentracion_p15 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 500))
nlsLM(concentracion_p15 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 500, Cl = 5))
nlsLM(concentracion_p15 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 4, Cl = 4))

nlsLM(concentracion_p17 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p17 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 500, Cl = 500))
nlsLM(concentracion_p17 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 500))
nlsLM(concentracion_p17 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 500, Cl = 5))
nlsLM(concentracion_p17 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 4, Cl = 4))

nlsLM(concentracion_p18 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p18 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 500, Cl = 500))
nlsLM(concentracion_p18 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 500))
nlsLM(concentracion_p18 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 500, Cl = 5))
nlsLM(concentracion_p18 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 4, Cl = 4))

nlsLM(concentracion_p17 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))
nlsLM(concentracion_p18 ~ modelo_ajuste_dfo_tac2(V, Cl), start = list(V = 6, Cl = 5))



simulacion_ambas_sol_dfo <- bind_rows(simulacion_tac1, simulacion_tac2) 
datos_dfo <- hplc_dfo_mean %>%
    filter(sol != "P")

grafico_dfo <- ggplot()+
  # geom_line(data = datos_dfo, aes(x=t_min, y=concentration, col = patient),size = 1, alpha = 0.5, linetype = "solid")+
  geom_boxplot(data = datos_dfo, aes(x=t_min, y=concentration, group = t_min),size = 1, alpha = 0.5, col = "#00000090")+
  # geom_boxplot(data = simulacion_ambas_sol_dfo, aes(x=tiempo_continuo, y=concentracion*1000, group = tiempo_continuo),size = 1, fill = "#0000003e", col = "#772929b8")+
  geom_line( data =  simulacion_ambas_sol_dfo, aes(x=tiempo_continuo, y=concentracion*1000, color = paciente),size = 0.7)+
  # geom_line(data = datos_dfo, aes(x=t_min, y=concentration, col = patient),size = 1, linetype = "dashed")+
  geom_point(data = datos_dfo, aes(x=t_min, y=concentration, col = patient),size = 4, alpha = 0.5)+
  scale_x_continuous(breaks = seq(0, 180, 30))+
  facet_grid(cols = vars(sol))+
  # scale_color_viridis_d()+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grafico_dfo
View(parametros_pk_DFO)

summary(parametros_pk_VC)
summary(parametros_pk_NAC)
summary(parametros_pk_DFO)
View(parametros_pk_VC)
View(parametros_pk_NAC)
View(parametros_pk_DFO)

densidad_vc_volumen <- density(parametros_pk_VC$V)
plot(densidad_vc_volumen)
densidad_vc_clearance <- density(parametros_pk_VC$Cl)
plot(densidad_vc_clearance)

densidad_nac_volumen <- density(parametros_pk_NAC$V)
plot(densidad_nac_volumen)
densidad_nac_clearance <- density(parametros_pk_NAC$Cl)
plot(densidad_nac_clearance)

densidad_dfo_volumen <- density(parametros_pk_DFO$V)
plot(densidad_dfo_volumen)
densidad_dfo_clearance <- density(parametros_pk_DFO$Cl)
plot(densidad_dfo_clearance)


grafico_vc
grafico_nac
grafico_dfo

# analisis de parámetros

caracteristicas_px <- tibble(read.csv("caracteristicas_px.csv", check.names = FALSE, fileEncoding = 'UTF-8-BOM'), .name_repair = "minimal")
todos_pk <- rename(parametros_pk_VC, V_vc = V, Cl_vc = Cl, Cph_vc = Cph)  %>%
    bind_cols(select(parametros_pk_NAC, V, Cl)  %>%
        rename(V_nac = V, Cl_nac = Cl)) %>%
    bind_cols(select(parametros_pk_DFO, V, Cl)  %>%
        rename(V_dfo = V, Cl_dfo = Cl))
write.csv(todos_pk, "00 parametros PK todos (DSM).csv", row.names  = FALSE)

caracteristicas_parametros <- bind_cols(caracteristicas_px, select(todos_pk, -c(paciente) ))
View(caracteristicas_parametros)

plot(select(caracteristicas_parametros, c(-c(paciente), -c(Cl_vc:Cl_dfo))))
plot(select(caracteristicas_parametros, -c(c(paciente, comorbid_cardiovasc, otras_comorbids, tratamiento), c(Cl_vc:Cl_dfo))))
plot(select(caracteristicas_parametros, -c(c(paciente, comorbid_cardiovasc, otras_comorbids, tratamiento, V_vc), c(V_nac:Cl_dfo))))

# install.packages("ggpubr")
library(ggpubr)
v_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  # stat_cor(label = paste0("R^2 = ", round(cor(x, y)^2, 2)), 
  #          label.x = 7, label.y = 16)+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = talla_cm, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = imc, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = edad, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = tbq, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(v_masa, v_talla, v_imc, v_edad, v_volumensang, v_tbq)


cl_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = talla_cm, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = imc, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = edad, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = tbq, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(cl_masa, cl_talla, cl_imc, cl_edad, cl_volumensang, cl_tbq)

# NAC
v_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  # stat_cor(label = paste0("R^2 = ", round(cor(x, y)^2, 2)), 
  #          label.x = 7, label.y = 16)+
  guides(title=NULL)+
  facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_masa
v_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = talla_cm, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_talla
v_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = imc, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_imc
v_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = edad, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = tbq, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(v_masa, v_talla, v_imc, v_edad, v_volumensang, v_tbq)

cl_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = talla_cm, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = imc, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = edad, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = tbq, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(cl_masa, cl_talla, cl_imc, cl_edad, cl_volumensang, cl_tbq)

# DFO
v_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  # stat_cor(label = paste0("R^2 = ", round(cor(x, y)^2, 2)), 
  #          label.x = 7, label.y = 16)+
  guides(title=NULL)+
  # facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
# v_masa
v_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = talla_cm, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  # facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
# v_talla
v_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = imc, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  # facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
# v_imc
v_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = edad, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = tbq, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(v_masa, v_talla, v_imc, v_edad, v_volumensang, v_tbq)

cl_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = talla_cm, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = imc, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = edad, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = tbq, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(cl_masa, cl_talla, cl_imc, cl_edad, cl_volumensang, cl_tbq)
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################

v_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  # stat_cor(label = paste0("R^2 = ", round(cor(x, y)^2, 2)), 
  #          label.x = 7, label.y = 16)+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = talla_cm, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = imc, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = edad, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = V_vc), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = tbq, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(v_masa, v_talla, v_imc, v_edad, v_volumensang, v_tbq)


cl_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = talla_cm, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = imc, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = edad, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = Cl_vc), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = tbq, y = V_vc), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(cl_masa, cl_talla, cl_imc, cl_edad, cl_volumensang, cl_tbq)

# NAC
v_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  # stat_cor(label = paste0("R^2 = ", round(cor(x, y)^2, 2)), 
  #          label.x = 7, label.y = 16)+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))

v_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = talla_cm, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))

v_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = imc, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))

v_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = edad, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = V_nac), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = tbq, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(v_masa, v_talla, v_imc, v_edad, v_volumensang, v_tbq)

cl_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = talla_cm, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = imc, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = edad, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = Cl_nac), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = tbq, y = V_nac), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(cl_masa, cl_talla, cl_imc, cl_edad, cl_volumensang, cl_tbq)

# DFO
v_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  # stat_cor(label = paste0("R^2 = ", round(cor(x, y)^2, 2)), 
  #          label.x = 7, label.y = 16)+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  # facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
# v_masa
v_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = talla_cm, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  # facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
# v_talla
v_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = imc, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  # facet_grid(cols = vars(tratamiento))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
# v_imc
v_volumensang <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = volumen_sanguineo_lbb, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = edad, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
v_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = V_dfo), size = 5, alpha = 0.5, col="#0a0aff")+
  geom_smooth(aes(x = tbq, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(v_masa, v_talla, v_imc, v_edad, v_volumensang, v_tbq)

cl_masa <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = masa_corporal_kg, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = masa_corporal_kg, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_talla <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = talla_cm, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = talla_cm, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_imc <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = imc, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = imc, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_volumensang <- ggplot(data = caracteristicas_parametros)+
  facet_grid(cols = vars(tratamiento))+
  geom_point(aes(x = volumen_sanguineo_lbb, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = volumen_sanguineo_lbb, y = V_dfo, group = tratamiento), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_edad <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = edad, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = edad, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
cl_tbq <- ggplot(data = caracteristicas_parametros)+
  geom_point(aes(x = tbq, y = Cl_dfo), size = 5, alpha = 0.5, col="#ff7c0a")+
  geom_smooth(aes(x = tbq, y = V_dfo), formula = y ~ x, method = "lm", se = TRUE, color = "#3e005b")+
  facet_grid(cols = vars(tratamiento))+
  guides(title=NULL)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8))
grid.arrange(cl_masa, cl_talla, cl_imc, cl_edad, cl_volumensang, cl_tbq)





