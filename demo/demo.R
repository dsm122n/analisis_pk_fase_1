# import libraries
install.packages("minpack.lm")
library(minpack.lm)
library(ggplot2)
# import data

asc <- read.csv("covariables_asc.csv", header = TRUE)[, -c(5:13)]
# only tac == tac1 patients
asc_filtered <- asc[asc$tac == "tac1" & asc$px == "p01", -c(2,4)]

# análisis VC
# tabla parámetros por paciente VC

# mg/ml * ml/min * 60min/h / 176.12g/mol  / 1000
# R0 para TAC 1
# 2475/250*6*60/176.12
# 20.2362

# R0 para TAC 2
# 2250/250*3*60/176.12
# 9.198274


# define function

modelo_ajuste_vc_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- c(30, 90)
  
  R <- 20.2362
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

# fit model

modelo_ajuste_vc_tac1_lm <- nlsLM(c_promedio/1000 ~ modelo_ajuste_vc_tac1(V, Cl), data = asc_filtered, start = list(V = 1, Cl = 1))

# animated LM algorithm demo showing the convergence steps and the final result and sum of squares  


modelo_simulacion_vc_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- c(30, 90)
  R <- 20.2362
  conc <- t

  conc[t <= t_cambio[1]] <- (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60))
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
V_vc <- summary(modelo_ajuste_vc_tac1_lm)$coefficients[1,1]
Cl_vc <- summary(modelo_ajuste_vc_tac1_lm)$coefficients[2,1]
simulacion_vc <- data.frame(tiempo = c(0:180), conc = modelo_simulacion_vc_tac1(V_vc, Cl_vc))

# plot
ggplot() +
    geom_point(data = asc_filtered, aes(x = Tiempo, y = c_promedio/1000), color = "blue")+
    geom_line(data = simulacion_vc, aes(x = tiempo, y = conc), color = "red")+
    labs(x = "Tiempo (min)", y = "Concentración (mmol/L)")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("vc_tac1.png", width = 10, height = 10, units = "cm")
# import libraries
install.packages("minpack.lm")
library(minpack.lm)
library(ggplot2)
# import data

asc <- read.csv("covariables_asc.csv", header = TRUE)[, -c(5:13)]
# only tac == tac1 patients
asc_filtered <- asc[asc$tac == "tac1" & asc$px == "p01", -c(2,4)]

# análisis VC
# tabla parámetros por paciente VC

# mg/ml * ml/min * 60min/h / 176.12g/mol  / 1000
# R0 para TAC 1
# 2475/250*6*60/176.12
# 20.2362

# R0 para TAC 2
# 2250/250*3*60/176.12
# 9.198274


# define function

modelo_ajuste_vc_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- c(30, 90)
  
  R <- 20.2362
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

# fit model

modelo_ajuste_vc_tac1_lm <- nlsLM(c_promedio/1000 ~ modelo_ajuste_vc_tac1(V, Cl), data = asc_filtered, start = list(V = 1, Cl = 1))

# animated LM algorithm demo showing the convergence steps and the final result and sum of squares  


# This function simulates the plasma concentration of a drug after a single dose
# of, say, 20 mg, administered by IV bolus injection, assuming a 1-compartment model
# with first-order elimination.
# V: volume of distribution, in L
# Cl: clearance, in L/h
# R: dose, in mg
# t: time, in min
# t_cambio: time of change in the infusion rate, in min
modelo_simulacion_vc_tac1 <- function(V, Cl, R, t, t_cambio) {
  k <- Cl / V
  t_cambio <- c(30, 90)
  conc <- t

  conc[t <= t_cambio[1]] <- (R/Cl) * (1 - exp(-k * t[t <= t_cambio[1]]/60))
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


V_vc <- summary(modelo_ajuste_vc_tac1_lm)$coefficients[1,1]
Cl_vc <- summary(modelo_ajuste_vc_tac1_lm)$coefficients[2,1]
simulacion_vc <- data.frame(tiempo = c(0:180), conc = modelo_simulacion_vc_tac1(V_vc, Cl_vc))

# plot
ggplot() +
    geom_point(data = asc_filtered, aes(x = Tiempo, y = c_promedio/1000), color = "blue")+
    geom_line(data = simulacion_vc, aes(x = tiempo, y = conc), color = "red")+
    labs(x = "Tiempo (min)", y = "Concentración (mmol/L)")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("vc_tac1.png", width = 10, height = 10, units = "cm")
