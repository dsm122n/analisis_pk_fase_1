install.packages("gridExtra")
install.packages("forcats")
install.packages("dplyr")
install.packages("tibble")
install.packages("ggplot2")
install.packages("minpack.lm")
install.packages("rsq") 
install.packages("DescTools")
while(TRUE){
  
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(gridExtra)
  library(forcats)
  library(httpgd)
  library(minpack.lm)
  library(rsq) 
  library(DescTools)
  break
}

# import data covariables
asc <- tibble(read.csv("output/covariables_asc.csv", header = TRUE, sep = ",")[,-13])
nac <- tibble(read.csv("output/covariables_nac.csv", header = TRUE, sep = ",")[,-13])
dfo <- tibble(read.csv("output/covariables_dfo.csv", header = TRUE, sep = ",")[,-13])

# ASC analysis

# area under the curve for each patient
asc_nca <- tibble(
  px = c("p01", "p03", "p05", "p06", "p08", "p09", "p10", "p11", "p14",
          "p15", "p17", "p18"),
  auc = c(0))

# group by patient (px), then calculate auc by patient
for (paciente in unique(asc$px)) {
  # filter by patient
  asc_px <- filter(asc, px == paciente)
  # calculate auc
  asc_nca$auc[asc_nca$px == paciente] <- AUC(asc_px$Tiempo/60, asc_px$c_promedio/1000) # convert to hours and mmol/L
}
# shapiro test for normality
shapiro.test(asc_nca$auc)

# calculate elimination constant (k) for each patient using linear regression at times 90-180 min
asc_nca <- mutate(asc_nca, k = 0)
for (paciente in unique(asc$px)) {
  # filter by patient
  asc_px <- filter(asc, px == paciente)
  # calculate k
  asc_nca$k[asc_nca$px == paciente] <- lm(formula = log(asc_px$c_promedio[asc_px$Tiempo >= 90]) ~ asc_px$Tiempo[asc_px$Tiempo >= 90])$coefficients[2] 
}
# convert k in micromol/L/min to mmol/L/h
asc_nca <- mutate(asc_nca, k = k * 60 / 1000)
# shapiro test for normality
shapiro.test(asc_nca$k)

# calculate elimination half-life (t1/2) for each patient
asc_nca <- mutate(asc_nca, t12 = 0)
for (paciente in unique(asc$px)) {
  # filter by patient
  asc_px <- filter(asc, px == paciente)
  # calculate t1/2
  asc_nca$t12[asc_nca$px == paciente] <- -log(2) / asc_nca$k[asc_nca$px == paciente]
}
# shapiro test for normality
shapiro.test(asc_nca$t12)

# calculate clearance (Cl) for each patient


# calculate 

# Cálculo velocidad infusión: 
# [mg/ml] * $ml/min$ * 60min/h / 176.12g/mol  / 1000

# R0 para TAC 1
# 2475/250*6*60/176.12
# 20.2362
# R1 para TAC 1 = R0/6, R2 para TAC 1 = 0
# dosis total

# R0 para TAC 2
# 2250/250*3*60/176.12
# 9.198274
# R1 para TAC 2 = 0

# ASC pk parameters table
(parametros_pk_VC <- tibble(
  px = c("P01", "P03", "P05", "P06", "P08", "P09", "P10", "P11", "P14", "P15", "P17", "P18"),
  V = 0, Cl= 0, Cph = 0
  )
)

# PK model
# model to adjust
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
# model to simulate
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

# group by patient (px), then fit c_promedio from tac1 by patient with nlsLM
for (paciente in unique(asc$tac)) {
  # filter tac1 c_promedio date
  asc <- filter(asc, asc$tac == paciente)
  
  # fit model
  fit <- nlsLM(
    conc ~ modelo_ajuste_vc_tac1(V, Cl, Cph),
    data = asc$c_promedio / 1000,
    start = list(V = 1, Cl = 1, Cph = 1),
    trace = TRUE,
    control = list(maxiter = 1000)
  )
  
  # save parameters
  parametros_pk_VC[parametros_pk_VC$px == paciente, 2:4] <- coef(fit)
  
  # plot
  p <- ggplot(asc, aes(x = time, y = conc)) +
    geom_point() +
    geom_line(aes(y = modelo_simulacion_vc_tac1(V, Cl, Cph)), color = "red") +
    labs(title = paste("Paciente", paciente)) +
    theme_bw()
  print(p)
}

for (i in unique(asc$px)) {
    # filter tac1 
    asc_filtered <- filter(asc, asc$tac == "tac1")
    fit <- nlsLM(
        conc ~ modelo_ajuste_vc_tac1(V, Cl, Cph),
        data = filter(asc_filtered, asc_filtered$px == i)$c_promedio / 1000,
        start = list(V = 1, Cl = 1, Cph = 1),
        trace = TRUE,
        control = list(maxiter = 1000)
        )
  
    # save parameters
    # parametros_pk_VC[i, 2:4] <- coef(fit)
  #p <- ggplot(filter(asc, asc$tac == tac1), aes(x = Tiempo, y = c_promedio)) +
  #      geom_point() +
  #      geom_line(aes(y = modelo_simulacion_vc_tac1(V, Cl, Cph)), color = "red") +
  #      labs(title = paste("Paciente", i)) +
  #      theme_bw()
  #  print(p)
}
