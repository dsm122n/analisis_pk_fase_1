# install.packages("gridExtra")
# install.packages("forcats")
# install.packages("dplyr")
# install.packages("tibble")
# install.packages("ggplot2")
# install.packages("minpack.lm")
# install.packages("rsq") 
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

# import data covariables
asc <- tibble(read.csv("output/covariables_asc_corrected.csv", header = TRUE, sep = ",")[,-13])
nac <- tibble(read.csv("output/covariables_nac_fixed.csv", header = TRUE, sep = ",")[,-13])
dfo <- tibble(read.csv("output/covariables_dfo_fixed.csv", header = TRUE, sep = ",")[,-13])
#
# ASC analysis

# Cálculo velocidad infusión: 
# [mg/ml] * $ml/min$ * 60min/h / 176.12g/mol  / 1000

# R0 para TAC 1
# 2475/250*6*60/176.12
# 20.2362
# R1 para TAC 1 = R0/6, R2 para TAC 1 = 0

# R0 para TAC 2
# 2250/250*3*60/176.12
# 9.198274
# R1 para TAC 2 = 0

# ASC pk parameters table
(pk_asc <- tibble(
  px = c("p01", "p03", "p05", "p06", "p08", "p09", "p10", "p11", "p14", "p15", "p17", "p18"),
  V = 0, Cl = 0
  )
)

# PK model
# model to adjust
modelo_ajuste_asc_tac1 <- function(V, Cl) {
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
# model to simulate
modelo_simulacion_asc_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- c(30, 90)
  R <- 20.2362
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


modelo_ajuste_asc_tac2 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- 90
  R <- 9.198274
  conc <- t
  conc[t <= t_cambio] <-  (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <-  
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}
modelo_simulacion_asc_tac2 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- 90
  R <- 9.198274
  
  conc <- t
  conc[t <= t_cambio] <-  (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}


# group by patient (px), then fit c_promedio from tac1 by patient with nlsLM

for (paciente in unique(asc$px[asc$tac == "tac1"])) {
  # fit model
  fit  <- nlsLM(
    c_promedio/1000 ~ modelo_ajuste_asc_tac1(V, Cl),
    data = asc[asc$px == paciente,],
    start = list(V = 1, Cl = 1)
  )
  # save parameters
  pk_asc$V[pk_asc$px == paciente]  <- summary(fit)$coefficients[1,1]
  pk_asc$Cl[pk_asc$px == paciente]  <- summary(fit)$coefficients[2,1]
}
for (paciente in unique(asc$px[asc$tac == "tac2"])) {
  # fit model
  fit  <- nlsLM(
    c_promedio/1000 ~ modelo_ajuste_asc_tac2(V, Cl),
    data = asc[asc$px == paciente,],
    start = list(V = 1, Cl = 1)
  )
  # save parameters
  pk_asc$V[pk_asc$px == paciente]  <- summary(fit)$coefficients[1,1]
  pk_asc$Cl[pk_asc$px == paciente]  <- summary(fit)$coefficients[2,1]
}

# graph simulated data and real date
# simulated data tibble
# tac1
asc_simulated_tac1 <- tibble (Tiempo = rep(c(0:180), 6),
                              tac = rep("tac1", 181*6),
                              px = rep(unique(asc$px[asc$tac == "tac1"]), each = 181),
                              c_simulated = c(0))
for(i in unique(asc$px[asc$tac == "tac1"])) {
  asc_simulated_tac1$c_simulated[asc_simulated_tac1$px == i] <- 1000*modelo_simulacion_asc_tac1(pk_asc$V[pk_asc$px == i], pk_asc$Cl[pk_asc$px == i])
}
# tac2
asc_simulated_tac2 <- tibble (Tiempo = rep(c(0:180), 6),
                              tac = rep("tac2", 181*6),
                              px = rep(unique(asc$px[asc$tac == "tac2"]), each = 181),
                              c_simulated = c(0))
for(i in unique(asc$px[asc$tac == "tac2"])) {
  asc_simulated_tac2$c_simulated[asc_simulated_tac2$px == i] <- 1000*modelo_simulacion_asc_tac2(pk_asc$V[pk_asc$px == i], pk_asc$Cl[pk_asc$px == i])
}
#bind rows
asc_simulated <- rbind(asc_simulated_tac1, asc_simulated_tac2)
# drop tac = p
asc_simulated <- asc_simulated[asc_simulated$tac != "p",]
asc_without_p <- asc[asc$tac != "p",]
# graph
ggplot()+
  geom_point(data = asc_without_p, aes(x = Tiempo, y = c_promedio, color = px))+
  geom_line(data = asc_simulated, aes(x = Tiempo, y = c_simulated, color = px))+
  facet_wrap(~tac)+
  scale_x_continuous(breaks = seq(0, 180, 30))+ 
  xlab("Time (min)")+
  ylab("Concentration (micromol/L)")+
  labs(title = "ASC")+
  theme_bw()
ggsave("output/comp_analysis/asc.png", width = 15, height = 10, units = "cm", dpi = 1000)
# summary tibble for pk_asc V and Cl

pk_asc_summary <- tibble(pk_parameter = c("V", "Cl"),
                         mean = c(mean(pk_asc$V), mean(pk_asc$Cl)),
                          sd = c(sd(pk_asc$V), sd(pk_asc$Cl)),
                          median = c(median(pk_asc$V), median(pk_asc$Cl)),
                          Q1 = c(quantile(pk_asc$V, 0.25), quantile(pk_asc$Cl, 0.25)),
                          Q3 = c(quantile(pk_asc$V, 0.75), quantile(pk_asc$Cl, 0.75)),
                          shapiro_pvalue = c(shapiro.test(pk_asc$V)$p.value, shapiro.test(pk_asc$Cl)$p.value))
View(pk_asc_summary)
write.csv(pk_asc_summary, "output/comp_analysis/pk_asc_summary.csv", row.names = FALSE)
write.csv(pk_asc, "output/comp_analysis/pk_asc.csv", row.names = FALSE)
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################



# nac analysis

# peso molecular NAC : 163.19 g/mol 

# mg/ml * ml/min * 60min/h /163.191 g/mol  / 1000
# R0 para TAC 1
# 2000/250*6*60/163.191 
# 17.64803
# R0 para TAC 2
# 4000/250*3*60/163.191
# 17.64803

# ASC pk parameters table
(pk_nac <- tibble(
  px = c("p01", "p03", "p05", "p06", "p08", "p09", "p10", "p11", "p14", "p15", "p17", "p18"),
  V = 0, Cl = 0
  )
)

# PK model
# model to adjust
modelo_ajuste_nac_tac1 <- function(V, Cl) {
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
# model to simulate
modelo_simulacion_nac_tac1 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0:180)
  t_cambio <- c(30, 90)
  R <- 20.2362
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


modelo_ajuste_nac_tac2 <- function(V, Cl) {
  k <- Cl / V
  t <- c(0, 15, 30, 60, 90, 120, 180)
  t_cambio <- 90
  R <- 9.198274
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
  R <- 9.198274
  
  conc <- t
  conc[t <= t_cambio] <-  (R/Cl) * (1 - exp(-k * t[t<=t_cambio]/60))

  conc[t > t_cambio] <- 
                  (R/Cl) * (1-exp(-k*(t_cambio/60))) * exp(-k*(t[t > t_cambio] - t_cambio)/60) +
                  (0/Cl) * (1-exp(-k*(t[t > t_cambio] - t_cambio)/60))  
  return(conc)
}


# group by patient (px), then fit c_promedio from tac1 by patient with nlsLM

for (paciente in unique(nac$px[nac$tac == "tac1"])) {
  # fit model
  fit  <- nlsLM(
    c_promedio/1000 ~ modelo_ajuste_nac_tac1(V, Cl),
    data = nac[nac$px == paciente,],
    start = list(V = 1, Cl = 1)
  )
  # save parameters
  pk_nac$V[pk_nac$px == paciente]  <- summary(fit)$coefficients[1,1]
  pk_nac$Cl[pk_nac$px == paciente]  <- summary(fit)$coefficients[2,1]
}
for (paciente in unique(nac$px[nac$tac == "tac2"])) {
  # fit model
  fit  <- nlsLM(
    c_promedio/1000 ~ modelo_ajuste_nac_tac2(V, Cl),
    data = nac[nac$px == paciente,],
    start = list(V = 1, Cl = 1)
  )
  # save parameters
  pk_nac$V[pk_nac$px == paciente]  <- summary(fit)$coefficients[1,1]
  pk_nac$Cl[pk_nac$px == paciente]  <- summary(fit)$coefficients[2,1]
}

# graph simulated data and real date
# simulated data tibble
# tac1
nac_simulated_tac1 <- tibble (Tiempo = rep(c(0:180), 6),
                              tac = rep("tac1", 181*6),
                              px = rep(unique(nac$px[nac$tac == "tac1"]), each = 181),
                              c_simulated = c(0))
for(i in unique(nac$px[nac$tac == "tac1"])) {
  nac_simulated_tac1$c_simulated[nac_simulated_tac1$px == i] <- 1000*modelo_simulacion_nac_tac1(pk_nac$V[pk_nac$px == i], pk_nac$Cl[pk_nac$px == i])
}
# tac2
nac_simulated_tac2 <- tibble (Tiempo = rep(c(0:180), 6),
                              tac = rep("tac2", 181*6),
                              px = rep(unique(nac$px[nac$tac == "tac2"]), each = 181),
                              c_simulated = c(0))
for(i in unique(nac$px[nac$tac == "tac2"])) {
  nac_simulated_tac2$c_simulated[nac_simulated_tac2$px == i] <- 1000*modelo_simulacion_nac_tac2(pk_nac$V[pk_nac$px == i], pk_nac$Cl[pk_nac$px == i])
}
#bind rows
nac_simulated <- rbind(nac_simulated_tac1, nac_simulated_tac2)
# drop tac = p
nac_simulated <- nac_simulated[nac_simulated$tac != "p",]
nac_without_p <- nac[nac$tac != "p",]
# graph
ggplot()+
  geom_point(data = nac_without_p, aes(x = Tiempo, y = c_promedio, color = px))+
  geom_line(data = nac_simulated, aes(x = Tiempo, y = c_simulated, color = px))+
  facet_wrap(~tac)+
  scale_x_continuous(breaks = seq(0, 180, 30))+ 
  xlab("Time (min)")+
  ylab("Concentration (micromol/L)")+
  labs(title = "NAC")+
  theme_bw()
ggsave("output/comp_analysis/nac.png", width = 15, height = 10, units = "cm", dpi = 1000)

# summary tibble for pk_nac V and Cl

pk_nac_summary <- tibble(pk_parameter = c("V", "Cl"),
                         mean = c(mean(pk_nac$V), mean(pk_nac$Cl)),
                          sd = c(sd(pk_nac$V), sd(pk_nac$Cl)),
                          median = c(median(pk_nac$V), median(pk_nac$Cl)),
                          Q1 = c(quantile(pk_nac$V, 0.25), quantile(pk_nac$Cl, 0.25)),
                          Q3 = c(quantile(pk_nac$V, 0.75), quantile(pk_nac$Cl, 0.75)),
                          shapiro_pvalue = c(shapiro.test(pk_nac$V)$p.value, shapiro.test(pk_nac$Cl)$p.value))
View(pk_nac_summary)
write.csv(pk_nac_summary, "output/comp_analysis/pk_nac_summary.csv", row.names = FALSE)
write.csv(pk_nac, "output/comp_analysis/pk_nac.csv", row.names = FALSE)

#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################



# dfo analysis

# peso molecular dfo : 163.19 g/mol 


# mg/ml * ml/min * 60min/h /656.8 g/mol  / 1000
# R0 para TAC 1
#1000/250*6*60/656.8
#2.192448
# R0 para TAC 2
#1600/250*3*60/656.8
#1.753959

(pk_dfo <- tibble(
  px = c("p01", "p03", "p05", "p06", "p08", "p09", "p10", "p11", "p14", "p15", "p17", "p18"),
  V = 0, Cl = 0
  )
)

# PK model
# model to adjust
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
# model to simulate
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


# group by patient (px), then fit c_promedio from tac1 by patient with nlsLM

for (paciente in unique(dfo$px[dfo$tac == "tac1"])) {
  # fit model
  fit  <- nlsLM(
    c_promedio/1000 ~ modelo_ajuste_dfo_tac1(V, Cl),
    data = dfo[dfo$px == paciente,],
    start = list(V = 1, Cl = 1)
  )
  # save parameters
  pk_dfo$V[pk_dfo$px == paciente]  <- summary(fit)$coefficients[1,1]
  pk_dfo$Cl[pk_dfo$px == paciente]  <- summary(fit)$coefficients[2,1]
}
for (paciente in unique(dfo$px[dfo$tac == "tac2"])) {
  # fit model
  fit  <- nlsLM(
    c_promedio/1000 ~ modelo_ajuste_dfo_tac2(V, Cl),
    data = dfo[dfo$px == paciente,],
    start = list(V = 1, Cl = 1)
  )
  # save parameters
  pk_dfo$V[pk_dfo$px == paciente]  <- summary(fit)$coefficients[1,1]
  pk_dfo$Cl[pk_dfo$px == paciente]  <- summary(fit)$coefficients[2,1]
}



# graph simulated data and real date
# simulated data tibble
# tac1
dfo_simulated_tac1 <- tibble (Tiempo = rep(c(0:180), 6),
                              tac = rep("tac1", 181*6),
                              px = rep(unique(dfo$px[dfo$tac == "tac1"]), each = 181),
                              c_simulated = c(0))
for(i in unique(dfo$px[dfo$tac == "tac1"])) {
  dfo_simulated_tac1$c_simulated[dfo_simulated_tac1$px == i] <- 1000*modelo_simulacion_dfo_tac1(pk_dfo$V[pk_dfo$px == i], pk_dfo$Cl[pk_dfo$px == i])
}
# tac2
dfo_simulated_tac2 <- tibble (Tiempo = rep(c(0:180), 6),
                              tac = rep("tac2", 181*6),
                              px = rep(unique(dfo$px[dfo$tac == "tac2"]), each = 181),
                              c_simulated = c(0))
for(i in unique(dfo$px[dfo$tac == "tac2"])) {
  dfo_simulated_tac2$c_simulated[dfo_simulated_tac2$px == i] <- 1000*modelo_simulacion_dfo_tac2(pk_dfo$V[pk_dfo$px == i], pk_dfo$Cl[pk_dfo$px == i])
}
#bind rows
dfo_simulated <- rbind(dfo_simulated_tac1, dfo_simulated_tac2)
# drop tac = p
dfo_simulated <- dfo_simulated[dfo_simulated$tac != "p",]
dfo_without_p <- dfo[dfo$tac != "p",]
# graph
ggplot()+
  geom_point(data = dfo_without_p, aes(x = Tiempo, y = c_promedio, color = px))+
  geom_line(data = dfo_simulated, aes(x = Tiempo, y = c_simulated, color = px))+
  facet_wrap(~tac)+
  scale_x_continuous(breaks = seq(0, 180, 30))+ 
  xlab("Time (min)")+
  ylab("Concentration (micromol/L)")+
  labs(title = "DFO")+
  theme_bw()
ggsave("output/comp_analysis/dfo.png", width = 15, height = 10, units = "cm", dpi = 1000)
# summary tibble for pk_dfo V and Cl

pk_dfo_summary <- tibble(pk_parameter = c("V", "Cl"),
                         mean = c(mean(pk_dfo$V), mean(pk_dfo$Cl)),
                          sd = c(sd(pk_dfo$V), sd(pk_dfo$Cl)),
                          median = c(median(pk_dfo$V), median(pk_dfo$Cl)),
                          Q1 = c(quantile(pk_dfo$V, 0.25), quantile(pk_dfo$Cl, 0.25)),
                          Q3 = c(quantile(pk_dfo$V, 0.75), quantile(pk_dfo$Cl, 0.75)),
                          shapiro_pvalue = c(shapiro.test(pk_dfo$V)$p.value, shapiro.test(pk_dfo$Cl)$p.value))
View(pk_dfo_summary)
write.csv(pk_dfo_summary, "output/comp_analysis/pk_dfo_summary.csv", row.names = FALSE)
write.csv(pk_dfo, "output/comp_analysis/pk_dfo.csv", row.names = FALSE)

