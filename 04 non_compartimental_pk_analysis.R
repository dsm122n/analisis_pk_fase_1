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
asc <- tibble(read.csv("output/covariables_asc_corrected.csv", header = TRUE, sep = ",")[,-13])
nac <- tibble(read.csv("output/covariables_nac_fixed.csv", header = TRUE, sep = ",")[,-13])
dfo <- tibble(read.csv("output/covariables_dfo_fixed.csv", header = TRUE, sep = ",")[,-13])

# ASC analysis

# area under the curve for each patient
asc_nca <- tibble(
  px = c("p01", "p03", "p05", "p06", "p08", "p09", "p10", "p11", "p14",
          "p15", "p17", "p18"),
  tac = c("tac1","tac1","tac1","tac1","tac1","tac1","tac2","tac2","tac2","tac2","tac2","tac2"),
  auc = c(0))
asc_nca
View(asc)
# group by patient (px), then calculate auc by patient
for (paciente in unique(asc$px)) {
  # filter by patient
  asc_px <- filter(asc, px == paciente)
  # calculate auc
  asc_nca$auc[asc_nca$px == paciente] <- AUC((asc_px$Tiempo)/60, (asc_px$c_promedio)/1000) # convert to hours and mmol/L
}
# shapiro test for normality
shapiro.test(asc_nca$auc) 
# p-value = 0.4062 > 0.05, therefore data is normal
asc_nca
View(asc)


# calculate elimination constant (k) for each patient using linear regression at times 90-180 min
asc_nca <- mutate(asc_nca, k = 0)
asc_milimolar <- mutate(asc, c_promedio = c_promedio / 1000) # convert to mmol/L
# log-linear concentration vs time plot
asc_log_linear <- tibble(Tiempo = asc_milimolar$Tiempo[(asc_milimolar$tac == "tac1" | asc_milimolar$tac == "tac2") & asc_milimolar$Tiempo >=90], 
                              px = asc_milimolar$px[(asc_milimolar$tac == "tac1" | asc_milimolar$tac == "tac2") & asc_milimolar$Tiempo >= 90],
                              log_c_promedio = log(asc_milimolar$c_promedio[(asc_milimolar$tac == "tac1" | asc_milimolar$tac == "tac2") & asc_milimolar$Tiempo >= 90]))
View(asc_log_linear)
ggplot() +
  geom_point(data = asc_log_linear, aes(x = Tiempo, y = log_c_promedio, col = asc_log_linear$px))+
  geom_line(data = asc_log_linear, aes(x = Tiempo, y = log_c_promedio, col = asc_log_linear$px))+
  theme_bw()

for (paciente in unique(asc_log_linear$px)) {
  # filter by patient
  asc_px <- filter(asc_log_linear, px == paciente)
  # calculate k
  asc_nca$k[asc_nca$px == paciente] <- lm(log_c_promedio ~ Tiempo, data = asc_px) %>%
    summary() %>%
    .$coefficients %>%
    .[2] * (-1)
}
View(asc_nca)
# convert k from micromol/L/min to mmol/L/h
asc_nca <- mutate(asc_nca, k = k * 60)
# shapiro test for normality
shapiro.test(asc_nca$k)


# calculate elimination half-life (t1/2) for each patient
asc_nca <- mutate(asc_nca, t12 = 0)
for (paciente in unique(asc$px)) {
  # filter by patient
  asc_px <- filter(asc, px == paciente)
  # calculate t1/2
  asc_nca$t12[asc_nca$px == paciente] <- log(2) / asc_nca$k[asc_nca$px == paciente]
}
# shapiro test for normality
shapiro.test(asc_nca$t12)

# dosis total en mg tac1
# AA 2475 NAC 2000 DFO 1000
# 2475/176.12 = 14.05 asc
# 2000/163.2 = 12.25 nac
# 1000/560.6 = 1.78 dfo
# dosis total en mg tac2
# AA 2250 NAC 4000 DFO 1600
# 2250/176.12 = 12.77 asc
# 4000/163.2 = 24.51 nac
# 1600/560.6 = 2.85 dfo

# add dosis total to asc_nca
asc_nca <- mutate(asc_nca, dosis_total = 0)
for (i in asc_nca$px) {
  if (i == "p01" | i == "p03" | i == "p05" | i == "p06" | i == "p08" | i == "p09") {
    asc_nca$dosis_total[asc_nca$px == i] <- 14.05
  } else {
    asc_nca$dosis_total[asc_nca$px == i] <- 12.77
  }
}

# calculate Clearance (Cl) for each patient
asc_nca <- mutate(asc_nca, Cl = 0)
for (paciente in unique(asc$px)) {
  # filter by patient
  asc_px <- filter(asc, px == paciente)
  # calculate Cl
  asc_nca$Cl[asc_nca$px == paciente] <- asc_nca$dosis_total[asc_nca$px == paciente] / asc_nca$auc[asc_nca$px == paciente]
}
shapiro.test(asc_nca$Cl)
# calculate Volume of distribution (V) for each patient
asc_nca <- mutate(asc_nca, V = 0)
for (paciente in unique(asc$px)) {
  # filter by patient
  asc_px <- filter(asc, px == paciente)
  # calculate V
  asc_nca$V[asc_nca$px == paciente] <- asc_nca$Cl[asc_nca$px == paciente] / asc_nca$k[asc_nca$px == paciente]
}
shapiro.test(asc_nca$V)
View(summary(asc_nca))
# export asc_nca summary table to csv

write.csv(asc_nca, "output/non_comp_analysis/asc_nca.csv", row.names = FALSE)

#shapiro test for normality for each parameter
shapiro_parameters <- tibble(
  parameter = c("auc", "k", "Cl", "V", "t12"),
  p_value = c(shapiro.test(asc_nca$auc)$p.value, shapiro.test(asc_nca$k)$p.value, shapiro.test(asc_nca$t12)$p.value, shapiro.test(asc_nca$Cl)$p.value, shapiro.test(asc_nca$V)$p.value)
)

asc_summary <- tibble(
  parameter = c("auc", "k", "t12", "Cl", "V"),
  mean = c(mean(asc_nca$auc), mean(asc_nca$k), mean(asc_nca$t12), mean(asc_nca$Cl), mean(asc_nca$V)),
  sd = c(sd(asc_nca$auc), sd(asc_nca$k), sd(asc_nca$t12), sd(asc_nca$Cl), sd(asc_nca$V)),
  median = c(median(asc_nca$auc), median(asc_nca$k), median(asc_nca$t12), median(asc_nca$Cl), median(asc_nca$V)),
  Q1 = c(quantile(asc_nca$auc, 0.25), quantile(asc_nca$k, 0.25), quantile(asc_nca$t12, 0.25), quantile(asc_nca$Cl, 0.25), quantile(asc_nca$V, 0.25)),
  Q3 = c(quantile(asc_nca$auc, 0.75), quantile(asc_nca$k, 0.75), quantile(asc_nca$t12, 0.75), quantile(asc_nca$Cl, 0.75), quantile(asc_nca$V, 0.75)),
  n = c(length(asc_nca$auc), length(asc_nca$k), length(asc_nca$t12), length(asc_nca$Cl), length(asc_nca$V)),
  shapiro_pvalue = c(shapiro.test(asc_nca$auc)$p.value, shapiro.test(asc_nca$k)$p.value, shapiro.test(asc_nca$t12)$p.value, shapiro.test(asc_nca$Cl)$p.value, shapiro.test(asc_nca$V)$p.value)
)
View(asc_summary)
write.csv(asc_summary, "output/non_comp_analysis/asc_nca_summary.csv", row.names = FALSE)
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################

# NAC analysis

# area under the curve for each patient
nac_nca <- tibble(
  px = c("p01", "p03", "p05", "p06", "p08", "p09", "p10", "p11", "p14",
          "p15", "p17", "p18"),
  tac = c("tac1","tac1","tac1","tac1","tac1","tac1","tac2","tac2","tac2","tac2","tac2","tac2"),
  auc = c(0))
nac_nca
View(nac)
# group by patient (px), then calculate auc by patient
for (paciente in unique(nac$px)) {
  # filter by patient
  nac_px <- filter(nac, px == paciente)
  # calculate auc
  nac_nca$auc[nac_nca$px == paciente] <- AUC((nac_px$Tiempo)/60, (nac_px$c_promedio)/1000) # convert to hours and mmol/L
}
# shapiro test for normality
shapiro.test(nac_nca$auc) 
# p-value = 0.4062 > 0.05, therefore data is normal
nac_nca
View(nac)


# calculate elimination constant (k) for each patient using linear regression at times 90-180 min
nac_nca <- mutate(nac_nca, k = 0)
nac_milimolar <- mutate(nac, c_promedio = c_promedio / 1000) # convert to mmol/L
# log-linear concentration vs time plot
View(nac)
nac_log_linear <- tibble(Tiempo = nac_milimolar$Tiempo[(nac_milimolar$tac == "tac1" | nac_milimolar$tac == "tac2") & nac_milimolar$Tiempo >=90], 
                              px = nac_milimolar$px[(nac_milimolar$tac == "tac1" | nac_milimolar$tac == "tac2") & nac_milimolar$Tiempo >= 90],
                              log_c_promedio = log(nac_milimolar$c_promedio[(nac_milimolar$tac == "tac1" | nac_milimolar$tac == "tac2") & nac_milimolar$Tiempo >= 90]))
View(nac_log_linear)
ggplot() +
  geom_point(data = nac_log_linear, aes(x = Tiempo, y = log_c_promedio, col = nac_log_linear$px))+
  geom_line(data = nac_log_linear, aes(x = Tiempo, y = log_c_promedio, col = nac_log_linear$px))+
  theme_bw()
ggplot() +
  geom_point(data = nac, aes(x = Tiempo, y = c_promedio, col = nac$px))+
  geom_line(data = nac, aes(x = Tiempo, y = c_promedio, col = nac$px))+
  theme_bw()

for (paciente in unique(nac_log_linear$px)) {
  # filter by patient
  nac_px <- filter(nac_log_linear, px == paciente)
  # calculate k
  nac_nca$k[nac_nca$px == paciente] <- lm(log_c_promedio ~ Tiempo, data = nac_px) %>%
    summary() %>%
    .$coefficients %>%
    .[2] * (-1)
}
View(nac_nca)
# convert k from 1/min to 1/h
nac_nca <- mutate(nac_nca, k = k * 60)
# shapiro test for normality
shapiro.test(nac_nca$k)


# calculate elimination half-life (t1/2) for each patient
nac_nca <- mutate(nac_nca, t12 = 0)
for (paciente in unique(nac$px)) {
  # filter by patient
  nac_px <- filter(nac, px == paciente)
  # calculate t1/2
  nac_nca$t12[nac_nca$px == paciente] <- log(2) / nac_nca$k[nac_nca$px == paciente]
}
# shapiro test for normality
shapiro.test(nac_nca$t12)

# dosis total en mg tac1
# AA 2475 NAC 2000 DFO 1000
# 2475/176.12 = 14.05 asc
# 2000/163.2 = 12.25 nac
# 1000/560.6 = 1.78 dfo
# dosis total en mg tac2
# AA 2250 NAC 4000 DFO 1600
# 2250/176.12 = 12.77 asc
# 4000/163.2 = 24.51 nac
# 1600/560.6 = 2.85 dfo

# add dosis total to nac_nca
nac_nca <- mutate(nac_nca, dosis_total = 0)
for (i in nac_nca$px) {
  if (i == "p01" | i == "p03" | i == "p05" | i == "p06" | i == "p08" | i == "p09") {
    nac_nca$dosis_total[nac_nca$px == i] <- 12.25
  } else {
    nac_nca$dosis_total[nac_nca$px == i] <- 24.51
  }
}

# calculate Clearance (Cl) for each patient
nac_nca <- mutate(nac_nca, Cl = 0)
for (paciente in unique(nac$px)) {
  # filter by patient
  nac_px <- filter(nac, px == paciente)
  # calculate Cl
  nac_nca$Cl[nac_nca$px == paciente] <- nac_nca$dosis_total[nac_nca$px == paciente] / nac_nca$auc[nac_nca$px == paciente]
}
shapiro.test(nac_nca$Cl)
# calculate Volume of distribution (V) for each patient
nac_nca <- mutate(nac_nca, V = 0)
for (paciente in unique(nac$px)) {
  # filter by patient
  nac_px <- filter(nac, px == paciente)
  # calculate V
  nac_nca$V[nac_nca$px == paciente] <- nac_nca$Cl[nac_nca$px == paciente] / nac_nca$k[nac_nca$px == paciente]
}
shapiro.test(nac_nca$V)
View(summary(nac_nca))
# export nac_nca summary table to csv
write.csv(summary(nac_nca), "output/non_comp_analysis/nac_nca_summary.csv", row.names = FALSE)
write.csv(nac_nca, "output/non_comp_analysis/nac_nca.csv", row.names = FALSE)

#shapiro test for normality for each parameter
(shapiro_parameters <- tibble(
  parameter = c("auc", "k", "Cl", "V", "t12"),
  p_value = c(shapiro.test(nac_nca$auc)$p.value, shapiro.test(nac_nca$k)$p.value, shapiro.test(nac_nca$t12)$p.value, shapiro.test(nac_nca$Cl)$p.value, shapiro.test(nac_nca$V)$p.value)
))

nac_summary <- tibble(
  parameter = c("auc", "k", "t12", "Cl", "V"),
  mean = c(mean(nac_nca$auc), mean(nac_nca$k), mean(nac_nca$t12), mean(nac_nca$Cl), mean(nac_nca$V)),
  sd = c(sd(nac_nca$auc), sd(nac_nca$k), sd(nac_nca$t12), sd(nac_nca$Cl), sd(nac_nca$V)),
  median = c(median(nac_nca$auc), median(nac_nca$k), median(nac_nca$t12), median(nac_nca$Cl), median(nac_nca$V)),
  Q1 = c(quantile(nac_nca$auc, 0.25), quantile(nac_nca$k, 0.25), quantile(nac_nca$t12, 0.25), quantile(nac_nca$Cl, 0.25), quantile(nac_nca$V, 0.25)),
  Q3 = c(quantile(nac_nca$auc, 0.75), quantile(nac_nca$k, 0.75), quantile(nac_nca$t12, 0.75), quantile(nac_nca$Cl, 0.75), quantile(nac_nca$V, 0.75)),
  n = c(length(nac_nca$auc), length(nac_nca$k), length(nac_nca$t12), length(nac_nca$Cl), length(nac_nca$V)),
  shapiro_pvalue = c(shapiro.test(nac_nca$auc)$p.value, shapiro.test(nac_nca$k)$p.value, shapiro.test(nac_nca$t12)$p.value, shapiro.test(nac_nca$Cl)$p.value, shapiro.test(nac_nca$V)$p.value)

)
View(nac_summary)
write.csv(nac_summary, "output/non_comp_analysis/nac_nca_summary.csv", row.names = FALSE)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################


# DFO analysis

# area under the curve for each patient
dfo_nca <- tibble(
  px = c("p01", "p03", "p05", "p06", "p08", "p09", "p10", "p11", "p14",
          "p15", "p17", "p18"),
  tac = c("tac1","tac1","tac1","tac1","tac1","tac1","tac2","tac2","tac2","tac2","tac2","tac2"),
  auc = c(0))
dfo_nca
View(dfo)
# group by patient (px), then calculate auc by patient
for (paciente in unique(dfo$px)) {
  # filter by patient
  dfo_px <- filter(dfo, px == paciente)
  # calculate auc
  dfo_nca$auc[dfo_nca$px == paciente] <- AUC((dfo_px$Tiempo)/60, (dfo_px$c_promedio)/1000, na.rm = TRUE) # convert to hours and mmol/L
}
# shapiro test for normality
shapiro.test(dfo_nca$auc) 


# calculate elimination constant (k) for each patient using linear regression at times 90-180 min
dfo_nca <- mutate(dfo_nca, k = 0)
dfo_milimolar <- mutate(dfo, c_promedio = c_promedio / 1000) # convert to mmol/L
# log-linear concentration vs time plot
View(dfo)
dfo_log_linear <- tibble(Tiempo = dfo_milimolar$Tiempo[(dfo_milimolar$tac == "tac1" | dfo_milimolar$tac == "tac2") & dfo_milimolar$Tiempo >=90], 
                              px = dfo_milimolar$px[(dfo_milimolar$tac == "tac1" | dfo_milimolar$tac == "tac2") & dfo_milimolar$Tiempo >= 90],
                              log_c_promedio = log(dfo_milimolar$c_promedio[(dfo_milimolar$tac == "tac1" | dfo_milimolar$tac == "tac2") & dfo_milimolar$Tiempo >= 90]))
View(dfo_log_linear)
ggplot() +
  geom_point(data = dfo_log_linear, aes(x = Tiempo, y = log_c_promedio, col = dfo_log_linear$px))+
  geom_line(data = dfo_log_linear, aes(x = Tiempo, y = log_c_promedio, col = dfo_log_linear$px))+
  theme_bw()
ggplot() +
  geom_point(data = dfo, aes(x = Tiempo, y = c_promedio, col = dfo$px))+
  geom_line(data = dfo, aes(x = Tiempo, y = c_promedio, col = dfo$px))+
  theme_bw()

for (paciente in unique(dfo_log_linear$px)) {
  # filter by patient
  dfo_px <- filter(dfo_log_linear, px == paciente)
  # calculate k
  dfo_nca$k[dfo_nca$px == paciente] <- lm(log_c_promedio ~ Tiempo, data = dfo_px) %>%
    summary() %>%
    .$coefficients %>%
    .[2] * (-1)
}
# convert k from 1/min to 1/h
dfo_nca <- mutate(dfo_nca, k = k * 60)
# shapiro test for normality
shapiro.test(dfo_nca$k)


# calculate elimination half-life (t1/2) for each patient
dfo_nca <- mutate(dfo_nca, t12 = 0)
for (paciente in unique(dfo$px)) {
  # filter by patient
  dfo_px <- filter(dfo, px == paciente)
  # calculate t1/2
  dfo_nca$t12[dfo_nca$px == paciente] <- log(2) / dfo_nca$k[dfo_nca$px == paciente]
}
# shapiro test for normality
shapiro.test(dfo_nca$t12)

# dosis total en mg tac1
# AA 2475 NAC 2000 DFO 1000
# 2475/176.12 = 14.05 asc
# 2000/163.2 = 12.25 nac
# 1000/560.6 = 1.78 dfo
# dosis total en mg tac2
# AA 2250 NAC 4000 DFO 1600
# 2250/176.12 = 12.77 asc
# 4000/163.2 = 24.51 nac
# 1600/560.6 = 2.85 dfo

# add dosis total to dfo_nca
dfo_nca <- mutate(dfo_nca, dosis_total = 0)
for (i in dfo_nca$px) {
  if (i == "p01" | i == "p03" | i == "p05" | i == "p06" | i == "p08" | i == "p09") {
    dfo_nca$dosis_total[dfo_nca$px == i] <- 1.78
  } else {
    dfo_nca$dosis_total[dfo_nca$px == i] <- 2.85
  }
}

# calculate Clearance (Cl) for each patient
dfo_nca <- mutate(dfo_nca, Cl = 0)
for (paciente in unique(dfo$px)) {
  # filter by patient
  dfo_px <- filter(dfo, px == paciente)
  # calculate Cl
  dfo_nca$Cl[dfo_nca$px == paciente] <- dfo_nca$dosis_total[dfo_nca$px == paciente] / dfo_nca$auc[dfo_nca$px == paciente]
}
shapiro.test(dfo_nca$Cl)
# calculate Volume of distribution (V) for each patient
dfo_nca <- mutate(dfo_nca, V = 0)
for (paciente in unique(dfo$px)) {
  # filter by patient
  dfo_px <- filter(dfo, px == paciente)
  # calculate V
  dfo_nca$V[dfo_nca$px == paciente] <- dfo_nca$Cl[dfo_nca$px == paciente] / dfo_nca$k[dfo_nca$px == paciente]
}
shapiro.test(dfo_nca$V)
View(summary(dfo_nca))

Desc(dfo_nca, stats = c("mean", "sd", "median", "iqr", "n"), verbose = "high")
summarise(dfo_nca)



# export dfo_nca summary table to csv
write.csv(dfo_nca, "output/non_comp_analysis/dfo_nca.csv", row.names = FALSE)

#shapiro test for normality for each parameter
(shapiro_parameters <- tibble(
  parameter = c("auc", "k", "Cl", "V", "t12"),
  p_value = c(shapiro.test(dfo_nca$auc)$p.value, shapiro.test(dfo_nca$k)$p.value, shapiro.test(dfo_nca$t12)$p.value, shapiro.test(dfo_nca$Cl)$p.value, shapiro.test(dfo_nca$V)$p.value)
))

dfo_summary <- tibble(
  parameter = c("auc", "k", "t12", "Cl", "V"),
  mean = c(mean(dfo_nca$auc), mean(dfo_nca$k), mean(dfo_nca$t12), mean(dfo_nca$Cl), mean(dfo_nca$V)),
  sd = c(sd(dfo_nca$auc), sd(dfo_nca$k), sd(dfo_nca$t12), sd(dfo_nca$Cl), sd(dfo_nca$V)),
  median = c(median(dfo_nca$auc), median(dfo_nca$k), median(dfo_nca$t12), median(dfo_nca$Cl), median(dfo_nca$V)),
  Q1 = c(quantile(dfo_nca$auc, 0.25), quantile(dfo_nca$k, 0.25), quantile(dfo_nca$t12, 0.25), quantile(dfo_nca$Cl, 0.25), quantile(dfo_nca$V, 0.25)),
  Q3 = c(quantile(dfo_nca$auc, 0.75), quantile(dfo_nca$k, 0.75), quantile(dfo_nca$t12, 0.75), quantile(dfo_nca$Cl, 0.75), quantile(dfo_nca$V, 0.75)),
  n = c(length(dfo_nca$auc), length(dfo_nca$k), length(dfo_nca$t12), length(dfo_nca$Cl), length(dfo_nca$V)),
  shapiro_pvalue = c(shapiro.test(dfo_nca$auc)$p.value, shapiro.test(dfo_nca$k)$p.value, shapiro.test(dfo_nca$t12)$p.value, shapiro.test(dfo_nca$Cl)$p.value, shapiro.test(dfo_nca$V)$p.value)
)
View(dfo_summary)
write.csv(dfo_nca, "output/non_comp_analysis/dfo_nca.csv", row.names = FALSE)
write.csv(dfo_summary, "output/non_comp_analysis/dfo_nca_summary.csv", row.names = FALSE)
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################



