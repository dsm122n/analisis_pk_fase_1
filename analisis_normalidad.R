library(tibble)
asc <- tibble(read.csv("output/covariables_asc.csv"))
dfo <- tibble(read.csv("output/covariables_dfo.csv"))
nac <- tibble(read.csv("output/covariables_nac.csv"))
shapiro_pvalue_asc <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_asc = c())
shapiro_pvalue_dfo <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_dfo = c())
shapiro_pvalue_nac <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_nac = c())
for (i in c(15, 30, 60, 90, 120, 180)) {
    shapiro_pvalue_asc$pvalue_asc[shapiro_pvalue_asc$tiempos == i] <- shapiro.test(na.omit(asc$c_promedio[asc$Tiempo == i]))$p.value
    shapiro_pvalue_dfo$pvalue_dfo[shapiro_pvalue_dfo$tiempos == i] <- shapiro.test(na.omit(dfo$c_promedio[dfo$Tiempo == i]))$p.value
    shapiro_pvalue_nac$pvalue_nac[shapiro_pvalue_nac$tiempos == i] <- shapiro.test(na.omit(nac$c_promedio[nac$Tiempo == i]))$p.value

}
#shapiro_pvalue_asc
todos <- cbind(shapiro_pvalue_asc, shapiro_pvalue_dfo, shapiro_pvalue_nac)
todos <- todos[, -c(3,5)]
#convert to tibble 
todos <- as.tibble(todos)
todos
#delete columns of time
todos
# mark the significant values
todos$asc <- ifelse(todos$pvalue_asc < 0.05, "sig", "no_sig")
todos$dfo <- ifelse(todos$pvalue_dfo < 0.05, "sig", "no_sig")
todos$nac <- ifelse(todos$pvalue_nac < 0.05, "sig", "no_sig")
todos
#reorder the columns
todos <- todos[, c(1, 2, 5, 3, 6, 4, 7)]
todos

#shapiro wilk test for normality for each time, filtered by tac column
(asc_tac1 <- asc[asc$tac == "tac1" & is.na(asc$c_promedio) == FALSE, ])
(dfo_tac1 <- dfo[dfo$tac == "tac1" & is.na(dfo$c_promedio) == FALSE, ])
(nac_tac1 <- nac[nac$tac == "tac1" & is.na(nac$c_promedio) == FALSE, ])
shapiro_pvalue_asc_tac1 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_asc = c())
shapiro_pvalue_dfo_tac1 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_dfo = c())
shapiro_pvalue_nac_tac1 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_nac = c())


for (i in c(15, 30, 60, 90, 120)) {
  if (sum(!is.na(asc_tac1$c_promedio[asc_tac1$Tiempo == i])) >= 3) {
    shapiro_pvalue_asc_tac1$pvalue_asc[shapiro_pvalue_asc_tac1$tiempos == i] <- shapiro.test(na.omit(asc_tac1$c_promedio[asc_tac1$Tiempo == i]))$p.value
  }
  
  if (sum(!is.na(dfo_tac1$c_promedio[dfo_tac1$Tiempo == i])) >= 3) {
    shapiro_pvalue_dfo_tac1$pvalue_dfo[shapiro_pvalue_dfo_tac1$tiempos == i] <- shapiro.test(na.omit(dfo_tac1$c_promedio[dfo_tac1$Tiempo == i]))$p.value
  }
  
  if (sum(!is.na(nac_tac1$c_promedio[nac_tac1$Tiempo == i])) >= 3) {
    shapiro_pvalue_nac_tac1$pvalue_nac[shapiro_pvalue_nac_tac1$tiempos == i] <- shapiro.test(na.omit(nac_tac1$c_promedio[nac_tac1$Tiempo == i]))$p.value
  }
}

todos_tac1 <- cbind(shapiro_pvalue_asc_tac1, shapiro_pvalue_dfo_tac1, shapiro_pvalue_nac_tac1)
todos_tac1 <- todos_tac1[, -c(3,5)]
#convert to tibble 
todos_tac1 <- as.tibble(todos_tac1)
todos_tac1
# mark the significant values
todos_tac1$asc <- ifelse(todos_tac1$pvalue_asc < 0.05, "sig", "no_sig")
todos_tac1$dfo <- ifelse(todos_tac1$pvalue_dfo < 0.05, "sig", "no_sig")
todos_tac1$nac <- ifelse(todos_tac1$pvalue_nac < 0.05, "sig", "no_sig")
todos_tac1
#reorder the columns
todos_tac1 <- todos_tac1[, c(1, 2, 5, 3, 6, 4, 7)]
todos_tac1
# transfer data to .csv file in output/normalidad folder
write.csv(todos_tac1, file = "output/normalidad/todos_tac1.csv", row.names = FALSE)

#shapiro wilk test for tac2
(asc_tac2 <- asc[asc$tac == "tac2" & is.na(asc$c_promedio) == FALSE, ])
(dfo_tac2 <- dfo[dfo$tac == "tac2" & is.na(dfo$c_promedio) == FALSE, ])
(nac_tac2 <- nac[nac$tac == "tac2" & is.na(nac$c_promedio) == FALSE, ])
shapiro_pvalue_asc_tac2 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_asc = c())
shapiro_pvalue_dfo_tac2 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_dfo = c())
shapiro_pvalue_nac_tac2 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_nac = c())

for (i in c(15, 30, 60, 90, 120)) {
  if (sum(!is.na(asc_tac2$c_promedio[asc_tac2$Tiempo == i])) >= 3) {
    shapiro_pvalue_asc_tac2$pvalue_asc[shapiro_pvalue_asc_tac2$tiempos == i] <- shapiro.test(na.omit(asc_tac2$c_promedio[asc_tac2$Tiempo == i]))$p.value
  }
  
  if (sum(!is.na(dfo_tac2$c_promedio[dfo_tac2$Tiempo == i])) >= 3) {
    shapiro_pvalue_dfo_tac2$pvalue_dfo[shapiro_pvalue_dfo_tac2$tiempos == i] <- shapiro.test(na.omit(dfo_tac2$c_promedio[dfo_tac2$Tiempo == i]))$p.value
  }
  
  if (sum(!is.na(nac_tac2$c_promedio[nac_tac2$Tiempo == i])) >= 3) {
    shapiro_pvalue_nac_tac2$pvalue_nac[shapiro_pvalue_nac_tac2$tiempos == i] <- shapiro.test(na.omit(nac_tac2$c_promedio[nac_tac2$Tiempo == i]))$p.value
  }
}

todos_tac2 <- cbind(shapiro_pvalue_asc_tac2, shapiro_pvalue_dfo_tac2, shapiro_pvalue_nac_tac2)
todos_tac2 <- todos_tac2[, -c(3,5)]
#convert to tibble
todos_tac2 <- as.tibble(todos_tac2)
todos_tac2
# mark the significant values
todos_tac2$asc <- ifelse(todos_tac2$pvalue_asc < 0.05, "sig", "no_sig")
todos_tac2$dfo <- ifelse(todos_tac2$pvalue_dfo < 0.05, "sig", "no_sig")
todos_tac2$nac <- ifelse(todos_tac2$pvalue_nac < 0.05, "sig", "no_sig")
todos_tac2
#reorder the columns
todos_tac2 <- todos_tac2[, c(1, 2, 5, 3, 6, 4, 7)]
todos_tac2
# transfer data to .csv file in output/normalidad folder
write.csv(todos_tac2, file = "output/normalidad/todos_tac2.csv", row.names = FALSE)

#qqplot for tac1 at each time point
# save all the plots in one pdf file
pdf("output/normalidad/qqplot_tac1_tac2.pdf")
#asc
asc_tac1_qqplot <- asc[asc$tac == "tac1" & is.na(asc$c_promedio) == FALSE, ]
asc_tac1_qqplot <- as.data.frame(asc_tac1_qqplot[order(asc_tac1_qqplot$Tiempo), ])
for (i in c(15, 30, 60, 90, 120, 180)) {
    asc_tac1_qqplot_tiempo <- asc_tac1_qqplot[asc_tac1_qqplot$Tiempo == i, "c_promedio"]
    qqnorm(asc_tac1_qqplot_tiempo, main = paste("asc_tac1_qqplot", i))
    qqline(asc_tac1_qqplot_tiempo)
  
}

#dfo
dfo_tac1_qqplot <- dfo[dfo$tac == "tac1" & is.na(dfo$c_promedio) == FALSE, ]
dfo_tac1_qqplot <- as.data.frame(dfo_tac1_qqplot[order(dfo_tac1_qqplot$Tiempo), ])
for (i in c(15, 30, 60, 90, 120, 180)) {
    dfo_tac1_qqplot_tiempo <- dfo_tac1_qqplot[dfo_tac1_qqplot$Tiempo == i, "c_promedio"]
    qqnorm(dfo_tac1_qqplot_tiempo, main = paste("dfo_tac1_qqplot", i))
    qqline(dfo_tac1_qqplot_tiempo)
  
}

#nac
nac_tac1_qqplot <- nac[nac$tac == "tac1" & is.na(nac$c_promedio) == FALSE, ]
nac_tac1_qqplot <- as.data.frame(nac_tac1_qqplot[order(nac_tac1_qqplot$Tiempo), ])
for (i in c(15, 30, 60, 90, 120, 180)) {
    nac_tac1_qqplot_tiempo <- nac_tac1_qqplot[nac_tac1_qqplot$Tiempo == i, "c_promedio"]
    qqnorm(nac_tac1_qqplot_tiempo, main = paste("nac_tac1_qqplot", i))
    qqline(nac_tac1_qqplot_tiempo)
  
}

#qqplot for tac2 at each time point
#asc
asc_tac2_qqplot <- asc[asc$tac == "tac2" & is.na(asc$c_promedio) == FALSE, ]
asc_tac2_qqplot <- as.data.frame(asc_tac2_qqplot[order(asc_tac2_qqplot$Tiempo), ])
for (i in c(15, 30, 60, 90, 120, 180)) {
    asc_tac2_qqplot_tiempo <- asc_tac2_qqplot[asc_tac2_qqplot$Tiempo == i, "c_promedio"]
    qqnorm(asc_tac2_qqplot_tiempo, main = paste("asc_tac2_qqplot", i))
    qqline(asc_tac2_qqplot_tiempo)
  
}

#dfo
dfo_tac2_qqplot <- dfo[dfo$tac == "tac2" & is.na(dfo$c_promedio) == FALSE, ]
dfo_tac2_qqplot <- as.data.frame(dfo_tac2_qqplot[order(dfo_tac2_qqplot$Tiempo), ])
for (i in c(15, 30, 60, 90, 120, 180)) {
    dfo_tac2_qqplot_tiempo <- dfo_tac2_qqplot[dfo_tac2_qqplot$Tiempo == i, "c_promedio"]
    qqnorm(dfo_tac2_qqplot_tiempo, main = paste("dfo_tac2_qqplot", i))
    qqline(dfo_tac2_qqplot_tiempo)
  
}

#nac
nac_tac2_qqplot <- nac[nac$tac == "tac2" & is.na(nac$c_promedio) == FALSE, ]
nac_tac2_qqplot <- as.data.frame(nac_tac2_qqplot[order(nac_tac2_qqplot$Tiempo), ])
for (i in c(15, 30, 60, 90, 120, 180)) {
    nac_tac2_qqplot_tiempo <- nac_tac2_qqplot[nac_tac2_qqplot$Tiempo == i, "c_promedio"]
    qqnorm(nac_tac2_qqplot_tiempo, main = paste("nac_tac2_qqplot", i))
    qqline(nac_tac2_qqplot_tiempo)
  
}

#all times pooled
#asc
asc_tac1_qqplot <- asc[asc$tac == "tac1" & is.na(asc$c_promedio) == FALSE, ]
asc_tac1_qqplot <- as.data.frame(asc_tac1_qqplot[order(asc_tac1_qqplot$Tiempo), ])
asc_tac1_qqplot_tiempo <- asc_tac1_qqplot[, "c_promedio"]
qqnorm(asc_tac1_qqplot_tiempo, main = "asc_tac1_qqplot")
qqline(asc_tac1_qqplot_tiempo)
#dfo
dfo_tac1_qqplot <- dfo[dfo$tac == "tac1" & is.na(dfo$c_promedio) == FALSE, ]
dfo_tac1_qqplot <- as.data.frame(dfo_tac1_qqplot[order(dfo_tac1_qqplot$Tiempo), ])
dfo_tac1_qqplot_tiempo <- dfo_tac1_qqplot[, "c_promedio"]
qqnorm(dfo_tac1_qqplot_tiempo, main = "dfo_tac1_qqplot")
qqline(dfo_tac1_qqplot_tiempo)
#nac
nac_tac1_qqplot <- nac[nac$tac == "tac1" & is.na(nac$c_promedio) == FALSE, ]
nac_tac1_qqplot <- as.data.frame(nac_tac1_qqplot[order(nac_tac1_qqplot$Tiempo), ])
nac_tac1_qqplot_tiempo <- nac_tac1_qqplot[, "c_promedio"]
qqnorm(nac_tac1_qqplot_tiempo, main = "nac_tac1_qqplot")
qqline(nac_tac1_qqplot_tiempo)

#tac2
#asc
asc_tac2_qqplot <- asc[asc$tac == "tac2" & is.na(asc$c_promedio) == FALSE, ]
asc_tac2_qqplot <- as.data.frame(asc_tac2_qqplot[order(asc_tac2_qqplot$Tiempo), ])
asc_tac2_qqplot_tiempo <- asc_tac2_qqplot[, "c_promedio"]
qqnorm(asc_tac2_qqplot_tiempo, main = "asc_tac2_qqplot")
qqline(asc_tac2_qqplot_tiempo)
#dfo
dfo_tac2_qqplot <- dfo[dfo$tac == "tac2" & is.na(dfo$c_promedio) == FALSE, ]
dfo_tac2_qqplot <- as.data.frame(dfo_tac2_qqplot[order(dfo_tac2_qqplot$Tiempo), ])
dfo_tac2_qqplot_tiempo <- dfo_tac2_qqplot[, "c_promedio"]
qqnorm(dfo_tac2_qqplot_tiempo, main = "dfo_tac2_qqplot")
qqline(dfo_tac2_qqplot_tiempo)
#nac
nac_tac2_qqplot <- nac[nac$tac == "tac2" & is.na(nac$c_promedio) == FALSE, ]
nac_tac2_qqplot <- as.data.frame(nac_tac2_qqplot[order(nac_tac2_qqplot$Tiempo), ])
nac_tac2_qqplot_tiempo <- nac_tac2_qqplot[, "c_promedio"]
qqnorm(nac_tac2_qqplot_tiempo, main = "nac_tac2_qqplot")
qqline(nac_tac2_qqplot_tiempo)


dev.off()