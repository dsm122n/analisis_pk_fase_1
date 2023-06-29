library(tibble)

#import data
asc <- tibble(read.csv("output/covariables_asc.csv"))
dfo <- tibble(read.csv("output/covariables_dfo.csv"))
nac <- tibble(read.csv("output/covariables_nac.csv"))

#shapiro wilk test for normality for each time TAC1 and TAC2 together
shapiro_pvalue_asc <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_asc = c(0))
shapiro_pvalue_dfo <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_dfo = c(0))
shapiro_pvalue_nac <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_nac = c(0))
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

#shapiro wilk test for normality for each time, filtered by "tac" column
# tac1
(asc_tac1 <- asc[asc$tac == "tac1" & is.na(asc$c_promedio) == FALSE, ])
(dfo_tac1 <- dfo[dfo$tac == "tac1" & is.na(dfo$c_promedio) == FALSE, ])
(nac_tac1 <- nac[nac$tac == "tac1" & is.na(nac$c_promedio) == FALSE, ])
shapiro_pvalue_asc_tac1 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_asc = c(0))
shapiro_pvalue_dfo_tac1 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_dfo = c(0))
shapiro_pvalue_nac_tac1 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_nac = c(0))


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
shapiro_pvalue_asc_tac2 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_asc = c(0))
shapiro_pvalue_dfo_tac2 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_dfo = c(0))
shapiro_pvalue_nac_tac2 <- tibble(tiempos = c(15, 30, 60, 90, 120), pvalue_nac = c(0))

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
# qqplots for tac1 at each time point
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


#shapiro test with standarized c_promedio over time
#tac1
#asc
asc_tac1_shapiro <- asc[asc$tac == "tac1" & is.na(asc$c_promedio) == FALSE, ]
asc_tac1_shapiro <- as.data.frame(asc_tac1_shapiro[order(asc_tac1_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
asc_tac1_shapiro_tiempo <- c()
for (i in c(
  0, 15, 30, 60, 90, 120, 180
  )) {
    asc_vector <- asc_tac1_shapiro$c_promedio[asc_tac1_shapiro$Tiempo == i]
    asc_tac1_shapiro_tiempo <- c(asc_tac1_shapiro_tiempo, scale(asc_vector))
  
}

hist(asc_tac1_shapiro_tiempo, main = "asc_tac1_shapiro")
#density plot
plot(density(asc_tac1_shapiro_tiempo), main = "asc_tac1_shapiro")
#print shapiro test in previous plot as text
p_value <- shapiro.test(asc_tac1_shapiro_tiempo)$p.value
text(0.5, 10, paste("p-value =", round(p_value, 4)), cex = 1.2)
#histogram
shapiro.test(asc_tac1_shapiro_tiempo)

#dfo
dfo_tac1_shapiro <- dfo[dfo$tac == "tac1" & is.na(dfo$c_promedio) == FALSE, ]
dfo_tac1_shapiro <- as.data.frame(dfo_tac1_shapiro[order(dfo_tac1_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
dfo_tac1_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    dfo_vector <- dfo_tac1_shapiro$c_promedio[dfo_tac1_shapiro$Tiempo == i]
    dfo_tac1_shapiro_tiempo <- c(dfo_tac1_shapiro_tiempo, scale(dfo_vector))
  
}
#density plot
plot(density(dfo_tac1_shapiro_tiempo), main = "dfo_tac1_shapiro")
#histogram
hist(dfo_tac1_shapiro_tiempo, main = "dfo_tac1_shapiro")
shapiro.test(dfo_tac1_shapiro_tiempo)

#nac
nac_tac1_shapiro <- nac[nac$tac == "tac1" & is.na(nac$c_promedio) == FALSE, ]
nac_tac1_shapiro <- as.data.frame(nac_tac1_shapiro[order(nac_tac1_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
nac_tac1_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    nac_vector <- nac_tac1_shapiro$c_promedio[nac_tac1_shapiro$Tiempo == i]
    nac_tac1_shapiro_tiempo <- c(nac_tac1_shapiro_tiempo, scale(nac_vector))
  
}
#density plot
plot(density(nac_tac1_shapiro_tiempo), main = "nac_tac1_shapiro")
#histogram
hist(nac_tac1_shapiro_tiempo, main = "nac_tac1_shapiro")
#shapiro test
shapiro.test(nac_tac1_shapiro_tiempo)


#tac2
#asc
asc_tac2_shapiro <- asc[asc$tac == "tac2" & is.na(asc$c_promedio) == FALSE, ]
asc_tac2_shapiro <- as.data.frame(asc_tac2_shapiro[order(asc_tac2_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
asc_tac2_shapiro_tiempo <- c()
for (i in c(
  0, 15, 30, 60, 90, 120, 180
  )) {
    asc_vector <- asc_tac2_shapiro$c_promedio[asc_tac2_shapiro$Tiempo == i]
    asc_tac2_shapiro_tiempo <- c(asc_tac2_shapiro_tiempo, scale(asc_vector))
  
}
#density plot
plot(density(asc_tac2_shapiro_tiempo), main = "asc_tac2_shapiro")
#histogram
hist(asc_tac2_shapiro_tiempo, main = "asc_tac2_shapiro")
shapiro.test(asc_tac2_shapiro_tiempo)

#dfo
dfo_tac2_shapiro <- dfo[dfo$tac == "tac2" & is.na(dfo$c_promedio) == FALSE, ]
dfo_tac2_shapiro <- as.data.frame(dfo_tac2_shapiro[order(dfo_tac2_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
dfo_tac2_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    dfo_vector <- dfo_tac2_shapiro$c_promedio[dfo_tac2_shapiro$Tiempo == i]
    dfo_tac2_shapiro_tiempo <- c(dfo_tac2_shapiro_tiempo, scale(dfo_vector))
  
}
#density plot
plot(density(dfo_tac2_shapiro_tiempo), main = "dfo_tac2_shapiro")
#histogram
hist(dfo_tac2_shapiro_tiempo, main = "dfo_tac2_shapiro")
shapiro.test(dfo_tac2_shapiro_tiempo)

#nac
nac_tac2_shapiro <- nac[nac$tac == "tac2" & is.na(nac$c_promedio) == FALSE, ]
nac_tac2_shapiro <- as.data.frame(nac_tac2_shapiro[order(nac_tac2_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
nac_tac2_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    nac_vector <- nac_tac2_shapiro$c_promedio[nac_tac2_shapiro$Tiempo == i]
    nac_tac2_shapiro_tiempo <- c(nac_tac2_shapiro_tiempo, scale(nac_vector))
  
}
#density plot
plot(density(nac_tac2_shapiro_tiempo), main = "nac_tac2_shapiro")
#histogram
hist(nac_tac2_shapiro_tiempo, main = "nac_tac2_shapiro")
#shapiro test
shapiro.test(nac_tac2_shapiro_tiempo)


######################################
######################################
######################################
######################################
######################################

# save the following plots in unique pdf
pdf("output/normalidad/density_shapiro_qqplots standarized data.pdf")

# shapiro test for both tac1 and tac2 pooled with standarized c_promedio values
#asc
asc_shapiro <- asc[asc$tac == "tac1" & is.na(asc$c_promedio) == FALSE, ]
asc_shapiro <- as.data.frame(asc_shapiro[order(asc_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
asc_shapiro_tiempo <- c()
for (i in c(
  0, 15, 30, 60, 90, 120, 180
  )) {
    asc_vector <- asc_shapiro$c_promedio[asc_shapiro$Tiempo == i]
    asc_shapiro_tiempo <- c(asc_shapiro_tiempo, scale(asc_vector))
  
}
# add tac2
asc_shapiro2 <- asc[asc$tac == "tac2" & is.na(asc$c_promedio) == FALSE, ]
asc_shapiro2 <- as.data.frame(asc_shapiro2[order(asc_shapiro2$Tiempo), ])
#Vector containing all standarized c_promedio values
asc_shapiro_tiempo2 <- c()

for (i in c(
  0, 15, 30, 60, 90, 120, 180
  )) {
    asc_vector2 <- asc_shapiro2$c_promedio[asc_shapiro2$Tiempo == i]
    asc_shapiro_tiempo2 <- c(asc_shapiro_tiempo2, scale(asc_vector2))
  
}
# join tac1 and tac2
asc_shapiro_tiempo <- c(asc_shapiro_tiempo, asc_shapiro_tiempo2)
(p_value <- shapiro.test(asc_shapiro_tiempo)$p.value)
#density plot
plot(density(asc_shapiro_tiempo), main = "asc density plot standarized_data (tac1 and tac2)")
# add shapiro_wilk p-value to plot
text(-1.2, 0.05, paste("Shapiro-Wilk p-value =", round(p_value, 4)), pos = 4, col = "#000e7b")
#histogram
hist(asc_shapiro_tiempo, main = "asc histogram standarized_data (tac1 and tac2)")
#shapiro test
shapiro.test(asc_shapiro_tiempo)
# qqplot with normal distribution, 45 degree line 
qqnorm(asc_shapiro_tiempo, main = "asc qqplot standarized_data (tac1 and tac2)", col = "#535001")
qqline(asc_shapiro_tiempo, col = "#535001", lwd = 2)
abline(a = 0, b = 1, col = "#000c7b", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000c7b")

# add 45° line named "normal distribution"

#dfo
dfo_shapiro <- dfo[dfo$tac == "tac1" & is.na(dfo$c_promedio) == FALSE, ]
dfo_shapiro <- as.data.frame(dfo_shapiro[order(dfo_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
dfo_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    dfo_vector <- dfo_shapiro$c_promedio[dfo_shapiro$Tiempo == i]
    dfo_shapiro_tiempo <- c(dfo_shapiro_tiempo, scale(dfo_vector))
  
}
# add tac2
dfo_shapiro2 <- dfo[dfo$tac == "tac2" & is.na(dfo$c_promedio) == FALSE, ]
dfo_shapiro2 <- as.data.frame(dfo_shapiro2[order(dfo_shapiro2$Tiempo), ])
#Vector containing all standarized c_promedio values
dfo_shapiro_tiempo2 <- c()

for (i in c(
  15, 30, 60, 90, 120
  )) {
    dfo_vector2 <- dfo_shapiro2$c_promedio[dfo_shapiro2$Tiempo == i]
    dfo_shapiro_tiempo2 <- c(dfo_shapiro_tiempo2, scale(dfo_vector2))
  
}
# join tac1 and tac2
dfo_shapiro_tiempo <- c(dfo_shapiro_tiempo, dfo_shapiro_tiempo2)
(p_value <- shapiro.test(dfo_shapiro_tiempo)$p.value)
#density plot
plot(density(dfo_shapiro_tiempo), main = "dfo density plot standarized_data (tac1 and tac2)")
# add shapiro_wilk p-value to plot
text(-1.2, 0.05, paste("Shapiro-Wilk p-value =", round(shapiro.test(dfo_shapiro_tiempo)$p.value, 4)), pos = 4, col = "#000e7b")

#histogram
hist(dfo_shapiro_tiempo, main = "dfo histogram standarized_data (tac1 and tac2)")
#shapiro test
shapiro.test(dfo_shapiro_tiempo)
# qqplot with normal distribution, 45 degree line 
qqnorm(dfo_shapiro_tiempo, main = "dfo qqplot standarized_data (tac1 and tac2)", col = "#535001")
qqline(dfo_shapiro_tiempo, col = "#535001", lwd = 2)
abline(a = 0, b = 1, col = "#000c7b", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000c7b")

#nac
nac_shapiro <- nac[nac$tac == "tac1" & is.na(nac$c_promedio) == FALSE, ]
nac_shapiro <- as.data.frame(nac_shapiro[order(nac_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
nac_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    nac_vector <- nac_shapiro$c_promedio[nac_shapiro$Tiempo == i]
    nac_shapiro_tiempo <- c(nac_shapiro_tiempo, scale(nac_vector))
  
}
# add tac2
nac_shapiro2 <- nac[nac$tac == "tac2" & is.na(nac$c_promedio) == FALSE, ]
nac_shapiro2 <- as.data.frame(nac_shapiro2[order(nac_shapiro2$Tiempo), ])
#Vector containing all standarized c_promedio values
nac_shapiro_tiempo2 <- c()

for (i in c(
  15, 30, 60, 90, 120
  )) {
    nac_vector2 <- nac_shapiro2$c_promedio[nac_shapiro2$Tiempo == i]
    nac_shapiro_tiempo2 <- c(nac_shapiro_tiempo2, scale(nac_vector2))
  
}
# join tac1 and tac2
nac_shapiro_tiempo <- c(nac_shapiro_tiempo, nac_shapiro_tiempo2)
p_value <- shapiro.test(nac_shapiro_tiempo)$p.value
#density plot
plot(density(nac_shapiro_tiempo), main = "nac density plot standarized_data (tac1 and tac2)")
# add shapiro_wilk p-value to plot
text(-1.2, 0.05, paste("Shapiro-Wilk p-value =", round(shapiro.test(nac_shapiro_tiempo)$p.value, 4)), pos = 4, col = "#000e7b")
#histogram
hist(nac_shapiro_tiempo, main = "nac histogram standarized_data (tac1 and tac2)")
#shapiro test
shapiro.test(nac_shapiro_tiempo)
# qqplot with normal distribution, 45 degree line
qqnorm(nac_shapiro_tiempo, main = "nac qqplot standarized_data (tac1 and tac2)", col = "#535001")
qqline(nac_shapiro_tiempo, col = "#535001", lwd = 2)
abline(a = 0, b = 1, col = "#000c7b", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000c7b")

# finish pdf and save in working directory
dev.off()
########################################################################
########################################################################
########################################################################
########################################################################
# repeat process for each tac

pdf("output/normalidad/tac1 vs tac2 density_shapiro_qqplots standarized data.pdf")

# shapiro test for tac1 standarized c_promedio values pooled
#asc
#tac1
asc_shapiro <- asc[asc$tac == "tac1" & is.na(asc$c_promedio) == FALSE, ]
asc_shapiro <- as.data.frame(asc_shapiro[order(asc_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
asc_shapiro_tiempo <- c()
for (i in c(
  0, 15, 30, 60, 90, 120, 180
  )) {
    asc_vector <- asc_shapiro$c_promedio[asc_shapiro$Tiempo == i]
    asc_shapiro_tiempo <- c(asc_shapiro_tiempo, scale(asc_vector))
  
}

(p_value_tac1_asc <- shapiro.test(asc_shapiro_tiempo)$p.value)

# tac2
asc_shapiro2 <- asc[asc$tac == "tac2" & is.na(asc$c_promedio) == FALSE, ]
asc_shapiro2 <- as.data.frame(asc_shapiro2[order(asc_shapiro2$Tiempo), ])
#Vector containing all standarized c_promedio values
asc_shapiro_tiempo2 <- c()
for (i in c(
  0, 15, 30, 60, 90, 120, 180
  )) {
    asc_vector2 <- asc_shapiro2$c_promedio[asc_shapiro2$Tiempo == i]
    asc_shapiro_tiempo2 <- c(asc_shapiro_tiempo2, scale(asc_vector2))
  
}

(p_value_tac2_asc <- shapiro.test(asc_shapiro_tiempo2)$p.value)

#density plot tac1 vs tac2, set y axis limits
plot(density(asc_shapiro_tiempo), ylim = c(0, 0.6), main = "asc density plot standarized_data (tac1 vs tac2)", , col = "#000e7b")
# add shapiro_wilk p-value to plot
text(-1.2, 0.05, paste("Shapiro-Wilk tac1 p-value =", round(p_value_tac1_asc, 4)), pos = 4, col = "#000e7b")
#density plot tac2
lines(density(asc_shapiro_tiempo2), col = "#620000")
# add shapiro_wilk p-value to plot
text(-1.2, 0.1, paste("Shapiro-Wilk tac2 p-value =", round(p_value_tac2_asc, 4)), pos = 4, col = "#620000")
# increment y axis limits


#histogram tac1 vs tac2
hist(asc_shapiro_tiempo, col = "#000e7b99", ylim = c(0, 16), main = "asc histogram standarized_data (tac1 vs tac2)")
# add tac2 histogram with transparency
hist(asc_shapiro_tiempo2, col = "#62000099", add = TRUE)


# qqplot tac1 with normal distribution, 45 degree line 
qqnorm(asc_shapiro_tiempo, main = "asc qqplot standarized_data (tac1)", col = "#000e7b99")
qqline(asc_shapiro_tiempo, col = "#000e7b", lwd = 2)
abline(a = 0, b = 1, col = "#000000", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000c7b")
#qqplot tac2 with normal distribution, 45 degree line
qqnorm(asc_shapiro_tiempo2, main = "asc qqplot standarized_data (tac2)", col = "#62000099", add = TRUE)
qqline(asc_shapiro_tiempo2, col = "#620000", lwd = 2)
abline(a = 0, b = 1, col = "#000000", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000000")


#dfo
#tac1
dfo_shapiro <- dfo[dfo$tac == "tac1" & is.na(dfo$c_promedio) == FALSE, ]
dfo_shapiro <- as.data.frame(dfo_shapiro[order(dfo_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
dfo_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    dfo_vector <- dfo_shapiro$c_promedio[dfo_shapiro$Tiempo == i]
    dfo_shapiro_tiempo <- c(dfo_shapiro_tiempo, scale(dfo_vector))
  
}
# add tac2
dfo_shapiro2 <- dfo[dfo$tac == "tac2" & is.na(dfo$c_promedio) == FALSE, ]
dfo_shapiro2 <- as.data.frame(dfo_shapiro2[order(dfo_shapiro2$Tiempo), ])
#Vector containing all standarized c_promedio values
dfo_shapiro_tiempo2 <- c()

for (i in c(
  15, 30, 60, 90, 120
  )) {
    dfo_vector2 <- dfo_shapiro2$c_promedio[dfo_shapiro2$Tiempo == i]
    dfo_shapiro_tiempo2 <- c(dfo_shapiro_tiempo2, scale(dfo_vector2))
  
}
# p-value tac1 and tac2
(p_value_tac1_dfo <- shapiro.test(dfo_shapiro_tiempo)$p.value)
(p_value_tac2_dfo <- shapiro.test(dfo_shapiro_tiempo2)$p.value)
#density plot
plot(density(dfo_shapiro_tiempo), main = "dfo density plot standarized_data (tac1 vs tac2)", col = "#000e7b")
# add shapiro_wilk p-value to plot
text(-1.2, 0.05, paste("Shapiro-Wilk tac1 p-value =", round(shapiro.test(dfo_shapiro_tiempo)$p.value, 4)), pos = 4, col = "#000e7b")
#add density plot tac2
lines(density(dfo_shapiro_tiempo2), col = "#620000")
# add shapiro_wilk p-value to plot
text(-1.2, 0.1, paste("Shapiro-Wilk tac2 p-value =", round(shapiro.test(dfo_shapiro_tiempo2)$p.value, 4)), pos = 4, col = "#620000")


#histogram
hist(dfo_shapiro_tiempo, main = "dfo histogram standarized_data (tac1 and tac2)", col = "#000e7b99")
# add tac2 histogram with transparency
hist(dfo_shapiro_tiempo2, col = "#62000099", add = TRUE)
#shapiro test
shapiro.test(dfo_shapiro_tiempo)
# qqplot with normal distribution, 45 degree line 
qqnorm(dfo_shapiro_tiempo, main = "dfo qqplot standarized_data (tac1)", col = "#000e7b99")
qqline(dfo_shapiro_tiempo, col = "#000e7b", lwd = 2)
abline(a = 0, b = 1, col = "#000000", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000000")
#qqplot tac2 with normal distribution, 45 degree line
qqnorm(dfo_shapiro_tiempo2, main = "dfo qqplot standarized_data (tac2)", col = "#62000099", add = TRUE)
qqline(dfo_shapiro_tiempo2, col = "#620000", lwd = 2)
abline(a = 0, b = 1, col = "#000000", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000000")


#nac
#tac1
nac_shapiro <- nac[nac$tac == "tac1" & is.na(nac$c_promedio) == FALSE, ]
nac_shapiro <- as.data.frame(nac_shapiro[order(nac_shapiro$Tiempo), ])
#Vector containing all standarized c_promedio values
nac_shapiro_tiempo <- c()
for (i in c(
  15, 30, 60, 90, 120
  )) {
    nac_vector <- nac_shapiro$c_promedio[nac_shapiro$Tiempo == i]
    nac_shapiro_tiempo <- c(nac_shapiro_tiempo, scale(nac_vector))
  
}
# add tac2
nac_shapiro2 <- nac[nac$tac == "tac2" & is.na(nac$c_promedio) == FALSE, ]
nac_shapiro2 <- as.data.frame(nac_shapiro2[order(nac_shapiro2$Tiempo), ])
#Vector containing all standarized c_promedio values
nac_shapiro_tiempo2 <- c()

for (i in c(
  15, 30, 60, 90, 120
  )) {
    nac_vector2 <- nac_shapiro2$c_promedio[nac_shapiro2$Tiempo == i]
    nac_shapiro_tiempo2 <- c(nac_shapiro_tiempo2, scale(nac_vector2))
  
}
# p-value tac1 and tac2
p_value_tac1_nac <- shapiro.test(nac_shapiro_tiempo)$p.value
p_value_tac2_nac <- shapiro.test(nac_shapiro_tiempo2)$p.value
#density plot tac1 and tac2
plot(density(nac_shapiro_tiempo), ylim = c(0, 0.65), main = "nac density plot standarized_data (tac1 and tac2)", col = "#000e7b")
# add shapiro_wilk p-value to plot
text(-1.2, 0.05, paste("Shapiro-Wilk tac1 p-value =", round(shapiro.test(nac_shapiro_tiempo)$p.value, 4)), pos = 4, col = "#000e7b")
#add density plot tac2
lines(density(nac_shapiro_tiempo2), col = "#620000")
# add shapiro_wilk p-value to plot
text(-1.2, 0.1, paste("Shapiro-Wilk tac2 p-value =", round(shapiro.test(nac_shapiro_tiempo2)$p.value, 4)), pos = 4, col = "#620000")

#histogram
hist(nac_shapiro_tiempo, col ="#000e7b99",  main = "nac histogram standarized_data (tac1 and tac2)")
# add tac2 histogram with transparency
hist(nac_shapiro_tiempo2, col = "#62000099", add = TRUE)
#shapiro test
shapiro.test(nac_shapiro_tiempo)
# qqplot with normal distribution, 45 degree line
qqnorm(nac_shapiro_tiempo, main = "nac qqplot standarized_data (tac1)", col = "#000e7b")
qqline(nac_shapiro_tiempo, col = "#000e7b", lwd = 2)
abline(a = 0, b = 1, col = "#000000", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000000")
#qqplot tac2 with normal distribution, 45 degree line
qqnorm(nac_shapiro_tiempo2, main = "nac qqplot standarized_data (tac2)", col = "#62000099", add = TRUE)
qqline(nac_shapiro_tiempo2, col = "#620000", lwd = 2)
abline(a = 0, b = 1, col = "#000000", lty = 2, lwd = 4)
text(-1, -1.25, "normal distribution", pos = 4, col = "#000000")



# finish pdf and save in working directory
# table with p-values tac1 and tac2 for each treatment for asc, dfo and nac
(p_value_table <- data.frame(
  row.names = c("asc", "dfo", "nac"),
  tac1 = c(p_value_tac1_asc, p_value_tac1_dfo, p_value_tac1_nac),
  tac2 = c(p_value_tac2_asc, p_value_tac2_dfo, p_value_tac2_nac)
))
# add p-value table to new page in pdf

library(gridExtra)
grid.table(p_value_table)

dev.off()
# install gridExtra package to save table in pdf
install.packages("gridExtra")

