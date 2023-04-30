library(tibble)
asc <- tibble(read.csv("output/covariables_asc.csv"))
dfo <- tibble(read.csv("output/covariables_dfo.csv"))
nac <- tibble(read.csv("output/covariables_nac.csv"))
shapiro_pvalue_asc <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_asc = c())
shapiro_pvalue_dfo <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_dfo = c())
shapiro_pvalue_nac <- tibble(tiempos = c(15, 30, 60, 90, 120, 180), pvalue_nac = c())
for (i in c(15, 30, 60, 90, 120, 180)) {
    # print(shapiro.test(asc$c_promedio[asc$Tiempo == i]))
    # print(shapiro.test(dfo$c_promedio[dfo$Tiempo == i]))
    # print(shapiro.test(nac$c_promedio[nac$Tiempo == i]))
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
