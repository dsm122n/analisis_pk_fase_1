library(tidyr)
library(dplyr)
library(ggplot2)

vc <- read.csv("datos_vc.csv", header = TRUE, sep = ",") %>%
        pivot_longer(cols = -c(Tiempo), values_to = "concentraciones", names_to = "paciente") %>%
            mutate(px = strtrim(paciente, 3)) %>%
                mutate(sample = substring(paciente, 5, 5))
vc
nac <- read.csv("datos_nac.csv", header = TRUE, sep = ",") %>%
        pivot_longer(cols = -c(Tiempo), values_to = "concentraciones", names_to = "paciente") %>%
            mutate(px = strtrim(paciente, 3)) %>%
                mutate(sample = substring(paciente, 5, 5))
nac
dfo <- read.csv("datos_dfo.csv", header = TRUE, sep = ",") %>%
        pivot_longer(cols = -c(Tiempo), values_to = "concentraciones", names_to = "paciente") %>%
            mutate(px = strtrim(paciente, 3)) %>%
                mutate(sample = substring(paciente, 5, 5))
dfo


vc_fco <- vc %>%
    mutate(fco = "vc")
nac_fco <- nac %>%  
    mutate(fco = "nac")
dfo_fco <- dfo %>%  
    mutate(fco = "dfo")
todos <- bind_rows(vc_fco, nac_fco, dfo_fco, .id = NULL)
ggplot() +
    stat_summary(data = todos, aes(x = Tiempo, y = concentraciones), fun = "mean", geom = "line", size = 1) +
    # stat_summary(data = todos, aes(x = Tiempo, y = concentraciones), fun = "sd", geom = "errorbar", size = 0.25) +
    geom_point(data = todos, aes(x = Tiempo, y = concentraciones, colour = sample), alpha = 0.3, size = 2) +
    facet_grid(cols = vars(px), rows = vars(fco), scales = "free") +
    theme_bw()
ggsave("todas_muestras_vc_nac_dfo.pdf", width = 70, height = 20, units = "cm")
 


# resumiendo datos

mean_asc <- asc %>%
    group_by(Tiempo, px) %>%
        summarise(c_promedio = mean(concentraciones, na.rm = TRUE))
ggplot()+
    geom_line(data = mean_asc, aes(x = Tiempo, y = c_promedio, colour = px)) +
    scale_x_continuous(breaks = seq(0, 180, 60),
                        minor_breaks = seq(0, 180, 15)) +
    theme_bw()

mean_nac <- nac %>%
    group_by(Tiempo, px) %>%
        summarise(c_promedio = mean(concentraciones, na.rm = TRUE))
ggplot()+
    geom_line(data = mean_nac, aes(x = Tiempo, y = c_promedio, colour = px)) +
    scale_x_continuous(breaks = seq(0, 180, 60),
                        minor_breaks = seq(0, 180, 15)) +
    theme_bw()

mean_dfo <- dfo %>%
    group_by(Tiempo, px) %>%
        summarise(c_promedio = mean(concentraciones, na.rm = TRUE))  
ggplot()+   
    geom_line(data = mean_dfo, aes(x = Tiempo, y = c_promedio, colour = px)) +
    scale_x_continuous(breaks = seq(0, 180, 60),
                        minor_breaks = seq(0, 180, 15)) +
    theme_bw()

covariables <- read.csv("covariables.csv", header = TRUE, sep = ",")

asc_covariables <- mean_asc %>%
    left_join(covariables, by = "px")
nac_covariables <- mean_nac %>%
    left_join(covariables, by = "px")
dfo_covariables <- mean_dfo %>%
    left_join(covariables, by = "px")
write.csv(asc_covariables, "covariables_asc.csv", row.names = FALSE)
write.csv(nac_covariables, "covariables_nac.csv", row.names = FALSE)
write.csv(dfo_covariables, "covariables_dfo.csv", row.names = FALSE)
