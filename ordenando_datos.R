library(tidyr)
library(dplyr)

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
library(ggplot2)
ggplot()+
    stat_summary(data = vc, aes(x = Tiempo, y = concentraciones), fun = "mean", geom = "line", size = 0.5) +
    geom_point(data = vc, aes(x = Tiempo, y = concentraciones, colour = sample), alpha = 0.3, size = 2) +
    facet_wrap(vars(px)) +
    theme_bw()
ggplot()+
    stat_summary(data = nac, aes(x = Tiempo, y = concentraciones), fun = "mean", geom = "line", size = 0.5) +
    geom_point(data = nac, aes(x = Tiempo, y = concentraciones, colour = sample), alpha = 0.3, size = 2) +
    facet_wrap(vars(px)) +
    theme_bw()
ggplot()+
    stat_summary(data = dfo, aes(x = Tiempo, y = concentraciones), fun = "mean", geom = "line", size = 0.5) +
    geom_point(data = dfo, aes(x = Tiempo, y = concentraciones, colour = sample), alpha = 0.3, size = 2) +
    facet_wrap(vars(px)) +
    theme_bw()

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
 