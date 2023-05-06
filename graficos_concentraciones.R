library(ggplot2)
library(dplyr)
library(tibble)

# Import data from covariables_asc.csv covariables_nac.csv and covariables_dfo.csv
asc <- tibble(read.csv("output/covariables_asc.csv", header = TRUE, sep = ","))
nac <- tibble(read.csv("output/covariables_nac.csv", header = TRUE, sep = ","))
dfo <- tibble(read.csv("output/covariables_dfo.csv", header = TRUE, sep = ","))

# Plot
# Plot mean and standard deviation for ASC c_promedio over Tiempo grouped by tac1, tac2 and p
asc_graph <- ggplot() +
    # add mean and standard deviation
    stat_summary(data = asc, aes(x = Tiempo, y = c_promedio, col = tac), fun.y = mean, geom = "line", size = 1) +
    stat_summary(data = asc, aes(x = Tiempo, y = c_promedio, col = tac, shape = tac), fun.y = mean, geom = "point", size = 3) +
    stat_summary(data = asc, aes(x = Tiempo, y = c_promedio, col = tac), fun.data = mean_se, fun.args=list(mult=1), geom = "errorbar", width = 6) +
    #stat_summary(data = asc, aes(x = Tiempo, y = c_promedio, col = tac), fun.data = mean_se, fun.args=list(mult=sqrt(length(unique(asc$px)))), geom = "errorbar", width = 6) +
    # if you want to add standard deviation adjusted to 95% confidence interval
    # stat_summary(data = asc, aes(x = Tiempo, y = c_promedio, col = tac), fun.data = mean_se, fun.args=list(mult=1.96), geom = "errorbar", width = 4) +
    # set colors
    scale_color_manual(values = c("#5c5c5c",  "#6c0000", "#04176b"), 
                        labels = c("Placebo", "CAT 1", "CAT 2")) +
    scale_shape_manual(values = c(17, 16, 15),
                        labels = c("Placebo", "CAT 1", "CAT 2")) +
    labs(x = "Time [min]", y = "Concentration [μM]", title = "Vitamin C") +
    # change font size
    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 12))
asc_graph


# Plot mean and standard deviation for NAC c_promedio over Tiempo grouped by tac1, tac2 and p
nac_graph <- ggplot() +
    # geom_point(data = nac, aes(x = Tiempo, y = c_promedio, col = tac)) +
    # add mean and standard deviation
    stat_summary(data = nac, aes(x = Tiempo, y = c_promedio, col = tac), fun.y = mean, geom = "line", size = 1) +
    stat_summary(data = nac, aes(x = Tiempo, y = c_promedio, col = tac, shape = tac), fun.y = mean, geom = "point", size = 3) +
    stat_summary(data = nac, aes(x = Tiempo, y = c_promedio, col = tac), fun.data = mean_se, fun.args=list(mult=1), geom = "errorbar", width = 6) +
    # set colors
    scale_color_manual(values = c("#6c0000", "#04176b"), labels = c("Placebo", "CAT 1", "CAT 2")) +
        scale_shape_manual(values = c(16, 15),
                        labels = c("Placebo", "CAT 1", "CAT 2")) +
    labs(x = "Time [min]", y = "Concentration [μM]", title = "N-Acetylcysteine") +
    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 12))
nac_graph

# Plot mean and standard deviation for DFO c_promedio over Tiempo grouped by tac1, tac2 and p
dfo_graph <- ggplot() +
    # geom_point(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac)) +
    # add mean and standard deviation
    stat_summary(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac), fun.y = mean, geom = "line", size = 1) +
    stat_summary(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac, shape = tac), fun.y = mean, geom = "point", size = 3) +
    stat_summary(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac), fun.data = mean_se, fun.args = list(mult=1), geom = "errorbar", width = 6) +
    # set colors
    scale_color_manual(values = c("#6c0000", "#04176b"), labels = c("CAT 1", "CAT 2")) +
    scale_shape_manual(values = c(16, 15),
                        labels = c("CAT 1", "CAT 2")) +
    labs(x = "Time [min]", y = "Concentration [μM]", title = "Deferoxamine") +
    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 12))
dfo_graph

# facet with the three graphs
# install.packages("gridExtra")
library(gridExtra)
all_plots <- grid.arrange(asc_graph, nac_graph, dfo_graph, ncol = 3)
all_plots
#change letters size
all_plots <- all_plots + theme(text = element_text(size=12))
all_plots
# save the grid.arrange plots
ggsave("output/concentraciones.png", all_plots, width = 12.5, height = 4, dpi = 1000)






ggplot()+
  geom_line(data = datos_nac, aes(x=t_min, y=concentration, col = patient),size = 1, alpha = 0.5, linetype = "solid")+
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