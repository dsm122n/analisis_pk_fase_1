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
    scale_x_continuous(breaks = seq(0, 180, 30)) +
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

#same graph but with boxplot instead of mean and standard deviation grouping by tac and time

asc_graph_boxplot <- ggplot() +
    geom_boxplot(data = asc, aes(x = Tiempo, y = c_promedio, group = interaction(Tiempo, tac), col = tac), size = 1) +
    # set colors
    scale_colour_manual(values = c("#5c5c5c",  "#6c0000", "#04176b"), 
                        labels = c("Placebo", "CAT 1", "CAT 2")) +
    labs(x = "Time [min]", y = "Concentration [μM]", title = "Vitamin C") +
    # change font size
    scale_x_continuous(breaks = seq(0, 180, 30)) +
    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 12))
asc_graph_boxplot

# Plot mean and standard deviation for NAC c_promedio over Tiempo grouped by tac1, tac2 and p
nac_graph <- ggplot() +
    # geom_point(data = nac, aes(x = Tiempo, y = c_promedio, col = tac)) +
    # add mean and standard deviation
    stat_summary(data = nac, aes(x = Tiempo, y = c_promedio, col = tac), fun.y = mean, geom = "line", size = 1) +
    stat_summary(data = nac, aes(x = Tiempo, y = c_promedio, col = tac, shape = tac), fun.y = mean, geom = "point", size = 3) +
    stat_summary(data = nac, aes(x = Tiempo, y = c_promedio, col = tac), fun.data = mean_se, fun.args=list(mult=1), geom = "errorbar", width = 6) +
    scale_x_continuous(breaks = seq(0, 180, 30)) +
    # set colors
    scale_color_manual(values = c("#6c0000", "#04176b"), labels = c("CAT 1", "CAT 2")) +
        scale_shape_manual(values = c(16, 15),
                        labels = c("CAT 1", "CAT 2")) +
    labs(x = "Time [min]", y = "Concentration [μM]", title = "N-Acetylcysteine") +
    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 12))
nac_graph

#box plot
nac_graph_boxplot <- ggplot() +
    geom_boxplot(data = nac, aes(x = Tiempo, y = c_promedio, group = interaction(Tiempo, tac), col = tac), size = 1) +
    # set colors
    scale_colour_manual(values = c("#6c0000", "#04176b"), 
                        labels = c("CAT 1", "CAT 2")) +
    labs(x = "Time [min]", y = "Concentration [μM]", title = "N-Acetylcysteine") +
    # change font size
    scale_x_continuous(breaks = seq(0, 180, 30)) +

    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 12))
nac_graph_boxplot

# Plot mean and standard deviation for DFO c_promedio over Tiempo grouped by tac1, tac2 and p
dfo_graph <- ggplot() +
    # geom_point(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac)) +
    # add mean and standard deviation
    stat_summary(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac), fun.y = mean, geom = "line", size = 1) +
    stat_summary(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac, shape = tac), fun.y = mean, geom = "point", size = 3) +
    stat_summary(data = dfo, aes(x = Tiempo, y = c_promedio, col = tac), fun.data = mean_se, fun.args = list(mult=1), geom = "errorbar", width = 6) +
    # set colors
    scale_x_continuous(breaks = seq(0, 180, 30)) +

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

#box plot
dfo_graph_boxplot <- ggplot() +
    geom_boxplot(data = dfo, aes(x = Tiempo, y = c_promedio, group = interaction(Tiempo, tac), col = tac), size = 1) +
    # set colors
    scale_colour_manual(values = c("#6c0000", "#04176b"), 
                        labels = c("CAT 1", "CAT 2")) +
    labs(x = "Time [min]", y = "Concentration [μM]", title = "Deferoxamine") +
    # change font size
    scale_x_continuous(breaks = seq(0, 180, 30)) +

    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size=8),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 12))
dfo_graph_boxplot

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

all_boxplots <- grid.arrange(asc_graph_boxplot, nac_graph_boxplot, dfo_graph_boxplot, ncol = 3)
all_boxplots_vert <- grid.arrange(asc_graph_boxplot, nac_graph_boxplot, dfo_graph_boxplot, nrow = 3)
ggplot(all_boxplots)

# save the grid.arrange plots
ggsave("output/concentraciones_boxplot.png", all_boxplots, width = 12.5, height = 4, dpi = 1000)
ggsave("output/concentraciones_boxplot_vert.png", all_boxplots_vert, width = 7.5, height = 10, dpi = 1000)

