


library(stringr)
library(rebus)
library(dplyr)
library(limma)
library(VennDiagram)
library(gplots)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(plotly)
library(tidyr)
library(corrgram)


# install.packages("extrafont")
# 
# library(extrafont)
# font_import()



library(ggExtra)
###################
#########################
#############################


#### Report files

plasma_PP <- read.csv(".../Report_file_R_plasma_PP.csv", stringsAsFactors = F)
cell_PP <- read.csv(".../Report_file_R_cell_PP.csv", stringsAsFactors = F)
plasma_cell_PP <- read.csv(".../Report_file_R_plasma_cell_PP.csv", stringsAsFactors = F)

plasma_TPP <- read.csv(".../Report_file_R_plasma_TPP.csv", stringsAsFactors = F)
cell_TPP <- read.csv(".../Report_file_R_cell_TPP.csv", stringsAsFactors = F)
plasma_cell_TPP <- read.csv(".../Report_file_R_plasma_cell_TPP.csv", stringsAsFactors = F)

plasma_all <- read.csv(".../Report_file_R_plasma_all.csv", stringsAsFactors = F)
cell_all <- read.csv(".../Report_file_R_cell_all.csv", stringsAsFactors = F)
plasma_cell_all <- read.csv(".../Report_file_R_plasma_cell_all.csv", stringsAsFactors = F)

plasma_all <- plasma_all %>%
  select(-ncol(plasma_all))
cell_all <- cell_all %>%
  select(-ncol(cell_all))
plasma_cell_all <- plasma_cell_all %>%
  select(-ncol(plasma_cell_all))


### Library files

lib_plasma_pp <- read.csv(".../Plasma_Library_PP.csv", stringsAsFactors = F)
lib_plasma_tpp <- read.csv(".../Plasma_Library_TPP.csv", stringsAsFactors = F)

lib_cell_pp <- read.csv(".../Cell_Calibrated_Library_PP.csv", stringsAsFactors = F)
lib_cell_tpp <- read.csv(".../Cell_Calibrated_Library_TPP.csv", stringsAsFactors = F)

lib_plasma_cell_pp <- read.csv(".../Plasma_Cell_Library_PP.csv", stringsAsFactors = F)
lib_plasma_cell_tpp <- read.csv(".../Plasma_Cell_Library_TPP.csv", stringsAsFactors = F)


lib_plasma_pp <- lib_plasma_pp[!duplicated(lib_plasma_pp$ModificationSequence),]
lib_plasma_tpp <- lib_plasma_tpp[!duplicated(lib_plasma_tpp$ModificationSequence),]


lib_cell_pp <- lib_cell_pp[!duplicated(lib_cell_pp$ModificationSequence),]
lib_cell_tpp <- lib_cell_tpp[!duplicated(lib_cell_tpp$ModificationSequence),]


lib_plasma_cell_pp <- lib_plasma_cell_pp[!duplicated(lib_plasma_cell_pp$ModificationSequence),]
lib_plasma_cell_tpp <- lib_plasma_cell_tpp[!duplicated(lib_plasma_cell_tpp$ModificationSequence),]


lib_plasma_pp$ModificationSequence <- str_replace_all(lib_plasma_pp$ModificationSequence, pattern = OPEN_BRACKET %R%
                                                        "CAM" %R%
                                                        CLOSE_BRACKET, 
                                                      replacement = OPEN_BRACKET %R%
                                                        "+57" %R%
                                                        CLOSE_BRACKET)

lib_plasma_tpp$ModificationSequence <- str_replace_all(lib_plasma_tpp$ModificationSequence, pattern = OPEN_BRACKET %R%
                                                         "Oxi" %R%
                                                         CLOSE_BRACKET, replacement = OPEN_BRACKET %R%
                                                         "+16" %R%
                                                         CLOSE_BRACKET)


lib_cell_pp$ModificationSequence <- str_replace_all(lib_cell_pp$ModificationSequence, pattern = OPEN_BRACKET %R%
                                                      "CAM" %R%
                                                      CLOSE_BRACKET, 
                                                    replacement = OPEN_BRACKET %R%
                                                      "+57" %R%
                                                      CLOSE_BRACKET)

lib_cell_tpp$ModificationSequence <- str_replace_all(lib_cell_tpp$ModificationSequence, pattern = OPEN_BRACKET %R%
                                                       "Oxi" %R%
                                                       CLOSE_BRACKET, replacement = OPEN_BRACKET %R%
                                                       "+16" %R%
                                                       CLOSE_BRACKET)


lib_plasma_cell_pp$ModificationSequence <- str_replace_all(lib_plasma_cell_pp$ModificationSequence, pattern = OPEN_BRACKET %R%
                                                             "CAM" %R%
                                                             CLOSE_BRACKET, 
                                                           replacement = OPEN_BRACKET %R%
                                                             "+57" %R%
                                                             CLOSE_BRACKET)

lib_plasma_cell_tpp$ModificationSequence <- str_replace_all(lib_plasma_cell_tpp$ModificationSequence, pattern = OPEN_BRACKET %R%
                                                              "Oxi" %R%
                                                              CLOSE_BRACKET, replacement = OPEN_BRACKET %R%
                                                              "+16" %R%
                                                              CLOSE_BRACKET)


### FREQUENCY DISTRIBUTION PLOTS

##labels = c("TPP", "PP", "TPP+PP"), values=c("#FFA500", "#ff3232", "#107dac")

plasma_PP2 <- plasma_PP[!duplicated(plasma_PP$Modified.Sequence),]
plasma_ppx <- aggregate(x = plasma_PP2$Modified.Sequence, by = list(Protein = plasma_PP2$Protein.Name), FUN = function(x) {length(x)})

p1 <- ggplot(data = plasma_ppx, aes(x = x, fill = )) +
  geom_histogram(binwidth = 1, fill = "#ff3232") +
  labs(title = "Plasma_PP", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#ff3232"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_TPP2 <- plasma_TPP[!duplicated(plasma_TPP$Modified.Sequence),]
plasma_tppx <- aggregate(x = plasma_TPP2$Modified.Sequence, by = list(Protein = plasma_TPP2$Protein.Name), FUN = function(x) {length(x)})

p2 <- ggplot(data = plasma_tppx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#FFA500") +
  labs(title = "Plasma_TPP", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#FFA500"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

cell_PP2 <- cell_PP[!duplicated(cell_PP$Modified.Sequence),]
cell_ppx <- aggregate(x = cell_PP2$Modified.Sequence, by = list(Protein = cell_PP2$Protein.Name), FUN = function(x) {length(x)})

p3 <- ggplot(data = cell_ppx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#ff3232") +
  labs(title = "Cell_PP", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#ff3232"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


cell_TPP2 <- cell_TPP[!duplicated(cell_TPP$Modified.Sequence),]
cell_tppx <- aggregate(x = cell_TPP2$Modified.Sequence, by = list(Protein = cell_TPP2$Protein.Name), FUN = function(x) {length(x)})

p4 <- ggplot(data = cell_tppx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#FFA500") +
  labs(title = "Cell_TPP", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#FFA500"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_PP2 <- plasma_cell_PP[!duplicated(plasma_cell_PP$Modified.Sequence),]
plasma_cell_ppx <- aggregate(x = plasma_cell_PP2$Modified.Sequence, by = list(Protein = plasma_cell_PP2$Protein.Name), FUN = function(x) {length(x)})

p5 <- ggplot(data = plasma_cell_ppx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#ff3232") +
  labs(title = "Plasma_Cell_PP", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#ff3232"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_TPP2 <- plasma_cell_TPP[!duplicated(plasma_cell_TPP$Modified.Sequence),]
plasma_cell_tppx <- aggregate(x = plasma_cell_TPP2$Modified.Sequence, by = list(Protein = plasma_cell_TPP2$Protein.Name), FUN = function(x) {length(x)})

p6 <- ggplot(data = plasma_cell_tppx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#FFA500") +
  labs(title = "Plasma_Cell_TPP", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#FFA500"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_all2 <- plasma_all[!duplicated(plasma_all$Modified.Sequence),]
plasma_allx <- aggregate(x = plasma_all2$Modified.Sequence, by = list(Protein = plasma_all2$Protein.Name), FUN = function(x) {length(x)})

p7 <- ggplot(data = plasma_allx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#107dac") +
  labs(title = "Plasma_All", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#107dac"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

cell_all2 <- cell_all[!duplicated(cell_all$Modified.Sequence),]
cell_allx <- aggregate(x = cell_all2$Modified.Sequence, by = list(Protein = cell_all2$Protein.Name), FUN = function(x) {length(x)})

p8 <- ggplot(data = cell_allx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#107dac") +
  labs(title = "Cell_All", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#107dac"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_all2 <- plasma_cell_all[!duplicated(plasma_cell_all$Modified.Sequence),]
plasma_cell_allx <- aggregate(x = plasma_cell_all2$Modified.Sequence, by = list(Protein = plasma_cell_all2$Protein.Name), FUN = function(x) {length(x)})

p9 <- ggplot(data = plasma_cell_allx, aes(x = x)) +
  geom_histogram(binwidth = 1, fill = "#107dac") +
  labs(title = "Plasma_Cell_All", x ="No. of Peptides", y = "Protein Count") +
  # scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values=c("#107dac"))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/frequency_distribution_DIA.tiff", units = "in", width = 10, height = 6, res = 300)
fdpp <- ggarrange(p1, p3, p5, ncol = 3, nrow = 1, legend = "right")
fdtpp <- ggarrange(p2, p4, p6, ncol = 3, nrow = 1, legend = "right")
fdall <- ggarrange(p7, p8, p9, ncol = 3, nrow = 1, legend = "right")
ggarrange(fdtpp, fdpp, fdall, ncol = 1, nrow = 3, legend = "right", labels = "AUTO")
# ggarrange(fdp1, fdp2, ncol = 1, nrow = 2, legend = "right", labels = "AUTO")

dev.off()


################################################################
##################################################################


####  DOTP VS MASS ERROR GRAPH


plasma_PP <- plasma_PP %>%
  mutate(Average.Library.Dot.Product = apply(m <- (plasma_PP[, str_which(colnames(plasma_PP), pattern = "Library.Dot.Product")]), 1, median))

plasma_PP <- plasma_PP %>%
  mutate(Average.Mass.Error = apply(m <- (plasma_PP[, str_which(colnames(plasma_PP), pattern = "Average.Mass.Error.PPM")]), 1, median))

plasma_PP <- plasma_PP %>%
  mutate("Library" = "aPlasma", "Type" = "bPP")

plasma_TPP <- plasma_TPP %>%
  mutate(Average.Library.Dot.Product = apply(m <- (plasma_TPP[, str_which(colnames(plasma_TPP), pattern = "Library.Dot.Product")]), 1, median))

plasma_TPP <- plasma_TPP %>%
  mutate(Average.Mass.Error = apply(m <- (plasma_TPP[, str_which(colnames(plasma_TPP), pattern = "Average.Mass.Error.PPM")]), 1, median))

plasma_TPP <- plasma_TPP %>%
  mutate("Library" = "aPlasma", "Type" = "bTPP")

plasma_all <- plasma_all %>%
  mutate(Average.Library.Dot.Product = apply(m <- (plasma_all[, str_which(colnames(plasma_all), pattern = "Library.Dot.Product")]), 1, median))

plasma_all <- plasma_all %>%
  mutate(Average.Mass.Error = apply(m <- (plasma_all[, str_which(colnames(plasma_all), pattern = "Average.Mass.Error.PPM")]), 1, median))

plasma_all <- plasma_all %>%
  mutate("Library" = "aPlasma", "Type" = "cAll")
##############

cell_PP <- cell_PP %>%
  mutate(Average.Library.Dot.Product = apply(m <- (cell_PP[, str_which(colnames(cell_PP), pattern = "Library.Dot.Product")]), 1, median))

cell_PP <- cell_PP %>%
  mutate(Average.Mass.Error = apply(m <- (cell_PP[, str_which(colnames(cell_PP), pattern = "Average.Mass.Error.PPM")]), 1, median))

cell_PP <- cell_PP %>%
  mutate("Library" = "bCell", "Type" = "bPP")

cell_TPP <- cell_TPP %>%
  mutate(Average.Library.Dot.Product = apply(m <- (cell_TPP[, str_which(colnames(cell_TPP), pattern = "Library.Dot.Product")]), 1, median))

cell_TPP <- cell_TPP %>%
  mutate(Average.Mass.Error = apply(m <- (cell_TPP[, str_which(colnames(cell_TPP), pattern = "Average.Mass.Error.PPM")]), 1, median))

cell_TPP <- cell_TPP %>%
  mutate("Library" = "bCell", "Type" = "bTPP")

cell_all <- cell_all %>%
  mutate(Average.Library.Dot.Product = apply(m <- (cell_all[, str_which(colnames(cell_all), pattern = "Library.Dot.Product")]), 1, median))

cell_all <- cell_all %>%
  mutate(Average.Mass.Error = apply(m <- (cell_all[, str_which(colnames(cell_all), pattern = "Average.Mass.Error.PPM")]), 1, median))

cell_all <- cell_all %>%
  mutate("Library" = "bCell", "Type" = "cAll")

###########

plasma_cell_PP <- plasma_cell_PP %>%
  mutate(Average.Library.Dot.Product = apply(m <- (plasma_cell_PP[, str_which(colnames(plasma_cell_PP), pattern = "Library.Dot.Product")]), 1, median))

plasma_cell_PP <- plasma_cell_PP %>%
  mutate(Average.Mass.Error = apply(m <- (plasma_cell_PP[, str_which(colnames(plasma_cell_PP), pattern = "Average.Mass.Error.PPM")]), 1, median))

plasma_cell_PP <- plasma_cell_PP %>%
  mutate("Library" = "cPlasma_Cell", "Type" = "bPP")

plasma_cell_TPP <- plasma_cell_TPP %>%
  mutate(Average.Library.Dot.Product = apply(m <- (plasma_cell_TPP[, str_which(colnames(plasma_cell_TPP), pattern = "Library.Dot.Product")]), 1, median))

plasma_cell_TPP <- plasma_cell_TPP %>%
  mutate(Average.Mass.Error = apply(m <- (plasma_cell_TPP[, str_which(colnames(plasma_cell_TPP), pattern = "Average.Mass.Error.PPM")]), 1, median))

plasma_cell_TPP <- plasma_cell_TPP %>%
  mutate("Library" = "cPlasma_Cell", "Type" = "bTPP")

plasma_cell_all <- plasma_cell_all %>%
  mutate(Average.Library.Dot.Product = apply(m <- (plasma_cell_all[, str_which(colnames(plasma_cell_all), pattern = "Library.Dot.Product")]), 1, median))

plasma_cell_all <- plasma_cell_all %>%
  mutate(Average.Mass.Error = apply(m <- (plasma_cell_all[, str_which(colnames(plasma_cell_all), pattern = "Average.Mass.Error.PPM")]), 1, median))

plasma_cell_all <- plasma_cell_all %>%
  mutate("Library" = "cPlasma_Cell", "Type" = "cAll")


dotp_mass <- rbind(rbind(rbind(rbind(rbind(rbind(rbind(rbind(plasma_TPP, cell_TPP), plasma_cell_TPP), plasma_PP), cell_PP), plasma_cell_PP),
                               plasma_all), cell_all), plasma_cell_all)

##labels = c("TPP", "PP", "TPP+PP"), values=c("#FFA500", "#ff3232", "#107dac")


dprod <- ggplot(data = dotp_mass, aes(y = Average.Library.Dot.Product, x= Library, fill = Type)) +
  geom_boxplot(alpha = 1,
               notch = TRUE, notchwidth = 0.8,
               outlier.colour = "red", outlier.fill = "red", outlier.size = 1)+
  # scale_fill_manual(values = c("#ff4945", "#75a3e7", "#75a3e7", "#2ac940")) +
  scale_fill_manual(labels = c("TPP", "PP", "TPP+PP"), values=c("#FFA500", "#ff3232", "#107dac")) +
  scale_x_discrete(labels= c("Plasma", "Cell", "Plasma_Cell")) +
  # geom_boxplot(width = 0.1, fill = "white") +
  # stat_summary(fun.data = "mean_sdl", mult = 1, geom = "crossbar", width = 0.2) +
  labs(x = "Libraries", y = "Dot Product") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "none")

merr <- ggplot(data = dotp_mass, aes(y = Average.Mass.Error, x= Library, fill = Type)) +
  geom_boxplot(alpha = 1,
               notch = TRUE, notchwidth = 0.8,
               outlier.colour = "red", outlier.fill = "red", outlier.size = 1)+
  # scale_fill_manual(values = c("#ff4945", "#75a3e7", "#75a3e7", "#2ac940")) +
  scale_fill_manual(labels = c("TPP", "PP", "TPP+PP"), values=c("#FFA500", "#ff3232", "#107dac")) +
  scale_x_discrete(labels= c("Plasma", "Cell", "Plasma_Cell")) +
  # geom_boxplot(width = 0.1, fill = "white") +
  # stat_summary(fun.data = "mean_sdl", mult = 1, geom = "crossbar", width = 0.2) +
  labs(x = "Libraries", y = "Mass Error (PPM)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "top")

tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/masserror_dotp.tiff", units = "in", width = 8, height = 5, res = 300)

dp_ms <- ggarrange(dprod, merr, common.legend = TRUE, ncol = 2, nrow = 1, legend = "top", labels = "AUTO")
dp_ms

dev.off()

#######################################################################################
####################################################################################################

#### CV plots for stagewise data
##### Intensity plots for stagewise data


plasma_cell_PP_AA <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_pp_A.csv")
plasma_cell_PP_BB <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_pp_B.csv")
plasma_cell_PP_CC <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_pp_C.csv")
plasma_cell_PP_DD <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_pp_D.csv")
plasma_cell_PP_EE <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_pp_E.csv")


plasma_cell_TPP_AA <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_tpp_A.csv")
plasma_cell_TPP_BB <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_tpp_B.csv")
plasma_cell_TPP_CC <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_tpp_C.csv")
plasma_cell_TPP_DD <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_tpp_D.csv")
plasma_cell_TPP_EE <- read.csv("E:/MQ-Charlie-Data_2018Nov20/All_libraries/ReportFiles_20190408/Report_File_MQ_plasma_cell_tpp_E.csv")

#### DotP Filter

plasma_cell_PP_AA <- plasma_cell_PP_AA %>%
  filter(stageA.1.Library.Dot.Product >= 0.70) %>%
  filter(stageA.2.Library.Dot.Product >= 0.70) %>%
  filter(stageA.3.Library.Dot.Product >= 0.70)
plasma_cell_PP_BB <- plasma_cell_PP_BB %>%
  filter(stageB.1.Library.Dot.Product >= 0.70) %>%
  filter(stageB.2.Library.Dot.Product >= 0.70) %>%
  filter(stageB.3.Library.Dot.Product >= 0.70)
plasma_cell_PP_CC <- plasma_cell_PP_CC %>%
  filter(stageC.1.Library.Dot.Product >= 0.70) %>%
  filter(stageC.2.Library.Dot.Product >= 0.70) %>%
  filter(stageC.3.Library.Dot.Product >= 0.70)
plasma_cell_PP_DD <- plasma_cell_PP_DD %>%
  filter(stageD.1.Library.Dot.Product >= 0.70) %>%
  filter(stageD.2.Library.Dot.Product >= 0.70) %>%
  filter(stageD.3.Library.Dot.Product >= 0.70)
plasma_cell_PP_EE <- plasma_cell_PP_EE %>%
  filter(stageE.1.Library.Dot.Product >= 0.70) %>%
  filter(stageE.2.Library.Dot.Product >= 0.70) %>%
  filter(stageE.3.Library.Dot.Product >= 0.70)

### Missed Cleavage filter

plasma_cell_PP_AA <- plasma_cell_PP_AA %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_PP_BB <- plasma_cell_PP_BB %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_PP_CC <- plasma_cell_PP_CC %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_PP_DD <- plasma_cell_PP_DD %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_PP_EE <- plasma_cell_PP_EE %>%
  filter(Missed.Cleavages <= 2)

### Precursor Charges filter

plasma_cell_PP_AA <- plasma_cell_PP_AA %>%
  filter(Precursor.Charge <= 4)
plasma_cell_PP_BB <- plasma_cell_PP_BB %>%
  filter(Precursor.Charge <= 4)
plasma_cell_PP_CC <- plasma_cell_PP_CC %>%
  filter(Precursor.Charge <= 4)
plasma_cell_PP_DD <- plasma_cell_PP_DD %>%
  filter(Precursor.Charge <= 4)
plasma_cell_PP_EE <- plasma_cell_PP_EE %>%
  filter(Precursor.Charge <= 4)

#### DotP Filter

plasma_cell_TPP_AA <- plasma_cell_TPP_AA %>%
  filter(stageA.1.Library.Dot.Product >= 0.70) %>%
  filter(stageA.2.Library.Dot.Product >= 0.70) %>%
  filter(stageA.3.Library.Dot.Product >= 0.70)
plasma_cell_TPP_BB <- plasma_cell_TPP_BB %>%
  filter(stageB.1.Library.Dot.Product >= 0.70) %>%
  filter(stageB.2.Library.Dot.Product >= 0.70) %>%
  filter(stageB.3.Library.Dot.Product >= 0.70)
plasma_cell_TPP_CC <- plasma_cell_TPP_CC %>%
  filter(stageC.1.Library.Dot.Product >= 0.70) %>%
  filter(stageC.2.Library.Dot.Product >= 0.70) %>%
  filter(stageC.3.Library.Dot.Product >= 0.70)
plasma_cell_TPP_DD <- plasma_cell_TPP_DD %>%
  filter(stageD.1.Library.Dot.Product >= 0.70) %>%
  filter(stageD.2.Library.Dot.Product >= 0.70) %>%
  filter(stageD.3.Library.Dot.Product >= 0.70)
plasma_cell_TPP_EE <- plasma_cell_TPP_EE %>%
  filter(stageE.1.Library.Dot.Product >= 0.70) %>%
  filter(stageE.2.Library.Dot.Product >= 0.70) %>%
  filter(stageE.3.Library.Dot.Product >= 0.70)

### Missed Cleavage filter

plasma_cell_TPP_AA <- plasma_cell_TPP_AA %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_TPP_BB <- plasma_cell_TPP_BB %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_TPP_CC <- plasma_cell_TPP_CC %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_TPP_DD <- plasma_cell_TPP_DD %>%
  filter(Missed.Cleavages <= 2)
plasma_cell_TPP_EE <- plasma_cell_TPP_EE %>%
  filter(Missed.Cleavages <= 2)

### Precursor Charges filter

plasma_cell_TPP_AA <- plasma_cell_TPP_AA %>%
  filter(Precursor.Charge <= 4)
plasma_cell_TPP_BB <- plasma_cell_TPP_BB %>%
  filter(Precursor.Charge <= 4)
plasma_cell_TPP_CC <- plasma_cell_TPP_CC %>%
  filter(Precursor.Charge <= 4)
plasma_cell_TPP_DD <- plasma_cell_TPP_DD %>%
  filter(Precursor.Charge <= 4)
plasma_cell_TPP_EE <- plasma_cell_TPP_EE %>%
  filter(Precursor.Charge <= 4)

plasma_cell_PP_AA <- plasma_cell_PP_AA[!duplicated(plasma_cell_PP_AA$Modified.Sequence),]
plasma_cell_TPP_AA <- plasma_cell_TPP_AA[!duplicated(plasma_cell_TPP_AA$Modified.Sequence),]

plasma_cell_PP_BB <- plasma_cell_PP_BB[!duplicated(plasma_cell_PP_BB$Modified.Sequence),]
plasma_cell_TPP_BB <- plasma_cell_TPP_BB[!duplicated(plasma_cell_TPP_BB$Modified.Sequence),]

plasma_cell_PP_CC <- plasma_cell_PP_CC[!duplicated(plasma_cell_PP_CC$Modified.Sequence),]
plasma_cell_TPP_CC <- plasma_cell_TPP_CC[!duplicated(plasma_cell_TPP_CC$Modified.Sequence),]

plasma_cell_PP_DD <- plasma_cell_PP_DD[!duplicated(plasma_cell_PP_DD$Modified.Sequence),]
plasma_cell_TPP_DD <- plasma_cell_TPP_DD[!duplicated(plasma_cell_TPP_DD$Modified.Sequence),]

plasma_cell_PP_EE <- plasma_cell_PP_EE[!duplicated(plasma_cell_PP_EE$Modified.Sequence),]
plasma_cell_TPP_EE <- plasma_cell_TPP_EE[!duplicated(plasma_cell_TPP_EE$Modified.Sequence),]


plasma_cell_PP_AA <- mutate(plasma_cell_PP_AA, "Stage" = "Stage A", "Type" = "bPP")
plasma_cell_TPP_AA <- mutate(plasma_cell_TPP_AA, "Stage" = "Stage A", "Type" = "aTPP")

plasma_cell_PP_BB <- mutate(plasma_cell_PP_BB, "Stage" = "Stage B", "Type" = "bPP")
plasma_cell_TPP_BB <- mutate(plasma_cell_TPP_BB, "Stage" = "Stage B", "Type" = "aTPP")

plasma_cell_PP_CC <- mutate(plasma_cell_PP_CC, "Stage" = "Stage C", "Type" = "bPP")
plasma_cell_TPP_CC <- mutate(plasma_cell_TPP_CC, "Stage" = "Stage C", "Type" = "aTPP")

plasma_cell_PP_DD <- mutate(plasma_cell_PP_DD, "Stage" = "Stage D", "Type" = "bPP")
plasma_cell_TPP_DD <- mutate(plasma_cell_TPP_DD, "Stage" = "Stage D", "Type" = "aTPP")

plasma_cell_PP_EE <- mutate(plasma_cell_PP_EE, "Stage" = "Stage E", "Type" = "bPP")
plasma_cell_TPP_EE <- mutate(plasma_cell_TPP_EE, "Stage" = "Stage E", "Type" = "aTPP")


column_names <- colnames(plasma_cell_PP_AA)
column_names <- str_replace_all(column_names, "stageA", "Rep")

colnames(plasma_cell_PP_AA) <- column_names
colnames(plasma_cell_PP_BB) <- column_names
colnames(plasma_cell_PP_CC) <- column_names
colnames(plasma_cell_PP_DD) <- column_names
colnames(plasma_cell_PP_EE) <- column_names

column_names <- colnames(plasma_cell_TPP_AA)
column_names <- str_replace_all(column_names, "stageA", "Rep")

colnames(plasma_cell_TPP_AA) <- column_names
colnames(plasma_cell_TPP_BB) <- column_names
colnames(plasma_cell_TPP_CC) <- column_names
colnames(plasma_cell_TPP_DD) <- column_names
colnames(plasma_cell_TPP_EE) <- column_names


plasma_cell_PP <- rbind(rbind(rbind(rbind(plasma_cell_PP_AA, plasma_cell_PP_BB), plasma_cell_PP_CC), plasma_cell_PP_DD), plasma_cell_PP_EE)
plasma_cell_TPP <- rbind(rbind(rbind(rbind(plasma_cell_TPP_AA, plasma_cell_TPP_BB), plasma_cell_TPP_CC), plasma_cell_TPP_DD), plasma_cell_TPP_EE)


plasma_cell_PP[, str_which(colnames(plasma_cell_PP), pattern = "Normalized.Area")] <- apply(m <- (plasma_cell_PP[, str_which(colnames(plasma_cell_PP), pattern = "Normalized.Area")]), 2, as.numeric)
plasma_cell_PP <- plasma_cell_PP %>%
  mutate("Average.Intensity" = apply(m <- (plasma_cell_PP[, str_which(colnames(plasma_cell_PP), pattern = "Normalized.Area")]), 1, mean))
plasma_cell_PP$Average.Intensity <- log(x = plasma_cell_PP$Average.Intensity)


plasma_cell_TPP[, str_which(colnames(plasma_cell_TPP), pattern = "Normalized.Area")] <- apply(m <- (plasma_cell_TPP[, str_which(colnames(plasma_cell_TPP), pattern = "Normalized.Area")]), 2, as.numeric)
plasma_cell_TPP <- plasma_cell_TPP %>%
  mutate("Average.Intensity" = apply(m <- (plasma_cell_TPP[, str_which(colnames(plasma_cell_TPP), pattern = "Normalized.Area")]), 1, mean))
plasma_cell_TPP$Average.Intensity <- log(x = plasma_cell_TPP$Average.Intensity)


plasma_cell_PP <- mutate(plasma_cell_PP, "CV.Normalized" = as.character(plasma_cell_PP$Cv.Total.Area.Normalized))
plasma_cell_PP$CV.Normalized <- str_replace(plasma_cell_PP$CV.Normalized, pattern = "%", replacement = "")
plasma_cell_PP$CV.Normalized <- round(as.numeric(plasma_cell_PP$CV.Normalized))

plasma_cell_TPP <- mutate(plasma_cell_TPP, "CV.Normalized" = as.character(plasma_cell_TPP$Cv.Total.Area.Normalized))
plasma_cell_TPP$CV.Normalized <- str_replace(plasma_cell_TPP$CV.Normalized, pattern = "%", replacement = "")
plasma_cell_TPP$CV.Normalized <- round(as.numeric(plasma_cell_TPP$CV.Normalized))

cvdat <- rbind(plasma_cell_TPP, plasma_cell_PP)
intdat <- rbind(plasma_cell_TPP, plasma_cell_PP)


tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/CV_plot_stages.tiff", units = "in", width = 9, height = 4, res = 300)
cv <- ggplot(data = cvdat, aes(y = CV.Normalized, x= Stage, fill = Type)) +
  geom_violin(alpha = 0.9, color = "black",
              trim = F)+
  # scale_fill_manual(values = c("#ff4945", "#75a3e7", "#75a3e7", "#2ac940")) +
  scale_fill_manual(labels = c("TPP", "PP"), values=c("#FFA500", "#ff3232")) +
  # scale_x_discrete(labels = c("Cattle", "Sheep", "Giraffe_C", "Giraffe_S")) +
  # geom_boxplot(width = 0.1, fill = "white") +
  # stat_summary(fun.data = "mean_sdl", mult = 1, geom = "crossbar", width = 0.2) +
  labs(x = "Replicates", y = "Normalized CVs (%)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")
cv
dev.off()


tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/Intensity_plot_stages.tiff", units = "in", width = 9, height = 4, res = 300)
int <- ggplot(data = intdat, aes(y = Average.Intensity, x= Stage, fill = Type)) +
  geom_boxplot(alpha = 1,
               notch = TRUE, notchwidth = 0.8,
               outlier.colour = "red", outlier.fill = "red", outlier.size = 1)+
  # scale_fill_manual(values = c("#ff4945", "#75a3e7", "#75a3e7", "#2ac940")) +
  scale_fill_manual(labels = c("TPP", "PP"), values=c("#FFA500", "#ff3232")) +
  # scale_x_discrete(labels = c("Cattle", "Sheep", "Giraffe_C", "Giraffe_S")) +
  # geom_boxplot(width = 0.1, fill = "white") +
  # stat_summary(fun.data = "mean_sdl", mult = 1, geom = "crossbar", width = 0.2) +
  labs(x = "Replicates", y = "Log Intensity") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")
int
dev.off()

tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/CV_intensity_plot_stages.tiff", units = "in", width = 7, height = 5.5, res = 300)

ggarrange(cv, int, ncol = 1, nrow= 2, common.legend = TRUE, legend = "top", labels = "AUTO", vjust = 1)

dev.off()


##################################################################################
########################################################

####################################################
############################

#### DDA DIA RT GRAPH PLUS DIA PP AND TPP GRAPH

plasma_pp_rt <- plasma_PP[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]
plasma_tpp_rt <- plasma_TPP[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]

cell_pp_rt <- cell_PP[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]
cell_tpp_rt <- cell_TPP[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]

plasma_cell_pp_rt <- plasma_cell_PP[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]
plasma_cell_tpp_rt <- plasma_cell_TPP[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]

plasma_all_rt <- plasma_all[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]
cell_all_rt <- cell_all[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]
plasma_cell_all_rt <- plasma_cell_all[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]

plasma_pp_rt <- mutate(plasma_pp_rt, "Library" = "aPlasma", "Type" = "bPP")
plasma_tpp_rt <- mutate(plasma_tpp_rt, "Library" = "aPlasma", "Type" = "aTPP")

cell_pp_rt <- mutate(cell_pp_rt, "Library" = "bCell", "Type" = "bPP")
cell_tpp_rt <- mutate(cell_pp_rt, "Library" = "bCell", "Type" = "aTPP")

plasma_cell_pp_rt <- mutate(plasma_cell_pp_rt, "Library" = "cPlasma-Cell", "Type" = "bPP")
plasma_cell_tpp_rt <- mutate(plasma_cell_pp_rt, "Library" = "cPlasma-Cell", "Type" = "aTPP")

plasma_all_rt <- mutate(plasma_all_rt, "Library" = "aPlasma", "Type" = "cAll")
cell_all_rt <- mutate(cell_all_rt, "Library" = "bCell", "Type" = "cAll")
plasma_cell_all_rt <- mutate(plasma_cell_all_rt, "Library" = "cPlasma-Cell", "Type" = "cAll")


comm_pep_lib_pp <- lib_plasma_cell_pp[!duplicated(lib_plasma_cell_pp$ModificationSequence), c("ModificationSequence", "Tr_recalibrated")]
comm_pep_lib1 <- intersect(plasma_cell_pp_rt$Modified.Sequence, comm_pep_lib_pp$ModificationSequence)
comm_pep_lib_pp <- comm_pep_lib_pp[comm_pep_lib_pp$ModificationSequence %in% comm_pep_lib1,]
plasma_cell_pp_rt <- plasma_cell_pp_rt[plasma_cell_pp_rt$Modified.Sequence %in% comm_pep_lib1,]
mergeRT_plasma_cell_pp <- merge(plasma_cell_pp_rt, comm_pep_lib_pp, by.x = 3, by.y = 1)

comm_pep_lib_tpp <- lib_plasma_cell_tpp[!duplicated(lib_plasma_cell_tpp$ModificationSequence), c("ModificationSequence", "Tr_recalibrated")]
comm_pep_lib2 <- intersect(plasma_cell_tpp_rt$Modified.Sequence, comm_pep_lib_tpp$ModificationSequence)
comm_pep_lib_tpp <- comm_pep_lib_tpp[comm_pep_lib_tpp$ModificationSequence %in% comm_pep_lib2,]
plasma_cell_tpp_rt <- plasma_cell_tpp_rt[plasma_cell_tpp_rt$Modified.Sequence %in% comm_pep_lib2,]
mergeRT_plasma_cell_tpp <- merge(plasma_cell_tpp_rt, comm_pep_lib_tpp, by.x = 3, by.y = 1)
##

comm_pep_lib_pp <- lib_plasma_pp[!duplicated(lib_plasma_pp$ModificationSequence), c("ModificationSequence", "Tr_recalibrated")]
comm_pep_lib1 <- intersect(plasma_pp_rt$Modified.Sequence, comm_pep_lib_pp$ModificationSequence)
comm_pep_lib_pp <- comm_pep_lib_pp[comm_pep_lib_pp$ModificationSequence %in% comm_pep_lib1,]
plasma_pp_rt <- plasma_pp_rt[plasma_pp_rt$Modified.Sequence %in% comm_pep_lib1,]
mergeRT_plasma_pp <- merge(plasma_pp_rt, comm_pep_lib_pp, by.x = 3, by.y = 1)

comm_pep_lib_tpp <- lib_plasma_tpp[!duplicated(lib_plasma_tpp$ModificationSequence), c("ModificationSequence", "Tr_recalibrated")]
comm_pep_lib2 <- intersect(plasma_tpp_rt$Modified.Sequence, comm_pep_lib_tpp$ModificationSequence)
comm_pep_lib_tpp <- comm_pep_lib_tpp[comm_pep_lib_tpp$ModificationSequence %in% comm_pep_lib2,]
plasma_tpp_rt <- plasma_tpp_rt[plasma_tpp_rt$Modified.Sequence %in% comm_pep_lib2,]
mergeRT_plasma_tpp <- merge(plasma_tpp_rt, comm_pep_lib_tpp, by.x = 3, by.y = 1)
###

comm_pep_lib_pp <- lib_cell_pp[!duplicated(lib_cell_pp$ModificationSequence), c("ModificationSequence", "Tr_recalibrated")]
comm_pep_lib1 <- intersect(cell_pp_rt$Modified.Sequence, comm_pep_lib_pp$ModificationSequence)
comm_pep_lib_pp <- comm_pep_lib_pp[comm_pep_lib_pp$ModificationSequence %in% comm_pep_lib1,]
cell_pp_rt <- cell_pp_rt[cell_pp_rt$Modified.Sequence %in% comm_pep_lib1,]
mergeRT_cell_pp <- merge(cell_pp_rt, comm_pep_lib_pp, by.x = 3, by.y = 1)

comm_pep_lib_tpp <- lib_cell_tpp[!duplicated(lib_plasma_cell_tpp$ModificationSequence), c("ModificationSequence", "Tr_recalibrated")]
comm_pep_lib2 <- intersect(cell_tpp_rt$Modified.Sequence, comm_pep_lib_tpp$ModificationSequence)
comm_pep_lib_tpp <- comm_pep_lib_tpp[comm_pep_lib_tpp$ModificationSequence %in% comm_pep_lib2,]
cell_tpp_rt <- cell_tpp_rt[cell_tpp_rt$Modified.Sequence %in% comm_pep_lib2,]
mergeRT_cell_tpp <- merge(cell_tpp_rt, comm_pep_lib_tpp, by.x = 3, by.y = 1)
##

mergeRT_cell_pp <- mergeRT_cell_pp[!duplicated(mergeRT_cell_pp$Modified.Sequence),]
mergeRT_cell_tpp <- mergeRT_cell_tpp[!duplicated(mergeRT_cell_tpp$Modified.Sequence),]

mergeRT_plasma_pp <- mergeRT_plasma_pp[!duplicated(mergeRT_plasma_pp$Modified.Sequence),]
mergeRT_plasma_tpp <- mergeRT_plasma_tpp[!duplicated(mergeRT_plasma_tpp$Modified.Sequence),]

mergeRT_plasma_cell_pp <- mergeRT_plasma_cell_pp[!duplicated(mergeRT_plasma_cell_pp$Modified.Sequence),]
mergeRT_plasma_cell_tpp <- mergeRT_plasma_cell_tpp[!duplicated(mergeRT_plasma_cell_tpp$Modified.Sequence),]


rtdat_pp <- rbind(rbind(mergeRT_plasma_pp, mergeRT_cell_pp),mergeRT_plasma_cell_pp)

rtdat_tpp <- rbind(rbind(mergeRT_plasma_tpp, mergeRT_cell_tpp),mergeRT_plasma_cell_tpp)


# scale_fill_manual(labels = c("TPP", "PP"), values=c("#FFA500", "#ff3232"))
### Put Correlation on graph  = 0.94

rt1 <- ggplot(mergeRT_cell_pp, aes(x = Average.Measured.Retention.Time, y = Tr_recalibrated)) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 68, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(mergeRT_cell_pp$Average.Measured.Retention.Time, mergeRT_cell_pp$Tr_recalibrated), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  # scale_color_manual(labels = c("PP"), values = c( "#ff3232")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("PP"), values = c("#ff3232"), aesthetics = ("fill"))+
  labs(title = "Cell_PP", x = "RT from Library (min)", y = " DIA RT (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))


rt2 <- ggplot(mergeRT_plasma_pp, aes(x = Average.Measured.Retention.Time, y = Tr_recalibrated)) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 75, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(mergeRT_plasma_pp$Average.Measured.Retention.Time, mergeRT_plasma_pp$Tr_recalibrated), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  scale_color_manual(labels = c("PP"), values = c( "#ff3232")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("PP"), values = c("#ff3232"), aesthetics = ("fill"))+
  labs(title = "Plasma_PP", x = "RT from Library (min)", y = " DIA RT (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 4)))

rt3 <- ggplot(mergeRT_plasma_cell_pp, aes(x = Average.Measured.Retention.Time, y = Tr_recalibrated)) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 75, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(mergeRT_plasma_cell_pp$Average.Measured.Retention.Time, mergeRT_plasma_cell_pp$Tr_recalibrated), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  scale_color_manual(labels = c("PP"), values = c("#ff3232")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("PP"), values = c("#ff3232"), aesthetics = ("fill"))+
  labs(title = "Plasma_Cell_PP", x = "RT from Library (min)", y = " DIA RT (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))

rt4 <- ggplot(mergeRT_cell_tpp, aes(x = Average.Measured.Retention.Time, y = Tr_recalibrated)) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 75, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(mergeRT_cell_tpp$Average.Measured.Retention.Time, mergeRT_cell_tpp$Tr_recalibrated), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  scale_color_manual(labels = c("TPP"), values = c("#FFA500")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("TPP"), values = c("#FFA500"), aesthetics = ("fill"))+
  labs(title = "Cell_TPP", x = "RT from Library (min)", y = " DIA RT (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))


rt5 <- ggplot(mergeRT_plasma_tpp, aes(x = Average.Measured.Retention.Time, y = Tr_recalibrated)) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 75, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(mergeRT_plasma_tpp$Average.Measured.Retention.Time, mergeRT_plasma_tpp$Tr_recalibrated), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  scale_color_manual(labels = c("TPP"), values = c("#FFA500")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("TPP"), values = c("#FFA500"), aesthetics = ("fill"))+
  labs(title = "Plasma_TPP", x = "RT from Library (min)", y = " DIA RT (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))

rt6 <- ggplot(mergeRT_plasma_cell_tpp, aes(x = Average.Measured.Retention.Time, y = Tr_recalibrated)) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 75, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(mergeRT_plasma_cell_tpp$Average.Measured.Retention.Time, mergeRT_plasma_cell_tpp$Tr_recalibrated), digits = 3)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  # scale_shape_manual(labels = c("Library"), values = 24) +
  scale_color_manual(labels = c("TPP"), values = c("#FFA500")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("TPP"), values = c("#FFA500"), aesthetics = ("fill"))+
  labs(title = "Plasma_Cell_TPP", x = "RT from Library (min)", y = " DIA RT (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))


rtp1 <- ggarrange(rt2, rt1, rt3, ncol = 3, nrow = 1, common.legend = T, legend = "right")
rtp2 <- ggarrange(rt5, rt4, rt6,  ncol = 3, nrow = 1, common.legend = T, legend = "right")




cell_comm_pep <- intersect(cell_PP$Peptide, cell_TPP$Peptide)
plasma_comm_pep <- intersect(plasma_PP$Peptide, plasma_TPP$Peptide)
plasma_cell_comm_pep <- intersect(plasma_cell_PP$Peptide, plasma_cell_TPP$Peptide)


cell_pp_comm <- cell_PP[cell_PP$Peptide %in% cell_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]

cell_tpp_comm <- cell_TPP[cell_TPP$Peptide %in% cell_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]

plasma_pp_comm <- plasma_PP[plasma_PP$Peptide %in% plasma_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]

plasma_tpp_comm <- plasma_TPP[plasma_TPP$Peptide %in% plasma_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]

plasma_cell_pp_comm <- plasma_cell_PP[plasma_cell_PP$Peptide %in% plasma_cell_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]

plasma_cell_tpp_comm <- plasma_cell_TPP[plasma_cell_TPP$Peptide %in% plasma_cell_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]

plasma_cor <- merge(plasma_pp_comm, plasma_tpp_comm, by = "Peptide", all = TRUE)
cell_cor <- merge(cell_pp_comm, cell_tpp_comm, by = "Peptide", all = TRUE)
plasma_cell_cor <- merge(plasma_cell_pp_comm, plasma_cell_tpp_comm, by = "Peptide", all= TRUE)

plasma_cor <- plasma_cor[!duplicated(plasma_cor$Peptide),]
cell_cor <- cell_cor[!duplicated(cell_cor$Peptide),]
plasma_cell_cor <- plasma_cell_cor[!duplicated(plasma_cell_cor$Peptide),]


cell_cor <- cell_cor %>%
  filter(Average.Measured.Retention.Time.x != "")
plasma_cor <- plasma_cor %>%
  filter(Average.Measured.Retention.Time.x != "")
plasma_cell_cor <- plasma_cell_cor %>%
  filter(Average.Measured.Retention.Time.x != "")

cell_cor <- cell_cor %>%
  filter(Average.Measured.Retention.Time.y != "")
plasma_cor <- plasma_cor %>%
  filter(Average.Measured.Retention.Time.y != "")
plasma_cell_cor <- plasma_cell_cor %>%
  filter(Average.Measured.Retention.Time.y != "")

rtcor1 <- ggplot(plasma_cor, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x = 25, y = 46.5, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cor$Average.Measured.Retention.Time.x, plasma_cor$Average.Measured.Retention.Time.y), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  # scale_x_continuous(limits = c(min(plasma_cor$Average.Measured.Retention.Time.x), max(plasma_cor$Average.Measured.Retention.Time.x)),
  #                    expand = c(0,0)) +
  # scale_y_continuous(limits = c(min(plasma_cor$Average.Measured.Retention.Time.y), max(plasma_cor$Average.Measured.Retention.Time.y)),
  #                    expand = c(0,0)) +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  # scale_color_manual(labels = c("TPP"), values = c("#FFA500")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("TPP"), values = c("#FFA500"), aesthetics = ("fill"))+
  labs(title = "Plasma", x = "DIA RT PP (min)", y = "DIA RT TPP (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


rtcor2 <- ggplot(cell_cor, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 55, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(cell_cor$Average.Measured.Retention.Time.x, cell_cor$Average.Measured.Retention.Time.y), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  # scale_x_continuous(limits = c(min(cell_cor$Average.Measured.Retention.Time.x), max(cell_cor$Average.Measured.Retention.Time.x)),
  #                    expand = c(0,0)) +
  # scale_y_continuous(limits = c(min(cell_cor$Average.Measured.Retention.Time.y), max(cell_cor$Average.Measured.Retention.Time.y)),
  #                    expand = c(0,0)) +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  # scale_color_manual(labels = c("TPP"), values = c("#FFA500")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("TPP"), values = c("#FFA500"), aesthetics = ("fill"))+
  labs(title = "Cell", x = "DIA RT PP (min)", y = "DIA RT TPP (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

rtcor3 <- ggplot(plasma_cell_cor, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x = 30, y = 55, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_cor$Average.Measured.Retention.Time.x, plasma_cell_cor$Average.Measured.Retention.Time.y), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  # scale_x_continuous(limits = c(min(plasma_cell_cor$Average.Measured.Retention.Time.x), max(plasma_cell_cor$Average.Measured.Retention.Time.x)),
  #                    expand = c(0,0)) +
  # scale_y_continuous(limits = c(min(plasma_cell_cor$Average.Measured.Retention.Time.y), max(plasma_cell_cor$Average.Measured.Retention.Time.y)),
  #                    expand = c(0,0)) +
  # scale_color_manual(labels = c("TPP"), values = c("#FFA500")) +
  # scale_color_manual(labels = c("TPP", "PP"), values = c("#FFA500","#ff3232")) +
  # scale_fill_manual(labels = c("TPP"), values = c("#FFA500"), aesthetics = ("fill"))+
  labs(title = "Plasma_Cell", x = "DIA RT PP (min)", y = "DIA RT TPP (min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


rt_cor <- ggarrange(rtcor1, rtcor2, rtcor3, ncol = 3, nrow = 1)



# plasma_all_rt <- mutate(plasma_all_rt, "Library" = "aPlasma", "Type" = "cAll")
# cell_all_rt <- mutate(cell_all_rt, "Library" = "bCell", "Type" = "cAll")
# plasma_cell_all_rt <- mutate(plasma_cell_all_rt, "Library" = "cPlasma-Cell", "Type" = "cAll")

pc_cell_all_comm_pep <- intersect(cell_all$Modified.Sequence, plasma_cell_all$Modified.Sequence)
pc_plasma_all_comm_pep <- intersect(plasma_all$Modified.Sequence, plasma_cell_all$Modified.Sequence)



cell_all_comm <- cell_all[cell_all$Modified.Sequence %in% pc_cell_all_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]
pc_cell_all_comm <- plasma_cell_all[plasma_cell_all$Modified.Sequence %in% pc_cell_all_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]

plasma_all_comm <- plasma_all[plasma_all$Modified.Sequence %in% pc_plasma_all_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]
pc_plasma_all_comm <- plasma_cell_all[plasma_cell_all$Modified.Sequence %in% pc_plasma_all_comm_pep, c("Protein.Name", "Peptide", "Modified.Sequence", "Average.Measured.Retention.Time")]


pc_plasma_all_cor <- merge(plasma_all_comm, pc_plasma_all_comm, by = "Peptide", all = TRUE)
pc_cell_all_cor <- merge(cell_all_comm, pc_cell_all_comm, by = "Peptide", all = TRUE)


pc_plasma_all_cor <- pc_plasma_all_cor[!duplicated(pc_plasma_all_cor$Peptide),]
pc_cell_all_cor <- pc_cell_all_cor[!duplicated(pc_cell_all_cor$Peptide),]



pc_plasma_all_cor <- pc_plasma_all_cor %>%
  filter(Average.Measured.Retention.Time.x != "")
pc_cell_all_cor <- pc_cell_all_cor %>%
  filter(Average.Measured.Retention.Time.x != "")

pc_plasma_all_cor <- pc_plasma_all_cor %>%
  filter(Average.Measured.Retention.Time.y != "")
pc_cell_all_cor <- pc_cell_all_cor %>%
  filter(Average.Measured.Retention.Time.y != "")

rt_pc_p <- ggplot(pc_plasma_all_cor, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x = 23, y = 60, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(pc_plasma_all_cor$Average.Measured.Retention.Time.x, pc_plasma_all_cor$Average.Measured.Retention.Time.y), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "Plasma - Plasma_Cell", x = "Plasma DIA RT(min)", y = "P_C DIA RT(min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

rt_pc_c <- ggplot(pc_cell_all_cor, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x = 23, y = 60, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(pc_cell_all_cor$Average.Measured.Retention.Time.x, pc_cell_all_cor$Average.Measured.Retention.Time.y), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  scale_x_continuous(breaks=seq(10, 100, 10)) +
  scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "Cell - Plasma_Cell", x = "Cell DIA RT(min)", y = "P_C DIA RT(min)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

rt_pc <- ggarrange(rt_pc_p, rt_pc_c, ncol = 2, nrow=1)

tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/RT_Correlations.tiff", units = "in", width = 9, height = 9, res = 300)

ggarrange(rt_cor, rtp1, rtp2, rt_pc, ncol = 1, nrow = 4, labels = "AUTO")

dev.off()


##################################################################################
########################################################

#### Intensity Correlation among technical replicates for all(PP+TPP) libraries


cell_all2 <- cell_all

cell_all2[, str_which(colnames(cell_all), pattern = "Normalized.Area")] <- apply(cell_all[, str_which(colnames(cell_all), pattern = "Normalized.Area")], 2, as.numeric)
cell_all2[, str_which(colnames(cell_all2), pattern = "Normalized.Area")] <- apply(cell_all2[, str_which(colnames(cell_all2), pattern = "Normalized.Area")], 2, log2)


plasma_all2 <- plasma_all

plasma_all2[, str_which(colnames(plasma_all), pattern = "Normalized.Area")] <- apply(plasma_all[, str_which(colnames(plasma_all), pattern = "Normalized.Area")], 2, as.numeric)
plasma_all2[, str_which(colnames(plasma_all2), pattern = "Normalized.Area")] <- apply(plasma_all2[, str_which(colnames(plasma_all2), pattern = "Normalized.Area")], 2, log)


plasma_cell_all2 <- plasma_cell_all

plasma_cell_all2[, str_which(colnames(plasma_cell_all), pattern = "Normalized.Area")] <- apply(plasma_cell_all[, str_which(colnames(plasma_cell_all), pattern = "Normalized.Area")], 2, as.numeric)
plasma_cell_all2[, str_which(colnames(plasma_cell_all2), pattern = "Normalized.Area")] <- apply(plasma_cell_all2[, str_which(colnames(plasma_cell_all2), pattern = "Normalized.Area")], 2, log)


plasma_cell_a12 <- ggplot(plasma_cell_all2, aes(x = stageA.1.Normalized.Area, y = stageA.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageA.1.Normalized.Area, plasma_cell_all2$stageA.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageA.1.Normalized.Area), max(plasma_cell_all2$stageA.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageA.2.Normalized.Area), max(plasma_cell_all2$stageA.2.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_a32 <- ggplot(plasma_cell_all2, aes(x = stageA.3.Normalized.Area, y = stageA.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageA.3.Normalized.Area, plasma_cell_all2$stageA.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageA.3.Normalized.Area), max(plasma_cell_all2$stageA.3.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageA.2.Normalized.Area), max(plasma_cell_all2$stageA.2.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R3", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_a13 <- ggplot(plasma_cell_all2, aes(x = stageA.1.Normalized.Area, y = stageA.3.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageA.1.Normalized.Area, plasma_cell_all2$stageA.3.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageA.1.Normalized.Area), max(plasma_cell_all2$stageA.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageA.3.Normalized.Area), max(plasma_cell_all2$stageA.3.Normalized.Area)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_b12 <- ggplot(plasma_cell_all2, aes(x = stageB.1.Normalized.Area, y = stageB.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageB.1.Normalized.Area, plasma_cell_all2$stageB.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageB.1.Normalized.Area), max(plasma_cell_all2$stageB.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageB.2.Normalized.Area), max(plasma_cell_all2$stageB.2.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_b32 <- ggplot(plasma_cell_all2, aes(x = stageB.3.Normalized.Area, y = stageB.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8.5, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageB.1.Normalized.Area, plasma_cell_all2$stageB.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageB.3.Normalized.Area), max(plasma_cell_all2$stageB.3.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageB.2.Normalized.Area), max(plasma_cell_all2$stageB.2.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R3", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_b13 <- ggplot(plasma_cell_all2, aes(x = stageB.1.Normalized.Area, y = stageB.3.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageB.1.Normalized.Area, plasma_cell_all2$stageB.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageB.1.Normalized.Area), max(plasma_cell_all2$stageB.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageB.3.Normalized.Area), max(plasma_cell_all2$stageB.3.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_c12 <- ggplot(plasma_cell_all2, aes(x = stageC.1.Normalized.Area, y = stageC.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageC.1.Normalized.Area, plasma_cell_all2$stageC.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageC.1.Normalized.Area), max(plasma_cell_all2$stageC.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageC.2.Normalized.Area), max(plasma_cell_all2$stageC.2.Normalized.Area)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_c32 <- ggplot(plasma_cell_all2, aes(x = stageC.3.Normalized.Area, y = stageC.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =9, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageC.1.Normalized.Area, plasma_cell_all2$stageC.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageC.3.Normalized.Area), max(plasma_cell_all2$stageC.3.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageC.2.Normalized.Area), max(plasma_cell_all2$stageC.2.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R3", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_c13 <- ggplot(plasma_cell_all2, aes(x = stageC.1.Normalized.Area, y = stageC.3.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageC.1.Normalized.Area, plasma_cell_all2$stageC.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageC.1.Normalized.Area), max(plasma_cell_all2$stageC.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageC.3.Normalized.Area), max(plasma_cell_all2$stageC.3.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")



plasma_cell_d12 <- ggplot(plasma_cell_all2, aes(x = stageD.1.Normalized.Area, y = stageD.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageD.1.Normalized.Area, plasma_cell_all2$stageD.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageD.1.Normalized.Area), max(plasma_cell_all2$stageD.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageD.2.Normalized.Area), max(plasma_cell_all2$stageD.2.Normalized.Area)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_d32 <- ggplot(plasma_cell_all2, aes(x = stageD.3.Normalized.Area, y = stageD.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =9, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageD.1.Normalized.Area, plasma_cell_all2$stageD.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageD.3.Normalized.Area), max(plasma_cell_all2$stageD.3.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageD.2.Normalized.Area), max(plasma_cell_all2$stageD.2.Normalized.Area)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Log Intensity R3", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_d13 <- ggplot(plasma_cell_all2, aes(x = stageD.1.Normalized.Area, y = stageD.3.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 14, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageD.1.Normalized.Area, plasma_cell_all2$stageD.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageD.1.Normalized.Area), max(plasma_cell_all2$stageD.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageD.3.Normalized.Area), max(plasma_cell_all2$stageD.3.Normalized.Area)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_e12 <- ggplot(plasma_cell_all2, aes(x = stageE.1.Normalized.Area, y = stageE.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 15, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageE.1.Normalized.Area, plasma_cell_all2$stageE.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageE.1.Normalized.Area), max(plasma_cell_all2$stageE.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageE.2.Normalized.Area), max(plasma_cell_all2$stageE.2.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_e32 <- ggplot(plasma_cell_all2, aes(x = stageE.3.Normalized.Area, y = stageE.2.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 15, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageE.1.Normalized.Area, plasma_cell_all2$stageE.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageE.3.Normalized.Area), max(plasma_cell_all2$stageE.3.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageE.2.Normalized.Area), max(plasma_cell_all2$stageE.2.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R3", y = "Log Intensity R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_e13 <- ggplot(plasma_cell_all2, aes(x = stageE.1.Normalized.Area, y = stageE.3.Normalized.Area)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =8, y = 15, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all2$stageE.1.Normalized.Area, plasma_cell_all2$stageE.2.Normalized.Area), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all2$stageE.1.Normalized.Area), max(plasma_cell_all2$stageE.1.Normalized.Area)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all2$stageE.3.Normalized.Area), max(plasma_cell_all2$stageE.3.Normalized.Area)),
                     expand = c(0,0)) +
  labs(title = "", x = "Log Intensity R1", y = "Log Intensity R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/Intensity_replicate_correlation.tiff", units = "in", width = 9, height = 13, res = 300)


sa <- ggarrange(plasma_cell_a12, plasma_cell_a32, plasma_cell_a13, ncol = 3, nrow = 1)
sb <- ggarrange(plasma_cell_b12, plasma_cell_b32, plasma_cell_b13, ncol = 3, nrow = 1)
sc <- ggarrange(plasma_cell_c12, plasma_cell_c32, plasma_cell_c13, ncol = 3, nrow = 1)
sd <- ggarrange(plasma_cell_d12, plasma_cell_d32, plasma_cell_d13, ncol = 3, nrow = 1)
se <- ggarrange(plasma_cell_e12, plasma_cell_e32, plasma_cell_e13, ncol = 3, nrow = 1)


ggarrange(sa, sb, sc, sd, se, ncol = 1, nrow = 5, labels = c("Stage A", "Stage B", "Stage C", "Stage D", "Stage E"))

dev.off()

##################################################################################
########################################################

#### Retention time among technical replicates


plasma_cell_a12 <- ggplot(plasma_cell_all, aes(x = stageA.1.Peptide.Retention.Time, y = stageA.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageA.1.Peptide.Retention.Time, plasma_cell_all$stageA.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageA.1.Peptide.Retention.Time), max(plasma_cell_all$stageA.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageA.2.Peptide.Retention.Time), max(plasma_cell_all$stageA.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_a32 <- ggplot(plasma_cell_all, aes(x = stageA.3.Peptide.Retention.Time, y = stageA.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageA.3.Peptide.Retention.Time, plasma_cell_all$stageA.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageA.3.Peptide.Retention.Time), max(plasma_cell_all$stageA.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageA.2.Peptide.Retention.Time), max(plasma_cell_all$stageA.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R3", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_a13 <- ggplot(plasma_cell_all, aes(x = stageA.1.Peptide.Retention.Time, y = stageA.3.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageA.1.Peptide.Retention.Time, plasma_cell_all$stageA.3.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageA.1.Peptide.Retention.Time), max(plasma_cell_all$stageA.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageA.3.Peptide.Retention.Time), max(plasma_cell_all$stageA.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_b12 <- ggplot(plasma_cell_all, aes(x = stageB.1.Peptide.Retention.Time, y = stageB.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageB.1.Peptide.Retention.Time, plasma_cell_all$stageB.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageB.1.Peptide.Retention.Time), max(plasma_cell_all$stageB.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageB.2.Peptide.Retention.Time), max(plasma_cell_all$stageB.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_b32 <- ggplot(plasma_cell_all, aes(x = stageB.3.Peptide.Retention.Time, y = stageB.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageB.1.Peptide.Retention.Time, plasma_cell_all$stageB.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageB.3.Peptide.Retention.Time), max(plasma_cell_all$stageB.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageB.2.Peptide.Retention.Time), max(plasma_cell_all$stageB.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R3", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_b13 <- ggplot(plasma_cell_all, aes(x = stageB.1.Peptide.Retention.Time, y = stageB.3.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageB.1.Peptide.Retention.Time, plasma_cell_all$stageB.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageB.1.Peptide.Retention.Time), max(plasma_cell_all$stageB.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageB.3.Peptide.Retention.Time), max(plasma_cell_all$stageB.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_c12 <- ggplot(plasma_cell_all, aes(x = stageC.1.Peptide.Retention.Time, y = stageC.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageC.1.Peptide.Retention.Time, plasma_cell_all$stageC.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageC.1.Peptide.Retention.Time), max(plasma_cell_all$stageC.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageC.2.Peptide.Retention.Time), max(plasma_cell_all$stageC.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_c32 <- ggplot(plasma_cell_all, aes(x = stageC.3.Peptide.Retention.Time, y = stageC.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageC.1.Peptide.Retention.Time, plasma_cell_all$stageC.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageC.3.Peptide.Retention.Time), max(plasma_cell_all$stageC.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageC.2.Peptide.Retention.Time), max(plasma_cell_all$stageC.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R3", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_c13 <- ggplot(plasma_cell_all, aes(x = stageC.1.Peptide.Retention.Time, y = stageC.3.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageC.1.Peptide.Retention.Time, plasma_cell_all$stageC.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageC.1.Peptide.Retention.Time), max(plasma_cell_all$stageC.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageC.3.Peptide.Retention.Time), max(plasma_cell_all$stageC.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")



plasma_cell_d12 <- ggplot(plasma_cell_all, aes(x = stageD.1.Peptide.Retention.Time, y = stageD.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageD.1.Peptide.Retention.Time, plasma_cell_all$stageD.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageD.1.Peptide.Retention.Time), max(plasma_cell_all$stageD.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageD.2.Peptide.Retention.Time), max(plasma_cell_all$stageD.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_d32 <- ggplot(plasma_cell_all, aes(x = stageD.3.Peptide.Retention.Time, y = stageD.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageD.1.Peptide.Retention.Time, plasma_cell_all$stageD.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageD.3.Peptide.Retention.Time), max(plasma_cell_all$stageD.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageD.2.Peptide.Retention.Time), max(plasma_cell_all$stageD.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Retention Time R3", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_d13 <- ggplot(plasma_cell_all, aes(x = stageD.1.Peptide.Retention.Time, y = stageD.3.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageD.1.Peptide.Retention.Time, plasma_cell_all$stageD.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageD.1.Peptide.Retention.Time), max(plasma_cell_all$stageD.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageD.3.Peptide.Retention.Time), max(plasma_cell_all$stageD.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  # scale_x_continuous(breaks=seq(10, 100, 10)) +
  # scale_y_continuous(breaks=seq(10, 100, 20)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_e12 <- ggplot(plasma_cell_all, aes(x = stageE.1.Peptide.Retention.Time, y = stageE.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageE.1.Peptide.Retention.Time, plasma_cell_all$stageE.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageE.1.Peptide.Retention.Time), max(plasma_cell_all$stageE.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageE.2.Peptide.Retention.Time), max(plasma_cell_all$stageE.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_e32 <- ggplot(plasma_cell_all, aes(x = stageE.3.Peptide.Retention.Time, y = stageE.2.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageE.1.Peptide.Retention.Time, plasma_cell_all$stageE.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageE.3.Peptide.Retention.Time), max(plasma_cell_all$stageE.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageE.2.Peptide.Retention.Time), max(plasma_cell_all$stageE.2.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R3", y = "Retention Time R2") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_e13 <- ggplot(plasma_cell_all, aes(x = stageE.1.Peptide.Retention.Time, y = stageE.3.Peptide.Retention.Time)) +
  geom_point(color = "Black", size = 1, alpha = 0.6) +
  annotate("text", x =30, y = 61, 
           label = deparse(bquote(italic(R)^2 ==. (format(
             cor(plasma_cell_all$stageE.1.Peptide.Retention.Time, plasma_cell_all$stageE.2.Peptide.Retention.Time), digits = 2)))),
           color = "black", parse = T)+
  geom_smooth(method = "auto", color = "black") +
  # ggtitle("RT Correlation") +
  scale_x_continuous(limits = c(min(plasma_cell_all$stageE.1.Peptide.Retention.Time), max(plasma_cell_all$stageE.1.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(min(plasma_cell_all$stageE.3.Peptide.Retention.Time), max(plasma_cell_all$stageE.3.Peptide.Retention.Time)),
                     expand = c(0,0)) +
  labs(title = "", x = "Retention Time R1", y = "Retention Time R3") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/retentionTime_replicate_correlation.tiff", units = "in", width = 9, height = 13, res = 300)


sa <- ggarrange(plasma_cell_a12, plasma_cell_a32, plasma_cell_a13, ncol = 3, nrow = 1)
sb <- ggarrange(plasma_cell_b12, plasma_cell_b32, plasma_cell_b13, ncol = 3, nrow = 1)
sc <- ggarrange(plasma_cell_c12, plasma_cell_c32, plasma_cell_c13, ncol = 3, nrow = 1)
sd <- ggarrange(plasma_cell_d12, plasma_cell_d32, plasma_cell_d13, ncol = 3, nrow = 1)
se <- ggarrange(plasma_cell_e12, plasma_cell_e32, plasma_cell_e13, ncol = 3, nrow = 1)


ggarrange(sa, sb, sc, sd, se, ncol = 1, nrow = 5, labels = c("Stage A", "Stage B", "Stage C", "Stage D", "Stage E"))

dev.off()

#####################################################################################################
#####################################################################################################

### AVERAGE RETENTION TIME COMPARISON BETWEEN STAGES FOR DIA FROM PLASMA_ALL, CELL_ALL, AND PLASMA_CELL_ALL LIBRARIES

plasma_all <- plasma_all %>%
  mutate("stageA.Average.Retention.Time" = apply(m <- (plasma_all[, c("stageA.1.Peptide.Retention.Time", "stageA.2.Peptide.Retention.Time", 
                                                                      "stageA.3.Peptide.Retention.Time")]), 1, median))
plasma_all <- plasma_all %>%
  mutate("stageB.Bverage.Retention.Time" = apply(m <- (plasma_all[, c("stageB.1.Peptide.Retention.Time", "stageB.2.Peptide.Retention.Time", 
                                                                      "stageB.3.Peptide.Retention.Time")]), 1, median))
plasma_all <- plasma_all %>%
  mutate("stageC.Cverage.Retention.Time" = apply(m <- (plasma_all[, c("stageC.1.Peptide.Retention.Time", "stageC.2.Peptide.Retention.Time", 
                                                                      "stageC.3.Peptide.Retention.Time")]), 1, median))
plasma_all <- plasma_all %>%
  mutate("stageD.Dverage.Retention.Time" = apply(m <- (plasma_all[, c("stageD.1.Peptide.Retention.Time", "stageD.2.Peptide.Retention.Time", 
                                                                      "stageD.3.Peptide.Retention.Time")]), 1, median))
plasma_all <- plasma_all %>%
  mutate("stageE.Everage.Retention.Time" = apply(m <- (plasma_all[, c("stageE.1.Peptide.Retention.Time", "stageE.2.Peptide.Retention.Time", 
                                                                      "stageE.3.Peptide.Retention.Time")]), 1, median))


cell_all <- cell_all %>%
  mutate("stageA.Average.Retention.Time" = apply(m <- (cell_all[, c("stageA.1.Peptide.Retention.Time", "stageA.2.Peptide.Retention.Time", 
                                                                    "stageA.3.Peptide.Retention.Time")]), 1, median))
cell_all <- cell_all %>%
  mutate("stageB.Bverage.Retention.Time" = apply(m <- (cell_all[, c("stageB.1.Peptide.Retention.Time", "stageB.2.Peptide.Retention.Time", 
                                                                    "stageB.3.Peptide.Retention.Time")]), 1, median))
cell_all <- cell_all %>%
  mutate("stageC.Cverage.Retention.Time" = apply(m <- (cell_all[, c("stageC.1.Peptide.Retention.Time", "stageC.2.Peptide.Retention.Time", 
                                                                    "stageC.3.Peptide.Retention.Time")]), 1, median))
cell_all <- cell_all %>%
  mutate("stageD.Dverage.Retention.Time" = apply(m <- (cell_all[, c("stageD.1.Peptide.Retention.Time", "stageD.2.Peptide.Retention.Time", 
                                                                    "stageD.3.Peptide.Retention.Time")]), 1, median))
cell_all <- cell_all %>%
  mutate("stageE.Everage.Retention.Time" = apply(m <- (cell_all[, c("stageE.1.Peptide.Retention.Time", "stageE.2.Peptide.Retention.Time", 
                                                                    "stageE.3.Peptide.Retention.Time")]), 1, median))


plasma_cell_all <- plasma_cell_all %>%
  mutate("stageA.Average.Retention.Time" = apply(m <- (plasma_cell_all[, c("stageA.1.Peptide.Retention.Time", "stageA.2.Peptide.Retention.Time", 
                                                                           "stageA.3.Peptide.Retention.Time")]), 1, median))
plasma_cell_all <- plasma_cell_all %>%
  mutate("stageB.Bverage.Retention.Time" = apply(m <- (plasma_cell_all[, c("stageB.1.Peptide.Retention.Time", "stageB.2.Peptide.Retention.Time", 
                                                                           "stageB.3.Peptide.Retention.Time")]), 1, median))
plasma_cell_all <- plasma_cell_all %>%
  mutate("stageC.Cverage.Retention.Time" = apply(m <- (plasma_cell_all[, c("stageC.1.Peptide.Retention.Time", "stageC.2.Peptide.Retention.Time", 
                                                                           "stageC.3.Peptide.Retention.Time")]), 1, median))
plasma_cell_all <- plasma_cell_all %>%
  mutate("stageD.Dverage.Retention.Time" = apply(m <- (plasma_cell_all[, c("stageD.1.Peptide.Retention.Time", "stageD.2.Peptide.Retention.Time", 
                                                                           "stageD.3.Peptide.Retention.Time")]), 1, median))
plasma_cell_all <- plasma_cell_all %>%
  mutate("stageE.Everage.Retention.Time" = apply(m <- (plasma_cell_all[, c("stageE.1.Peptide.Retention.Time", "stageE.2.Peptide.Retention.Time", 
                                                                           "stageE.3.Peptide.Retention.Time")]), 1, median))


plasma_all <- plasma_all %>%
  select(Modified.Sequence, 102:106)
colnames(plasma_all) <- c("Peptides", "StageA", "StageB", "StageC", "StageD", "StageE")

cell_all <- cell_all %>%
  select(Modified.Sequence, 102:106)
colnames(cell_all) <- c("Peptides", "StageA", "StageB", "StageC", "StageD", "StageE")

plasma_cell_all <- plasma_cell_all %>%
  select(Modified.Sequence, 102:106)
colnames(plasma_cell_all) <- c("Peptides", "StageA", "StageB", "StageC", "StageD", "StageE")


p_all <- corrgram(plasma_all, lower.panel=panel.pts, upper.panel=panel.conf, order=FALSE, text.panel=panel.txt, label.pos = c(0.4, 0.85),
                  cex.labels = 3, diag.panel = panel.density,  cor.method = "pearson")
mtext("DIA Retention Time (min)", side=1, cex=2, line = -2, outer=TRUE, xpd=NA)
mtext("DIA Retention Time (min)", side=2, cex=2, line = -2, outer=TRUE, xpd=NA)
grid.echo()
P1 <- grid.grab()
P1 <- as_ggplot(P1)

c_all <- corrgram(cell_all, lower.panel=panel.pts, upper.panel=panel.conf, order=FALSE, text.panel=panel.txt, label.pos = c(0.4, 0.85),
                  cex.labels = 3, diag.panel = panel.density,  cor.method = "pearson")
mtext("DIA Retention Time (min)", side=1, cex=2, line = -2, outer=TRUE, xpd=NA)
mtext("DIA Retention Time (min)", side=2, cex=2, line = -2, outer=TRUE, xpd=NA)
grid.echo()
P2 <- grid.grab()
P2 <- as_ggplot(P2)

pc_all <- corrgram(plasma_cell_all, lower.panel=panel.pts, upper.panel=panel.conf, order=FALSE, text.panel=panel.txt, label.pos = c(0.4, 0.85),
                   cex.labels = 3, diag.panel = panel.density,  cor.method = "pearson")
mtext("DIA Retention Time (min)", side=1, cex=2, line = -2, outer=TRUE, xpd=NA)
mtext("DIA Retention Time (min)", side=2, cex=2, line = -2, outer=TRUE, xpd=NA)
grid.echo()
P3 <- grid.grab()
P3 <- as_ggplot(P3)


tiff("E:/MQ-Charlie-Data_2018Nov20/All_libraries/Figures/Plasma_Cell_Plasma-Cell_Stage_RT_cor.tiff", units = "in", width = 18, height = 14, res = 300)
ggarrange(P1, P2, P3, nrow = 2, ncol = 2, labels = c("A. Plasma", "B. Cell", "C. Plasma_Cell"), font.label = list(size = 25))


dev.off()



#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################