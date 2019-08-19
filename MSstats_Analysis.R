
library(ggplot2)
library(ggpubr)
library(dplyr)
library(VennDiagram)
library(stringr)
library(MSstats)
library(rebus)

###########
### Libraries

lib_plasma_pp <- read.csv(".../Plasma_Library_PP.csv", stringsAsFactors = F)
lib_plasma_tpp <- read.csv(".../Plasma_Library_TPP.csv", stringsAsFactors = F)

lib_cell_pp <- read.csv(".../Cell_Calibrated_Library_PP.csv", stringsAsFactors = F)
lib_cell_tpp <- read.csv(".../Cell_Calibrated_Library_TPP.csv", stringsAsFactors = F)

lib_plasma_cell_pp <- read.csv(".../Plasma_Cell_Library_PP.csv", stringsAsFactors = F)
lib_plasma_cell_tpp <- read.csv(".../Plasma_Cell_Library_TPP.csv", stringsAsFactors = F)

##################
###########


##### Reports MSstats Stages A and E _ Plasma_Cell


plasma_cell_PP_ms <- read.csv(".../MSstats_R_plasma_cell_pp.csv", stringsAsFactors = F)

plasma_cell_TPP_ms <- read.csv(".../MSstats_R_plasma_cell_tpp.csv", stringsAsFactors = F)

### Combining PP and TPP

protein.names <- plasma_cell_PP_ms[, 1:2]
protein.names <- protein.names %>%
  mutate("UniprotID" = "")

for(i in 1:length(protein.names$Protein.Name))
{
  uniprotID <- lib_plasma_cell_pp[lib_plasma_cell_pp$ProteinName == protein.names$Protein.Name[i], "UniprotID"]
  protein.names$UniprotID[i] <- toString(uniprotID[1])
}

plasma_cell_PP_ms$Protein.Name <- protein.names$UniprotID


###
protein.names <- plasma_cell_TPP_ms[, 1:2]
protein.names <- protein.names %>%
  mutate("UniprotID" = "")

for(i in 1:length(protein.names$Protein.Name))
{
  uniprotID <- lib_plasma_cell_tpp[lib_plasma_cell_tpp$ProteinName == protein.names$Protein.Name[i], "UniprotID"]
  protein.names$UniprotID[i] <- toString(uniprotID[1])
}

plasma_cell_TPP_ms$Protein.Name <- protein.names$UniprotID

comm_pep <- intersect(plasma_cell_PP_ms$Peptide.Modified.Sequence, plasma_cell_TPP_ms$Peptide.Modified.Sequence)

unique_PP <- plasma_cell_PP_ms[!(plasma_cell_PP_ms$Peptide.Modified.Sequence %in% comm_pep), ]
unique_TPP <- plasma_cell_TPP_ms[!(plasma_cell_TPP_ms$Peptide.Modified.Sequence %in% comm_pep), ]

plasma_cell_all_ms <- rbind(plasma_cell_PP_ms, unique_TPP)

write.csv(plasma_cell_all_ms, ".../MSstats_R_plasma_cell_all.csv", row.names = FALSE)

quant <- SkylinetoMSstatsFormat(plasma_cell_all_ms, removeProtein_with1Feature = TRUE, filter_with_Qvalue = F, useUniquePeptide = T)

sky_data_all <- dataProcess(quant, normalization = "None", summaryMethod = "TMP", cutoffCensored = "minFeature",
                            censoredInt = 0, MBimpute = TRUE, maxQuantileforCensored = 0.999)


levels(sky_data_all$ProcessedData$GROUP_ORIGINAL)


comparison1 <- matrix(c(-1, 1, 0, 0, 0), nrow = 1)
comparison2 <- matrix(c(-1, 0, 1, 0, 0), nrow = 1)
comparison3 <- matrix(c(-1, 0, 0, 1, 0), nrow = 1)
comparison4 <- matrix(c(-1, 0, 0, 0, 1), nrow = 1)


comparison <- rbind(comparison1, comparison2, comparison3, comparison4)
row.names(comparison) <- c("Stage A - Healthy", "Stage B - Healthy", "Stage C - Healthy", "Stage D - Healthy")


sky_data_all.comparisons <- groupComparison(contrast.matrix = comparison, data = sky_data_all)

View(sky_data_all.comparisons$ComparisonResult)


SignificantProteins_pc_all <- with(sky_data_all.comparisons, ComparisonResult[ComparisonResult$pvalue <= 0.05, ])


Proteins_Plasma_Cell <- sky_data_all.comparisons$ComparisonResult

Proteins_Plasma_Cell <- Proteins_Plasma_Cell %>%
  filter(adj.pvalue != 0.000)

Proteins_Plasma_Cell <- Proteins_Plasma_Cell %>%
  mutate("Dysregulation" = "")

Proteins_Plasma_Cell[Proteins_Plasma_Cell$log2FC >= 0, "Dysregulation"] <- "cNo regulation"
Proteins_Plasma_Cell[Proteins_Plasma_Cell$log2FC <= 0, "Dysregulation"] <- "cNo regulation"
Proteins_Plasma_Cell[Proteins_Plasma_Cell$log2FC >= 1.5, "Dysregulation"] <- "aUp-regulated"
Proteins_Plasma_Cell[Proteins_Plasma_Cell$log2FC <= -1.5, "Dysregulation"] <- "bDown-regulated"
Proteins_Plasma_Cell[Proteins_Plasma_Cell$pvalue >= 0.05, "Dysregulation"] <- "cNo regulation"

Proteins_Plasma_Cell_sAB <- Proteins_Plasma_Cell %>%
  filter(Label == "Stage A - Healthy" | Label == "Stage B - Healthy")

Proteins_Plasma_Cell_sAB <- Proteins_Plasma_Cell_sAB %>%
  mutate("UniprotID" = str_replace(Protein , pattern =  START %R%
                                     "sp" %R% ANY_CHAR, replacement = "")) %>%
  mutate("UniprotID" = str_extract(UniprotID , pattern =  START %R%
                                     one_or_more(WRD)))

for(i in 1:length(Proteins_Plasma_Cell_sAB$Dysregulation))
{
  if(Proteins_Plasma_Cell_ae_sAB$Dysregulation[i] == "cNo regulation")
    Proteins_Plasma_Cell_ae_sAB$UniprotID[i] <-  " "
}

tiff(".../Figures/Plasma_cell_StageAB_volc.tiff", units = "in", width = 5, height = 3.5, res = 300)

pc_plot <- ggplot(Proteins_Plasma_Cell_sAB, aes(x = log2FC, y = -log(pvalue, base = 10), color = Dysregulation)) +
  geom_point(alpha = 1, size = 0.7)+
  annotate("text", x = -2, y = -log(0.01, base = 10),
           label = "p-value < 0.05", color = "black", size = 2.5) +
  facet_wrap(.~Label) +
  scale_x_continuous(breaks=seq(-4, 4, 1)) +
  geom_hline(aes(yintercept=-log(0.05, base = 10)), linetype="dashed", color = "black", size = 0.3) +
  geom_vline(aes(xintercept= 1.5), linetype="dashed", color = "black", size = 0.3) +
  geom_vline(aes(xintercept= -1.5), linetype="dashed", color = "black", size = 0.3) +
  scale_color_manual(labels = c("Up-regulated", "Down-regulated", "No regulation"), values=c("#ff4d4d", "#569cdd", "darkgray")) +
  labs(subtitle = "Plasma_Cell", x = "Log2 Fold Change", y = "-Log10 (pvalue)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold", vjust = 0.8),
        axis.title.y = element_text(size = 10, face = "bold", vjust = 0.8),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

pc_plot

dev.off()

###################################################################################################################################
##################################################################################################################


###############################
##############################################
###########################

#### post analysis

library(tidyr)


data_pc <- Proteins_Plasma_Cell_ae %>%
  filter(pvalue <= 0.05)
data_pc <- select(data_pc, Protein, Label, log2FC, Dysregulation)

data2_pc <- data_pc %>% 
  group_by(Label) %>% 
  spread(Label, log2FC) 

write.csv(data2_pc, ".../allstageAnalysisPC_07.csv", row.names = F)
###



###############################
##############################################
###########################


###############################
##############################################
###########################
