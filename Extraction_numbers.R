
##

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
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(plotrix)
library(cowplot)
library(Hmisc)
library(tidyr)
library(rtiff)
library(grid)
library(ggplotify)
library(gridExtra)


###########
##################
###########

#### Libraries

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


#######
#######
### Reports

plasma_PP <- read.csv(".../Report_file_R_plasma_PP.csv", stringsAsFactors = F)
cell_PP <- read.csv(".../Report_file_R_cell_PP.csv", stringsAsFactors = F)
plasma_cell_PP <- read.csv(".../Report_file_R_plasma_cell_pp.csv", stringsAsFactors = F)


plasma_TPP <- read.csv(".../Report_file_R_plasma_TPP.csv", stringsAsFactors = F)
cell_TPP <- read.csv(".../Report_file_R_cell_TPP.csv", stringsAsFactors = F)
plasma_cell_TPP <- read.csv(".../Report_file_R_plasma_cell_tpp.csv", stringsAsFactors = F)


##### Removing NAs and any Q_val = >0.01

plasma_PP <- plasma_PP %>%
  filter_all(all_vars(. != "#N/A"))

plasma_PP <- plasma_PP %>%
  filter_at(vars(ends_with("Value")), any_vars(. <=0.01))

cell_PP <- cell_PP %>%
  filter_all(all_vars(. != "#N/A"))

cell_PP <- cell_PP %>%
  filter_at(vars(ends_with("Value")), any_vars(. <=0.01))

plasma_cell_PP <- plasma_cell_PP %>%
  filter_all(all_vars(. != "#N/A"))

plasma_cell_PP <- plasma_cell_PP %>%
  filter_at(vars(ends_with("Value")), any_vars(. <=0.01))


plasma_TPP <- plasma_TPP %>%
  filter_all(all_vars(. != "#N/A"))

plasma_TPP <- plasma_TPP %>%
  filter_at(vars(ends_with("Value")), any_vars(. <=0.01))

cell_TPP <- cell_TPP %>%
  filter_all(all_vars(. != "#N/A"))

cell_TPP <- cell_TPP %>%
  filter_at(vars(ends_with("Value")), any_vars(. <=0.01))

plasma_cell_TPP <- plasma_cell_TPP %>%
  filter_all(all_vars(. != "#N/A"))

plasma_cell_TPP <- plasma_cell_TPP %>%
  filter_at(vars(ends_with("Value")), any_vars(. <=0.01))


#######################
####

###########################
#####################################

## Number of quantified proteins and peptides
### Combining results of TPP and PP

plasma_PP <- plasma_PP[!duplicated(plasma_PP$Modified.Sequence),]
plasma_TPP <- plasma_TPP[!duplicated(plasma_TPP$Modified.Sequence),]

check_PP <- merge(plasma_PP, lib_plasma_pp, by.x = 1, by.y = 4, all.x = TRUE)
check_PP2 <- check_PP[check_PP$Modified.Sequence %in% plasma_PP$Modified.Sequence, c(1:101, 108)]
plasma_PP <- check_PP2[!duplicated(check_PP2$Modified.Sequence),]

check_TPP <- merge(plasma_TPP, lib_plasma_tpp, by.x = 1, by.y = 4, all.x = TRUE)
check_TPP2 <- check_TPP[check_TPP$Modified.Sequence %in% plasma_TPP$Modified.Sequence, c(1:101, 108)]
plasma_TPP <- check_TPP2[!duplicated(check_TPP2$Modified.Sequence),]

comm_pep <- intersect(plasma_PP$Modified.Sequence, plasma_TPP$Modified.Sequence)

unique_PP <- plasma_PP[!(plasma_PP$Modified.Sequence %in% comm_pep), ]
unique_TPP <- plasma_TPP[!(plasma_TPP$Modified.Sequence %in% comm_pep), ]

plasma_all <- rbind(plasma_PP, unique_TPP)

write.csv(plasma_all, ".../Report_file_R_plasma_all.csv", row.names = F)

###########################
#####################################


cell_PP <- cell_PP[!duplicated(cell_PP$Modified.Sequence),]
cell_TPP <- cell_TPP[!duplicated(cell_TPP$Modified.Sequence),]

check_PP <- merge(cell_PP, lib_cell_pp, by.x = 1, by.y = 4, all.x = TRUE)
check_PP2 <- check_PP[check_PP$Modified.Sequence %in% cell_PP$Modified.Sequence, c(1:101, 108)]
cell_PP <- check_PP2[!duplicated(check_PP2$Modified.Sequence),]

check_TPP <- merge(cell_TPP, lib_cell_tpp, by.x = 1, by.y = 4, all.x = TRUE)
check_TPP2 <- check_TPP[check_TPP$Modified.Sequence %in% cell_TPP$Modified.Sequence, c(1:101, 108)]
cell_TPP <- check_TPP2[!duplicated(check_TPP2$Modified.Sequence),]

comm_pep <- intersect(cell_PP$Modified.Sequence, cell_TPP$Modified.Sequence)

unique_PP <- cell_PP[!(cell_PP$Modified.Sequence %in% comm_pep), ]
unique_TPP <- cell_TPP[!(cell_TPP$Modified.Sequence %in% comm_pep), ]

cell_all <- rbind(cell_PP, unique_TPP)

write.csv(cell_all, ".../Report_file_R_cell_all.csv", row.names = F)


###########################
#####################################


plasma_cell_PP <- plasma_cell_PP[!duplicated(plasma_cell_PP$Modified.Sequence),]
plasma_cell_TPP <- plasma_cell_TPP[!duplicated(plasma_cell_TPP$Modified.Sequence),]

check_PP <- merge(plasma_cell_PP, lib_plasma_cell_pp, by.x = 1, by.y = 4, all.x = TRUE)
check_PP2 <- check_PP[check_PP$Modified.Sequence %in% plasma_cell_PP$Modified.Sequence, c(1:101, 108)]
plasma_cell_PP <- check_PP2[!duplicated(check_PP2$Modified.Sequence),]

check_TPP <- merge(plasma_cell_TPP, lib_plasma_cell_tpp, by.x = 1, by.y = 4, all.x = TRUE)
check_TPP2 <- check_TPP[check_TPP$Modified.Sequence %in% plasma_cell_TPP$Modified.Sequence, c(1:101, 108)]
plasma_cell_TPP <- check_TPP2[!duplicated(check_TPP2$Modified.Sequence),]

comm_pep <- intersect(plasma_cell_PP$Modified.Sequence, plasma_cell_TPP$Modified.Sequence)

unique_PP <- plasma_cell_PP[!(plasma_cell_PP$Modified.Sequence %in% comm_pep), ]
unique_TPP <- plasma_cell_TPP[!(plasma_cell_TPP$Modified.Sequence %in% comm_pep), ]

plasma_cell_all <- rbind(plasma_cell_PP, unique_TPP)

write.csv(plasma_cell_all, ".../Report_file_R_plasma_cell_all.csv", row.names = F)

##############################################
#######################################

#####
###########

plasma_all_prot <- plasma_all[!duplicated(plasma_all$UniprotID),"UniprotID"]
cell_all_prot <- cell_all[!duplicated(cell_all$UniprotID), "UniprotID"]
plasma_cell_all_prot <- plasma_cell_all[!duplicated(plasma_cell_all$UniprotID), "UniprotID"]


plasma_all_pep <- plasma_all[!duplicated(plasma_all$Modified.Sequence), "Modified.Sequence"]
cell_all_pep <- cell_all[!duplicated(cell_all$Modified.Sequence), "Modified.Sequence"]
plasma_cell_all_pep <- plasma_cell_all[!duplicated(plasma_cell_all$Modified.Sequence), "Modified.Sequence"]

##2pep

plasma_all_prot <- plasma_all[!duplicated(plasma_all$UniprotID),"UniprotID"]


plasma_all_prot2 <- aggregate(x = plasma_all$Modified.Sequence, by = list(Protein = plasma_all$UniprotID), FUN = function(x) {length(x)})
plasma_all_prot <- plasma_all_prot2[plasma_all_prot2$x >1,]

length(plasma_all_prot$Protein)
sum(plasma_all_prot$x)

plasma_all_nn <- plasma_all %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

plasma_all_prot2 <- aggregate(x = plasma_all_nn$Modified.Sequence, by = list(Protein = plasma_all_nn$UniprotID), FUN = function(x) {length(x)})
plasma_all_prot <- plasma_all_prot2[plasma_all_prot2$x >1,]

length(plasma_all_prot$Protein)
sum(plasma_all_prot$x)
##

cell_all_prot <- cell_all[!duplicated(cell_all$Protein.Name),"UniprotID"]


cell_all_prot2 <- aggregate(x = cell_all$Modified.Sequence, by = list(Protein = cell_all$UniprotID), FUN = function(x) {length(x)})
cell_all_prot <- cell_all_prot2[cell_all_prot2$x >1,]

length(cell_all_prot$Protein)
sum(cell_all_prot$x)

cell_all_nn <- cell_all %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

cell_all_prot2 <- aggregate(x = cell_all_nn$Modified.Sequence, by = list(Protein = cell_all_nn$UniprotID), FUN = function(x) {length(x)})
cell_all_prot <- cell_all_prot2[cell_all_prot2$x >1,]

length(cell_all_prot$Protein)
sum(cell_all_prot$x)
##


plasma_cell_all_prot <- plasma_cell_all[!duplicated(plasma_cell_all$Protein.Name),"UniprotID"]


plasma_cell_all_prot2 <- aggregate(x = plasma_cell_all$Modified.Sequence, by = list(Protein = plasma_cell_all$UniprotID), FUN = function(x) {length(x)})
plasma_cell_all_prot <- plasma_cell_all_prot2[plasma_cell_all_prot2$x >1,]

length(plasma_cell_all_prot$Protein)
sum(plasma_cell_all_prot$x)

plasma_cell_all_nn <- plasma_cell_all %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

plasma_cell_all_prot2 <- aggregate(x = plasma_cell_all_nn$Modified.Sequence, by = list(Protein = plasma_cell_all_nn$UniprotID), FUN = function(x) {length(x)})
plasma_cell_all_prot <- plasma_cell_all_prot2[plasma_cell_all_prot2$x >1,]

length(plasma_cell_all_prot$Protein)
sum(plasma_cell_all_prot$x)

#####
##########

plasma_PP_prot <- plasma_PP[!duplicated(plasma_PP$UniprotID),"UniprotID"]


plasma_PP_prot2 <- aggregate(x = plasma_PP$Modified.Sequence, by = list(Protein = plasma_PP$UniprotID), FUN = function(x) {length(x)})
plasma_PP_prot <- plasma_PP_prot2[plasma_PP_prot2$x >1,]

length(plasma_PP_prot$Protein)
sum(plasma_PP_prot$x)

plasma_PP_nn <- plasma_PP %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

plasma_PP_prot2 <- aggregate(x = plasma_PP_nn$Modified.Sequence, by = list(Protein = plasma_PP_nn$UniprotID), FUN = function(x) {length(x)})
plasma_PP_prot <- plasma_PP_prot2[plasma_PP_prot2$x >1,]

length(plasma_PP_prot$Protein)
sum(plasma_PP_prot$x)
##

cell_PP_prot <- cell_PP[!duplicated(cell_PP$UniprotID), "UniprotID"]

cell_PP_prot2 <- aggregate(x = cell_PP$Modified.Sequence, by = list(Protein = cell_PP$UniprotID), FUN = function(x) {length(x)})
cell_PP_prot <- cell_PP_prot2[cell_PP_prot2$x >1, ]

length(cell_PP_prot$Protein)
sum(cell_PP_prot$x)

cell_PP_nn <- cell_PP %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

cell_PP_prot2 <- aggregate(x = cell_PP_nn$Modified.Sequence, by = list(Protein = cell_PP_nn$UniprotID), FUN = function(x) {length(x)})
cell_PP_prot <- cell_PP_prot2[cell_PP_prot2$x >1, ]

length(cell_PP_prot$Protein)
sum(cell_PP_prot$x)
##


plasma_cell_PP_prot <- plasma_cell_PP[!duplicated(plasma_cell_PP$UniprotID), "UniprotID"]

plasma_cell_PP_prot2 <- aggregate(x = plasma_cell_PP$Modified.Sequence, by = list(Protein = plasma_cell_PP$UniprotID), FUN = function(x) {length(x)})
plasma_cell_PP_prot <- plasma_cell_PP_prot2[plasma_cell_PP_prot2$x >1, ]

length(plasma_cell_PP_prot$Protein)
sum(plasma_cell_PP_prot$x)

plasma_cell_PP_nn <- plasma_cell_PP %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

plasma_cell_PP_prot2 <- aggregate(x = plasma_cell_PP_nn$Modified.Sequence, by = list(Protein = plasma_cell_PP_nn$UniprotID), FUN = function(x) {length(x)})
plasma_cell_PP_prot <- plasma_cell_PP_prot2[plasma_cell_PP_prot2$x >1, ]

length(plasma_cell_PP_prot$Protein)
sum(plasma_cell_PP_prot$x)


plasma_PP_pep <- plasma_PP[!duplicated(plasma_PP$Modified.Sequence), "Modified.Sequence"]
cell_PP_pep <- cell_PP[!duplicated(cell_PP$Modified.Sequence), "Modified.Sequence"]
plasma_cell_PP_pep <- plasma_cell_PP[!duplicated(plasma_cell_PP$Modified.Sequence), "Modified.Sequence"]

#######
###########


plasma_TPP_prot <- plasma_TPP[!duplicated(plasma_TPP$UniprotID),"UniprotID"]

plasma_TPP_prot2 <- aggregate(x = plasma_TPP$Modified.Sequence, by = list(Protein = plasma_TPP$UniprotID), FUN = function(x) {length(x)})
plasma_TPP_prot <- plasma_TPP_prot2[plasma_TPP_prot2$x >1,]

length(plasma_TPP_prot$Protein)
sum(plasma_TPP_prot$x)

plasma_TPP_nn <- plasma_TPP %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

plasma_TPP_prot2 <- aggregate(x = plasma_TPP_nn$Modified.Sequence, by = list(Protein = plasma_TPP_nn$UniprotID), FUN = function(x) {length(x)})
plasma_TPP_prot <- plasma_TPP_prot2[plasma_TPP_prot2$x >1,]

length(plasma_TPP_prot$Protein)
sum(plasma_TPP_prot$x)
##


cell_TPP_prot <- cell_TPP[!duplicated(cell_TPP$UniprotID), "UniprotID"]

cell_TPP_prot2 <- aggregate(x = cell_TPP$Modified.Sequence, by = list(Protein = cell_TPP$UniprotID), FUN = function(x) {length(x)})
cell_TPP_prot <- cell_TPP_prot2[cell_TPP_prot2$x >1, ]

length(cell_TPP_prot$Protein)
sum(cell_TPP_prot$x)

cell_TPP_nn <- cell_TPP %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

cell_TPP_prot2 <- aggregate(x = cell_TPP_nn$Modified.Sequence, by = list(Protein = cell_TPP_nn$UniprotID), FUN = function(x) {length(x)})
cell_TPP_prot <- cell_TPP_prot2[cell_TPP_prot2$x >1, ]

length(cell_TPP_prot$Protein)
sum(cell_TPP_prot$x)
##

plasma_cell_TPP_prot <- plasma_cell_TPP[!duplicated(plasma_cell_TPP$UniprotID), "UniprotID"]

plasma_cell_TPP_prot2 <- aggregate(x = plasma_cell_TPP$Modified.Sequence, by = list(Protein = plasma_cell_TPP$UniprotID), FUN = function(x) {length(x)})
plasma_cell_TPP_prot <- plasma_cell_TPP_prot2[plasma_cell_TPP_prot2$x >1, ]

length(plasma_cell_TPP_prot$Protein)
sum(plasma_cell_TPP_prot$x)


plasma_cell_TPP_nn <- plasma_cell_TPP %>%
  filter(Missed.Cleavages == 0) %>%
  filter(!str_detect(Modified.Sequence, pattern = OPEN_BRACKET))

plasma_cell_TPP_prot2 <- aggregate(x = plasma_cell_TPP_nn$Modified.Sequence, by = list(Protein = plasma_cell_TPP_nn$UniprotID), FUN = function(x) {length(x)})
plasma_cell_TPP_prot <- plasma_cell_TPP_prot2[plasma_cell_TPP_prot2$x >1, ]

length(plasma_cell_TPP_prot$Protein)
sum(plasma_cell_TPP_prot$x)


plasma_TPP_pep <- plasma_TPP[!duplicated(plasma_TPP$Modified.Sequence), "Modified.Sequence"]
cell_TPP_pep <- cell_TPP[!duplicated(cell_TPP$Modified.Sequence), "Modified.Sequence"]
plasma_cell_TPP_pep <- plasma_cell_TPP[!duplicated(plasma_cell_TPP$Modified.Sequence), "Modified.Sequence"]


extractions_data_actual <- data.frame("Libraries" = c("aPlasma", "bCell", "cPlas-Cell", "aPlasma", "bCell", "cPlas-Cell", "aPlasma", "bCell", "cPlas-Cell"), 
                                      "Proteins" = c(length(plasma_TPP_prot), length(cell_TPP_prot), length(plasma_cell_TPP_prot),
                                                     length(plasma_PP_prot), length(cell_PP_prot), length(plasma_cell_PP_prot),
                                                     length(plasma_all_prot), length(cell_all_prot), length(plasma_cell_all_prot)),
                                      "Peptides" = c(length(plasma_TPP_pep), length(cell_TPP_pep), length(plasma_cell_TPP_pep),
                                                     length(plasma_PP_pep), length(cell_PP_pep), length(plasma_cell_PP_pep),
                                                     length(plasma_all_pep), length(cell_all_pep), length(plasma_cell_all_pep)),
                                      "Type" = c("aTPP", "aTPP", "aTPP", "bPP", "bPP", "bPP", "cPP+TPP", "cPP+TPP", "cPP+TPP"))


extractions_data <- data.frame("Libraries" = c("aPlasma", "bCell", "cPlas-Cell", "aPlasma", "bCell", "cPlas-Cell", "aPlasma", "bCell", "cPlas-Cell"), 
                               "Proteins" = c(length(plasma_TPP_prot)/100, length(cell_TPP_prot)/100, length(plasma_cell_TPP_prot)/100,
                                              length(plasma_PP_prot)/100, length(cell_PP_prot)/100, length(plasma_cell_PP_prot)/100,
                                              length(plasma_all_prot)/100, length(cell_all_prot)/100, length(plasma_cell_all_prot)/100),
                               "Peptides" = c(length(plasma_TPP_pep)/100, length(cell_TPP_pep)/100, length(plasma_cell_TPP_pep)/100,
                                              length(plasma_PP_pep)/100, length(cell_PP_pep)/100, length(plasma_cell_PP_pep)/100,
                                              length(plasma_all_pep)/100, length(cell_all_pep)/100, length(plasma_cell_all_pep)/100),
                               "Type" = c("aTPP", "aTPP", "aTPP", "bPP", "bPP", "bPP", "cPP+TPP", "cPP+TPP", "cPP+TPP"))

tiff(".../Figures/extraction_numbers.tiff", units = "in", width = 8, height = 5, res = 300)

prot <- ggplot(data=extractions_data, aes(x = Libraries, y = Proteins, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.9, position=position_dodge()) +
  # geom_text(aes(label=Proteins), vjust=0, color="black", fontface = "bold",
  #           position = position_dodge(0.9), size=5) +
  scale_y_continuous(breaks=seq(0, 28, 4))+
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(labels = c("TPP", "PP", "TPP+PP"), values=c("#FFA500", "#ff3232", "#107dac"))+
  scale_x_discrete(labels = c("Plasma", "Cell", "Plasma-Cell")) +
  labs(y = deparse(bquote(Proteins(10^2)))) +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right") 

pep <- ggplot(data=extractions_data, aes(x = Libraries, y = Peptides, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.9, position=position_dodge()) +
  # geom_text(aes(label=Peptides), vjust=0, color="black", fontface = "bold",
  #           position = position_dodge(0.9), size=5) +
  scale_y_continuous(breaks=seq(0, 52, 4))+
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(labels = c("TPP", "PP", "TPP+PP"), values=c("#FFA500", "#ff3232", "#107dac"))+
  scale_x_discrete(labels = c("Plasma", "Cell", "Plasma-Cell")) +
  labs(y = deparse(bquote(Peptides(10^2))))+
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        # plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
        # text = element_text(family = "Tahoma"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right") 

ggarrange(prot, pep, ncol = 2, nrow = 1, common.legend = TRUE, labels = "AUTO", font.label = list(size = 18))

dev.off()

##### Extractions Venn diagram PP and TPP

plasma_pep <- venn.diagram(list(plasma_TPP_pep, plasma_PP_pep),
                           category.names = c("TPP_Peptides" , "PP_Peptides"),
                           resolution = 300,
                           height = 2000,
                           width = 2000,
                           cat.default.pos = "outer",
                           lwd = 1,
                           fill = c("#FFA500", "#ff3232"),
                           cat.pos = c(-170, 170),
                           alpha = c(0.7, 0.7), 
                           cex = 2,
                           cat.fontface = 4,
                           lty =2, 
                           cat.fontfamily = "sans",
                           fontfamily = "sans",
                           filename = NULL)

plasma_pro <- venn.diagram(list(plasma_TPP_prot, plasma_PP_prot),
                           category.names = c("TPP_Proteins" , "PP_Proteins"),
                           resolution = 300,
                           height = 2000,
                           width = 2000,
                           cat.default.pos = "outer",
                           lwd = 1,
                           fill = c("#FFA500", "#ff3232"),
                           cat.pos = c(-170, 170),
                           alpha = c(0.7, 0.7), 
                           cex = 2,
                           cat.fontface = 4,
                           lty =2, 
                           cat.fontfamily = "sans",
                           fontfamily = "sans",
                           filename = NULL)

## Cell PP and Cell TPP

cell_pep <- venn.diagram(list(cell_TPP_pep, cell_PP_pep),
                         category.names = c("TPP_Peptides" , "PP_Peptides"),
                         resolution = 300,
                         height = 2000,
                         width = 2000,
                         cat.default.pos = "outer",
                         lwd = 1,
                         fill = c("#FFA500", "#ff3232"),
                         cat.pos = c(-170, 170),
                         alpha = c(0.7, 0.7), 
                         cex = 2,
                         cat.fontface = 4,
                         lty =2, 
                         cat.fontfamily = "sans",
                         fontfamily = "sans",
                         filename = NULL)

cell_pro <- venn.diagram(list(cell_TPP_prot, cell_PP_prot),
                         category.names = c("TPP_Proteins" , "PP_Proteins"),
                         resolution = 300,
                         height = 2000,
                         width = 2000,
                         cat.default.pos = "outer",
                         lwd = 1,
                         fill = c("#FFA500", "#ff3232"),
                         cat.pos = c(-170, 170),
                         alpha = c(0.7, 0.7), 
                         cex = 2,
                         cat.fontface = 4,
                         lty =2, 
                         cat.fontfamily = "sans",
                         fontfamily = "sans",
                         filename = NULL)

## Plasma-Cell PP and Plasma-Cell TPP

plasma_cell_pep <- venn.diagram(list(plasma_cell_TPP_pep, plasma_cell_PP_pep),
                                category.names = c("TPP_Peptides" , "PP_Peptides"),
                                resolution = 300,
                                height = 2000,
                                width = 2000,
                                cat.default.pos = "outer",
                                lwd = 1,
                                fill = c("#FFA500", "#ff3232"),
                                cat.pos = c(-170, 170),
                                alpha = c(0.7, 0.7), 
                                cex = 2,
                                cat.fontface = 4,
                                lty =2, 
                                cat.fontfamily = "sans",
                                fontfamily = "sans",
                                filename = NULL)

plasma_cell_pro <- venn.diagram(list(plasma_cell_TPP_prot, plasma_cell_PP_prot),
                                category.names = c("TPP_Proteins" , "PP_Proteins"),
                                resolution = 300,
                                height = 2000,
                                width = 2000,
                                cat.default.pos = "outer",
                                lwd = 1,
                                fill = c("#FFA500", "#ff3232"),
                                cat.pos = c(-170, 170),
                                alpha = c(0.7, 0.7), 
                                cex = 2,
                                cat.fontface = 4,
                                lty =2, 
                                cat.fontfamily = "sans",
                                fontfamily = "sans",
                                filename = NULL)


plasma_pro2 <- gTree(children=plasma_pro)
plasma_pep2 <- gTree(children=plasma_pep)
cell_pro2 <- gTree(children=cell_pro)
cell_pep2 <- gTree(children=cell_pep)
plasma_cell_pro2 <- gTree(children=plasma_cell_pro)
plasma_cell_pep2 <- gTree(children=plasma_cell_pep)


plasma_pro2 <- as_ggplot(plasma_pro2)

plasma_pep2 <- as_ggplot(plasma_pep2)

cell_pro2 <- as_ggplot(cell_pro2)

cell_pep2 <- as_ggplot(cell_pep2)

plasma_cell_pro2 <- as_ggplot(plasma_cell_pro2)

plasma_cell_pep2 <- as_ggplot(plasma_cell_pep2)


venn1 <- ggarrange(plasma_pro2, cell_pro2, plasma_cell_pro2, ncol = 3, nrow = 1, labels = c("Plasma", "Cell", "Plasma_Cell"),
                   font.label = list(size = 14))
venn2 <- ggarrange(plasma_pep2, cell_pep2, plasma_cell_pep2, ncol = 3, nrow = 1, labels = c("Plasma", "Cell", "Plasma_Cell"),
                   font.label = list(size = 14))

tiff(".../Figures/extractions_venn_diagram.tiff", units = "in", width = 10, height = 7, res = 300)

ggarrange(venn1, venn2, ncol = 1, nrow = 2, labels = "AUTO", font.label = list(size = 26), 
          hjust = 0, vjust = 1)

dev.off()

#### Extraction Venn Diagram Plasma, Cell and Plasma_Cell

##### Venn diagram

pc_p_pp_pep <- venn.diagram(list(plasma_cell_PP_pep, plasma_PP_pep),
                            category.names = c("Plasma_Cell_Peptides" , "Plasma_Peptides"),
                            resolution = 300,
                            height = 2000,
                            width = 2000,
                            cat.default.pos = "outer",
                            lwd = 1,
                            fill = c("#5db8e0", "#158abd"),
                            cat.pos = c(-170, 170),
                            alpha = c(0.7, 0.7), 
                            cex = 2,
                            cat.fontface = 4,
                            lty =2, 
                            cat.fontfamily = "sans",
                            fontfamily = "sans",
                            filename = NULL)

pc_p_pp_prot <- venn.diagram(list(plasma_cell_PP_prot, plasma_PP_prot),
                             category.names = c("Plasma_Cell_Proteins" , "Plasma_Proteins"),
                             resolution = 300,
                             height = 2000,
                             width = 2000,
                             cat.default.pos = "outer",
                             lwd = 1,
                             fill = c("#5db8e0", "#158abd"),
                             cat.pos = c(-170, 170),
                             alpha = c(0.7, 0.7), 
                             cex = 2,
                             cat.fontface = 4,
                             lty =2, 
                             cat.fontfamily = "sans",
                             fontfamily = "sans",
                             filename = NULL)


pc_c_pp_pep <- venn.diagram(list(plasma_cell_PP_pep, cell_PP_pep),
                            category.names = c("Plasma_Cell_Peptides" , "Cell_Peptides"),
                            resolution = 300,
                            height = 2000,
                            width = 2000,
                            cat.default.pos = "outer",
                            lwd = 1,
                            fill = c("#5db8e0", "#158abd"),
                            cat.pos = c(-170, 165),
                            alpha = c(0.7, 0.7), 
                            cex = 2,
                            cat.fontface = 4,
                            lty =2, 
                            cat.fontfamily = "sans",
                            fontfamily = "sans",
                            filename = NULL)

pc_c_pp_prot <- venn.diagram(list(plasma_cell_PP_prot, cell_PP_prot),
                             category.names = c("Plasma_Cell_Proteins" , "Cell_Proteins"),
                             resolution = 300,
                             height = 2000,
                             width = 2000,
                             cat.default.pos = "outer",
                             lwd = 1,
                             fill = c("#5db8e0", "#158abd"),
                             cat.pos = c(-170, 165),
                             alpha = c(0.7, 0.7), 
                             cex = 2,
                             cat.fontface = 4,
                             lty =2, 
                             cat.fontfamily = "sans",
                             fontfamily = "sans",
                             filename = NULL)

## Plasma-Cell PP and Plasma-Cell TPP

pc_p_tpp_pep <- venn.diagram(list(plasma_cell_TPP_pep, plasma_TPP_pep),
                             category.names = c("Plasma_Cell_Peptides" , "Plasma_Peptides"),
                             resolution = 300,
                             height = 2000,
                             width = 2000,
                             cat.default.pos = "outer",
                             lwd = 1,
                             fill = c("#5db8e0", "#158abd"),
                             cat.pos = c(-170, 170),
                             alpha = c(0.7, 0.7), 
                             cex = 2,
                             cat.fontface = 4,
                             lty =2, 
                             cat.fontfamily = "sans",
                             fontfamily = "sans",
                             filename = NULL)

pc_p_tpp_prot <- venn.diagram(list(plasma_cell_TPP_prot, plasma_TPP_prot),
                              category.names = c("Plasma_Cell_Proteins" , "Plasma_Proteins"),
                              resolution = 300,
                              height = 2000,
                              width = 2000,
                              cat.default.pos = "outer",
                              lwd = 1,
                              fill = c("#5db8e0", "#158abd"),
                              cat.pos = c(-170, 170),
                              alpha = c(0.7, 0.7), 
                              cex = 2,
                              cat.fontface = 4,
                              lty =2, 
                              cat.fontfamily = "sans",
                              fontfamily = "sans",
                              filename = NULL)


pc_c_tpp_pep <- venn.diagram(list(plasma_cell_TPP_pep, cell_TPP_pep),
                             category.names = c("Plasma_Cell_Peptides" , "Cell_Peptides"),
                             resolution = 300,
                             height = 2000,
                             width = 2000,
                             cat.default.pos = "outer",
                             lwd = 1,
                             fill = c("#5db8e0", "#158abd"),
                             cat.pos = c(-170, 165),
                             alpha = c(0.7, 0.7), 
                             cex = 2,
                             cat.fontface = 4,
                             lty =2, 
                             cat.fontfamily = "sans",
                             fontfamily = "sans",
                             filename = NULL)

pc_c_tpp_prot <- venn.diagram(list(plasma_cell_TPP_prot, cell_TPP_prot),
                              category.names = c("Plasma_Cell_Proteins" , "Cell_Proteins"),
                              resolution = 300,
                              height = 2000,
                              width = 2000,
                              cat.default.pos = "outer",
                              lwd = 1,
                              fill = c("#5db8e0", "#158abd"),
                              cat.pos = c(-170, 155),
                              alpha = c(0.7, 0.7), 
                              cex = 2,
                              cat.fontface = 4,
                              lty =2, 
                              cat.fontfamily = "sans",
                              fontfamily = "sans",
                              filename = NULL)


pc_p_pp_pep2 <- gTree(children=pc_p_pp_pep)
pc_p_pp_prot2 <- gTree(children=pc_p_pp_prot)
pc_c_pp_pep2 <- gTree(children=pc_c_pp_pep)
pc_c_pp_prot2 <- gTree(children=pc_c_pp_prot)
pc_p_tpp_pep2 <- gTree(children=pc_p_tpp_pep)
pc_p_tpp_prot2 <- gTree(children=pc_p_tpp_prot)
pc_c_tpp_pep2 <- gTree(children=pc_c_tpp_pep)
pc_c_tpp_prot2 <- gTree(children=pc_c_tpp_prot)


pc_p_pp_pep2 <- as_ggplot(pc_p_pp_pep2)

pc_p_pp_prot2 <- as_ggplot(pc_p_pp_prot2)

pc_c_pp_pep2 <- as_ggplot(pc_c_pp_pep2)

pc_c_pp_prot2 <- as_ggplot(pc_c_pp_prot2)

pc_p_tpp_pep2 <- as_ggplot(pc_p_tpp_pep2)

pc_p_tpp_prot2 <- as_ggplot(pc_p_tpp_prot2)

pc_c_tpp_pep2 <- as_ggplot(pc_c_tpp_pep2)

pc_c_tpp_prot2 <- as_ggplot(pc_c_tpp_prot2)


venn1 <- ggarrange(pc_p_pp_pep2, pc_p_pp_prot2, ncol = 2, nrow = 1, labels = c(),
                   font.label = list(size = 14))
venn2 <- ggarrange(pc_c_pp_pep2, pc_c_pp_prot2, ncol = 2, nrow = 1, labels = c(),
                   font.label = list(size = 14))
venn3 <- ggarrange(pc_p_tpp_pep2, pc_p_tpp_prot2, ncol = 2, nrow = 1, labels = c(),
                   font.label = list(size = 14))
venn4 <- ggarrange(pc_c_tpp_pep2, pc_c_tpp_prot2, ncol = 2, nrow = 1, labels = c(),
                   font.label = list(size = 14))

tiff(".../Figures/plasma_cell_PPcompare_venn_diagram.tiff", units = "in", width = 10, height = 8, res = 300)

ggarrange(venn1, venn2, ncol = 1, nrow = 2, labels = c("A. PP", "B. PP"), font.label = list(size = 26), 
          hjust = 0, vjust = 1) +labs(title = "PP")

dev.off()

tiff(".../Figures/plasma_cell_TPPcompare_venn_diagram.tiff", units = "in", width = 10, height = 8, res = 300)

ggarrange(venn3, venn4, ncol = 1, nrow = 2, labels = c("A. TPP", "B. TPP"), font.label = list(size = 26), 
          hjust = 0, vjust = 1)

dev.off()


#############################################################################################
######################################################################################

######################################################################################
#############################################################################################


##### Number of identified proteins and peptides in each of the libraries


lib_plasma_pp <- read.csv(".../Plasma_Library_PP.csv", stringsAsFactors = F)
lib_plasma_tpp <- read.csv(".../Plasma_Library_TPP.csv", stringsAsFactors = F)

lib_cell_pp <- read.csv(".../Cell_Calibrated_Library_PP.csv", stringsAsFactors = F)
lib_cell_tpp <- read.csv(".../Cell_Calibrated_Library_TPP.csv", stringsAsFactors = F)

lib_plasma_cell_pp <- read.csv(".../Plasma_Cell_Library_PP.csv", stringsAsFactors = F)
lib_plasma_cell_tpp <- read.csv(".../Plasma_Cell_Library_TPP.csv", stringsAsFactors = F)


plasma_lib_PP_prot <- lib_plasma_pp[!duplicated(lib_plasma_pp$ProteinName),"UniprotID"]
cell_lib_PP_prot <- lib_cell_pp[!duplicated(lib_cell_pp$ProteinName), "UniprotID"]
plasma_cell_lib_PP_prot <- lib_plasma_cell_pp[!duplicated(lib_plasma_cell_pp$ProteinName), "UniprotID"]


plasma_lib_PP_pep <- lib_plasma_pp[!duplicated(lib_plasma_pp$ModificationSequence), "ModificationSequence"]
cell_lib_PP_pep <- lib_cell_pp[!duplicated(lib_cell_pp$ModificationSequence), "ModificationSequence"]
plasma_cell_lib_PP_pep <- lib_plasma_cell_pp[!duplicated(lib_plasma_cell_pp$ModificationSequence), "ModificationSequence"]


plasma_lib_TPP_prot <- lib_plasma_tpp[!duplicated(lib_plasma_tpp$ProteinName), "UniprotID"]
cell_lib_TPP_prot <- lib_cell_tpp[!duplicated(lib_cell_tpp$ProteinName), "UniprotID"]
plasma_cell_lib_TPP_prot <- lib_plasma_cell_tpp[!duplicated(lib_plasma_cell_tpp$ProteinName), "UniprotID"]


plasma_lib_TPP_pep <- lib_plasma_tpp[!duplicated(lib_plasma_tpp$ModificationSequence), "ModificationSequence"]
cell_lib_TPP_pep <- lib_cell_tpp[!duplicated(lib_cell_tpp$ModificationSequence), "ModificationSequence"]
plasma_cell_lib_TPP_pep <- lib_plasma_cell_tpp[!duplicated(lib_plasma_cell_tpp$ModificationSequence), "ModificationSequence"]

extractions_lib_data_actual <- data.frame("Libraries" = c("aPlasma", "bCell", "cPlas-Cell", "aPlasma", "bCell", "cPlas-Cell"), 
                                          "Proteins" = c(length(plasma_lib_TPP_prot) , length(cell_lib_TPP_prot) , length(plasma_cell_lib_TPP_prot) ,
                                                         length(plasma_lib_PP_prot) , length(cell_lib_PP_prot) , length(plasma_cell_lib_PP_prot) ),
                                          "Peptides" = c(length(plasma_lib_TPP_pep) , length(cell_lib_TPP_pep) , length(plasma_cell_lib_TPP_pep) ,
                                                         length(plasma_lib_PP_pep) , length(cell_lib_PP_pep) , length(plasma_cell_lib_PP_pep) ),
                                          "Type" = c("aTPP", "aTPP", "aTPP", "bPP", "bPP", "bPP"))

extractions_lib_data <- data.frame("Libraries" = c("aPlasma", "bCell", "cPlas-Cell", "aPlasma", "bCell", "cPlas-Cell"), 
                                   "Proteins" = c(length(plasma_lib_TPP_prot)/1000, length(cell_lib_TPP_prot)/1000, length(plasma_cell_lib_TPP_prot)/1000,
                                                  length(plasma_lib_PP_prot)/1000, length(cell_lib_PP_prot)/1000, length(plasma_cell_lib_PP_prot)/1000),
                                   "Peptides" = c(length(plasma_lib_TPP_pep)/1000, length(cell_lib_TPP_pep)/1000, length(plasma_cell_lib_TPP_pep)/1000,
                                                  length(plasma_lib_PP_pep)/1000, length(cell_lib_PP_pep)/1000, length(plasma_cell_lib_PP_pep)/1000),
                                   "Type" = c("aTPP", "aTPP", "aTPP", "bPP", "bPP", "bPP"))


tiff(".../Figures/libraries_numbers.tiff", units = "in", width = 6, height = 4, res = 300)

prot <- ggplot(data=extractions_lib_data, aes(x = Libraries, y = Proteins, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.9, position=position_dodge()) +
  # geom_text(aes(label=Proteins), vjust=0, color="black", fontface = "bold",
  # position = position_dodge(0.9), size=5) +
  scale_y_continuous(breaks=seq(0, 6, 1))+
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(labels = c("TPP", "PP"), values=c("#FFA500", "#ff3232"))+
  scale_x_discrete(labels = c("Plasma", "Cell", "Plasma-Cell"))+
  labs(y = deparse(bquote(Proteins(10^3)))) +
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
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

pep <- ggplot(data=extractions_lib_data, aes(x = Libraries, y = Peptides, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.9, position=position_dodge()) +
  # geom_text(aes(label=Peptides), vjust=0, color="black", fontface = "bold",
  #           position = position_dodge(0.9), size=5) +
  scale_y_continuous(breaks=seq(0, 45, 10))+
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(labels = c("TPP", "PP"), values=c("#FFA500", "#ff3232"))+
  scale_x_discrete(labels = c("Plasma", "Cell", "Plasma-Cell"))+
  labs(y = deparse(bquote(Peptides(10^3)))) +
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
        legend.title = element_text(size = 14, face = "bold"),
        legend.justification = "center",
        legend.position = "right") 

ggarrange(prot, pep, ncol = 2, nrow = 1, common.legend = TRUE, labels = "AUTO", font.label = list(size = 14))

dev.off()

#############################################################################################
######################################################################################

#### Libraries Venn Diagram

## Plasma PP and Plasma TPP

plasma_pep <- venn.diagram(list(plasma_lib_PP_pep, plasma_lib_TPP_pep),
                           category.names = c("PP_Peptides" , "TPP_Peptides"),
                           resolution = 300,
                           height = 2000,
                           width = 2000,
                           cat.default.pos = "outer",
                           lwd = 1,
                           fill = c("#FFA500", "#ff3232"),
                           cat.pos = c(-170, 170),
                           alpha = c(0.7, 0.7), 
                           cex = 3,
                           cat.cex = 2,
                           cat.fontface = 4,
                           lty =2, 
                           cat.fontfamily = "sans",
                           fontfamily = "sans",
                           filename = NULL)

plasma_pro <- venn.diagram(list(plasma_lib_PP_prot, plasma_lib_TPP_prot),
                           category.names = c("PP_Proteins" , "TPP_Proteins"),
                           resolution = 300,
                           height = 1000,
                           width = 1000,
                           cat.default.pos = "outer",
                           lwd = 1,
                           fill = c("#FFA500", "#ff3232"),
                           cat.pos = c(-170, 170),
                           alpha = c(0.7, 0.7), 
                           cex = 3,
                           cat.cex = 2,
                           cat.fontface = 4,
                           lty =2, 
                           cat.fontfamily = "sans",
                           fontfamily = "sans",
                           filename = NULL)

## Cell PP and Cell TPP

cell_pep <- venn.diagram(list(cell_lib_PP_pep, cell_lib_TPP_pep),
                         category.names = c("PP_Peptides" , "TPP_Peptides"),
                         resolution = 300,
                         height = 1000,
                         width = 1000,
                         cat.default.pos = "outer",
                         lwd = 1,
                         fill = c("#FFA500", "#ff3232"),
                         cat.pos = c(-170, 170),
                         alpha = c(0.7, 0.7), 
                         cex = 3,
                         cat.cex = 2,
                         cat.fontface = 4,
                         lty =2, 
                         cat.fontfamily = "sans",
                         fontfamily = "sans",
                         filename = NULL)

cell_pro <- venn.diagram(list(cell_lib_PP_prot, cell_lib_TPP_prot),
                         category.names = c("PP_Proteins" , "TPP_Proteins"),
                         resolution = 300,
                         height = 1000,
                         width = 1000,
                         cat.default.pos = "outer",
                         lwd = 1,
                         fill = c("#FFA500", "#ff3232"),
                         cat.pos = c(-170, 170),
                         alpha = c(0.7, 0.7), 
                         cex = 3,
                         cat.cex = 2,
                         cat.fontface = 4,
                         lty =2, 
                         cat.fontfamily = "sans",
                         fontfamily = "sans",
                         filename = NULL)

## Plasma-Cell PP and Plasma-Cell TPP

plasma_cell_pep <- venn.diagram(list(plasma_cell_lib_PP_pep, plasma_cell_lib_TPP_pep),
                                category.names = c("PP_Peptides" , "TPP_Peptides"),
                                resolution = 300,
                                height = 1000,
                                width = 1000,
                                cat.default.pos = "outer",
                                lwd = 1,
                                fill = c("#FFA500", "#ff3232"),
                                cat.pos = c(-170, 170),
                                alpha = c(0.7, 0.7), 
                                cex = 3,
                                cat.cex = 2,
                                cat.fontface = 4,
                                lty =2, 
                                cat.fontfamily = "sans",
                                fontfamily = "sans",
                                filename = NULL)

plasma_cell_pro <- venn.diagram(list(plasma_cell_lib_PP_prot, plasma_cell_lib_TPP_prot),
                                category.names = c("PP_Proteins" , "TPP_Proteins"),
                                resolution = 300,
                                height = 1000,
                                width = 1000,
                                cat.default.pos = "outer",
                                lwd = 1,
                                fill = c("#FFA500", "#ff3232"),
                                cat.pos = c(-170, 170),
                                alpha = c(0.7, 0.7), 
                                cex = 3,
                                cat.cex = 2,
                                cat.fontface = 4,
                                lty =2, 
                                cat.fontfamily = "sans",
                                fontfamily = "sans",
                                filename = NULL)


plasma_pro2 <- gTree(children=plasma_pro)
plasma_pep2 <- gTree(children=plasma_pep)
cell_pro2 <- gTree(children=cell_pro)
cell_pep2 <- gTree(children=cell_pep)
plasma_cell_pro2 <- gTree(children=plasma_cell_pro)
plasma_cell_pep2 <- gTree(children=plasma_cell_pep)


plasma_pro2 <- as_ggplot(plasma_pro2)

plasma_pep2 <- as_ggplot(plasma_pep2)

cell_pro2 <- as_ggplot(cell_pro2)

cell_pep2 <- as_ggplot(cell_pep2)

plasma_cell_pro2 <- as_ggplot(plasma_cell_pro2)

plasma_cell_pep2 <- as_ggplot(plasma_cell_pep2)


venn1 <- ggarrange(plasma_pro2, cell_pro2, plasma_cell_pro2, ncol = 3, nrow = 1, labels = c("Plasma", "Cell", "Plasma_Cell"),
                   font.label = list(size = 22), hjust = c(-0.5, 0.5, 0.5))
venn2 <- ggarrange(plasma_pep2, cell_pep2, plasma_cell_pep2, ncol = 3, nrow = 1, labels = c("Plasma", "Cell", "Plasma_Cell"),
                   font.label = list(size = 22), hjust = c(-0.5, 0.5, 0.5))

tiff(".../Figures/libraries_venn_diagram.tiff", units = "in", width = 23, height = 15, res = 300)

ggarrange(venn1, venn2, ncol = 1, nrow = 2, labels = "AUTO", font.label = list(size = 26), 
          hjust = 0, vjust = 1)

dev.off()



################################################################################################33
#########################################################################################333
#########################################################################################333

### Libraries frequency distribution

plasma_lib_PP_prot <- lib_plasma_pp[!duplicated(lib_plasma_pp$ModificationSequence),]
plasma_ppx <- aggregate(x = plasma_lib_PP_prot$ModificationSequence, by = list(Protein = plasma_lib_PP_prot$UniprotID), FUN = function(x) {length(x)})

plasma_ppx <- plasma_ppx %>%
  filter(x <= 75)

plasma_2pep <- plasma_ppx %>%
  filter(x >= 2)

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
        axis.text.x = element_text(colour = "black", size = 15, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_lib_TPP_prot <- lib_plasma_tpp[!duplicated(lib_plasma_tpp$ModificationSequence),]
plasma_tppx <- aggregate(x = plasma_lib_TPP_prot$ModificationSequence, by = list(Protein = plasma_lib_TPP_prot$UniprotID), FUN = function(x) {length(x)})

plasma_tppx <- plasma_tppx %>%
  filter(x <= 75)

plasma_2pep <- plasma_tppx %>%
  filter(x >= 2)

p2 <- ggplot(data = plasma_tppx, aes(x = x, fill = )) +
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
        axis.text.x = element_text(colour = "black", size = 15, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

cell_lib_PP_prot <- lib_cell_pp[!duplicated(lib_cell_pp$ModificationSequence),]
cell_ppx <- aggregate(x = cell_lib_PP_prot$ModificationSequence, by = list(Protein = cell_lib_PP_prot$UniprotID), FUN = function(x) {length(x)})

cell_ppx <- cell_ppx %>%
  filter(x <= 75)

cell_2pep <- cell_ppx %>%
  filter(x >= 2)

p3 <- ggplot(data = cell_ppx, aes(x = x, fill = )) +
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
        axis.text.x = element_text(colour = "black", size = 15, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

cell_lib_TPP_prot <- lib_cell_tpp[!duplicated(lib_cell_tpp$ModificationSequence),]
cell_tppx <- aggregate(x = cell_lib_TPP_prot$ModificationSequence, by = list(Protein = cell_lib_TPP_prot$UniprotID), FUN = function(x) {length(x)})

cell_tppx <- cell_tppx %>%
  filter(x <= 75)

cell_2pep <- cell_tppx %>%
  filter(x >= 2)

p4 <- ggplot(data = cell_tppx, aes(x = x, fill = )) +
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
        axis.text.x = element_text(colour = "black", size = 15, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


plasma_cell_lib_PP_prot <- lib_plasma_cell_pp[!duplicated(lib_plasma_cell_pp$ModificationSequence),]
plasma_cell_ppx <- aggregate(x = plasma_cell_lib_PP_prot$ModificationSequence, by = list(Protein = plasma_cell_lib_PP_prot$UniprotID), FUN = function(x) {length(x)})

plasma_cell_ppx <- plasma_cell_ppx %>%
  filter(x <= 75)

plasma_cell_2pep <- plasma_cell_ppx %>%
  filter(x >= 2)

p5 <- ggplot(data = plasma_cell_ppx, aes(x = x, fill = )) +
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
        axis.text.x = element_text(colour = "black", size = 15, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")

plasma_cell_lib_TPP_prot <- lib_plasma_cell_tpp[!duplicated(lib_plasma_cell_tpp$ModificationSequence),]
plasma_cell_tppx <- aggregate(x = plasma_cell_lib_TPP_prot$ModificationSequence, by = list(Protein = plasma_cell_lib_TPP_prot$UniprotID), FUN = function(x) {length(x)})

plasma_cell_tppx <- plasma_cell_tppx %>%
  filter(x <= 75)

plasma_cell_2pep <- plasma_cell_tppx %>%
  filter(x >= 2)

p6 <- ggplot(data = plasma_cell_tppx, aes(x = x, fill = )) +
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
        axis.text.x = element_text(colour = "black", size = 15, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold", vjust = 0.5),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.justification = "center",
        legend.position = "right")


tiff(".../Figures/frequency_distribution_Lib.tiff", units = "in", width = 12, height = 7, res = 300)
fdpp <- ggarrange(p1, p3, p5, ncol = 3, nrow = 1, legend = "right")
fdtpp <- ggarrange(p2, p4, p6, ncol = 3, nrow = 1, legend = "right")
ggarrange(fdtpp, fdpp, ncol = 1, nrow = 2, legend = "right", labels = "AUTO", font.label = list(size = 20))

dev.off()

#############################################################################################
######################################################################################
#############################################################################################

