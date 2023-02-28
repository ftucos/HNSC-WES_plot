wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)

source("src/reorder_patients.R")

drop.na <- function(vector) {vector[!is.na(vector)]}

# configuration


# to do
# [] distinguish HPV from non HPV
# [] keep only Oral cavity and tongue



# load data ------------------------------------
# https://gdc.cancer.gov/about-data/publications/pancanatlas
nurd <- c("CHD3", "CHD4", "CHD5", "HDAC1", "HDAC2", "MTA1", "MTA2", "MTA3", "MBD2", "GATAD2A", "GATAD2B", "PWWP2A", "PWWP2B", "MBD2", "MBD3", "RBBP4", "RBBP7", "CDK2AP1")
HN_driver <- c("CDKN2A", "FAT1", "TP53", "CASP8", "NOTCH1", "PIK3CA", "AJUBA", "KMT2D", "KRT5", "TGFBR2", "CTCF", "HLA-A", "HRAS", "EPHA2", "H1-5", "HRNR", "NECAB1", "PSIP1", "FOSL2")

# limited to the initiial nature cohort
HPV <- read_excel("data/nature14129-s2/1.2.xlsx") %>%
  select(sample = Barcode, Tumor_Site, Final_HPV_Status)


site <- data.table::fread("data/hnsc_tcga/data_bcr_clinical_data_patient.txt", skip=4) %>%
  select(sample=PATIENT_ID, Site=PRIMARY_SITE_PATIENT)

df.raw <- data.table::fread("data/manual_cBioportal_download/CNA_and_MUT-HNSC_drivers_and_NURD.tsv")

CNA <- df.raw %>%
  filter(track_type %in% c("CNA")) %>%
  gather(key="sample", value="CNA", -c(1:2)) %>%
  select(-track_type)

mut <- df.raw %>%
    filter(track_type %in% c("MUTATIONS")) %>%
  gather(key="sample", value="SNV", -c(1:2)) %>%
  select(-track_type)

# Preprocessing ----------------------------

df <- full_join(CNA, mut) %>%
  mutate(group = ifelse(track_name %in% nurd, "NuRD Complex", "Recurrent mutations in OSCC"),
         CNA = ifelse(CNA == "", NA, CNA),
         # remove het loss
         # CNA = ifelse(CNA %in% c("Shallow Deletion", "hetloss_rec"), NA, CNA),
         # homdel_rec is the same as Deep Deletion but for genes for wich it's recognised as recurrent (Putative Driver) event.
         CNA = recode(CNA,   "amp_rec" = "Amplification","homdel_rec" =  "Deep Deletion", "Deep Deletion" = "Deep Deletion", "Shallow Deletion" = "Shallow Deletion", "hetloss_rec" = "Shallow Deletion"),
         SNV = recode(SNV, "splice_rec" = "Splice Mutation (putative driver)"),
         SNV = factor(SNV, levels= c("Truncating mutation (putative driver)", "Inframe Mutation (putative driver)", "Splice Mutation (putative driver)", "Missense Mutation (putative driver)", "Truncating mutation (putative passenger)", "Inframe Mutation (putative passenger)", "Missense Mutation (putative passenger)")),
         SNV_simp = str_remove(SNV, "\\s\\(.*\\)"),
         #mutations = paste0(CNA, SNV_simp, collapse = ","),
         ) %>%
  left_join(site) %>%
  rowwise() %>%
  mutate(mutations = list(c(as.character(CNA), as.character(SNV_simp)) %>% drop.na() %>% paste0(collapse=";")))

# plot configuration -----------------------

col <- c("Deep Deletion" = "steelblue", "Amplification" = "firebrick", "Shallow Deletion" = "steelblue",
         "Truncating mutation" = "#2D3E50", "Inframe Mutation" = "#8F44AD", "Splice Mutation" = "#F49C12", "Missense Mutation" = "#25AF60")

h_marg <- 0.15
v_marg <- 3


# PLOT: HN driver with hetloss  --------------------------
mat.hn <- df %>%
  filter(track_name %in% HN_driver) %>%
  select(track_name, sample, mutations) %>%
  spread(key="sample", value="mutations") %>%
  column_to_rownames("track_name") %>%
  as.matrix()

pdf("results/HN_drivers-TCGA-with_hetloss.pdf", width=10, height=nrow(mat.hn)/3)
oncoplot <- oncoPrint(mat.hn, column_order = reorder_patients(mat.hn), pct_digits = 1,
          pct_side = "right", row_names_side = "left",
          alter_fun = list(
            background = alter_graphic("rect", fill = "#F1F1F1", horiz_margin = unit(0, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Amplification` = alter_graphic("rect", fill = col["Amplification"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], alpha=0.4, horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Truncating mutation` = alter_graphic("rect", height = 0.5, fill = col["Truncating mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Inframe Mutation` = alter_graphic("rect", height = 0.5, fill = col["Inframe Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Splice Mutation` = alter_graphic("rect", height = 0.5, fill = col["Splice Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Missense Mutation` = alter_graphic("rect", height = 0.5, fill = col["Missense Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt"))
          ))
draw(oncoplot)
dev.off()

# PLOT: HN driver without hetloss  --------------------------
mat.hn <- df %>%
  mutate(mutations = str_remove(mutations, "Shallow Deletion;?"))%>%
  filter(track_name %in% HN_driver) %>%
  select(track_name, sample, mutations) %>%
  spread(key="sample", value="mutations") %>%
  column_to_rownames("track_name") %>%
  as.matrix()

pdf("results/HN_drivers-TCGA.pdf", width=10, height=nrow(mat.hn)/3)
oncoplot <- oncoPrint(mat.hn, column_order = reorder_patients(mat.hn), pct_digits = 1,
                      pct_side = "right", row_names_side = "left",
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "#F1F1F1", horiz_margin = unit(0, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Amplification` = alter_graphic("rect", fill = col["Amplification"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], alpha=0.4, horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Truncating mutation` = alter_graphic("rect", height = 0.5, fill = col["Truncating mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Inframe Mutation` = alter_graphic("rect", height = 0.5, fill = col["Inframe Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Splice Mutation` = alter_graphic("rect", height = 0.5, fill = col["Splice Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Missense Mutation` = alter_graphic("rect", height = 0.5, fill = col["Missense Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt"))
                      ))
draw(oncoplot)
dev.off()


# PLOT: NURD driver with hetloss  --------------------------
mat.nurd <- df %>%
  filter(track_name %in% nurd) %>%
  select(track_name, sample, mutations) %>%
  spread(key="sample", value="mutations") %>%
  column_to_rownames("track_name") %>%
  as.matrix()

pdf("results/NuRD-TCGA-with_hetloss.pdf", width=10, height=nrow(mat.nurd)/3)
oncoplot <- oncoPrint(mat.nurd, column_order = reorder_patients(mat.nurd), pct_digits = 1,
                      pct_side = "right", row_names_side = "left",
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "#F1F1F1", horiz_margin = unit(0, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Amplification` = alter_graphic("rect", fill = col["Amplification"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], alpha=0.4, horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Truncating mutation` = alter_graphic("rect", height = 0.5, fill = col["Truncating mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Inframe Mutation` = alter_graphic("rect", height = 0.5, fill = col["Inframe Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Splice Mutation` = alter_graphic("rect", height = 0.5, fill = col["Splice Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Missense Mutation` = alter_graphic("rect", height = 0.5, fill = col["Missense Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt"))
                      ))
draw(oncoplot)
dev.off()

# PLOT: NURD driver without hetloss  --------------------------
mat.nurd <- df %>%
  mutate(mutations = str_remove(mutations, "Shallow Deletion;?"))%>%
  filter(track_name %in% nurd) %>%
  select(track_name, sample, mutations) %>%
  spread(key="sample", value="mutations") %>%
  column_to_rownames("track_name") %>%
  as.matrix()

pdf("results/NuRD-TCGA.pdf", width=10, height=nrow(mat.nurd)/3)
oncoplot <- oncoPrint(mat.nurd, column_order = reorder_patients(mat.nurd),
                      pct_side = "right", row_names_side = "left", pct_digits = 1,
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "#F1F1F1", horiz_margin = unit(0, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Amplification` = alter_graphic("rect", fill = col["Amplification"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], alpha=0.4, horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Truncating mutation` = alter_graphic("rect", height = 0.5, fill = col["Truncating mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Inframe Mutation` = alter_graphic("rect", height = 0.5, fill = col["Inframe Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Splice Mutation` = alter_graphic("rect", height = 0.5, fill = col["Splice Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
                        `Missense Mutation` = alter_graphic("rect", height = 0.5, fill = col["Missense Mutation"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt"))
                      ))
draw(oncoplot)
dev.off()

