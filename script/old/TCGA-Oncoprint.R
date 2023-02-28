wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)

source("script/src/reorder_patient.R")

drop.na <- function(vector) {vector[!is.na(vector)]}

# configuration


# to do
# [] distinguish HPV from non HPV
# [] keep only Oral cavity and tongue



# load data ------------------------------------
nurd <- c("CHD3", "CHD4", "CHD5", "HDAC1", "HDAC2", "MTA1", "MTA2", "MTA3", "GATAD2A", "GATAD2B", "PWWP2A", "PWWP2B", "MBD2", "MBD3", "RBBP4", "RBBP7", "CDK2AP1")
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
         SNV = recode(SNV, "splice_rec" = "Splice site (putative driver)", "splice" = "Splice site"),
         SNV_simp = str_remove(SNV, "\\s\\(.*\\)") %>% str_remove(regex("\\s?mutation\\s?", ignore_case = T)),
         #mutations = paste0(CNA, SNV_simp, collapse = ","),
         ) %>%
  left_join(site) %>%
  rowwise() %>%
  mutate(mutations = list(c(as.character(CNA), as.character(SNV_simp)) %>% drop.na() %>% paste0(collapse=";")))

# plot configuration -----------------------

col <- c("Deep Deletion" = "#1C5FA0", "Amplification" = "#A51428", "Shallow Deletion" = "#89C0D9",
         "Truncating" = "#2B2B2B", "Inframe" = "#3598DB", "Splice site" = "#FFCD00", "Missense" = "#25AF60")

# PLOT: HN driver with hetloss  --------------------------
mat.hn <- df %>%
  filter(track_name %in% HN_driver) %>%
  select(track_name, sample, mutations) %>%
  spread(key="sample", value="mutations") %>%
  column_to_rownames("track_name") %>%
  as.matrix()

h_pct = 0.7
w_pct = 0.7

pdf("results/HN_drivers-TCGA-with_hetloss.pdf", width=10, height=nrow(mat.hn)/3)
oncoplot <- oncoPrint(mat.hn, column_order = reorder_patients(mat.hn), 
          pct_side = "right", row_names_side = "left",pct_digits = 1,
          alter_fun = list(
            background = alter_graphic("rect", fill = "gray94", height = h_pct, width = 1),
            `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], height = h_pct, width = w_pct),
            `Amplification` = alter_graphic("rect", fill = col["Amplification"], height = h_pct, width = w_pct),
            `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], height = h_pct, width = w_pct),
            `Truncating` = alter_graphic("rect", fill = col["Truncating"], height = h_pct/2, width = w_pct),
            `Inframe` = alter_graphic("rect", fill = col["Inframe"],  height = h_pct/2, width = w_pct),
            `Splice site` = alter_graphic("rect", fill = col["Splice site"],   height = h_pct/2, width = w_pct),
            `Missense` = alter_graphic("rect", fill = col["Missense"],  height = h_pct/2, width = w_pct)
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
oncoplot <- oncoPrint(mat.hn, column_order = reorder_patients(mat.hn), 
                      pct_side = "right", row_names_side = "left",pct_digits = 1,
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "gray94", height = h_pct, width = 1),
                        `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], height = h_pct, width = w_pct),
                        `Amplification` = alter_graphic("rect", fill = col["Amplification"], height = h_pct, width = w_pct),
                        `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], height = h_pct, width = w_pct),
                        `Truncating` = alter_graphic("rect", fill = col["Truncating"], height = h_pct/2, width = w_pct),
                        `Inframe` = alter_graphic("rect", fill = col["Inframe"],  height = h_pct/2, width = w_pct),
                        `Splice site` = alter_graphic("rect", fill = col["Splice site"],   height = h_pct/2, width = w_pct),
                        `Missense` = alter_graphic("rect", fill = col["Missense"],  height = h_pct/2, width = w_pct)
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
oncoplot <- oncoPrint(mat.nurd, column_order = reorder_patients(mat.nurd), 
                      pct_side = "right", row_names_side = "left",pct_digits = 1,
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "gray94", height = h_pct, width = 1),
                        `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], height = h_pct, width = w_pct),
                        `Amplification` = alter_graphic("rect", fill = col["Amplification"], height = h_pct, width = w_pct),
                        `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], height = h_pct, width = w_pct),
                        `Truncating` = alter_graphic("rect", fill = col["Truncating"], height = h_pct/2, width = w_pct),
                        `Inframe` = alter_graphic("rect", fill = col["Inframe"],  height = h_pct/2, width = w_pct),
                        `Splice site` = alter_graphic("rect", fill = col["Splice site"],   height = h_pct/2, width = w_pct),
                        `Missense` = alter_graphic("rect", fill = col["Missense"],  height = h_pct/2, width = w_pct)
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
                      pct_side = "right", row_names_side = "left",pct_digits = 1,
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "gray94", height = h_pct, width = 1),
                        `Deep Deletion` = alter_graphic("rect", fill = col["Deep Deletion"], height = h_pct, width = w_pct),
                        `Amplification` = alter_graphic("rect", fill = col["Amplification"], height = h_pct, width = w_pct),
                        `Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], height = h_pct, width = w_pct),
                        `Truncating` = alter_graphic("rect", fill = col["Truncating"], height = h_pct/2, width = w_pct),
                        `Inframe` = alter_graphic("rect", fill = col["Inframe"],  height = h_pct/2, width = w_pct),
                        `Splice site` = alter_graphic("rect", fill = col["Splice site"],   height = h_pct/2, width = w_pct),
                        `Missense` = alter_graphic("rect", fill = col["Missense"],  height = h_pct/2, width = w_pct)
                      ))
draw(oncoplot)
dev.off()

