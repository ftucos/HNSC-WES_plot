wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)

source("script/src/reorder_patient.R")

drop.na <- function(vector) {vector[!is.na(vector)]}


# Load raw data --------------------
# https://gdc.cancer.gov/about-data/publications/pancanatlas
# nurd <- c("CHD3", "CHD4", "CHD5", "HDAC1", "HDAC2", "MTA1", "MTA2", "MTA3", "MBD2", "GATAD2A", "GATAD2B", "PWWP2A", "PWWP2B", "MBD2", "MBD3", "RBBP4", "RBBP7", "CDK2AP1")
# HN_driver <- c("CDKN2A", "FAT1", "TP53", "CASP8", "NOTCH1", "PIK3CA", "AJUBA", "KMT2D", "KRT5", "TGFBR2", "CTCF", "HLA-A", "HRAS", "EPHA2", "H1-5", "HRNR", "NECAB1", "PSIP1", "FOSL2")

nurd <- c("CDK2AP1")
HN_driver <- c("CDKN2A", "FAT1", "TP53", "CASP8", "NOTCH1")



SNV.raw <- data.table::fread("data/hnsc_tcga/data_mutations_mskcc.txt") %>%
  # remove duplicates due to multiple matched normal
  select(Hugo_Symbol, Consequence, Variant_Classification, Variant_Type, Tumor_Sample_Barcode, Start_Position) %>%
  distinct() %>%
  filter(!Variant_Classification %in% c("Silent", "RNA", "Intron",  "3'Flank", "3'UTR", "5'Flank", "5'UTR", "Splice_Region")) 

SNV_priority <- c("Nonsense", "Start Loss", "Frameshift", "Stop Loss", "Inframe InDel", "Missense", "Splice Site")

SNV <- SNV.raw %>%
  filter(Hugo_Symbol %in% c(nurd, HN_driver)) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, SNV = Variant_Classification) %>%
  mutate(SNV = recode(SNV,
                           "Missense_Mutation" = "Missense",
                           "Splice_Site" = "Splice Site",
                           "Frame_Shift_Del" = "Frameshift",
                           "Frame_Shift_Ins" = "Frameshift",
                           "In_Frame_Del" = "Inframe InDel",
                           "In_Frame_Ins" = "Inframe InDel",
                           "Nonstop_Mutation" = "Stop Loss",
                           "Nonsense_Mutation" = "Nonsense",
                           "Translation_Start_Site" = "Start Loss",
                           ) %>% factor(levels=SNV_priority))

# keep only higher priority mutation when concurrent are present
SNV.1 <- SNV %>% group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  arrange(SNV) %>%
  slice_head(n = 1)

# validate prioritizzation
# SNV.2 <-  SNV %>% group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarize(types = list(as.character(SNV)), count = n())
# 
# x <- left_join(SNV.1, SNV.2)

# number of patients
SNV.raw$Tumor_Sample_Barcode %>% unique() %>% length()
sample_id.SNV = SNV.raw$Tumor_Sample_Barcode %>% unique()

# CNA
CNA.raw <- data.table::fread("data/hnsc_tcga/data_CNA.txt") 

CNA.raw %>% colnames(-c(1:2)) %>% unique() %>% length()
sample_id.CNA = CNA.raw %>% colnames(-c(1:2)) %>% unique()

CNA <- CNA.raw %>%
  filter(Hugo_Symbol %in% c(nurd, HN_driver)) %>%
  select(-Entrez_Gene_Id) %>%
  gather(key="Tumor_Sample_Barcode", value = "CNA", -1) %>%
  mutate(CNA = recode(CNA, '-2' = "Deep Deletion",
                           '2' = "Amplification")) %>%
  filter(!is.na(CNA))

# should be 504 !!!!!!
complete_patients <- intersect(sample_id.SNV, sample_id.CNA)

# viral status --------------

viral <- read_tsv("data/HPV_status/viral.tsv")

viral.HNSC <- viral %>%
  filter(Study == "HNSC") %>%
  mutate(Tumor_Sample_Barcode = str_remove(SampleBarcode, "[A-Z]$"),
# cutoffs: NRPM Thresholds HPV: 10, EBV: 5, HBV; 5.  from The Immune Landscape of Cancer
# https://www-sciencedirect-com.pros1.lib.unimi.it/science/article/pii/S1074761318301213
  HPV_status = ifelse(HPV > 10, "Positive", "Negative"), 
  EBV_status = ifelse(EBV > 5, "Positive", "Negative"),
  HBV_status = ifelse(HBV > 5, "Positive", "Negative"))

# 73 patients
HPV_positive <- viral.HNSC %>% filter(HPV_status == "Positive") %>% pull("Tumor_Sample_Barcode")
write_tsv(as.data.frame(HPV_positive), "processed/HPV_positive_HNSC_to_exclude.tsv",col_names = F)
EBV_positive <- viral.HNSC %>% filter(EBV_status == "Positive") %>% pull("Tumor_Sample_Barcode")

# n=438
complete_patients.no_HPV <- complete_patients[!complete_patients %in% HPV_positive]

# merge CNA and SNVs----------------------
df <- full_join(CNA , SNV.1) %>%
  filter(Tumor_Sample_Barcode %in% complete_patients.no_HPV) %>%
  mutate(group = ifelse(Hugo_Symbol %in% nurd, "NuRD Complex", "Recurrent mutations in OSCC")) %>%
  rowwise() %>%
  mutate(mutation = list(c(as.character(CNA), as.character(SNV)) %>% drop.na() %>% paste0(collapse=";")) %>% unlist()) %>%
  select(-CNA, -SNV)

mutations.mat.partial <- df %>%
  select(-group) %>%
  spread(key="Tumor_Sample_Barcode", value="mutation") %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

# add missing patients with no mutations
no_mut_patients <- complete_patients.no_HPV[!complete_patients.no_HPV %in% colnames(mutations.mat.partial)]
no_mut_matrix <- matrix(nrow = nrow(mutations.mat.partial), ncol = length(no_mut_patients),
                        dimnames = list(rownames(mutations.mat.partial), no_mut_patients))

mutations.mat <- cbind(mutations.mat.partial, no_mut_matrix)
  
  
# replace NA with ""
mutations.mat[is.na(mutations.mat)] <- ""


# plot configuration -----------------------

col <- c("Deep Deletion" = "#1C5FA0", "Amplification" = "#A51428", "Shallow Deletion" = "#89C0D9",
         "Nonsense" = "#2B2B2B", "Frameshift" = "#755EC5", "Inframe InDel" = "#3598DB", "Splice Site" = "#FFCD00", "Missense" = "#25AF60", "Start Loss" = "#EE717A", "Stop Loss" = "green")


h_marg <- 0.15
v_marg <- 3


# PLOT: HN driver with hetloss  --------------------------
# mat.hn <- df %>%
#   filter(Hugo_Symbol %in% HN_driver) %>%
#   select(Hugo_Symbol, sample, mutations) %>%
#   spread(key="sample", value="mutations") %>%
#   column_to_rownames("track_name") %>%
#   as.matrix()

pdf("results/Top_HN_drivers_and_CDK2AP1-TCGA.pdf", width=10, height=nrow(mutations.mat)/3)
oncoplot <- oncoPrint(mutations.mat, column_order = reorder_patients(mutations.mat), pct_digits = 1,
          pct_side = "right", row_names_side = "left",
          alter_fun = list(
            background = alter_graphic("rect", fill = "#F1F1F1", horiz_margin = unit(0, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Deep Deletion` = alter_graphic("rect", alpha = 0.7, fill = col["Deep Deletion"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Amplification` = alter_graphic("rect", alpha = 0.7, fill = col["Amplification"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            #`Shallow Deletion` = alter_graphic("rect", fill = col["Shallow Deletion"], alpha=0.4, horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Splice Site` = alter_graphic("rect", height = 0.5, fill = col["Splice Site"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Missense` = alter_graphic("rect", height = 0.5, fill = col["Missense"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Inframe InDel` = alter_graphic("rect", height = 0.5, fill = col["Inframe InDel"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Stop Loss` = alter_graphic("rect", height = 0.5, fill = col["Stop Loss"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Start Loss` = alter_graphic("rect", height = 0.5, fill = col["Start Loss"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Frameshift` = alter_graphic("rect", height = 0.5, fill = col["Frameshift"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Nonsense` = alter_graphic("rect", height = 0.5, fill = col["Nonsense"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt"))
            ))
draw(oncoplot)
dev.off()
