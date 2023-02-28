wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(tidyverse)
library(ComplexHeatmap)

source("script/src/reorder_patient.R")

drop.na <- function(vector) {vector[!is.na(vector)]}


# Load data --------------------
# https://gdc.cancer.gov/about-data/publications/pancanatlas
# nurd <- c("CHD3", "CHD4", "CHD5", "HDAC1", "HDAC2", "MTA1", "MTA2", "MTA3", "MBD2", "GATAD2A", "GATAD2B", "PWWP2A", "PWWP2B", "MBD2", "MBD3", "RBBP4", "RBBP7", "CDK2AP1")
HN_driver <- c("CDKN2A", "FAT1", "TP53", "CASP8", "NOTCH1", "PIK3CA", "AJUBA", "KMT2D", "KRT5", "TGFBR2", "CTCF", "HLA-A", "HRAS", "EPHA2", "H1-5", "HRNR", "NECAB1", "PSIP1", "FOSL2")
# top mutated including non significant driver
# HN_driver <- c("CDKN2A","FAT1","TP53","CASP8","NOTCH1","PIK3CA","AJUBA","KMT2D","KRT5","TGFBR2","CTCF","HLA-A","HRAS","EPHA2","HIST1H1B","HRNR","NECAB1","PSIP1","FOSL2","NFE2L2","FBXW7","RAC1","SURF2","HIST1H4E","HNRNPL","TANC1","NSD1","OR52E2","CXorf48","TARP","LCP1","KRTAP13-2","EGFR","C14orf179","SAT1","SATB2","ERC2","STC1","ZNF750","PTEN","EPDR1","MUDENG","C6","ITPA","KCNT2","FLG","IL1RAPL1","RB1","ADCY4","ANP32C","ZNF451","TCEAL3","VPS33B","CUL3","FGF23","PEX2","SEMA5A","KIAA1529","UBA52","FCRL4","B2M","CDH3","NPY5R")


SNV.raw <-  data.table::fread("processed/Cell_Lines-mutect2-Annovar-filtered.tsv") %>%
  rename(Hugo_Symbol = Gene.refGene,
         Variant_Classification = ExonicFunc.simplified)

SNV_priority <- c("Nonsense", "Start Loss", "Frameshift", "Stop Loss", "Inframe InDel", "Missense", "Splice Site")

SNV <- SNV.raw %>%
  filter(Hugo_Symbol %in% c( HN_driver)) %>%
  select(Hugo_Symbol, Sample,  SNV = Variant_Classification) %>%
  mutate(SNV = SNV %>% factor(levels=SNV_priority))

# keep only higher priority mutation when concurrent are present
SNV.1 <- SNV %>% group_by(Sample, Hugo_Symbol) %>%
  arrange(SNV) %>%
  slice_head(n = 1)

# validate prioritizzation
# SNV.2 <-  SNV %>% group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarize(types = list(as.character(SNV)), count = n())
# 
# x <- left_join(SNV.1, SNV.2)

# number of patients
SNV.raw$Sample %>% unique() %>% length()
sample_id.SNV = SNV.raw$Sample %>% unique()

# CNA -------------------------------------
CNA.raw <- data.table::fread("processed/Cell_Lines-CNVKit.tsv") %>%
  rename(Sample = sample,
         Hugo_Symbol = gene,
         CNA = status)

CNA.raw$Sample %>% unique() %>% length()
sample_id.CNA = CNA.raw$Sample %>% unique()

CNA <- CNA.raw %>%
  filter(Hugo_Symbol %in% c(HN_driver),
         CNA != "Normal") %>%
  select(Sample, Hugo_Symbol, CNA)

# should be 504 !!!!!!
complete_samples <- intersect(sample_id.SNV, sample_id.CNA)

# merge CNA and SNVs----------------------
df <- full_join(CNA , SNV.1) %>%
  rowwise() %>%
  mutate(mutation = list(c(as.character(CNA), as.character(SNV)) %>% drop.na() %>% paste0(collapse=";")) %>% unlist()) %>%
  select(-CNA, -SNV) %>%
  ungroup()

# add missing genes with no event
df <- df %>%
  plyr::rbind.fill(data.frame(Hugo_Symbol = setdiff(c(HN_driver), unique(df$Hugo_Symbol))))

mutations.mat.partial <- df %>%
  spread(key="Sample", value="mutation") %>%
  column_to_rownames("Hugo_Symbol") %>%
  select(-`<NA>`) %>%
  as.matrix()

# add missing patients with no mutations
no_mut_samples <- complete_samples[!complete_samples %in% colnames(mutations.mat.partial)]
no_mut_matrix <- matrix(nrow = nrow(mutations.mat.partial), ncol = length(no_mut_samples),
                        dimnames = list(rownames(mutations.mat.partial), no_mut_samples))

mutations.mat <- cbind(mutations.mat.partial, no_mut_matrix)
# replace NA with ""
mutations.mat[is.na(mutations.mat)] <- ""


# plot configuration -----------------------

col <- c("Loss" = "#1C5FA0", "Gain" = "#A51428",
         "Nonsense" = "#2B2B2B", "Frameshift" = "#755EC5", "Inframe InDel" = "#3598DB", "Splice Site" = "#FFCD00", "Missense" = "#25AF60", "Start Loss" = "#EE717A", "Stop Loss" = "green")


h_marg <- 2
v_marg <- 3



# PLOT: HN driver with hetloss  --------------------------
# for longer plot and more mutations set width to 3.5 (in place of 3.19)
pdf("results/HN_driver-HNSCC_cell_lines.pdf", width=3.19, height=nrow(mutations.mat)/3)
oncoplot <- oncoPrint(mutations.mat, column_order = c("CA1", "LM", "LUC4", "SCC4", "SCC9", "SCC15", "SCC25"), pct_digits = 1,
          pct_side = "right", row_names_side = "left", show_column_names = TRUE,
          alter_fun = list(
            background = alter_graphic("rect", fill = "#F1F1F1", horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Loss` = alter_graphic("rect", alpha = 0.7, fill = col["Loss"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
            `Gain` = alter_graphic("rect", alpha = 0.7, fill = col["Gain"], horiz_margin = unit(h_marg, "pt"), vertical_margin = unit(v_marg, "pt")),
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
