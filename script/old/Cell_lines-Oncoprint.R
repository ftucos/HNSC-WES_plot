wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)

source("script/src/reorder_cell_lines.R")

drop.na <- function(vector) {vector[!is.na(vector)]}

# show all genes
# column name order from TCGA



# load data ------------------------------------
nurd <- c("CHD3", "CHD4", "CHD5", "HDAC1", "HDAC2", "MTA1", "MTA2", "MTA3", "GATAD2A", "GATAD2B", "PWWP2A", "PWWP2B", "MBD2", "MBD3", "RBBP4", "RBBP7", "CDK2AP1")
# defined by TCGA freequency of mutation
nurd_order <- c("CHD4", "CHD3", "RBBP7", "MTA1", "MTA2", "CHD5", "MBD2", "GATAD2B", "HDAC2", "PWWP2B", "RBBP4", "MBD3", "MTA3", "GATAD2A", "HDAC1", "CDK2AP1", "PWWP2A")
HN_driver <- c("CDKN2A", "FAT1", "TP53", "CASP8", "NOTCH1", "PIK3CA", "AJUBA", "KMT2D", "KRT5", "TGFBR2", "CTCF", "HLA-A", "HRAS", "EPHA2", "H1-5", "HRNR", "NECAB1", "PSIP1", "FOSL2")
# defined by TCGA freequency of mutation
HN_driver_order <- c("TP53", "CDKN2A", "PIK3CA", "FAT1", "NOTCH1", "KMT2D", "CASP8", "AJUBA", "HRAS", "TGFBR2", "HLA-A","HRNR", "NECAB1", "EPHA2", "PSIP1","CTCF", "FOSL2", "H1-5","KRT5")
# mutations <- data.table::fread("/Volumes/TucosHDD/Bioinformatics/archive/Whole_Exome_Sequencing/WES_annotated_ANNOVAR/output/filtered_mutations-larger_freq.csv")
mutations <- data.table::fread("/Volumes/TucosHDD/Bioinformatics/archive/Whole_Exome_Sequencing/WES_annotated_ANNOVAR/output/filtered_mutations.csv")

mutation.1 <- mutations %>%
  select(track_name = Gene, sample = Sample, SNV = Func) %>% 
  mutate(SNV = recode_factor(SNV,
    "frameshift" = "Frameshift",
    "inframe indel" = "Inframe",
    "missense" = "Missense",
    "splice site" = "Splice site",
    "stop gain" = "Nonsense",
    "start loss" = "Other non-synonymous mutations",
    "stop loss" = "Other non-synonymous mutations"
    )) %>%
  group_by(track_name, sample) %>%
  summarize(SNV = paste0(SNV, collapse=";")) %>%
  # if multiple SNV, add a _2 to recognized tham and plot as an overlaying triangle
  mutate(SNV = ifelse(str_detect(SNV, ";"), paste0(SNV, "_2"), SNV))

# CNV
CNV <- read.csv('/Volumes/TucosHDD/Bioinformatics/archive/Whole_Exome_Sequencing/WES_annotated_ANNOVAR/data/merged_CNV.csv', stringsAsFactors=FALSE)
CNV.1 <- CNV %>%
  filter(BF > 20) %>%
  filter(is.na(Conrad.hg19)) %>%
  filter(!(type == 'deletion' & reads.ratio >= 1)) %>% # almost 80 genes from a couple of regions
  mutate(Func = recode(type, 'duplication'='copy gain', 'deletion'='copy loss'),
         Func = ifelse(.$reads.observed < 2 & .$reads.expected > 30, 'large deletion', Func)) %>%
  filter(standardizedValue <= -1 | standardizedValue >= 1)

CNV.2 <- CNV.1 %>%
  select(track_name = GeneSymbol, sample, CNA = Func) %>%
  mutate(CNA = recode(CNA, "copy gain" = "Amplification",
               "copy loss" = "Shallow Deletion",
               "large deletion" = "Deep Deletion"))
  

# Preprocessing ----------------------------
df <- full_join(CNV.2, mutation.1) %>%
  mutate(group = ifelse(track_name %in% nurd, "NuRD Complex", "Recurrent mutations in OSCC")) %>%
  rowwise() %>%
  mutate(mutations = list(c(as.character(CNA), as.character(SNV)) %>% drop.na() %>% paste0(collapse=";")))

# plot configuration -----------------------
col <- c("Deep Deletion" = "#1C5FA0", "Amplification" = "#A51428", "Shallow Deletion" = "#89C0D9",
         "Nonsense" = "#2B2B2B", "Frameshift" = "#755EC5", "Inframe" = "#3598DB", "Splice site" = "#FFCD00", "Missense" = "#25AF60", "Other non-synonymous mutations" = "#EE717A")
         
# PLOT: HN driver with hetloss  --------------------------
empty_hn = data.frame(track_name = HN_driver)

mat.hn <- df %>%
  filter(track_name %in% HN_driver) %>%
  select(track_name, sample, mutations)%>%
  pivot_wider(names_from="sample", values_from = "mutations") %>%
  full_join(empty_hn) %>%
  column_to_rownames("track_name") %>%
  replace(.=="NULL", "") %>%
  # add eventual missing columns
  add_column(CA1 = "", LM = "", LUC4 = "", SCC4 = "", SCC9 = "", SCC15 = "", SCC25 = "") %>%
  select(-ends_with(".1")) %>%
  as.matrix()

h_pct = 0.7
w_pct = 0.7

pdf("results/HN_drivers-cell_lines-with_hetloss.pdf", width=3.5, height=nrow(mat.hn)/3+0.5)
oncoplot <- oncoPrint(mat.hn, column_order = c("CA1", "LM", "LUC4", "SCC4", "SCC9", "SCC15", "SCC25" ), 
                      show_column_names=T, row_order = HN_driver_order, 
          pct_side = "right", row_names_side = "left", pct_digits = 1,
          alter_fun = list(
            background = alter_graphic("rect", fill = "gray94", height = h_pct, width = h_pct),
            `Deep Deletion` = alter_graphic("rect", col = "gray94", fill = col["Deep Deletion"],  height = h_pct, width = w_pct),
            `Amplification` = alter_graphic("rect", col = "gray94", fill = col["Amplification"], height = h_pct, width = w_pct),
            `Shallow Deletion` = alter_graphic("rect", col = "gray94", fill = col["Shallow Deletion"],  height = h_pct, width = w_pct),
            `Nonsense` = alter_graphic("rect", col = "gray94",  fill = col["Nonsense"],  height = h_pct/2, width = w_pct),
            `Frameshift` = alter_graphic("rect", col = "gray94", fill = col["Frameshift"], height = h_pct/2, width = w_pct),
            `Other non-synonymous mutations` = alter_graphic("rect", col = "gray94", fill = col["Other non-synonymous mutations"], height = h_pct/2, width = w_pct),
            `Inframe` = alter_graphic("rect", col = "gray94", fill = col["Inframe"], height = h_pct/2, width = w_pct),
            `Splice site` = alter_graphic("rect", col = "gray94", fill = col["Splice site"],  height = h_pct/2, width = w_pct),
            `Missense` = alter_graphic("rect", col = "gray94", fill = col["Missense"],  height = h_pct/2, width = w_pct),
            `Nonsense_2` = function(x, y, w, h) {
              grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                           y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                           gp = gpar(col = "gray94", fill = col["Nonsense"]))
            }, 
            `Frameshift_2` = function(x, y, w, h) {
              grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                           y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                           gp = gpar(col =" gray94", fill = col["Frameshift"]))
            }, 
            `Other non-synonymous mutations_2` = function(x, y, w, h) {
              grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                           y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                           gp = gpar( col = "gray94", fill = col["Other non-synonymous mutations"]))
            }, 
            `Splice site_2` = function(x, y, w, h) {
              grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                           y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                           gp = gpar(col = "gray94", fill = col["Splice site"]))
            }, 
            `Missense_2` = function(x, y, w, h) {
              grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                           y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                           gp = gpar(col = "gray94", fill = col["Missense"]))
            }, 
            
            `Inframe_2` = function(x, y, w, h) {
              grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                           y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                           gp = gpar(col = "gray94", fill = col["Inframe"]))
            } 
          ))
draw(oncoplot)
dev.off()

# PLOT: HN driver without hetloss  --------------------------
mat.hn <- df %>%
  mutate(mutations = str_remove(mutations, "Shallow Deletion;?"))%>%
  filter(track_name %in% HN_driver) %>%
  select(track_name, sample, mutations) %>%
  pivot_wider(names_from="sample", values_from = "mutations") %>%
  full_join(empty_hn) %>%
  column_to_rownames("track_name") %>%
  replace(.=="NULL", "") %>%
  # add eventual missing columns
  add_column(CA1 = "", LM = "", LUC4 = "", SCC4 = "", SCC9 = "", SCC15 = "", SCC25 = "") %>%
  select(-ends_with(".1")) %>%
  as.matrix()

pdf("results/HN_drivers-cell_line.pdf", width=3.5, height=nrow(mat.hn)/3+0.5)
oncoplot <- oncoPrint(mat.hn, column_order = c("CA1", "LM", "LUC4", "SCC4", "SCC9", "SCC15", "SCC25" ), 
                      show_column_names=T, row_order = HN_driver_order,
                      pct_side = "right", row_names_side = "left", pct_digits = 1,
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "gray94", height = h_pct, width = h_pct),
                        `Deep Deletion` = alter_graphic("rect", col = "gray94", fill = col["Deep Deletion"],  height = h_pct, width = w_pct),
                        `Amplification` = alter_graphic("rect", col = "gray94", fill = col["Amplification"], height = h_pct, width = w_pct),
                        `Shallow Deletion` = alter_graphic("rect", col = "gray94", fill = col["Shallow Deletion"],  height = h_pct, width = w_pct),
                        `Nonsense` = alter_graphic("rect", col = "gray94",  fill = col["Nonsense"],  height = h_pct/2, width = w_pct),
                        `Frameshift` = alter_graphic("rect", col = "gray94", fill = col["Frameshift"], height = h_pct/2, width = w_pct),
                        `Other non-synonymous mutations` = alter_graphic("rect", col = "gray94", fill = col["Other non-synonymous mutations"], height = h_pct/2, width = w_pct),
                        `Inframe` = alter_graphic("rect", col = "gray94", fill = col["Inframe"], height = h_pct/2, width = w_pct),
                        `Splice site` = alter_graphic("rect", col = "gray94", fill = col["Splice site"],  height = h_pct/2, width = w_pct),
                        `Missense` = alter_graphic("rect", col = "gray94", fill = col["Missense"],  height = h_pct/2, width = w_pct),
                        `Nonsense_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Nonsense"]))
                        }, 
                        `Frameshift_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col =" gray94", fill = col["Frameshift"]))
                        }, 
                        `Other non-synonymous mutations_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Other non-synonymous mutations"]))
                        }, 
                        `Splice site_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Splice site"]))
                        }, 
                        `Missense_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Missense"]))
                        }, 
                        
                        `Inframe_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Inframe"]))
                        }
                      ))
draw(oncoplot)
dev.off()


# PLOT: NURD driver with hetloss  --------------------------
empty_nurd = data.frame(track_name = nurd)

mat.nurd <- df %>%
  filter(track_name %in% nurd) %>%
  select(track_name, sample, mutations) %>%
  pivot_wider(names_from="sample", values_from = "mutations") %>%
  full_join(empty_nurd) %>%
  column_to_rownames("track_name") %>%
  replace(.=="NULL", "") %>%
  # add eventual missing columns
  add_column(CA1 = "", LM = "", LUC4 = "", SCC4 = "", SCC9 = "", SCC15 = "", SCC25 = "") %>%
  select(-ends_with(".1")) %>%
  as.matrix()

pdf("results/NuRD-cell_line-with_hetloss.pdf", width=3.3, height=nrow(mat.nurd)/3+0.5)
oncoplot <- oncoPrint(mat.nurd, column_order = c("CA1", "LM", "LUC4", "SCC4", "SCC9", "SCC15", "SCC25" ), 
                      show_column_names=T, row_order = nurd_order,
                      pct_side = "right", row_names_side = "left", pct_digits = 1,
                      alter_fun = list(
                        background = alter_graphic("rect", fill = "gray94", height = h_pct, width = h_pct),
                        `Deep Deletion` = alter_graphic("rect", col = "gray94", fill = col["Deep Deletion"],  height = h_pct, width = w_pct),
                        `Amplification` = alter_graphic("rect", col = "gray94", fill = col["Amplification"], height = h_pct, width = w_pct),
                        `Shallow Deletion` = alter_graphic("rect", col = "gray94", fill = col["Shallow Deletion"],  height = h_pct, width = w_pct),
                        `Nonsense` = alter_graphic("rect", col = "gray94",  fill = col["Nonsense"],  height = h_pct/2, width = w_pct),
                        `Frameshift` = alter_graphic("rect", col = "gray94", fill = col["Frameshift"], height = h_pct/2, width = w_pct),
                        `Other non-synonymous mutations` = alter_graphic("rect", col = "gray94", fill = col["Other non-synonymous mutations"], height = h_pct/2, width = w_pct),
                        `Inframe` = alter_graphic("rect", col = "gray94", fill = col["Inframe"], height = h_pct/2, width = w_pct),
                        `Splice site` = alter_graphic("rect", col = "gray94", fill = col["Splice site"],  height = h_pct/2, width = w_pct),
                        `Missense` = alter_graphic("rect", col = "gray94", fill = col["Missense"],  height = h_pct/2, width = w_pct),
                        `Nonsense_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Nonsense"]))
                        }, 
                        `Frameshift_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col =" gray94", fill = col["Frameshift"]))
                        }, 
                        `Other non-synonymous mutations_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Other non-synonymous mutations"]))
                        }, 
                        `Splice site_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Splice site"]))
                        }, 
                        `Missense_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Missense"]))
                        }, 
                        
                        `Inframe_2` = function(x, y, w, h) {
                          grid.polygon(x = c(x+0.5*w*w_pct, x+0.5*w*w_pct, x-0.5*w*w_pct),
                                       y = c(y+0.5*h*h_pct/2, y-0.5*h*h_pct/2, y+0.5*h*h_pct/2),
                                       gp = gpar( col = "gray94", fill = col["Inframe"]))
                        }
                      ))
draw(oncoplot)
dev.off()

