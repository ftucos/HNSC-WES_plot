wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
#library(ComplexHeatmap)

# SCC9, 4, 15 are triploid, SCC25 are diploid

# source("script/src/reorder_patient.R")

# Note: Control Freec only reports CNA and not non altered regions

drop.na <- function(vector) {vector[!is.na(vector)]}

# CNVKIT --------------
path <- "/Volumes/TucosHDD/Bioinformatics/data/Fodde_Lab/WES_reprocessing-OSCC/result/variant_calling/cnvkit"
# no differences between call.cns and cns
files <- list.files(path = path, recursive = T, pattern = "*.call.cns", full.names = T)

parse_CNVkit <- function(path){
  df <- data.table::fread(path)
  df.1 <- df %>%
    rowwise() %>%
    mutate(gene = list(str_split(gene, ",") %>%
             unlist() %>%
             {.[str_detect(., "ref\\|")] %>% str_remove("ref\\|")} %>% 
             {.[!str_detect(.,"N[A-Z]_")]})) %>%
    unnest("gene") %>%
    distinct()
  
  df.2 <- df.1 %>%
    group_by(chromosome, gene) %>%
    # when change in CN happens across a gene region, the minimum number of cn is kept as reference status for that gene
    arrange(cn) %>%
    slice_min(cn, n=1) %>%
    select(chromosome, gene, cn, log2, depth, p_ttest, weight) %>%
    mutate(sample = str_extract(path, "(?<=cnvkit/)[^/]+"))
  
  return(df.2)
}

df.cnvkit <- map(files, parse_CNVkit) %>%
  do.call(what=rbind) %>%
  select(sample, everything())  

# gene_cp <- df.cnvkit %>% group_by(gene) %>% summarize(status = median(cn))

# ControlFreec ---------------------
path <- "/Volumes/TucosHDD/Bioinformatics/data/Fodde_Lab/WES_reprocessing-OSCC/custom_results/controlfreec"
# no differences between call.cns and cns
files <- list.files(path = path, recursive = T, pattern = "*.annotated.p.value.txt", full.names = T)

parse_ControlFreec <- function(path){
  df <- data.table::fread(path, col.names = c("chromosome", "start", "end", "other", "gene")) %>%
    mutate(gene = str_extract_all(gene, '(?<=gene_name=")[^"]+')) %>% 
    unnest(gene) %>%
    distinct() %>%
    mutate(other = str_split(other, "\\|")) %>%
    unnest_wider(other, names_sep = "_") %>%
    rename(cn = other_1, status = other_2, genotype = other_3, uncertainty = other_4, WilcoxonRankSumTestPvalue = other_5, KolmogorovSmirn = other_6) %>%
    select(-status, -genotype)
  
  df.1 <- df %>%
    group_by(chromosome, gene) %>%
    arrange(cn) %>%
    slice_min(cn, n=1) %>%
    # when change in CN happens across a gene region, the minimum number of cn is kept as reference status for that gene
    select(chromosome, gene, cn, WilcoxonRankSumTestPvalue, KolmogorovSmirn, uncertainty) %>%
    mutate(cn = as.numeric(cn),
           WilcoxonRankSumTestPvalue = as.numeric(WilcoxonRankSumTestPvalue),
           KolmogorovSmirn = as.numeric(KolmogorovSmirn),
           uncertainty = as.numeric(uncertainty),
           sample = str_extract(path, "(?<=controlfreec/)[^\\.]+"))
  
  return(df.1)
}
  
df.controlfreec <- map(files, parse_ControlFreec) %>%
  do.call(what=rbind) %>%
  select(sample, everything())

# gene_cp.2 <- df.controlfreec %>% group_by(gene) %>% summarize(status = median(cn))

# merge --------------
df <- full_join(df.cnvkit %>% rename(CN.cnvkit = cn, log2.cnvkit = log2, depth.cnvk = depth, p_ttest.cnvkit = p_ttest, weight.cnvkit = weight),
                df.controlfreec %>% rename(CN.controlfreec = cn, WilcoxonRankSumTestPvalue.controlfreec = WilcoxonRankSumTestPvalue, KolmogorovSmirn.comtrolfreec = KolmogorovSmirn, uncertainty.controlfreec = uncertainty)) %>%
  ungroup() %>%
  mutate(CN.cnvkit = ifelse(CN.cnvkit >= 10, 10, CN.cnvkit),
         CN.controlfreec = ifelse(CN.controlfreec >= 20, 20, CN.controlfreec),
         #  missing genes with ploidy = 2
         CN.cnvkit = replace_na(CN.cnvkit, 2),
         CN.controlfreec = replace_na(CN.controlfreec, 2)
         )

# compare CCLE ---------------
sccs_ref <- data.table::fread("data/CCLE/Copy_Number_(Absolute)_subsetted.csv") %>%
  select(-depmap_id, -lineage_1, -lineage_2, -lineage_3, -lineage_4, ) %>%
  column_to_rownames("cell_line_display_name") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  gather(key="sample", value="CN.CCLE", -1) %>%
  rowwise() %>%
  # scale by 1 (min is set to 0) nearly triploid cell lines
  mutate(CN.CCLE_as_diploid = ifelse(sample %in% c("SCC4", "SCC9", "SCC15"), max(0, CN.CCLE-1), CN.CCLE)) %>%
  ungroup() %>%
  # categorize CCLE CNA stuats
  mutate(CN_status.CCLE = case_when(
    CN.CCLE == 0 ~ "Deep Deletion",
    CN.CCLE == 1 ~ "Shallow Deletion",
    CN.CCLE >= 5 ~ "Amplification",
    CN.CCLE >= 4 ~ "Gain",
    !is.na(CN.CCLE) ~ "Normal"
  ) %>% factor(levels = c("Deep Deletion", "Shallow Deletion", "Normal", "Gain", "Amplification")),
  CN_status_3cl.CCLE = recode_factor(CN_status.CCLE,
                              "Shallow Deletion" = "Normal",
                              "Gain" = "Normal") %>%
    factor(levels =  c("Deep Deletion", "Normal", "Amplification")),
  CN_status.CCLE_as_diploid = case_when(
    CN.CCLE_as_diploid == 0 ~ "Deep Deletion",
    CN.CCLE_as_diploid == 1 ~ "Shallow Deletion",
    CN.CCLE_as_diploid >= 5  ~ "Amplification",
    CN.CCLE_as_diploid >= 3 ~ "Gain",
    !is.na(CN.CCLE_as_diploid) ~ "Normal"
  ) %>% factor(levels = c("Deep Deletion", "Shallow Deletion", "Normal", "Gain", "Amplification")),
  CN_status_3cl.CCLE_as_diploid = recode_factor(CN_status.CCLE_as_diploid,
                              "Shallow Deletion" = "Normal",
                              "Gain" = "Normal") %>%
  factor(levels =  c("Deep Deletion", "Normal", "Amplification"))
)

compare <- inner_join(sccs_ref, df) %>%
  mutate(CN.mean = (CN.cnvkit + CN.controlfreec)/2,
         CN.min_dev = case_when(
           is.na(CN.cnvkit) & !is.na(CN.controlfreec) ~ CN.controlfreec,
           !is.na(CN.cnvkit) & is.na(CN.controlfreec) ~ CN.cnvkit,
           abs(CN.cnvkit - 2) < abs(CN.controlfreec - 2) ~ CN.cnvkit,
           abs(CN.cnvkit - 2) > abs(CN.controlfreec - 2) ~ CN.controlfreec,
           TRUE ~ (CN.cnvkit + CN.controlfreec)/2
         ))

cor.test(compare$CN.CCLE, compare$CN.controlfreec, method = "spearman")
cor.test(compare$CN.CCLE, compare$CN.cnvkit, method = "spearman")
cor.test(compare$CN.CCLE, compare$log2.cnvkit, method = "spearman")
cor.test(compare$CN.CCLE, compare$CN.mean, method = "spearman")
cor.test(compare$CN.CCLE, compare$CN.min_dev, method = "spearman")

cor.test(compare$CN.CCLE_as_diploid, compare$CN.controlfreec, method = "spearman")
cor.test(compare$CN.CCLE_as_diploid, compare$CN.cnvkit, method = "spearman")
cor.test(compare$CN.CCLE_as_diploid, compare$log2.cnvkit, method = "spearman")
cor.test(compare$CN.CCLE_as_diploid, compare$CN.mean, method = "spearman")
cor.test(compare$CN.CCLE_as_diploid, compare$CN.min_dev, method = "spearman")


# best correlation within CN.CCLE_as_diploid and CN.mean

# categorized with mean combined score ---------------

ggplot(compare, aes(y = CN.controlfreec, x=CN.cnvkit, color = CN_status.CCLE_as_diploid)) +
  geom_jitter(alpha = 0.5, size = 0.02) +
  scale_color_manual(values = c("black", "blue", "green", "orange", "red"), name="CCLE CN\n(real. to normal ploidy)") +
  theme_bw() +
  ggh4x::force_panelsizes(rows=unit(12, "cm"),cols = unit(12, "cm")) +
  xlab("CN (CNVKit)") + ylab("CN (Control Freec)")

ggplot(compare %>% mutate(log2.cnvkit = ifelse(log2.cnvkit < -3, -3, log2.cnvkit)), aes(y = CN.controlfreec, x=log2.cnvkit, color = CN_status.CCLE_as_diploid)) +
  geom_jitter(alpha = 0.5, size = 0.02) +
  scale_color_manual(values = c("black", "blue", "green", "orange", "red"), name="CCLE CN\n(real. to normal ploidy)") +
  theme_bw() +
  ggh4x::force_panelsizes(rows=unit(12, "cm"),cols = unit(12, "cm")) +
  xlab("log2 (CNVKit)") + ylab("CN (Control Freec)")

ggsave("results/exploratory_plots/log2_CNVKIT_vs_ControlFreec-SCCs_CCLE_control.png", plot=last_plot(), width = 20, height = 18, units="cm")  


min_p = min(compare$p_ttest[compare$p_ttest != 0], na.rm = T)

ggplot(compare, aes(x = CN.CCLE_as_diploid, y=CN.cnvkit, color = -log10(p_ttest.cnvkit + min_p))) +
  geom_jitter(alpha = 0.5, size = 0.05) +
  scale_color_viridis_c()

ggplot(compare, aes(x = CN.cnvkit, y=CN.controlfreec, color = CN.CCLE_as_diploid)) +
  geom_jitter(alpha = 0.5, size = 0.05) +
  scale_color_viridis_c()


# Export -----------------
export <- df.cnvkit %>%
  mutate(status = case_when(
    cn >= 5 ~ "Amplification",
    cn >=3 ~ "Gain",
    cn == 2 ~ "Normal",
    cn == 1 ~ "Shallow Deletion",
    cn == 0 ~ "Deep Deletion"
  ))

# Export simlplified -----------------
compare.1 <- compare %>%
  mutate(neglog_p.cnvkit = -log10(p_ttest.cnvkit+min(df.cnvkit$p_ttest %>% .[.!=0])),
         # make NS gray
         neglog_p.cnvkit = ifelse(p_ttest.cnvkit > 0.05, NA, neglog_p.cnvkit),
         trimmed_log2.cnvkit = ifelse(log2.cnvkit < -4, -4, log2.cnvkit))

ggplot(compare.1, aes(x=CN.CCLE, y=trimmed_log2.cnvkit, color=neglog_p.cnvkit)) + 
  geom_hline(yintercept = c(-0.5, 0.5), color = "red", linetype='dotted') + 
  geom_jitter(height = 0.1, alpha=0.4, size=0.2) +
  theme_bw() +
  scale_color_viridis_c(na.value = "red") +
  xlab("CCLE Absolute CN\n(SNP Array)")+ylab("Log2 copy ratio\n(CNVkit)")+
  ggh4x::force_panelsizes(rows=unit(12, "cm"),cols = unit(12, "cm"))

ggsave("results/exploratory_plots/log2_CNVKIT_vs_SCCs_CCLE_control.png", plot=last_plot(), width = 20, height = 18, units="cm")  


ggplot(compare.1, aes(x = trimmed_log2.cnvkit, y=neglog_p.cnvkit, color = CN_status.CCLE_as_diploid)) +
  geom_jitter(alpha = 0.5, size = 0.02) +
  scale_color_manual(values = c("black", "blue", "green", "orange", "red"), name="CCLE CN\n(real. to normal ploidy)") +
  theme_bw() +
  ggh4x::force_panelsizes(rows=unit(12, "cm"),cols = unit(12, "cm")) +
  xlab("log2 (CNVKit)") 

export.simmplified <- df.cnvkit %>%
  mutate(status = case_when(
    log2 < -0.5 & p_ttest < 0.05 ~ "Loss",
    log2 > 0.5 & p_ttest < 0.05 ~ "Gain",
    TRUE ~ "Normal"
  )) 

write_tsv(export.simmplified, "processed/Cell_Lines-CNVKit.tsv")
