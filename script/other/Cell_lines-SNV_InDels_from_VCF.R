wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)


# source("script/src/reorder_patient.R")

drop.na <- function(vector) {vector[!is.na(vector)]}

# CNVKIT --------------
path <- "/Volumes/TucosHDD/Bioinformatics/data/Fodde_Lab/WES_reprocessing-OSCC/result/annotation_vcf/mutect2/"
files <- list.files(path = path, recursive = T, pattern = "ann.vcf.gz$", full.names = T)

vcf <- data.table::fread(cmd = paste0("gunzip -c ", files[1], "| awk '/^#?[^#]/ { print $0 }' - "))

vcf.1 <- vcf %>% filter(FILTER == "PASS") %>%
  rename(STATS = CA1_CA1) 

# define new column names
nm1 <- c(setdiff(names(vcf.1), 'STATS'),
         c('GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'FAD', 'PGT', 'PID', 'PS', 'SB'))

vcf.2 <- vcf.1 %>%
  # convert GT:AD:AF:DP:F1R2:F2R1:FAD[:]SB formatted rows to GT:AD:AF:DP:F1R2:F2R1:FAD[:PGT:PID:PS:]SB
  mutate(STATS = ifelse(FORMAT == "GT:AD:AF:DP:F1R2:F2R1:FAD:SB",
                        str_replace(STATS, pattern=":(?=[^:]+$)", replacement ="::::"),
                        STATS)) %>%
  mutate(STATS = str_split(STATS, pattern = ":")) %>%
# split the STATS column and rename fiels
  unnest_wider(STATS, names_repair = ~ nm1, names_sep = "") %>%
  select(-FORMAT) 

nm2 <- c(setdiff(names(vcf.2), 'INFO'),
         c('Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'UNIPROT_ISOFORM', 'GENE_PHENO', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'FREQS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'TRANSCRIPTION_FACTORS'))
# split the INFO column 
vcf.3 <- vcf.2 %>%
  relocate(INFO, .after = last_col()) %>%
  mutate(INFO = str_split(INFO, pattern = "|")) %>%
  # split the INFO column and rename fiels
  unnest_wider(INFO, names_repair = ~ nm2, names_sep = "") %>%
  select(-INFO) 

  



