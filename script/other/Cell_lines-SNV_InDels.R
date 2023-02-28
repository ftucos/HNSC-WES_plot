wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)


# source("script/src/reorder_patient.R")

drop.na <- function(vector) {vector[!is.na(vector)]}

# CNVKIT --------------
path <- "/Volumes/TucosHDD/Bioinformatics/data/Fodde_Lab/WES_reprocessing-OSCC/result/annotation"
files <- list.files(path = path, recursive = T, pattern = "*.mutect2.filtered_snpEff_VEP.ann.tab", full.names = T)

parse_variants <- function(path) {
  df <- data.table::fread(cmd = paste0("awk '/^#?[^#]/ { print $0 }' ", path), na.strings = "-",check.names = T) %>%
    mutate(Sample = str_extract(path, "(?<=annotation\\/)[^/]+")) %>%
    # keep only columns relevant for further filtering
    select(Sample, SYMBOL, Uploaded_variation = `X.Uploaded_variation`, Location, Allele, Gene, Feature, Feature_type, Consequence, Amino_acids, Codons, CLIN_SIG, Existing_variation,  IMPACT, VARIANT_CLASS, BIOTYPE, SIFT, PolyPhen, DOMAINS, AF, gnomAD_AF, ExAC_AF, ExAC_nonTCGA_AF, MAX_AF, FREQS, GERP_RS_rankscore = `GERP.._RS_rankscore`)
  return(df)
}

df <- map(files, parse_variants) %>%
  do.call(what=rbind)

df.1 <- df %>%
  # filter protein coding only 
  filter(BIOTYPE == "protein_coding") %>%
  # remove low impact
  filter(IMPACT != "LOW")

# a <- as.data.frame(table(df.1$Consequence, df.1$IMPACT)) %>% filter(Freq != 0)  

df.2 <- df.1 %>%
  filter(str_detect(Consequence, "(frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_gained|stop_lost|inframe_deletion|inframe_insertion|protein_altering_variant|missense_variant)")) %>%
  mutate(Cosmic_ID = ifelse(!is.na(Existing_variation), str_extract_all(Existing_variation, "COSV[^,]+"), NA)) %>%
  unnest(Cosmic_ID, keep_empty = TRUE) %>%
  mutate(Cosmic_ID = ifelse(Cosmic_ID == "", NA, Cosmic_ID))

# add cosmic fequency for consensus mutations
cosmic_freq <- data.table::fread("/Volumes/TucosHDD/Bioinformatics/resources/cosmic/cmc_export.tsv.gz") %>%
  select(Cosmic_ID = GENOMIC_MUTATION_ID, COSMIC_SAMPLE_MUTATED, MUTATION_SIGNIFICANCE_TIER)

df.3 <- df.2 %>%
  left_join(cosmic_freq) %>%
  group_by(across(all_of(colnames(df.2)[colnames(df.2) != "Cosmic_ID"]))) %>%
  summarize(Cosmic_ID = unique(Cosmic_ID) %>% drop.na() %>% paste0(collapse=","),
            MUTATION_SIGNIFICANCE_TIER = unique(MUTATION_SIGNIFICANCE_TIER) %>% drop.na() %>% paste0(collapse = ","),
            COSMIC_SAMPLE_MUTATED = COSMIC_SAMPLE_MUTATED %>% sum(na.rm = TRUE)) 

recurrent_in_samples <- df.3 %>% group_by(Uploaded_variation, SYMBOL) %>% summarize(count_in_samples = n())

df.4 <- df.3 %>%
  left_join(recurrent_in_samples)

plot(log10(df.4$COSMIC_SAMPLE_MUTATED), log10(df.4$MAX_AF))

# conflicting mutations
a <- df.4 %>% filter(COSMIC_SAMPLE_MUTATED > 50, MAX_AF > 0.5)

# x <- data.table::fread(cmd = paste0("awk '/^#?[^#]/ { print $0 }' ", "/Volumes/TucosHDD/Bioinformatics/data/Fodde_Lab/WES_reprocessing-OSCC/result/variant_calling/mutect2/CA1/CA1.mutect2.filtered.vcf.gz"), na.strings = "-",check.names = T) 
df.filtered <- df.4 %>%
  ungroup() %>%
  # filter MAF < 1% or missing
  filter(MAX_AF < 0.01 %>% replace_na(TRUE))
  
