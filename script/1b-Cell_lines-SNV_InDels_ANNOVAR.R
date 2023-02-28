wd <- "/Users/tucos/Documents/bioinformatic_projects/HNSC_WES_plot"
setwd(wd)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)

### Note: if changing filters, consider checking again ExoniFunc field to retain (splicing related?)
# source("script/src/reorder_patient.R")

drop.na <- function(vector) {vector[!is.na(vector)]}

# CNVKIT --------------
path <- "/Volumes/TucosHDD/Bioinformatics/data/Fodde_Lab/WES_reprocessing-OSCC/custom_results/ANNOVAR/"
files <- list.files(path = path, recursive = T, pattern = "multianno.txt$", full.names = T)

parse_annovar <- function(path) {
df <- data.table::fread(path, na.strings = ".", check.names = T) %>%
  mutate(Sample = str_extract(path, "[0-9A-Za-z]+(?=\\.mutect2\\.hg38_multianno\\.)")) %>%
  select(Sample, Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, cytoBand, rmsk, Interpro_domain, X1000g2015aug_all, esp6500siv2_all, gnomAD_genome_ALL, gnomad312_AF, CLNSIG, avsnp150, SIFT_pred, SIFT4G_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, LRT_pred, MutationTaster_pred, MutationAssessor_pred, FATHMM_pred, PrimateAI_pred, `GERP++_RS` = GERP.._RS, cosmic97_coding, ICGC_Id, ICGC_Occurrence) %>%
  mutate(X1000g2015aug_all = replace_na(X1000g2015aug_all, 0),
         esp6500siv2_all = replace_na(esp6500siv2_all, 0),
         gnomAD_genome_ALL = replace_na(gnomAD_genome_ALL, 0),
         gnomad312_AF = replace_na(gnomad312_AF, 0))

df.1 <- df %>%
  mutate(cosmic97 = str_extract(cosmic97_coding, "COSV[^;]+"),
         cosmic97_occurrence = str_extract_all(cosmic97_coding, "[0-9]+(?=\\()"),
         ICGC_Occurrence = str_extract_all(ICGC_Occurrence, "(?<=(^|,)[^\\|]{1,10}\\|)[0-9]+"),
         ) %>%
  rowwise() %>%
  mutate(cosmic97_occurrence = cosmic97_occurrence %>% unlist() %>% as.numeric() %>% sum(na.rm=T),
         ICGC_Occurrence = ICGC_Occurrence %>% unlist() %>% as.numeric() %>% sum(na.rm=T)) %>%
  select(-cosmic97_coding) %>%
  ungroup()

return(df.1)
}

df<- map(files, parse_annovar) %>%
  do.call(what=rbind)

recurrence <- df %>%
  distinct() %>%
  group_by(Chr, Start, End, Ref, Alt) %>% 
  summarize(reccurrence_in_samples = n())

df.1 <- df %>%
  left_join(recurrence) %>%
  # report splicing effect also in the exonic func 
  mutate(ExonicFunc.refGene = ifelse(Func.refGene == "splicing", "splicing", ExonicFunc.refGene)) %>%
  filter(reccurrence_in_samples < 5,
         X1000g2015aug_all < 0.01, esp6500siv2_all < 0.01, gnomAD_genome_ALL < 0.01, gnomad312_AF < 0.01,
         Func.refGene %in% c("exonic", "exonic;splicing", "splicing"),
         ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion",
                                    "nonframeshift deletion", "nonframeshift insertion",
                                    "nonframeshift substitution", "nonsynonymous SNV", "startloss", "stopgain", "stoploss",
                                    "splicing"),
         # In-frame deletions/insertions in a repetitive region without a known function
         !(ExonicFunc.refGene %in% c("nonframeshift insertion", "nonframeshift substitution") & !is.na(rmsk))) %>%
  mutate(ExonicFunc.simplified = recode(ExonicFunc.refGene,
                                        "frameshift deletion" = "Frameshift",
                                        "frameshift insertion" = "Frameshift",
                                        "nonframeshift deletion" = "Inframe InDel",
                                        "nonframeshift insertion" = "Inframe InDel",
                                        "nonframeshift substitution" = "Inframe InDel",
                                        "nonsynonymous SNV" = "Missense",
                                        "startloss" = "Start Loss",
                                        "stopgain" = "Nonsense",
                                        "stoploss" = "Stop Loss",
                                        "splicing" = "Splice Site"
                                        ))

a <- df %>% filter(Gene.refGene == "TP53")

write.table(df.1, "processed/Cell_Lines-mutect2-Annovar-filtered.tsv",row.names = F, sep = "\t", quote = F, na = "")
  
# compare SCCS from depmap --------------
mutations <- data.table::fread("data/CCLE/SCC4-Cosmic_Mutations.tsv") %>% mutate(sample = "SCC4") %>%
  rbind(data.table::fread("data/CCLE/SCC9-Cosmic_Mutations.tsv") %>% mutate(sample = "SCC9")) %>%
  rbind(data.table::fread("data/CCLE/SCC15-Cosmic_Mutations.tsv") %>% mutate(sample = "SCC15")) %>%
  rbind(data.table::fread("data/CCLE/SCC25-Cosmic_Mutations.tsv") %>% mutate(sample = "SCC25")) %>%
  filter(!str_detect(Gene, "ENST[0-9]{8}")) %>%
  filter(Type != "Substitution - coding silent")

mutations.min <- mutations %>%
  select(Gene.refGene = Gene, Sample = sample, Cosmic.Type = Type,  Census_tier_1 = `Census Tier 1`, Cosmic.cds_Mutation = `CDS Mutation`, Cosmic.Validated = Validated, Position) %>%
  mutate(approx_start = str_extract(Position, "(?<=:)[0-9]+") %>% as.numeric(),
         approx_start = round(approx_start/10, 0)*10)

compare <- df %>% filter(Sample %in% mutations.min$Sample) %>%
  mutate(approx_start = round(Start/10, 0)*10) %>%
  select(Sample, Chr, Start, approx_start, End, Ref, Alt, Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, rmsk, Interpro_domain, gnomad312_AF, CLNSIG, avsnp150, cosmic97, cosmic97_occurrence, ICGC_Occurrence) %>%
  full_join(mutations.min, .)
