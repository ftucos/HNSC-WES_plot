layer1 <- c("Deep Deletion", "Amplification", "Shallow Deletion", "Gain")
layer2 <-  c("Nonsense", "Frameshift", "Start Loss", "Stop Loss", "Inframe InDel", "Missense", "Splice Site")
# layer3 <- "Shallow Deletion"

reorder_patients <- function(mat) {
    # ordered base on the frequenciy of patients altered and not total number of mutations
  gene_order = rowSums(mat == "") %>% sort(decreasing = F) %>% names()
  
  column_order <- data.frame(gene_name = rep(gene_order, each=2), # each=3
                             layer = rep(c("layer1", "layer2"), times=length(gene_order))) %>% # "layer3"
    mutate(col_name = paste0(layer, "_", gene_name)) %>%
    pull("col_name")
  
  mat.order <- mat %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(key="patient", value="mutation", -1) %>%
    mutate(mutation = str_split(mutation, ";"),
           gene  = factor(gene, levels = gene_order)) %>%
    rowwise() %>%
    mutate(layer1 = list(mutation[mutation %in% layer1]),
           layer2 = list(mutation[mutation %in% layer2]),
#           layer3 = list(mutation[mutation %in% layer3]),
           # remove character0()
           layer1 = ifelse(length(layer1) == 0, NA, layer1),
           layer2 = ifelse(length(layer2) == 0, NA, layer2),
 #          layer3 = ifelse(length(layer3) == 0, NA, layer3),
           
           # prioritize mutation order
           layer1 = factor(layer1, levels = layer1),
           layer2 = factor(layer2, levels = layer2),
  #         layer3 = factor(layer3, levels = layer3)
           
    ) %>% 
    select(patient, gene, layer1, layer2) %>% # "layer3"
    pivot_wider(
      id_cols = patient,
      names_from = gene,
      values_from = c("layer1", "layer2")) %>% # "layer3"
    arrange(.[column_order]) %>% 
    pull("patient")
  
  return(mat.order)
}


  
