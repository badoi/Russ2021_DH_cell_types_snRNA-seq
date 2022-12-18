library(Seurat)
library(SeuratDisk)
library(harmony)

library(tidyverse)
library(cowplot)
library(ArchR)
library(here)

library(SingleCellExperiment)
library(BiocParallel)
library(scran)

library(data.table)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/tables'
dir.create(here(DATADIR), showWarnings = F)

## set to be parallel over 28 cores with 500Gb
plan("multicore", workers = 12)
BPPARAM = MulticoreParam(12)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 150 * 1024^3)

########################################################
## 1) load the object with the dorsal horn neuron clusters
obj_neuron = LoadH5Seurat(here('data/tidy_data/rdas/Russ_et_al_dorsal_horn_neurons_countsplit.h5seurat'), 
                          assay = 'countsplit')
sce_neuron = as.SingleCellExperiment(obj_neuron)
clusters = unique(sce_neuron$final_cluster_assignment) %>% sort()
alpha = 0.05

## choose the pairwise comparisions to do the DEG tests
cluster_grpr = paste0('Excit.', c(12:17))
cluster_exclude = clusters[! clusters %in% cluster_grpr]
clusters = unique(sce_neuron$final_cluster_assignment) %>% sort()
pairwise = expand.grid(clusters,clusters) %>% filter(Var1 != Var2) %>% 
  filter(Var1 %in% cluster_grpr)

test.wilcox <- pairwiseWilcox( 
  sce_neuron, groups = sce_neuron$final_cluster_assignment, 
  block = sce_neuron$dataset, direction = "up", 
  assay.type = 'logcounts', BPPARAM = BPPARAM)

test.binom <- pairwiseBinom( 
  sce_neuron, groups = sce_neuron$final_cluster_assignment, 
  block = sce_neuron$dataset, direction = "up", 
  assay.type = 'logcounts', BPPARAM = BPPARAM)

test.ttest <- pairwiseTTests( 
  sce_neuron, groups = sce_neuron$final_cluster_assignment, 
  block = sce_neuron$dataset, direction = "up",
  assay.type = 'logcounts', BPPARAM = BPPARAM)

## find the consistent Grpr markers across wilcoxon rank sum comparisons
i1 = which(test.wilcox$pairs$first %in% cluster_grpr)
grpr.wilcox.list = combineMarkers(
  de.lists = test.wilcox$statistics[i1], pairs = test.wilcox$pairs[i1,],
  pval.field = "p.value", effect.field = "AUC", pval.type = "some", 
  BPPARAM = BPPARAM)
grpr.wilcox.df = grpr.wilcox.list %>% as.list() %>% lapply(as.data.frame) %>% 
  lapply(rownames_to_column, 'gene') %>% 
  rbindlist(idcol = 'final_cluster_assignment',  fill=TRUE) %>% 
  dplyr::select(-starts_with('AUC')) %>% filter(FDR < alpha) %>% 
  mutate(test = 'wilcox') %>% dplyr::rename('summary' = 'summary.AUC')

## find the consistent Grpr markers across binomialcomparisons
i2 = which(test.binom$pairs$first %in% cluster_grpr)
grpr.binom.list = combineMarkers(
  de.lists = test.binom$statistics[i2], pairs = test.binom$pairs[i2,],
  pval.field = "p.value", effect.field = "logFC", pval.type = "some", 
  BPPARAM = BPPARAM)
grpr.binom.df = grpr.binom.list %>% as.list() %>% lapply(as.data.frame) %>% 
  lapply(rownames_to_column, 'gene') %>% 
  rbindlist(idcol = 'final_cluster_assignment',  fill=TRUE) %>% 
  dplyr::select(-starts_with('logFC')) %>% filter(FDR < alpha) %>% 
  mutate(test = 'binom') %>% dplyr::rename('summary' = 'summary.logFC')

## find the consistent Grpr markers across ttests rank sum comparisons
i3 = which(test.ttest$pairs$first %in% cluster_grpr)
grpr.ttest.list = combineMarkers(
  de.lists = test.ttest$statistics[i3], pairs = test.ttest$pairs[i3,],
  pval.field = "p.value", effect.field = "logFC", pval.type = "some", 
  BPPARAM = BPPARAM)
grpr.ttest.df = grpr.ttest.list %>% as.list() %>% lapply(as.data.frame) %>% 
  lapply(rownames_to_column, 'gene') %>% 
  rbindlist(idcol = 'final_cluster_assignment',  fill=TRUE) %>% 
  dplyr::select(-starts_with('logFC')) %>% filter(FDR < alpha) %>% 
  mutate(test = 'ttest') %>% dplyr::rename('summary' = 'summary.logFC')

## combine across tests
grpr.markers.df = bind_rows(grpr.wilcox.df, grpr.binom.df, grpr.ttest.df) %>% 
  pivot_wider(id_cols = c(final_cluster_assignment:gene), 
              names_from = test, values_from = c(p.value:summary)) %>% 
  mutate(meta_score = -(log(p.value_wilcox) + log(p.value_binom) +log(`p.value_ttest`))/3) %>% 
  filter(!is.na(meta_score)) %>% arrange(dplyr::desc(meta_score))
grpr.markers.df$meta_score
table(grpr.markers.df$final_cluster_assignment)

table_fn = here(DATADIR, 'Russ_et_al_dorsal_horn_Grpr+_subcluster_markerGenes.xlsx')
grpr.markers.df %>% split(f = .$final_cluster_assignment) %>% writexl::write_xlsx(table_fn)
grpr.markers.df %>% saveRDS(here('data/tidy_data/rdas', 
                                 'Russ_et_al_dorsal_horn_Grpr+_subcluster_markerGenes.rds'))
