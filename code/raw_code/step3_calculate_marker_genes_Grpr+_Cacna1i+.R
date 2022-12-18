library(Seurat)
library(SeuratDisk)

library(tidyverse)
library(cowplot)
library(ArchR)
library(here)

library(Rmagic)
library(phateR)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/tables'
dir.create(here(DATADIR), showWarnings = F)

## set to be parallel over 28 cores with 500Gb
plan("multicore", workers = 12)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 150 * 1024^3)

########################################################
## 1) load the object with the dorsal horn neuron clusters
obj_neuron = LoadH5Seurat(here('data/tidy_data/rdas/Russ_et_al_dorsal_horn_neurons.h5seurat'))
obj_neuron$group = ifelse(obj_neuron$family %in% c('Rreb1', 'Sox5'), 'Grpr', 'non-Grpr')
Idents(obj_neuron) ='group'

#############################################################################
## 2) compute the marker genes that are consistently differentially regulated 
## in Grpr+ neurons between datasets
markers_grpr = FindConservedMarkers( obj_neuron, ident.1 = 'Grpr',ident.2 = 'non-Grpr',
  grouping.var = 'dataset', assay = "RNA", slot = "data", min.cells.group = 3,
  ## lower the min.pct and logfc.threshold or else lowly expressed genes won't be tested
  verbose = TRUE, logfc.threshold = .01, min.pct = .01) %>% 
  rownames_to_column('gene')

# check Grpr should be in table and that it's upregulated
markers_grpr %>% filter(gene == 'Grpr') 

## rename the columns to be interpretable
markers_df = markers_grpr %>%
  rename_all(~ str_replace_all(., 'pct.1', 'pct.Exp.Grpr') %>% 
               str_replace_all('pct.2', 'pct.Exp.non.Grpr'))

table_fn = here(DATADIR, 'Russ_et_al_dorsal_horn_GrprVsNonGrpr_markerGenes.xlsx')
markers_df %>% writexl::write_xlsx(table_fn)

