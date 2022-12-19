library(Seurat)
library(SeuratDisk)

library(cowplot)
library(ArchR)
library(here)
library(ggpubr)
library(tidyverse)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

PLOTDIR='figures/exploratory'
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)

############################################################
## 1) load the object with the dorsal horn neuron clusters
obj_neuron = LoadH5Seurat(here('data/tidy_data/rdas/Russ_et_al_dorsal_horn_neurons.h5seurat'), 
                          assay = 'RNA')
cols = paletteDiscrete(values = unique(obj_neuron$family))
cols2 = paletteDiscrete(values = unique(obj_neuron$final_cluster_assignment))
keep_dataset = c('Haring', 'Zeisel', 'Sathyamurthy')
cluster_grpr = paste0('Excit.', c(12:16))

## subset this down to the relevant ages and groupings
obj_neuron = obj_neuron[, obj_neuron$dataset %in% keep_dataset]
obj_exc = obj_neuron[, obj_neuron$coarse_cell_types == 'Excit']
obj_inh = obj_neuron[, obj_neuron$coarse_cell_types == 'Inhib']

table(obj_neuron$dataset)
obj_grpr = obj_exc[, obj_exc$final_cluster_assignment %in% cluster_grpr]
with(obj_grpr[[]], table(dataset, final_cluster_assignment))

###############################################################################
## 2) plot the expression of GRPR subcluster marker genes cross all DH neurons
grpr_df1 = obj_grpr[[]] %>% 
  group_by(dataset, final_cluster_assignment) %>% 
  summarize(num = n()) %>% ungroup() %>% group_by(dataset) %>% 
  summarize(final_cluster_assignment = final_cluster_assignment, 
            num = num, prop = num / sum(num))

grpr_df2 = obj_grpr[[]] %>% 
  group_by(final_cluster_assignment) %>% 
  summarize(dataset = 'All', num = n()) %>% 
  ungroup() %>% 
  mutate(num = num, prop = num / sum(num))

grpr_df = bind_rows(grpr_df1, grpr_df2) %>% 
  mutate(percent = prop * 100)

###############################################################
## 3) plot the marker gene expression across Grpr+ sub clusters

## put the 2 dot plots aligned vertically
pdf(here(PLOTDIR, 'plots', 'barPlot_grpr_subcluster_proportion.pdf'), 
    width = 8,  height =4)
ggplot(grpr_df, aes(x = final_cluster_assignment, y = percent)) + 
  geom_bar(stat = 'identity', aes(fill =final_cluster_assignment)) + 
  geom_text(aes(y = percent + 2, label=round(percent) )) + 
  scale_fill_brewer(palette = 'Set1') + 
  facet_wrap(dataset~.) + theme_classic() + 
  theme(legend.position = 'none')
dev.off()


