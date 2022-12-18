library(Seurat)
library(SeuratDisk)

library(tidyverse)
library(cowplot)
library(ArchR)
library(here)
library(ggpubr)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

PLOTDIR='figures/exploratory'
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)

## load the object with the dorsal horn neuron clusters
obj_neuron = LoadH5Seurat(here('data/tidy_data/rdas/Russ_et_al_dorsal_horn_neurons.h5seurat'))
cols = paletteDiscrete(values = unique(obj_neuron$family))
cols2 = paletteDiscrete(values = unique(obj_neuron$final_cluster_assignment))

cluster_grpr = paste0('Excit.', c(12:17))

obj_exc = obj_neuron[, obj_neuron$coarse_cell_types == 'Excit']
obj_inh = obj_neuron[, obj_neuron$coarse_cell_types == 'Inhib']
obj_grpr = obj_exc[, obj_exc$final_cluster_assignment %in% cluster_grpr]
grpr.markers.df = readRDS(here('data/tidy_data/rdas', 
                                 'Russ_et_al_dorsal_horn_Grpr+_subcluster_markerGenes.rds'))

###############################################################################
## 1) plot the expression of GRPR subcluster marker genes cross all DH neurons
features = grpr.markers.df %>% group_by(final_cluster_assignment) %>% 
  top_n(10, meta_score) %>% ungroup() %>% 
  top_n(40, summary_wilcox) %>% 
  arrange(final_cluster_assignment) %>% pull(gene) %>% unique()

p1 = DotPlot(obj_exc, features = features, assay = "RNA", col.min = 0,
             group.by = 'final_cluster_assignment', scale = T, scale.by = 'radius') + RotatedAxis() + 
  theme(legend.position = 'none', axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) + scale_y_discrete(limits = rev) + 
  ylab('Excitatory cluster')
p2 = DotPlot(obj_inh, features = features, assay = "RNA", col.min = 0,
             group.by = 'final_cluster_assignment', scale = T, scale.by = 'radius') + RotatedAxis()+ 
  theme(legend.position = 'bottom') + scale_y_discrete(limits = rev) +
  ylab('Inhibitory cluster') + xlab('Gene') +
  guides(color = guide_colourbar(title.position="top", 'Average Expression'),
         size = guide_legend(title.position="top",'Percent Expressed' ))

## put the 2 dot plots aligned vertically
pdf(here(PLOTDIR, 'plots', 'dotPlot_dh_grpr_subcluster_markers.pdf'), width = 8, height =12)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 16, 
          ncol = 1, align = "v", rel_heights = c(1, 1))
dev.off()



###############################################################
## 2) plot the marker gene expression across Grpr+ sub clusters
p3 = DotPlot(obj_grpr, features = features, assay = "RNA", col.min = 0,
             group.by = 'final_cluster_assignment', scale = T, scale.by = 'radius') + RotatedAxis()+ 
  theme(legend.position = 'right') + scale_y_discrete(limits = rev) +
  ylab('Grpr sub-cluster') + xlab('Gene')

## put the 2 dot plots aligned vertically
pdf(here(PLOTDIR, 'plots', 'dotPlot_dh_grpr_subcluster_markers.grpr_subset.pdf'), width = 8, height =4)
p3
dev.off()


