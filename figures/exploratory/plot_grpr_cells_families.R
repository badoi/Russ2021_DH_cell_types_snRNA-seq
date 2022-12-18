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
obj_neuron$ageGroup = case_when(
  obj_neuron$dataset %in% c('Baek', 'Hayashi', 'Rosenberg') ~ 'Post-natal', 
  obj_neuron$dataset %in% c('Haring', 'Zeisel') ~ 'Juvenile', 
  obj_neuron$dataset %in% c('Sathyamurthy') ~ 'Adult'
) %>% factor(c('Post-natal', 'Juvenile',  'Adult'))

obj_exc = obj_neuron[, obj_neuron$coarse_cell_types == 'Excit']
obj_inh = obj_neuron[, obj_neuron$coarse_cell_types == 'Inhib']

#############################################################################
## 2) make plot to see which cell type has most Grpr expression, by families
p1 = VlnPlot(obj_exc,group.by = 'family', features = 'Grpr', cols = cols,
             log = T, assay = "MAGIC_RNA", pt.size = 0) + 
  xlab('Excitatory neuron cluster') + ylim(c(NA, 1.2)) +
  ylab('Relative Normalized Expression') + 
  theme(legend.position = 'none')
p2 = VlnPlot(obj_inh,group.by = 'family', features = 'Grpr', cols = cols,
             log = T, assay = "MAGIC_RNA", pt.size = 0) + 
  xlab('Inhibitory neuron cluster') + ylim(c(NA, 1.2)) +
  ylab('Relative Normalized Expression') + 
  theme(legend.position = 'none')


## put into 2 panel grid
pdf(here(PLOTDIR, 'plots', 'violinPlot_dh_neurons_families_grpr.pdf'), height = 4, width =7)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 16, align = "h",
          rel_widths = c(1.2, 1))
dev.off()

####################################################################################
## 3) make plot to see which cell type has most Grpr expression, by final clusters

### plot the final clusters
p1 = VlnPlot(obj_exc,group.by = 'final_cluster_assignment', features = 'Grpr', cols = cols2,
             log = T, assay = "MAGIC_RNA", pt.size = 0) + 
  xlab('Excitatory neuron cluster') + ylim(c(NA, 1.2)) +
  ylab('Relative Normalized Expression') + 
  theme(legend.position = 'none')
p2 = VlnPlot(obj_inh,group.by = 'final_cluster_assignment', features = 'Grpr', cols = cols2,
             log = T, assay = "MAGIC_RNA", pt.size = 0) + 
  xlab('Inhibitory neuron cluster') + ylim(c(NA, 1.2)) +
  ylab('Relative Normalized Expression') + 
  theme(legend.position = 'none')


## put into 2 panel grid
pdf(here(PLOTDIR, 'plots', 'violinPlot_dh_neurons_finalClusters_grpr.pdf'), height = 4, width =12)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 16, align = "h",
          rel_widths = c(1.2, 1))
dev.off()


#################################################
## 3) plot Grpr across developmental age groups
p3 = VlnPlot(obj_exc, group.by = 'family', split.by = 'ageGroup',
             features = 'Grpr', cols = cols,
             log = T, assay = "MAGIC_RNA", pt.size = 0) + 
  xlab('Excitatory neuron cluster') + ylim(c(NA, 1.2)) +
  ylab('Rel. Norm. Expression') + 
  theme(axis.text.x = element_text(angle = 0))
p4 = VlnPlot(obj_inh, group.by = 'family',split.by = 'ageGroup', 
             features = 'Grpr', cols = cols,
             log = T, assay = "MAGIC_RNA", pt.size = 0) + 
  xlab('Inhibitory neuron cluster') + ylim(c(NA, 1.2)) +
  ylab('Rel. Norm. Expression') + 
  theme(axis.text.x = element_text(angle = 0))

## put into 2 panel grid
pdf(here(PLOTDIR, 'plots', 'violinPlot_dh_neurons_grpr_ageGroup.pdf'), height = 6, width =7)
plot_grid(p3, p4, labels = c('A', 'B'), label_size = 16, ncol = 1, align = "v")
dev.off()


###############################################################################
## 4) plot the expression of T-type volt gated calcium change genes in dot plot
features = c('Grpr','Cacna1i', 'Cacna1g', 'Cacna1h')

p5 = DotPlot(obj_exc, features = features, assay = "RNA", col.min = 0,
             group.by = 'family', scale = T, scale.by = 'radius') + RotatedAxis() + 
  theme(legend.position = 'none', axis.text.x = element_blank(),
        axis.ticks.x=element_blank(), axis.title.x=element_blank()) + 
  ylab('Excitatory cluster')
p6 = DotPlot(obj_inh, features = features, assay = "RNA", col.min = 0,
             group.by = 'family', scale = T, scale.by = 'radius') + RotatedAxis()+ 
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 0)) + 
  ylab('Inhibitory cluster') + xlab('Gene') +
  guides(color = guide_colourbar(title.position="top", 'Average Expression'),
         size = guide_legend(title.position="top",'Percent Expressed' ))

## put the 2 dot plots aligned vertically
pdf(here(PLOTDIR, 'plots', 'dotPlot_dh_neurons_grpr_cacna1.pdf'), height = 6, width =7)
plot_grid(p5, p6, labels = c('A', 'B'), label_size = 16, 
          ncol = 1, align = "v", rel_heights = c(1, 1.4))
dev.off()



