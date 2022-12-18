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

PLOTDIR='figures/exploratory'
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)

############################################################
## 0) load the object with the dorsal horn neuron familys
obj_neuron = LoadH5Seurat(here('data/tidy_data/rdas/Russ_et_al_dorsal_horn_neurons.h5seurat'))
cols = paletteDiscrete(values = unique(obj_neuron$family))

## Grpr is excitatory, only look inside these neurons
obj_exc = obj_neuron[, obj_neuron$coarse_cell_types == 'Excit']

############################################################
## 1) scatter plot of Grpr vs. t-type calcium channel genes 
scatter_theme = theme(legend.position = 'none') + theme_bw(base_size = 8)

## make sure to use the imputed gene counts (MAGIC_RNA or MAGIC_SCT)
DefaultAssay(obj_exc) = 'MAGIC_RNA'; Idents(obj_exc) = 'family'

## gather the genes for all the data along w/ the Russ et al families 
df = data.frame(family = Idents(obj_exc), 
                Grpr = FetchData(obj_exc,"Grpr"),
                Cacna1i = FetchData(obj_exc,"Cacna1i"), 
                Cacna1g = FetchData(obj_exc,"Cacna1g"), 
                Cacna1h = FetchData(obj_exc,"Cacna1h")) %>% 
  pivot_longer(-c(family, Grpr))

## put the 2 dot plots aligned vertically
pdf(here(PLOTDIR, 'plots', 'scatterPlot_dh_neurons_grpr_cacna_isoforms.pdf'), height = 6, width =11)
ggplot(df, aes(x= Grpr , y= value)) + geom_point(size=1, aes(col = family)) +
  facet_grid(name ~ family) +
  stat_cor(method = "pearson", label.y.npc="top", label.x.npc = "left",  
           label.sep = "\n", r.digits = 2, p.digits = 2)  + 
  scale_color_manual(values = cols, guide = 'none') + scatter_theme + 
  ylab('Relative expression level, t-Type calcium channel gene') + 
  xlab('Relative expression level, Grpr')
dev.off()


