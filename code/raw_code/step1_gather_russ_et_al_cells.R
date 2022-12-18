library(Seurat)
library(SeuratDisk)
library(harmony)
library(SingleCellExperiment)
library(tidyverse)
library(here)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

####################################################
## 1) find the raw mouse spinal cord snRNA-seq data
dir.create(here('data/tidy_data/rdas'), showWarnings = F)
h5Seurat_file = here('data/tidy_data/rdas/final_cluster_assignment.h5seurat')
if(!file.exists(h5Seurat_file)){
  obj = readRDS('/projects/pfenninggroup/singleCell/MacaqueMouse_SealDorsalHorn_snRNA-Seq/mmData/final_cluster_assignment.rds') %>% UpdateSeuratObject()
  SaveH5Seurat(obj, h5Seurat_file, overwrite = T)
}

## clean up the cell types labels 
obj <- LoadH5Seurat(h5Seurat_file)
head(obj[[]])

## drop the bad cells
obj = obj[, !obj$final_cluster_assignment %in% c('Junk', 'Doublets')]
obj = obj %>% RunPCA() %>% RunUMAP(dims = 1:30)

## coarse cell types aren't all correct
obj$coarse_cell_types = obj$final_cluster_assignment %>% 
  as.character() %>% make.names() %>% ss('\\.')
obj$final_cluster_assignment = obj$final_cluster_assignment %>% 
  as.character() %>% str_trim() %>% make.names()
table(obj$final_cluster_assignment)
table(obj$coarse_cell_types)

## add the mouse families
families = read_csv('/projects/pfenninggroup/singleCell/MacaqueMouse_SealDorsalHorn_snRNA-Seq/mmData/cluster_families.csv') %>%  mutate_all(make.names) %>% deframe()
families[families == 'CC.'] = 'CC'

table(obj$final_cluster_assignment[!obj$final_cluster_assignment %in% names(families)])
obj$family = families[obj$final_cluster_assignment ]
table(obj$family, obj$coarse_cell_types)

SaveH5Seurat(obj, h5Seurat_file, overwrite = T)

#################################
## 3) subset to just the DH cells
h5Seurat_file2 = here('data/tidy_data/rdas/Russ_et_al_dorsal_horn_cells.h5seurat')
obj_dh = obj[, !obj$family %in% c('VI', 'VE', 'MI', 'ME', 'MN', 'Lmx1b.ME')]
obj_dh = obj_dh %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
SaveH5Seurat(obj_dh, h5Seurat_file2, overwrite = T)


###################################
## 4) subset to just the DH neurons
h5Seurat_file3 = here('data/tidy_data/rdas/Russ_et_al_dorsal_horn_neurons.h5seurat')
obj_dhn = obj_dh[, obj_dh$coarse_cell_types %in% c('Excit', 'Inhib')]
obj_dhn = obj_dhn %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

table(obj_dhn$family, obj_dhn$coarse_cell_types)
table(obj_dhn$family, obj_dhn$dataset)
head(obj_dhn[[]])

## use SCTransform normalized counts across these datasets
obj_dhn[["percent.rb"]] <- PercentageFeatureSet(object = obj_dhn, pattern = "^Rp[sl]")
obj_dhn <- SCTransform(obj_dhn, method = "glmGamPoi", 
                       vars.to.regress = c('percent.rb','percent.mt'))


## smooth out the genes with imputation across full data
DefaultAssay(obj_dhn) <- "SCT"
obj_dhn <- magic(obj_dhn, genes='all_genes')
DefaultAssay(obj_dhn) <- "MAGIC_SCT"

## smooth out the genes with imputation across full data
DefaultAssay(obj_dhn) <- "RNA"
obj_dhn <- magic(obj_dhn, genes='all_genes')
DefaultAssay(obj_dhn) <- "MAGIC_RNA"

SaveH5Seurat(obj_dhn, h5Seurat_file3, overwrite = T)







