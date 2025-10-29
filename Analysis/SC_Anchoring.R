rm(list = ls())
# BiocManager::install('harmony')
library(harmony)
library(Seurat)
library(tidyverse)
library(data.table)
library(future)
library(SingleR)
library(scuttle)
# BiocManager::install('pheatmap')
library(pheatmap)
library(scran)
## Loading Individual Datasets

set.seed(9703)

gc()

# plan("multiprocess", workers = 32)
# plan()

GSE134355 <- readRDS('./GSE134355_filtered.rds')

GSE182416 <- readRDS('./GSE182416_filtered.rds') %>%
  as.SingleCellExperiment() %>%
  logNormCounts()

GSE191288 <- readRDS('./GSE191288_filtered.rds')

GSE134355$Study <- c(rep('GSE134355', ncol(GSE134355)))
GSE182416$Study <- c(rep('GSE182416', ncol(GSE182416)))
GSE191288$Study <- c(rep('GSE191288', ncol(GSE191288)))

## Trying SingleR annotation based on GSE182416 dataset's cell type annotation.
head(GSE182416)

GSE182416$celltype_global %>% head

M.Data <- merge(x = GSE134355, 
                y = GSE191288, 
                add.cell.id = c('GSE134355', 'GSE191288')) %>%
  as.SingleCellExperiment() %>%
  logNormCounts()

str(M.Data)
head(M.Data)

Pred.Thyr <- SingleR(test = M.Data,
                     ref = GSE182416,
                     labels = GSE182416$celltype_global,
                     de.method = 'wilcox')

table(Pred.Thyr$labels)

my.heat <- plotScoreHeatmap(Pred.Thyr)
par(cex.lab = 0.5, cex.axis = 0.5, cex.main = 0.6, cex.sub = 0.5)


png(filename = 'AnnotationHeatmap.png', width = 5, height = 5, units = 'in', res = 300)
print(my.heat)
dev.off()


plotDeltaDistribution(Pred.Thyr, ncol = 3)


summary(is.na(Pred.Thyr$pruned.labels))

# Good quality: 10,159  Poor quality: 44

all.markers <- metadata(Pred.Thyr)$de.genes
head(all.markers)

M.Data$labels <- Pred.Thyr$labels
M.Data$pruned.labels <- Pred.Thyr$pruned.labels

library(scater)

plotHeatmap(M.Data, order_columns_by = 'labels',
            features = unique(unlist(all.markers$`T/NK cells`)))


Test.Data <- M.Data[, !is.na(colData(M.Data)$pruned.labels)] %>%
  as.Seurat
dim(Test.Data)

Test.Data$pruned.labels %>% unique

head(Test.Data)
colData(Test.Data) %>% head

GSE182416 <- GSE182416 %>%
  as.Seurat

GSE182416$labels <- GSE182416$celltype_global
GSE182416$pruned.labels <- GSE182416$celltype_global

head(GSE182416)

Comb.Data <- merge(x = GSE182416,
                   y = Test.Data)

dim(Comb.Data)

head(Comb.Data)

gc()

Comb.Data <- SCTransform(Comb.Data,
                         vars.to.regress = c('percent.mt', 'Study'))

gc()

Thyroid_data <- Comb.Data %>%
  RunPCA(.)

ElbowPlot(Thyroid_data)

DimHeatmap(Thyroid_data, dims = 1:6, balanced = T)

# Only after saving the file in the end of this page
# Thyroid_data <- read_rds('./Thyroid_Harmony.rds')

Thyroid_data <- Thyroid_data %>%
  RunUMAP(., dims = 1:6, reduction = 'pca') %>%
  FindNeighbors(., reduction = 'pca') %>%
  FindClusters(., resolution = 0.4)

gc()

gc()

celltype_colors <- c(
  "B cells" = "#1F77B4",            # Blue
  "Endothelial cells" = "#FF7F0E",  # Orange
  "Epithelial cells" = "#2CA02C",   # Green
  "Fibroblasts" = "#D62728",        # Red
  "Myeloid cells" = "#9467BD",      # Purple
  "Proliferating cells" = "#BCBD22",# Olive
  "SMCs/Pericytes" = "#17FDAA",     # Cyan
  "T/NK cells" = "#7F7F7F"          # Gray
)


p1 <- DimPlot(Thyroid_data, group.by='labels', reduction='umap',
              cols = celltype_colors) +
  ggtitle('Seurat Clusters\n(Before Correction)') +
  labs(x = 'UMAP_1', y = 'UMAP_2')

Thyroid_data <- Thyroid_data %>%
  RunUMAP(., dims = 1:6, reduction = 'pca', reduction.name = 'umap.pca') # <-- reduction.name 지정

p2 <- DimPlot(Thyroid_data, group.by='Study', reduction='umap.pca') + 
  ggtitle('Sequencing Batch\n(Before Correction)') +
  labs(x = 'UMAP_1', y = 'UMAP_2') +
  theme(legend.position = 'none')

p1 | p2

Thyroid_data <- Thyroid_data %>%
  RunHarmony(., 'Study', assay.use = 'SCT')

ElbowPlot(Thyroid_data, reduction = 'harmony')

Thyroid_data <- Thyroid_data %>%
  RunUMAP(., reduction = 'harmony', dims = 1:11) %>%
  FindNeighbors(., reduction = 'harmony') %>%
  FindClusters(., resolution = 0.4)


p3 <- DimPlot(Thyroid_data, group.by='labels', reduction='umap',
              cols = celltype_colors) +
  ggtitle('Cell Annotation\n(After Correction)') +
  labs(x = 'UMAP_1', y = 'UMAP_2')

p4 <- DimPlot(Thyroid_data, group.by='Study', reduction='umap') +
  ggtitle('Sequencing Batch\n(After Correction)') +
  labs(x = 'UMAP_1', y = 'UMAP_2')


ggsave(p3, filename = './Revised_Cell_Annotation_plot.png', width = 7, height = 6, dpi = 300, units = 'in')

p1 | p3

ggsave('Compare_Clust.png', width = 14, height = 6, dpi = 300, units = 'in')

p2 | p4

ggsave('Correction_effect.png', width = 13, height = 6, dpi = 300, units = 'in')

saveRDS(Thyroid_data, file = './Thyroid_Harmony.rds')

gc()

head(Thyroid_data)

Thyroid_data






