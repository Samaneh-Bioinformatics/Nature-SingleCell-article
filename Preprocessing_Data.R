devtools::install_github("satijalab/seurat-data")
devtools::install_github("mojaveazure/seurat-disk", force = T)
install.packages("Seurat")

library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Convert data to h5seurat to load in R 
Convert('covid_portal_210320_with_raw.h5ad'
        ,assay = 'RNA',dest = 'main.h5seurat')

# load seurat object :
so = LoadH5Seurat('main.h5seurat')