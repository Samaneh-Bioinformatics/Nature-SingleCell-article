# Fig 3
#####
# Fig. 3a
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

# load seurat object :
so = LoadH5Seurat('main.h5seurat')
c1 <- c(c("CD4.Naive", "CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Prolif"
          , "CD4.Th1", "CD4.Th2", "CD4.Th17", "CD4.Tfh",
          "Treg", "CD8.Naive", "CD8.Activated", "CD8.Prolif",
          "CD8.CM", "CD8.TE", "CD8.EM", "gdT", "MAIT", "NKT"))

Tcell <- so[,so@meta.data$full_clustering %in% c1]
DefaultAssay(Tcell) <- "RNA"
Tcell <- ScaleData(Tcell, verbose = FALSE)
Tcell <- FindVariableFeatures(Tcell,selection.method = 'vst', nfeatures = 3000)
Tcell <- RunPCA(Tcell, npcs = 50, verbose = FALSE)
Tcell <- RunUMAP(Tcell, reduction = "pca_harmony", dims = 1:50)

library(RColorBrewer)
DimPlot(Tcell,reduction = "umap", group.by = "full_clustering",
        cols = c(brewer.pal(9, "Set3"),brewer.pal(8, "Dark2")))

# Density
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
tcell.umap.merge <- as.data.frame(Tcell@reductions[["umap"]]@cell.embeddings)
all(rownames(tcell.umap.merge) == row.names(Tcell@meta.data))
tcell.umap.merge$Sub.Annotation <- Tcell@meta.data$full_clustering
tcell.umap.merge$CellID <- row.names(tcell.umap.merge)

dens.umap.list <- list()
loop.cells <- intersect(c1, unique(as.character(tcell.umap.merge$Sub.Annotation)))
for(i in seq_along(loop.cells)){
  i.c <- loop.cells[i]
  # change this to 1% of all cells
  i.n <- ceiling(sum(tcell.umap.merge$Sub.Annotation %in% i.c)/100)
  if(i.n < 100){
    i.n <- 100
  }
  i.dens <- get_density(tcell.umap.merge$UMAP_1[tcell.umap.merge$Sub.Annotation %in% i.c],
                        tcell.umap.merge$UMAP_2[tcell.umap.merge$Sub.Annotation %in% i.c], n=i.n)
  
  dens.umap.list[[i.c]] <- data.frame("Sub.Annotation"=i.c,
                                      "CellID"=tcell.umap.merge$CellID[tcell.umap.merge$Sub.Annotation %in% i.c],
                                      "Dens"=i.dens)
}

umap.dens <- do.call(rbind.data.frame, dens.umap.list)
umap.dens.merge <- merge(tcell.umap.merge, umap.dens, by=c('Sub.Annotation', 'CellID'))

library(ggsci)
library(ggplot2)
library(scattermore)
library(viridis)
library(cowplot)
ggplot(umap.dens.merge,
       aes(x=UMAP_1, y=UMAP_2)) +
  geom_scattermore(data=umap.dens.merge[, c("UMAP_1", "UMAP_2")],
                   colour='grey80', alpha=0.2) +
  geom_scattermore(data=umap.dens.merge[umap.dens.merge$Sub.Annotation %in% c("CD4.Tfh"), ],
                   aes(colour=Dens)) +
  scale_colour_viridis() +
  theme_cowplot() +
  theme(legend.key.size = unit(0.2, "cm"),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6)) +
  guides(colour=guide_colourbar(title="CD4.Tfh density", override.aes = list(size=3)))






#####
# Fig. 3b
DotPlot(object = Tcell, group.by = 'full_clustering',assay = 'RNA',
        features = c("CD3G", "CD3E", "FCGR3A", "CD8A", "CD8B",
                     "CD4", "CD80", "CD86", "CXCL10",
                     "VEGFA","GATA3", "PECAM1", "RORC",
                     "CCR7", "CCR2", "CCR5", "CXCR3", "CD36", 
                     "FOXP3", "MS4A1", "LYZ", "GZMB", "CD68", 
                     "FCGR2A", "FCGR2B", "FCGR3B", "TPO", "HBB",
                     "NCR1", "NCR2", "F8", "HIF1A", "CTLA4", "CD28","CXCR1",
                     "PDCD1", "LAG3", "TYMS", "PITX1", "TIGIT",
                     "CXCL13", "BCL6","CD44", "IKZF2","CD27","F2","F5", 
                     "IKZF4","ITGA2B", "ITGAX", "TOX", "FYN", "NCAM1", 
                     "LRRC32","SELL", "CXCR5", "CSF1", "IL7R", "CD40LG",
                     "MYC", "MKI67","IFNG", "IL6", "THRB", "CD19", "CD38"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

# Antibody protein
ADT=Tcell[(Tcell@assays[["RNA"]]@meta.features[["feature_types"]] =='Antibody Capture'),]
adt_assay = GetAssay(ADT)
Tcell[['adt']] = adt_assay
DefaultAssay(Tcell) <- "adt"

DotPlot(object = Tcell, group.by = 'full_clustering',assay = 'adt',
        features = c("AB-CD3", "AB-CD4", "AB-CCR7","AB-CD45RA",
                     "AB-CD45RO","AB-CD27", "AB-CD28","AB-CD62L",
                     "AB-CD25","AB-CTLA4","AB-CXCR5", "AB-CD40LG",
                     "AB-ICOS","AB-CXCR3","AB-CD8","AB-CD274","AB-TCR-Vg9",
                     "AB-TCR-Vg2","AB-TCR-Va7.2","AB-TCR-Va24-Ja18",
                     "AB-CD56","AB-CD16"),col.min = 0,col.max = 1)+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_colour_gradient2(low = "#FFFFFF", high = "red")+
  geom_point(aes(size=pct.exp), shape = 21, colour="grey35", stroke=0.5)

#####
# Fig.3c
# load seurat object :

cell <- c("IL2","IL1A","IL1B","IL10","IL22","TNF","INFG",
          "LTF","IL4","IL5","IL13","IL21","IL17A","IL17F",
          "IL12A","IL12B")

DotPlot(object = Tcell, group.by = 'full_clustering',assay = 'RNA',
        features = cell, col.min = 0,col.max = 1)+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  geom_point(aes(size=pct.exp), shape = 21, colour="grey35", stroke=0.5)
