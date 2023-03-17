# Fig 2
#####
# Fig. 2a
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

# load seurat object :
so = LoadH5Seurat('main.h5seurat')

#Fig. 2a(left)
c1 <- c('DC1', 'DC2', 'DC3', 'pDC', 'ASDC', 'DC_prolif',
        'CD83_CD14_mono', 'CD14_mono', 'CD16_mono', 'C1_CD16_mono',
        'CD4.Prolif','Mono_prolif')
Mey <- so[,so@meta.data$full_clustering %in% c1]
DotPlot(object = Mey, group.by = 'full_clustering',assay = 'RNA',
        features = c('CLEC9A', 'CADM1', 'CLEC10A','CD1C', 'CD14', 'VCAN',
                     'CCR7', 'LAMP3', 'AXL', 'SIGLEC6','LILRA4', 'ITM2C',
                     'GZMB','IL1B', 'IER3', 'LDLR', 'CD83','S100A12',
                     'CSF3R', 'FCGR3A', 'MS4A7', 'LILRB1', 'CSF1R', 
                     'CDKN1C','C1QA', 'C1QB', 'C1QC', 'CCR1',
                     'MARCO', 'MKI67', 'TOP2A'), col.min = 0)+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_colour_gradient2(low = "#FFFFFF", high = "#073A7B")+
  geom_point(aes(size=pct.exp), shape = 21, colour="grey35", stroke=0.5)

# creat DotPlot for all cells of BALF
Bso = LoadH5Seurat('Bso2.h5Seurat')
DotPlot(object = Bso, group.by = 'full_clustering',assay = 'RNA',
        features = c('CLEC9A', 'CADM1', 'CLEC10A','CD1C', 'CD14', 'VCAN',
                     'CCR7', 'LAMP3', 'AXL', 'SIGLEC6','LILRA4', 'ITM2C',
                     'GZMB','IL1B', 'IER3', 'LDLR', 'CD83','S100A12',
                     'CSF3R', 'FCGR3A', 'MS4A7', 'LILRB1', 'CSF1R', 
                     'CDKN1C','C1QA', 'C1QB', 'C1QC', 'CCR1',
                     'MARCO', 'MKI67', 'TOP2A'), col.min = 0)+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_colour_gradient2(low = "#FFFFFF", high = "#073A7B")+
  geom_point(aes(size=pct.exp), shape = 21, colour="grey35", stroke=0.5)

#Fig. 2a(left)
DotPlot(object = Mey, group.by = 'full_clustering',assay = 'RNA',
        features = c('AB-CD141', 'AB-CLEC9A', 'AB-KIT','AB-BTLA',
                     'AB-CD1C', 'AB-CD101', 'AB-FcERIa',
                     'AB-CD5','AB-CD123', 'AB-CD45RA','AB-CD304',
                     'AB-CD14', 'AB-CD99', 'AB-CD64','AB-CX3CR1',  
                     'AB-CR1', 'AB-ITGAM','AB-CD16', 'AB-C5AR1'),
        col.min = -0.25, col.max = 0.75)+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_colour_gradient2(low = "#FFFFFF", high = "#A00602")+
  geom_point(aes(size=pct.exp), shape = 21, colour="grey35", stroke=0.5)

######
#Fig. 2b
c1 <- c('DC1', 'DC2', 'DC3', 'pDC', 'ASDC', 'DC_prolif',
        'CD83_CD14_mono', 'CD14_mono', 'CD16_mono', 'C1_CD16_mono',
        'CD4.Prolif','Mono_prolif')
Mey <- so[,so@meta.data$full_clustering %in% c1]
Status <- Mey@meta.data$Status_on_day_collection_summary


Cell_type <- factor(Mey@meta.data$full_clustering, levels = c1)
data <- data.frame(Status,Cell_type)
data <- data[!data$Status %in% c("Non_covid","LPS_90mins",
                                 "LPS_10hours"),]
c3 <- factor(data$Status, 
             levels = c("Healthy","Asymptomatic","Mild","Moderate","Severe","Critical"))

data$Status <- c3
S1 <- unique(data$Status)
D1 <- as.data.frame(prop.table(table(data[data$Status == S1[1],]$Cell_type)))
D1$status <- S1[1]
for (i in 2:length(S1)) {
  D2 <- as.data.frame(prop.table(table(data[data$Status == S1[i],]$Cell_type)))
  D2$status <- S1[i]
  D1 <- rbind(D1,D2)
}
colnames(D1)[1:2] <- c("Cell_type","Cell_proportions")

c2 <- c(brewer.pal(n = 12, name = "Set3"))

# Stacked
ggplot(D1, aes(fill = Cell_type, y = Cell_proportions, x = status)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = c2) +
  coord_flip()+ theme_bw()+
  theme(axis.text = element_text(size = 15, family = "Serif", face = "bold"))

#####
#Fig. 2d

# prepare BALF data
Bso = LoadH5Seurat('Bso2.h5Seurat')
BAL_mac <- Bso[,Bso@meta.data$full_clustering == 'Macrophages']
BAL_mac <- NormalizeData(BAL_mac, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
BAL_mac <- FindVariableFeatures(BAL_mac,selection.method = "vst", 
                                nfeatures = 2000)
BAL_mac <- ScaleData(BAL_mac, features = rownames(BAL_mac))

# prepare BPMC data
so = LoadH5Seurat('test.h5seurat')
cell1 <- c("CD83_CD14_mono","CD14_mono","CD16_mono",
           "C1_CD16_mono")
CD_data <- so[,so@meta.data$full_clustering %in% cell1]
CD_data <- NormalizeData(CD_data, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
CD_data <- FindVariableFeatures(CD_data,selection.method = "vst", 
                                nfeatures = 2000)
CD_data <- ScaleData(CD_data, features = rownames(CD_data))

# Create Combined Seurat object
data_combined <- merge(CD_data, y = BAL_mac, project = "PBMC_MALF",
                       merge.data = TRUE)
data_combined@meta.data$full_clustering <- replace(data_combined@meta.data$full_clustering,
                                                   is.na(data_combined@meta.data$full_clustering),"BAL_mac")

# Check Genes in Python Code and figure 2.d
gene_names = c('NFKBIA','FOSB','IL1B','S100A8','S100A9',
               'CD74','B2M','FTL','C1QA','C1QB','C1QC',
               'IL16','CCL2','CCL8','CCL7','CCL4L2','CCL3',
               'CCL4','CCL18','CXCL8','CXCL16', 'CXCL10',
               'KLF6', 'VIM', 'CD14', 
               'HLA-DPA1', 'HLA-DPB1', 'FCGR3A','CTSC', 
               'CTSL','CCR1','RSAD2','TNFSF10')

which(rownames(BAL_mac@assays$RNA) %in% gene_names)
which(rownames(CD_data@assays$RNA) %in% gene_names)

# Creat data for Heatmap (First appraoach)
DoHeatmap(object = CD_data, group.by = 'full_clustering',
          features = gene_names)


# Second approach
# Just for limitation in Hardware
CD_data@meta.data$full_clustering <- factor(CD_data@meta.data$full_clustering,
                                            levels = c("CD14_mono","CD83_CD14_mono","CD16_mono","C1_CD16_mono"))

NUM <- c(which(CD_data@meta.data$full_clustering == "CD14_mono"))
NUM <- c(NUM,which(CD_data@meta.data$full_clustering == "CD83_CD14_mono"))
NUM <- c(NUM,which(CD_data@meta.data$full_clustering == "CD16_mono"))
NUM <- c(NUM,which(CD_data@meta.data$full_clustering == "C1_CD16_mono"))

data_heat <- as.data.frame(as.matrix(CD_data@assays[["RNA"]]@counts[gene_names,NUM]))
col_ann <- as.character(CD_data@meta.data$full_clustering[NUM])

library(pheatmap)
pheatmap(data_heat, annotation_col = col_ann,width = 0.01,
         height = 0.5)
#cannot allocate vector of size 74.7 Gb

#####
# Fig. 2f

data <- read.delim("fig2f.txt")
data <- filter(data, full_clustering %in% c("bal_DC2","bal_Mac","CD1_CD16_mono",
                                            "CD14_mono","CD16_mono","CD83_CD14_mono",
                                            "DC1","DC2","DC3"),
               Status_on_day_collection_summary %in% c("bal_healthy","bal_mild",
                                                       "bal_severe","Healthy",
                                                       "Mild","Severe"))
M1 <- aggregate(TNF_score ~ full_clustering + Status_on_day_collection_summary,
                data, mean)
M2 <- aggregate(IFN_score ~ full_clustering + Status_on_day_collection_summary,
                data, mean)
M3 <- aggregate(JAK_score ~ full_clustering + Status_on_day_collection_summary,
                data, mean)

M1 <- arrange(M1, full_clustering, Status_on_day_collection_summary)
M2 <- arrange(M2, full_clustering, Status_on_day_collection_summary)
M3 <- arrange(M3, full_clustering, Status_on_day_collection_summary)
M <- cbind(M1,M2[,3]) %>%
  cbind(M3[,3])
colnames(M)[4:5] <- c(colnames(M2)[3],colnames(M3)[3]) 

bal_DC2 <- M[M$full_clustering == "bal_DC2",]
bal_Mac <- M[M$full_clustering == "bal_Mac",]
CD1_CD16_mono <- M[M$full_clustering == "CD1_CD16_mono",]
CD14_mono <- M[M$full_clustering == "CD14_mono",]
CD16_mono <- M[M$full_clustering == "CD16_mono",]
CD83_CD14_mono <- M[M$full_clustering == "CD83_CD14_mono",]
DC1 <- M[M$full_clustering == "DC1",]
DC2 <- M[M$full_clustering == "DC2",]
DC3 <- M[M$full_clustering == "DC3",]
bal_DC2 <- bal_DC2[,-1]
bal_Mac <- bal_Mac[,-1]
CD1_CD16_mono <- CD1_CD16_mono[,-1]
CD14_mono <- CD14_mono[,-1]
CD16_mono <- CD16_mono[,-1]
CD83_CD14_mono <- CD83_CD14_mono[,-1]
DC1 <- DC1[,-1]
DC2 <- DC2[,-1]
DC3 <- DC3[,-1]
rownames(bal_DC2) <- NULL
rownames(bal_Mac) <- NULL
rownames(CD1_CD16_mono) <- NULL
rownames(CD14_mono) <- NULL
rownames(CD16_mono) <- NULL
rownames(CD83_CD14_mono) <- NULL
rownames(DC1) <- NULL
rownames(DC2) <- NULL
rownames(DC3) <- NULL

bal_DC2 <- tibble::column_to_rownames(bal_DC2, "Status_on_day_collection_summary")
bal_Mac <- tibble::column_to_rownames(bal_Mac, "Status_on_day_collection_summary")
CD1_CD16_mono <- tibble::column_to_rownames(CD1_CD16_mono, "Status_on_day_collection_summary")
CD14_mono <- tibble::column_to_rownames(CD14_mono, "Status_on_day_collection_summary")
CD16_mono <- tibble::column_to_rownames(CD16_mono, "Status_on_day_collection_summary")
CD83_CD14_mono <- tibble::column_to_rownames(CD83_CD14_mono, "Status_on_day_collection_summary")
DC1 <- tibble::column_to_rownames(DC1, "Status_on_day_collection_summary")
DC2 <- tibble::column_to_rownames(DC2, "Status_on_day_collection_summary")
DC3 <- tibble::column_to_rownames(DC3, "Status_on_day_collection_summary")

bal_DC2 <- t(bal_DC2)
bal_Mac <- t(bal_Mac)
CD1_CD16_mono <- t(CD1_CD16_mono)
CD14_mono <- t(CD14_mono)
CD16_mono <- t(CD16_mono)
CD83_CD14_mono <- t(CD83_CD14_mono)
DC1 <- t(DC1)
DC2 <- t(DC2)
DC3 <- t(DC3)
colnames(bal_DC2) <- paste("bal_DC2",colnames(bal_DC2), sep = "_")
colnames(bal_Mac) <- paste("bal_Mac",colnames(bal_Mac), sep = "_")
colnames(CD1_CD16_mono) <- paste("CD1_CD16_mono",colnames(CD1_CD16_mono), sep = "_")
colnames(CD14_mono) <- paste("CD14_mono",colnames(CD14_mono), sep = "_")
colnames(CD16_mono) <- paste("CD16_mono",colnames(CD16_mono), sep = "_")
colnames(CD83_CD14_mono) <- paste("CD83_CD14_mono",colnames(CD83_CD14_mono), sep = "_")
colnames(DC1) <- paste("DC1",colnames(DC1), sep = "_")
colnames(DC2) <- paste("DC2",colnames(DC2), sep = "_")
colnames(DC3) <- paste("DC3",colnames(DC3), sep = "_")
l <- list(CD1_CD16_mono,CD14_mono,CD16_mono,CD83_CD14_mono,bal_Mac,
          DC1,DC2,DC3,bal_DC2)
M4 <- Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
row.names(M4) <- NULL
M4 <- tibble::column_to_rownames(M4,"rn")

col1 <- colorRampPalette(c("#AED3F0","white"))(57)
col2 <- colorRampPalette(c("white","#B53A37"))(270)
Heatmap(M4,col = c(col1,col2), width = unit(20, "cm"),
        height = unit(2.5, "cm"),
        cluster_rows = F, cluster_columns = F,
        name = "Enrichment score")


#####
#Fig. 2h  

P1 <- so[,so@meta.data$full_clustering == "Platelets"]
P1 <- P1[,!P1@meta.data$Status_on_day_collection_summary %in% c("LPS_90mins",
                                                                "Non_covid",
                                                                "LPS_10hours")]
unique(P1@meta.data$Status_on_day_collection_summary)
P1@meta.data$Status_on_day_collection_summary <- factor(P1@meta.data$Status_on_day_collection_summary,
                                                        levels = c("Healthy","Asymptomatic","Mild","Moderate",
                                                                   "Severe","Critical"))


DotPlot(object = P1, assay = 'RNA', group.by = 'Status_on_day_collection_summary',
        features = c('TIMP1', 'GP9', 'MPIG6B','PF4', 'CLEC1B', 'RAP1B',
                     'SRGN', 'PPBP', 'STXBP2', 'PFN1'), col.min = 0, col.max = 4)+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_colour_gradient2(low = "#FFFFFF", high = "#073A7B")+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

#####
# Fig. 2i (Part1)
# load Seurat object :
so = LoadH5Seurat('test.h5seurat')

#Creat new Seurat object based on HSC cells
c1 <- c("HSC_CD38pos","HSC_MK", "HSC_prolif", "HSC_erythroid",
        "HSC_CD38neg", "HSC_myeloid")

HSP <- so[,so@meta.data$full_clustering %in% c1]
HSP@meta.data$full_clustering <- factor(HSP@meta.data$full_clustering, 
                                        levels = c("HSC_CD38pos","HSC_CD38neg",
                                                   "HSC_erythroid","HSC_prolif",
                                                   "HSC_myeloid","HSC_MK","CD38pos_HSPC",
                                                   "CD38neg_HSPC","Erythroid_prog",
                                                   "proliferating_prog","Myeloid_prog",
                                                   "MK_prog"))

# Replace features' names based on article
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_CD38pos",
                                         "CD38pos_HSPC")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_CD38neg",
                                         "CD38neg_HSPC")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_MK",
                                         "MK_prog")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_prolif",
                                         "proliferating_prog")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_erythroid",
                                         "Erythroid_prog")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_myeloid",
                                         "Myeloid_prog")

HSP@meta.data$full_clustering <- factor(HSP@meta.data$full_clustering, 
                                        levels = c("CD38pos_HSPC","CD38neg_HSPC",
                                                   "Erythroid_prog","proliferating_prog",
                                                   "Myeloid_prog","MK_prog"))

# Explore data
HSP[["MTpercent"]] <- PercentageFeatureSet(HSP, pattern = "^MT-")
VlnPlot(HSP,features = c("nFeature_RNA","nCount_RNA","MTpercent"),
        ncol = 3)

# Select valid subset, Normalization and Scaling
HSP <- subset(HSP,subset = nFeature_RNA > 200 &  nFeature_RNA < 5000 &
                MTpercent < 10)
HSP <- NormalizeData(HSP, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
HSP <- FindVariableFeatures(HSP,selection.method = "vst", nfeatures = 2000)
HSP <- ScaleData(HSP, features = rownames(HSP))


# Dimention reduction
HSP <- RunPCA(HSP, features = VariableFeatures(HSP))
HSP <- JackStraw(HSP, num.replicate = 100)
HSP <- ScoreJackStraw(HSP, dims = 1:20)
JackStrawPlot(HSP, dims = 1:20)
ElbowPlot(HSP)

# Clustering
HSP <- FindNeighbors(HSP, dim = 1:18)
HSP <- FindClusters(HSP, resolution = 0.5)
Clus <- as.data.frame(table(HSP@meta.data$full_clustering,
                            HSP@meta.data$seurat_clusters))
Clus <- Clus[Clus$Freq > 0,]

# tSNE
HSP <- RunTSNE(HSP, dims = 1:18)
DimPlot(HSP,reduction = "tsne", pt.size = 1.5)

# UMAP
#Fig. 2i (part1)
HSP <- RunUMAP(HSP,dims = 1:18)
DimPlot(HSP,reduction = "umap", pt.size = 2,group.by = "full_clustering",
        cols = c("steelblue3","orange","red","violet","green",
                             "brown"))
                             
# Fig. 2i (Part2)
#Creat new Seurat object based on HSC cells
c1 <- c("HSC_CD38pos","HSC_MK", "HSC_prolif", "HSC_erythroid",
        "HSC_CD38neg", "HSC_myeloid")

HSP <- so[,so@meta.data$full_clustering %in% c1]
HSP@meta.data$full_clustering <- factor(HSP@meta.data$full_clustering, 
                                        levels = c("HSC_CD38pos","HSC_CD38neg",
                                                   "HSC_erythroid","HSC_prolif",
                                                   "HSC_myeloid","HSC_MK","CD38pos_HSPC",
                                                   "CD38neg_HSPC","Erythroid_prog",
                                                   "proliferating_prog","Myeloid_prog",
                                                   "MK_prog"))

# Replace features' names based on article
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_CD38pos",
                                         "CD38pos_HSPC")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_CD38neg",
                                         "CD38neg_HSPC")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_MK",
                                         "MK_prog")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_prolif",
                                         "proliferating_prog")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_erythroid",
                                         "Erythroid_prog")
HSP@meta.data$full_clustering <- replace(HSP@meta.data$full_clustering,
                                         HSP@meta.data$full_clustering == "HSC_myeloid",
                                         "Myeloid_prog")

HSP@meta.data$full_clustering <- factor(HSP@meta.data$full_clustering, 
                                        levels = c("CD38pos_HSPC","CD38neg_HSPC",
                                                   "Erythroid_prog","proliferating_prog",
                                                   "Myeloid_prog","MK_prog"))

DotPlot(HSP, assay = 'RNA', group.by = 'full_clustering',
        features = c('CD34', 'AVP', 'CD38','GATA1', 'MPO','PF4'))+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_colour_gradient2(low = "#FFFF00", high = "#073A7B")+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)