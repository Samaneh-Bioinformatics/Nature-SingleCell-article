# Fig 4
#####
# Fig. 4a
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

# load seurat object :
so = LoadH5Seurat('main.h5seurat')
cell <- c("B_naive","Plasma_cell_IgG","B_non-switched_memory","Plasma_cell_IgM",
          "Plasmablast","B_exhausted","B_immature","B_switched_memory",
          "Plasma_cell_IgA")
Bcells <- so[,so@meta.data$full_clustering %in% cell]
Bcells <- ScaleData(Bcells, verbose = FALSE)
Bcells <- FindVariableFeatures(Bcells,selection.method = 'vst', nfeatures = 3000)

# Dimention reduction
Bcells <- RunPCA(Bcells, npcs = 50, verbose = FALSE)
Bcells <- JackStraw(Bcells, num.replicate = 100)
Bcells <- ScoreJackStraw(Bcells, dims = 1:20)
JackStrawPlot(Bcells, dims = 1:20)
ElbowPlot(Bcells)

# tSNE
Bcells <- RunTSNE(Bcells, dims = 1:30)
DimPlot(Bcells,reduction = "tsne", pt.size = 1.5, 
        cols = c("steelblue3","orange","red","violet","green",
                             "brown","blue","yellow2","olivedrab"),
                             group.by = "full_clustering")

library(RColorBrewer)
DimPlot(Bcells,reduction = "tsne", pt.size = 1.5, 
        cols = c(brewer.pal(12, "Set3"),brewer.pal(8, "Dark2")),
        group.by = "seurat_clusters")

# UMAP
Bcells <- RunUMAP(Bcells, reduction = "pca_harmony", dims = 1:50)
DimPlot(Bcells,reduction = "umap", pt.size = 2,group.by = "full_clustering",
        cols = c("steelblue3","orange","red","violet","green",
                             "brown","blue","yellow2","olivedrab"))
                             library(RColorBrewer)
DimPlot(Bcells,reduction = "umap", pt.size = 1.5, 
        cols = c(brewer.pal(9, "Set3"),brewer.pal(8, "Dark2")),
        group.by = "seurat_clusters")

######
# Fig.4b

DotPlot(object = Bcells, group.by = 'full_clustering',assay = 'RNA',
        features = c("MS4A1","CD19","CD40","CD69","CD22","FCER2","CD24",
                     "CR2","MME","MKI67","IGHD","IGHM","IGHA1","IGHG1",
                     "MZB1","CD27","TNFRSF13B","CD38","SDC1"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_colour_gradient2(low = "#FFFFFF", high = "blue4")+
  geom_point(aes(size=pct.exp), shape = 21, colour="gray35", stroke=0.5)

#####
# Fig.4c
Bcells_in <- Bcells[,Bcells@meta.data$Resample %in% c("Initial")]
table(Bcells@meta.data$Status_on_day_collection_summary)


# Replace Status_on_day_collection_summary' names based on article
Bcells_in@meta.data$Status_on_day_collection_summary <- 
  factor(Bcells_in@meta.data$Status_on_day_collection_summary,
         levels =c("Healthy","Asymptomatic","Mild","Moderate","Severe",
                   "Critical","Non_covid","LPS","LPS_90mins","LPS_10hours"))

Bcells_in@meta.data$Status_on_day_collection_summary <- 
  replace(Bcells_in@meta.data$Status_on_day_collection_summary,
          Bcells_in@meta.data$Status_on_day_collection_summary %in%
            c("LPS_10hours","LPS_90mins"),"LPS")
Bcells_in@meta.data$Status_on_day_collection_summary <- 
  factor(Bcells_in@meta.data$Status_on_day_collection_summary,
         levels =c("Healthy","Asymptomatic","Mild","Moderate","Severe",
                   "Critical","Non_covid","LPS"))

library(ggplot2)
library(RColorBrewer)

Status <- Bcells_in@meta.data$Status_on_day_collection_summary
Cell_type <- as.character(Bcells_in@meta.data$full_clustering)
data <- data.frame(Status,Cell_type)
S1 <- unique(data$Status)
D1 <- as.data.frame(prop.table(table(data[data$Status == S1[1],]$Cell_type)))
D1$status <- S1[1]
for (i in 2:length(S1)) {
  D2 <- as.data.frame(prop.table(table(data[data$Status == S1[i],]$Cell_type)))
  D2$status <- S1[i]
  D1 <- rbind(D1,D2)
}
colnames(D1)[1:2] <- c("Cell_type","Cell_proportions")
D1 <- D1[D1$Cell_proportions > 0,]

c2 <- c(brewer.pal(n = 9, name = "Dark2"), brewer.pal(n = 10, name = "Set3"))

# Stacked
ggplot(D1, aes(fill = Cell_type, y = Cell_proportions, x = status)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = c2) +theme_bw()+
  theme(axis.text = element_text(size = 15, family = "Serif", face = "bold"))

#####
#Fig.4d
df <- Bcells_in@meta.data %>% 
  group_by(sample_id, Status_on_day_collection_summary, full_clustering) %>%
  summarise(total = n()) %>% mutate(freq = total / sum(total))
df

df1 <- df[df$full_clustering %in% c("Plasma_cell_IgA",
                                    "Plasma_cell_IgG",
                                    "Plasma_cell_IgM",
                                    "Plasmablast"),]

df1 <- df1[!df1$Status_on_day_collection_summary %in% c("Non_covid","LPS"),]
df1$Status_on_day_collection_summary <- 
  factor(df1$Status_on_day_collection_summary, levels = c("Healthy",
                                                          "Asymptomatic",
                                                          "Mild","Moderate","Severe",
                                                          "Critical"))
df1$sample_id <- as.character(df1$sample_id)
length(unique(df1$sample_id))

df_IgA <- df1[df1$full_clustering == "Plasma_cell_IgA",]
df_IgG <- df1[df1$full_clustering == "Plasma_cell_IgG",]
df_IgM <- df1[df1$full_clustering == "Plasma_cell_IgM",]
df_Plas <- df1[df1$full_clustering == "Plasmablast",]

pdf(file = "results/Fig.4d.pdf", width = 5.5, height = 3.5)
boxplot(freq ~ Status_on_day_collection_summary, 
        data = df_IgA, col = "white", main="Plasma_cell_IgA")
# Points
stripchart(freq ~ Status_on_day_collection_summary,
           data = df_IgA,
           method = "jitter",
           pch = 19,
           col = 2:4,
           vertical = TRUE,
           add = TRUE)

boxplot(freq ~ Status_on_day_collection_summary, 
        data = df_IgG, col = "white", main="Plasma_cell_IgG")
# Points
stripchart(freq ~ Status_on_day_collection_summary,
           data = df_IgG,
           method = "jitter",
           pch = 19,
           col = 2:4,
           vertical = TRUE,
           add = TRUE)

boxplot(freq ~ Status_on_day_collection_summary, 
        data = df_IgM, col = "white", main="Plasma_cell_IgM")
# Points
stripchart(freq ~ Status_on_day_collection_summary,
           data = df_IgM,
           method = "jitter",
           pch = 19,
           col = 2:4,
           vertical = TRUE,
           add = TRUE)

boxplot(freq ~ Status_on_day_collection_summary, 
        data = df_Plas, col = "white", main="Plasmablast")
# Points
stripchart(freq ~ Status_on_day_collection_summary,
           data = df_Plas,
           method = "jitter",
           pch = 19,
           col = 2:4,
           vertical = TRUE,
           add = TRUE)

dev.off()
