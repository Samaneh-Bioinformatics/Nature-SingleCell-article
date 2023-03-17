# Fig 1
#####
# Fig. 1b
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(RColorBrewer)

# load seurat object :
so = LoadH5Seurat('main.h5seurat')

# Umap
DimPlot(so,reduction = "umap", group.by = "initial_clustering",label = T,
        cols = c(brewer.pal(10, "Set3"),brewer.pal(9, "Dark2")))

#####
# Fig. 1c
so@meta.data$cluster1 <- so@meta.data$initial_clustering
so@meta.data$cluster1 <- factor(so@meta.data$cluster1,
                                levels = c(as.character(unique(so@meta.data$cluster1))
                                           ,'CD14_mono','CD16_mono','HSPC'))
for (i in seq(nrow(so@meta.data))) {
  if (so@meta.data$full_clustering[i] == "CD14_mono") 
  {so@meta.data$cluster1[i] <- "CD14_mono"}
  if (so@meta.data$full_clustering[i] == "CD16_mono") 
  {so@meta.data$cluster1[i] <- "CD16_mono"}
  if (so@meta.data$cluster1[i] == "HSC") 
  {so@meta.data$cluster1[i] <- "HSPC"}
}

so@meta.data$cluster1 <- factor(so@meta.data$cluster1,
                                levels = c(as.character(unique(so@meta.data$cluster1))))
Status <- so@meta.data$Status_on_day_collection_summary
Cell_type <- so@meta.data$cluster1
data <- data.frame(Status,Cell_type)
data <- data[data$Status != "Non_covid",]

data$Status <- factor(data$Status,levels =c("Healthy","Asymptomatic","Mild",
                                            "Moderate","Severe","Critical"
                                            ,"LPS","LPS_90mins","LPS_10hours"))

S1 <- unique(data$Status)
D1 <- as.data.frame(prop.table(table(data[data$Status == S1[1],]$Cell_type)))
D1$status <- S1[1]
for (i in 2:length(S1)) {
  D2 <- as.data.frame(prop.table(table(data[data$Status == S1[i],]$Cell_type)))
  D2$status <- S1[i]
  D1 <- rbind(D1,D2)
}
colnames(D1)[1:2] <- c("Cell_type","Cell_proportions")

c2 <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 12, name = "Set3"))

# Stacked
ggplot(D1, aes(fill = Cell_type, y = Cell_proportions, x = status)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = c2) +
  coord_flip()+ theme_bw()+
  theme(axis.text = element_text(size = 15, family = "Serif", face = "bold"))

######
# Fig. 1d
GO <- read.delim("Human_GO-0034340.txt")
Genes <- intersect(GO$Symbol_synonym, rownames(so))

