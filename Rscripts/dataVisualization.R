######
# Script : Visualizacion grafica de los resultados de DEG
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 27/02/2024
# Description: El siguiente script nos permite realiza el Analisis de Terminos GO
# a partir de los datos provenientes del Analisis de DEG
# Primero correr el script "DEG_analysis.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: 
#       - dds_Times_vs_control.RData (dds), 
#       - vst_Times_vs_control.RData (vsdata) 
#       - archivos de salida de DEG en formato CSV (res_15t, res_30t, res_4t) 
#   - Output: Volcano plot y Heatmap
#######

# --- Load packages ----------
library(dplyr)
library(pheatmap)
library(ggplot2)

# --- Load data -----
# Cargar archivos
figdir <- '/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/figures/'

#Cargar variable "dds", proveniente del script "DEG_analysis.R"
load("/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/dds_treatment_vs_control.RData")

#Cargar variable "vsdata", proveniente del script "DEG_analysis.R" 
load("/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/vst_treatment_vs_control.RData") 

#Cargar variable "res_30t", proveniente del script "DEG_analysis.R"
#load("/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/results/DE_30min_vs_control.csv") 
load("/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/DE_treatment_vs_control.csv") 

# ---- volcano plot ----
df <- as.data.frame(res_t)
# padj 0.05 y log2FoldChange de 2
df <- df %>% 
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -(2) & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))

# visualizacion
png(file = paste0(figdir, "VolcanoPlot_treatment_vs_control.png"))

ggplot(df, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 0.7) +
  labs(title = "treatment vs control") +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-adj")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = -(2), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5)

dev.off()

# --- Heatmap (vsd) -----
topGenes <- head(order(res_t$padj), 20) # Obtener el nombre de los 20 genes con p value mas significativo

png(file = paste0(figdir, "Heatmap_vsd_topgenes.png"))
pheatmap(assay(vsdata)[topGenes,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE)
dev.off()

# --- Heatmap  (por contrastes) (log2 Fold Change) -----
betas <- coef(dds)
colnames(betas)

mat <- betas[topGenes, -c(1,2)] # crear la matriz con el topgene de genes

# Filtro de 3 log2foldchange
thr <- 1
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

# Almacenar la grafica
png(file = paste0(figdir, "Heatmap_log2FoldChage_topgenes.png"))
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
dev.off()

# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
