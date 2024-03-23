######
# Script : Analisis de expresion diferencial
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 27/02/2024
# Description: El siguiente script nos permite realiza el Analisis de expresion Diferencial
# a partir de los datos provenientes del alineamiento de STAR a R,
# Primero correr el script "load_data_inR.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: Cargar la variable raw_counts.RData que contiene la matriz de cuentas y la metadata
#   - Output: DEG
#######

# --- Load packages ----------
library(DESeq2)

# --- Load data -----
# Cargar archivos
outdir <- "/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/"
figdir <- '/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/figures/'

#Cargar variable "counts", proveniente del script "load_data_inR.R"
load("/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/counts/raw_counts.RData") 
samples <- metadata$sample_id # Extraer los nombres de los Transcriptomas
metadata$type <- as.factor(metadata$type) # convertir a factor

# --- DEG ----
counts <- counts[which(rowSums(counts) > 10),] #Seleccionamos genes con mas de 10 cuentas

# Convertir al formato dds
dds <- DESeqDataSetFromMatrix(countData =  counts, 
            colData = metadata, design = ~type) #Se hace un DESeqDataSet para realizar un analisis

dim(dds) # checar las dimensiones

##  -- Asignar la referencia y generar contrastes -----
# Las comparaciones se realizan por pares
#Si no se indica de manera explicita que se va a comparara, lo va a tomar de manera alfabetica, 
# en este caso se indica que control es la referencia, 
dds$type <- relevel(dds$type, ref = "CONTROL") 

## --- Obtener archivo dds ----

dds <- DESeq(dds)

# Obtener la lista de coeficientes o contrastes
resultsNames(dds)

# Guardar la salida del diseno
save(metadata, dds, file = paste0(outdir, 'dds_treatment_vs_control.RData'))

## --- Normalizacion de los datos ---------
# Opcion 1. log2(n + 1)
ntd <- normTransform(dds)

# Opcion 2. regularized logarithm or rlog
# Normalizacion de las cuentas por logaritmo y podrias hacer el analisis usando este objeto en lugar del dds
ddslog <- rlog(dds, blind = F) 

# Opcion 3. vsd
# Estima la tendencia de dispersion de los datos y calcula la varianza, hace una normalizacion de las 
# cuentas con respecto al tamaño de la libreria
vsdata <- vst(dds, blind = F) 

## --- Deteccion de batch effect ----

# Almacenar la grafica
png(file = paste0(figdir, "PCA_rlog.png"))
plt <- plotPCA(ddslog, intgroup = "type")
print(plt)
dev.off()

# Almacenar la grafica
png(file = paste0(figdir, "PCA_vsd.png"))
plt <- plotPCA(vsdata, intgroup = "type")
print(plt)
dev.off()

# Guardar la salida del diseno (vsdata)
save(metadata, vsdata, file = paste0(outdir, 'vst_treatment_vs_control.RData'))

# En la grafica de las primeras dos componentes principales son notorias las diferencias 
# entre tipos de muestras con respecto a las componente principales que capturan su varianza, 
# cada componente principal representa una combinacion lineal de las variables (en este caso genes) 
# que explican la mayor cantidad de varianza en nuestros datos (las cuentas).


## ---- Obtener informacion del contraste 1 ----
# results(dds, contrast=c("condition","treated","untreated"))
res_t <- results(dds, name = "type_TREATMENT_vs_CONTROL")
res_t

summary(res_t)

# Guardar los resultados
write.csv(res_t, file=paste0(outdir, 'DE_treatment_vs_control.csv'))

```
For the visualization:
```{r, eval = FALSE}
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
