######
# Script : Importar datos de cuentas en R
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 27/02/2024
# Description: El siguiente script nos permite importar los datos provenientes del alineamiento de STAR a R,
# para el posterior analisis de Expresion diferencial con DESEq2.
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: metadata.csv, cuentas de STAR (Terminacion ReadsPerGene.out.tab)
#   - Output: Matriz de cuentas (CSV y RData)
#######

# --- Load data -----
# Cargar archivos
#indir <- "/mnt/Guanina/bioinfo24/data/Clase_RNASeq2024/STAR_output"
indir <- "/mnt/Guanina/bioinfo24/Equipo3/Proyecto/STAR_output/"
outdir <- "/mnt/Guanina/bioinfo24/Equipo3/Proyecto/results/"

# Opcion B -  sin movernos de carpeta
files <- dir(indir, pattern = "ReadsPerGene.out.tab")

# crear matriz de cuentas
counts <- c() # esta sera la matriz
for(i in seq_along(files)){
  x <- read.table(file = files[i], sep = "\t", header = F, as.is = T)
  # as.is para no convertir tipo de datos
  counts <- cbind(counts, x[,2])
}

# Cargar Metadatos 
metadata <- read.csv("/mnt/Guanina/bioinfo24/Equipo3/Proyecto/metadata.csv", header = F)
# Renombrar columnas en la metadata
colnames(metadata) <- c("sample_id", "type")
# Convertir a formato dataframe
counts <- as.data.frame(counts)
rownames(counts) <- x[,1] # Renombrar las filas con el nombre de los genes
colnames(counts) <- sub("_ReadsPerGene.out.tab", "", files)

# Eliminar las 4 primeras filas
counts <- counts[-c(1:4),]

# Almacenar metadata y matriz de cuentas
save(metadata, counts, file = paste0(outdir, "counts/raw_counts.RData"))
write.csv(counts, file = paste0(outdir,"counts/raw_counts.csv"))

# Guardar informacion de ejecucion
sessionInfo()
```
*The R scripts were provided in the RNA-seq course by Evelia Coss, we just modified the pathways and file names 

Then, it is necessary normalize the data to avoid technical and biological biases and visualize it:   
```{r, eval = FALSE}
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
