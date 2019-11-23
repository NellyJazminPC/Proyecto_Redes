#############################################
#############################################
#### Minimum Spanning Network  ( M S N ) #### 
#### Análisis de Red de Expansión Mínima ####
#############################################

#Limpiar ambiente de trabajo:
rm(list = ls(all=TRUE))

# Establecer directorio de trabajo:
setwd("/media/nel_pc/n311_pc/Campos_Project_003_Analysis/data_metadata_5/CGUM_79/Figuras/")

#Librerías:
library(vcfR)
library(utils)
library(adegenet)
library(poppr)
library(devtools)
library(plotly)

# Importar matriz de SNPs (archivo VCF: Variant Call Format)
CGUM_79_1.VCF <- read.vcfR("/media/nel_pc/n311_pc/Campos_Project_003_Analysis/data_metadata_5/CGUM_79/var.79.2.1.sorted.vcf")
CGUM_79_1.VCF

# Importar base de datos global de los individuos muestreados:
pop.data <- read.csv("/media/nel_pc/n311_pc/Campos_Project_003_Analysis/data_metadata_5/Qmacdougalli_79ind_stacks.csv")
head (pop.data)
#Verificar que ambos archivos tengan la misma cantidad de individuos:
all(colnames(CGUM_79_1.VCF@gt)[-1] == pop.data$NOM_IND)

#Convertir a genlight:
cgum_79_1_vcf <- vcfR2genlight(CGUM_79_1.VCF)

#Establecer ploidia:
ploidy(cgum_79_1_vcf) <- 2
#Verificar la ploidia:
cgum_79_1_vcf@ploidy

#Establecer los sitios de muestreo como poblaciones artificiales:
pop(cgum_79_1_vcf) <- pop.data$SITIO

# Establecer un vector con los colores para cada sitio de muestreo:

cols <- c("#7570B3", "#075277","#00B1E8","#1FC944",
          "#E6AB02", "#E7298A","#E07E34", "#F15858")

##########################################################  
#Analisis de Red de Expansión Mínima:

#Calcular las distancias genéticas
cgum_79_1.dist <- bitwise.dist(cgum_79_1_vcf)

#MSN
cgum_79_1.msn <- poppr.msn(cgum_79_1_vcf, cgum_79_1.dist, showplot = FALSE, include.ties = T)

#Definir los nombres y tamaños de los nodos
node.size <- rep(2, times = nInd(cgum_79_1_vcf))
names(node.size) <- indNames(cgum_79_1_vcf)
vertex.attributes(cgum_79_1.msn$graph)$size <- node.size

#Graficar MSN
set.seed(12345)
plot_poppr_msn(cgum_79_1_vcf, cgum_79_1.msn, palette = cols, gadj = 500)

#Modo interactivo para el análisis MSN
imsn()

