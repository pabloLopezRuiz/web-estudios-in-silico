library(DESeq2)
library(ggplot2)
library(tidyverse)
library(PCAtools)
library(heatmaply)
library(mixOmics)
setwd(".")

# Cargo el summarized experiment
load("web_esb/GSE124799/gse.RData")

# Construyo el dataframe para hacer los barplots
coldata <- colData(g) %>%
  as.data.frame() %>%
  mutate(gender = factor(gender, levels = c("male", "female"))) %>%
  mutate(group = factor(group, levels = c("Sed", "Ex", "ExSed"))) %>%
  dplyr::select(gender, group) 

save(coldata, file = "web_esb/data/coldata.RData")

# Distribución de sexos
table(coldata$gender)
ggplot(coldata, aes(x = gender)) +
    scale_x_discrete(drop = FALSE) +
    geom_bar(width = 0.3, fill = "salmon") +
    theme_minimal() 

# Distribución de grupos experimentales
table(colData(g)$group)
ggplot(coldata, aes(x = group, fill = group)) +
    geom_bar(width = 0.3) +
    theme_minimal()

# Distribución de conteos
limma::plotDensities(
  object = log2(assay(g) + 1), 
  main = "Densidad de conteos (log2 + 1)",
  legend = F,
)
legend("topright", legend = colnames(g), col = 1:ncol(g), lty = 1)

# Mapeo de los grupos a los individuos
meta_df <- coldata %>%
  dplyr::select(group) %>%
  rownames_to_column("ind")

# Boxplot previo a la normalización con el logaritmo de los conteos
counts_long <- as.data.frame(assay(g)) %>% 
  stack() %>%
  rename(counts=1, ind=2) %>%
  left_join(meta_df, by = "ind")

bxpl <- ggplot(counts_long, aes(x = ind, y = log2(counts + 1), fill = group)) +
    geom_boxplot() +
    ggtitle("Distribución de conteos (log2 + 1) por muestra") + 
    ylab("Counts") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

save(bxpl, file="web_esb/data/boxplot_prenorm.RData")


# Filtro por genes poco detectados
keep <- rowSums(edgeR::cpm(assay(g))> 1) >= 3
g_sel <- g[keep, ]

# Estabelzco el grupo sedentario como control
g_sel$group <- relevel(g_sel$group, ref = "Sed")

# DEA
ddsSE <- DESeqDataSet(g_sel, design = ~group)
dsNorm <- estimateSizeFactors(ddsSE, type="ratio")
dds <- DESeq(dsNorm)

# Matriz de diseño
dmtx <- as.data.frame(model.matrix(design(dds), colData(dds))) %>%
  rename(Sed=1, Ex=2, ExSed=3)
save(dmtx, file = "web_esb/data/dmtx.RData")

# Transformo la varianza para estudiar la correlación entre muestras
vst_dds <- vst(dds, blind = TRUE) 

cor_matrix <- cor(assay(vst_dds), method = "spearman")
save(cor_matrix, file = "web_esb/data/corr_matrix.RData")


# Correlation
heatmaply(
  cor_matrix,
  ColSideColors = coldata$group,
  colors = colorRampPalette(c("navy","white","red"))(50),
  plot_method = "plotly"
)

# Boxplot tras normalizar
norm_counts <- log2(counts(dds, normalized = TRUE) + 1) %>%
  stack() %>%
  as.data.frame() %>%
  rename(gene=1, ind=2, counts=3) %>%
  left_join(meta_df, by = "ind")

bxpl <- ggplot(norm_counts, aes(x = ind, y = counts, fill = group)) +
  geom_boxplot() +
  ggtitle("Distribución de conteos (log2 + 1) por muestra") + 
  ylab("Counts") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

save(bxpl, file="web_esb/data/boxplot_norm.RData")

# PCA
pca.res <- PCAtools::pca(mat=assay(vst_dds), metadata=colData(g),  scale=TRUE)

# Elimino los atributos que no voy a usar para aligerar el fichero
pca.res$sdev = NULL
pca.res$rotated = NULL
pca.res$xvars = NULL
pca.res$yvars = NULL
pca.res$metadata = NULL

save(pca.res, file="web_esb/data/pca_res.RData")

# plot varianza explicada
screeplot(pca.res)

# Para el plot del PCA, prefiero como se ve el generado por mixOmics
X <- mixOmics::pca(t(assay(vst_dds)), ncomp = 5)
save(X, file = "data/mixomics_pca.RData")

# Plot de los loadings
plotLoadings(X, comp=1, ndisplay = 50)

# PCA plot
plotIndiv(X, group = colData(g)$group, pch = "circle", legend = TRUE)



# Guardo los distintos contrastes
SvE <- results(dds, contrast = c("group", "Sed", "Ex"))
SvE <- na.omit(SvE)
save(SvE, file = "web_esb/GSE124799/DE_Sed_vs_Ex.RData")

SvEs <- results(dds, contrast = c("group", "Sed", "ExSed"))
SvEs <- na.omit(SvEs)
save(SvEs, file = "web_esb/GSE124799/DE_Sed_vs_ExSed.RData")

EvEs <- results(dds, contrast = c("group", "Ex", "ExSed"))
EvEs <- na.omit(EvEs)
save(EvEs, file = "web_esb/GSE124799/DE_Ex_vs_ExSed.RData")

