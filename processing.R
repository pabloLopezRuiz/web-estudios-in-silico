library(GEOquery)
library(biomaRt)
library(SummarizedExperiment)

# coldata
gse <- getGEO("GSE124799")
gse <- gse$GSE124799_series_matrix.txt.gz
coldata <- pData(gse)

# expdata
expdata <- experimentData(gse)

# counts
filename <- rownames(getGEOSuppFiles("GSE124799")[1])
mtx <- read.csv(filename, sep = "\t")

rowmtx <- mtx[, 1]
mtx <- mtx[, -1]
rownames(mtx) <- rowmtx

# Es lo mismo
colnames(mtx)
coldata$title

# Prefiero tener el GEO accession como nombre de las columnas
rownames(coldata)
colnames(mtx) <- coldata$geo_accession

# Busco anotaciones adicionales para los genes
ensembl <- useMart(biomart = "ensembl", "mmusculus_gene_ensembl")
mygenes <- rownames(mtx)
map_genes <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "start_position", "end_position", "chromosome_name"),
    filters = "ensembl_gene_id",
    values = mygenes,
    mart = ensembl
)
map_genes <- map_genes[!duplicated(map_genes$ensembl_gene_id), ]
rownames(map_genes) <- map_genes$ensembl_gene_id

# Filtro aquellos genes que no tienen
mtx_fil <- mtx[rownames(map_genes), ]

# map_genes[map_genes$ensembl_gene_id %in% rownames(mtx_fil), ]

g <- SummarizedExperiment(
    assays = list(counts = as.matrix(mtx_fil)),
    colData = coldata,
    rowData = map_genes,
    metadata = list(experimentData = expdata)
)


colData(g)$title <- make.names(colData(g)$title)

colData(g)$group <- sapply(strsplit(colData(g)$title, "\\."), "[", 1)
colData(g)$group <- as.factor(colData(g)$group)
colData(g)$group <- relevel(colData(g)$group, ref = "Sed")

# Cambio algunos nombres del summarized experiment para facilitar el análisis
colnames(colData(g)) <- make.names(colnames(colData(g)))
colnames(colData(g))[45] <- "gender"
colData(g)$gender <- as.factor(colData(g)$gender) 

save(g, file = "web_esb/GSE124799/gse.RData")
