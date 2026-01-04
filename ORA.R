library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
library(STRINGdb)
library(biomaRt)
library(tidyverse)
library(visNetwork)

# Cargo los resultados del DEA
load("web_esb/GSE124799/DE_Sed_vs_Ex.RData")
load("web_esb/GSE124799/DE_Sed_vs_ExSed.RData")
load("web_esb/GSE124799/DE_Ex_vs_ExSed.RData")

# Establezco los umbrales
pval = 0.05
lfc = 3

# Filtro los genes por p valor y logFoldChange en ambos sentidos
# Lo hago todo por triplicado para ver los tres contrastes
top_sve <- SvE[SvE$padj < pval & SvE$log2FoldChange > lfc, ]
top_genes_sve <- rownames(top_sve)
top_eves <- EvEs[EvEs$padj < pval & EvEs$log2FoldChange > lfc, ]
top_genes_eves <- rownames(top_eves)
top_sves <- SvEs[SvEs$padj < pval & SvEs$log2FoldChange > lfc, ]
top_genes_sves <- rownames(top_sves)

bot_sve <- SvE[SvE$padj < pval & SvE$log2FoldChange < -lfc, ]
bot_genes_sve <- rownames(bot_sve)
bot_eves <- EvEs[EvEs$padj < pval & EvEs$log2FoldChange < -lfc, ]
bot_genes_eves <- rownames(bot_eves)
bot_sves <- SvEs[SvEs$padj < pval & SvEs$log2FoldChange < -lfc, ]
bot_genes_sves <- rownames(bot_sves)

# ORA
# No filtro las rutas (pval = 1), eso lo haré después
oraTop_sve <- enrichGO(gene        = top_genes_sve,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 500)
save(oraTop_sve, file = "data/oraTop_SvE.RData")

oraTop_eves <- enrichGO(gene        = top_genes_eves,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 500)
save(oraTop_eves, file = "data/oraTop_EvEs.RData")

oraTop_sves <- enrichGO(gene        = top_genes_sves,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 500)
save(oraTop_sves, file = "web_esb/data/oraTop_SvEs.RData")

oraBot_sve <- enrichGO(gene        = bot_genes_sve,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 500)
save(oraBot_sve, file = "web_esb/data/oraBot_SvE.RData")

oraBot_eves <- enrichGO(gene        = bot_genes_eves,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 500)
save(oraBot_eves, file = "web_esb/data/oraBot_EvEs.RData")

oraBot_sves <- enrichGO(gene        = bot_genes_sves,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 500)
save(oraBot_sves, file = "web_esb/data/oraBot_SvEs.RData")


# Barplots
barplot(oraTop_sve, showCategory = 5)
barplot(oraTop_eves, showCategory = 5)
barplot(oraTop_sves, showCategory = 5)

# Para el GSEA ordeno los genes significativos por logFoldChange
sve_ord <- SvE %>%
  as.data.frame() %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(log2FoldChange) %>%
  rownames_to_column() %>%
  deframe()

eves_ord <- EvEs %>%
  as.data.frame() %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(log2FoldChange) %>%
  rownames_to_column() %>%
  deframe()

sves_ord <- SvEs %>%
  as.data.frame() %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(log2FoldChange) %>%
  rownames_to_column() %>%
  deframe()


# GSEA, vuelvo a no filtrar, lo haré luego
gsea_go_sve <- gseGO(geneList  = sve_ord,
                 ont       = "BP",
                 OrgDb     = org.Mm.eg.db, 
                 keyType   = 'ENSEMBL',
                 verbose   = T,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1
)
save(gsea_go_sve, file = "web_esb/data/GSEA_GO_SvE.RData")

gsea_go_eves <- gseGO(geneList  = eves_ord,
                 ont       = "BP",
                 OrgDb     = org.Mm.eg.db, 
                 keyType   = 'ENSEMBL',
                 verbose   = T,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1
)
save(gsea_go_eves, file = "web_esb/data/GSEA_GO_EvEs.RData")

gsea_go_sves <- gseGO(geneList  = sves_ord,
                 ont       = "BP",
                 OrgDb     = org.Mm.eg.db, 
                 keyType   = 'ENSEMBL',
                 verbose   = T,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1
)
save(gsea_go_sves, file = "web_esb/data/GSEA_GO_SvEs.RData")

# Enrichplot
enrichplot::dotplot(gsea_go_sve, showCategory = 3) 
enrichplot::dotplot(gsea_go_eves, showCategory = 3) 
enrichplot::dotplot(gsea_go_sves, showCategory = 3) 

# GSEAplot, este gráfico finalmente no lo he usado
gseaplot2(gsea_go_sves, geneSetID = 2, title = gsea_go_sves$Description[1])

# Repito con KEGG
# hago un biomart para conseguir los entrez
all_genes <- unique(c(names(sve_ord), names(sves_ord), names(eves_ord)))

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
conversion <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters    = "ensembl_gene_id",
  values     = all_genes,
  mart       = mart
)

# Uno las notaciones
sve_kegg <- sve_ord %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(conversion, by= "ensembl_gene_id") %>%
  drop_na(entrezgene_id) %>%
  distinct(entrezgene_id, .keep_all = T) %>%
  dplyr::select("entrezgene_id", ".") %>%
  deframe()

eves_kegg <- eves_ord %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(conversion, by= "ensembl_gene_id") %>%
  drop_na(entrezgene_id) %>%
  distinct(entrezgene_id, .keep_all = T) %>%
  dplyr::select("entrezgene_id", ".") %>%
  deframe()

sves_kegg <- sves_ord %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(conversion, by= "ensembl_gene_id") %>%
  drop_na(entrezgene_id) %>%
  distinct(entrezgene_id, .keep_all = T) %>%
  dplyr::select("entrezgene_id", ".") %>%
  deframe()

gsea_kegg_sve = gseKEGG(geneList     = sve_kegg,
                    organism     = 'mmu',
                    keyType      = "ncbi-geneid",
                    pvalueCutoff = 0.5,
                    verbose      = TRUE)
save(gsea_kegg_sve, file = "web_esb/data/GSEA_KEGG_SvE.RData")

gsea_kegg_eves = gseKEGG(geneList     = eves_kegg,
                    organism     = 'mmu',
                    keyType      = "ncbi-geneid",
                    pvalueCutoff = 0.5,
                    verbose      = TRUE)
save(gsea_kegg_eves, file = "web_esb/data/GSEA_KEGG_EvEs.RData")

gsea_kegg_sves = gseKEGG(geneList     = sves_kegg,
                    organism     = 'mmu',
                    keyType      = "ncbi-geneid",
                    pvalueCutoff = 0.5,
                    verbose      = TRUE)
save(gsea_kegg_sves, file = "web_esb/data/GSEA_KEGG_SvEs.RData")

enrichplot::dotplot(gsea_kegg_sve, showCategory = 3)
enrichplot::dotplot(gsea_kegg_eves, showCategory = 3)
enrichplot::dotplot(gsea_kegg_sves, showCategory = 3)

# STRING
string_db <- STRINGdb$new(
  version = "12",
  species = 10090,
  score_threshold = 400
)

top_sve <- SvE %>% 
  as.data.frame() %>%
  dplyr::select("log2FoldChange", "padj") %>%
  rownames_to_column()
colnames(top_sve) <- c("gene", "logFC", "P.Value")

top_eves <- EvEs %>% 
  as.data.frame() %>%
  dplyr::select("log2FoldChange", "padj") %>%
  rownames_to_column()
colnames(top_eves) <- c("gene", "logFC", "P.Value")

top_sves <- SvEs %>% 
  as.data.frame() %>%
  dplyr::select("log2FoldChange", "padj") %>%
  rownames_to_column()
colnames(top_sves) <- c("gene", "logFC", "P.Value")


mapped_sve <- string_db$map(top_sve, "gene", removeUnmappedRows = TRUE)
mapped_eves <- string_db$map(top_eves, "gene", removeUnmappedRows = TRUE)
mapped_sves <- string_db$map(top_sves, "gene", removeUnmappedRows = TRUE)

save(mapped_sve, file="web_esb/data/string_SvE.RData")
save(mapped_eves, file="web_esb/data/string_EvEs.RData")
save(mapped_sves, file="web_esb/data/string_SvEs.RData")

# Visualizo la Red
hits_sve <- mapped_sve$STRING_id[1:50]
hits_eves <- mapped_eves$STRING_id[1:50]
hits_sves <- mapped_sves$STRING_id[1:50]

edges_sve <- string_db$get_interactions(hits_sve)
edges_eves <- string_db$get_interactions(hits_eves)
edges_sves <- string_db$get_interactions(hits_sves)

nodes_sve <- data.frame(
  id = unique(c(edges_sve$from, edges_sve$to)),
  label = unique(c(edges_sve$from, edges_sve$to))
)
nodes_eves <- data.frame(
  id = unique(c(edges_eves$from, edges_eves$to)),
  label = unique(c(edges_eves$from, edges_eves$to))
)
nodes_sves <- data.frame(
  id = unique(c(edges_sves$from, edges_sves$to)),
  label = unique(c(edges_sves$from, edges_sves$to))
)

nwk_sve <- visNetwork(nodes_sve, edges_sve) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = FALSE)
nwk_eves <- visNetwork(nodes_eves, edges_eves) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = FALSE)
nwk_sves <- visNetwork(nodes_sves, edges_sves) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = FALSE)



# Para facilitarme el trabajo, aglomero los resultados por contraste

sve_fun <- list(
  ORA_Top  = oraTop_sve,
  ORA_Bot  = oraBot_sve,
  GSEA_GO  = gsea_go_sve,
  GSEA_KEGG = gsea_kegg_sve,
  Network    = nwk_sve
)

eves_fun <- list(
  ORA_Top   = oraTop_eves,
  ORA_Bot   = oraBot_eves,
  GSEA_GO   = gsea_go_eves,
  GSEA_KEGG = gsea_kegg_eves,
  Network    = nwk_eves
)

sves_fun <- list(
  ORA_Top   = oraTop_sves,
  ORA_Bot   = oraBot_sves,
  GSEA_GO   = gsea_go_sves,
  GSEA_KEGG = gsea_kegg_sves,
  Network    = nwk_sves
)

save(sve_fun, file="web_esb/data/sve_fun.RData")
save(eves_fun, file="web_esb/data/eves_fun.RData")
save(sves_fun, file="web_esb/data/sves_fun.RData")
