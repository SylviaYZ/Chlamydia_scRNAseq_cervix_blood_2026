# reactome pathways

for (folder in c("All_INFvsUNINF/")) dir.create(folder, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(ggrepel)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
packageVersion("Seurat")
# Load data
cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

cervix.integrated$major_celltypes[which(cervix.integrated$major_celltypes%in% c("CD14_Mono","CD16_Mono"))] <- "Monocyte"

Idents(cervix.integrated) <- "major_celltypes"
# Enrichment anlaysis based on Reactome database
celltypes <- c("CD4_T", "CD8_T", "MAIT", "Monocyte", "B", "NK")

for (ct in celltypes){
  All.DE <- read.csv(paste0("All_INFvsUNINF/",ct,"_INFvsUNINF_All_Wilcox.csv"))

  All.DE <-  All.DE  %>%
    filter(p_val_adj < 0.05)

  DE.Entrez.SYMBOL <- bitr(All.DE$X, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% mutate(Mapping = "SYMBOL", USE = 1*!duplicated(SYMBOL))

  print("DE genes don't have Entrez ID:")

  print(sort(All.DE$X[!All.DE$X %in% DE.Entrez.SYMBOL$SYMBOL]))

  enrichedP <- enrichPathway(gene = DE.Entrez.SYMBOL$ENTREZID,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                readable = TRUE)

  output.df <- as.data.frame(enrichedP@result) %>% dplyr::select(ID, Description,GeneRatio,BgRatio, p.adjust,geneID) %>% mutate(geneID_upFC = "", geneID_downFC="", geneID_alphabet = "")

  for (row_i in 1:nrow(output.df)){
    genes_i <- unlist(strsplit(output.df$geneID[row_i], split = "/"))
    gene_i_log2FC <- All.DE$avg_log2FC[which(All.DE$X %in% genes_i)]
    names(gene_i_log2FC) <- genes_i

    # Separate into up-regulated and down-regulated genes
    upregulated <- genes_i[gene_i_log2FC[genes_i] > 1]
    downregulated <- genes_i[gene_i_log2FC[genes_i] < 1]

    # Order each group by fold change
    upregulated <- upregulated[order(gene_i_log2FC[upregulated], decreasing = TRUE)]
    downregulated <- downregulated[order(gene_i_log2FC[downregulated])]

    output.df[row_i, "geneID_upFC"] <- paste(upregulated,collapse="/")
    output.df[row_i, "geneID_downFC"] <- paste(downregulated,collapse="/")
    output.df[row_i, "geneID_alphabet"] <- paste(sort(genes_i),collapse="/")
  }

   xlsx::write.xlsx(output.df %>% dplyr::select(-c(geneID)) ,
            file = paste0("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_DE_INFvsUNINF/Reactome_pathways/",ct,"_Reactome_All_DE.xlsx"))

  print(
    barplot(enrichedP, showCategory=20, x = "GeneRatio") + ggtitle(paste0("REACTOME for ",ct))
  )
  print(
    dotplot(enrichedP, showCategory=20, x = "GeneRatio") + ggtitle(paste0("REACTOME for ",ct))
  )

    print( cnetplot(enrichedP, foldChange=All.DE[na.omit(match(All.DE$X, enrichedP@gene2Symbol)), "avg_log2FC"]) + ggtitle(paste0("REACTOME for ",ct)))

  enrichedP <- pairwise_termsim(enrichedP)

   print( emapplot(enrichedP, showCategory = 20) + ggtitle(paste0("REACTOME for ",ct)))

   print( emapplot(enrichedP, showCategory = 10) + ggtitle(paste0("REACTOME for ",ct)))

}

