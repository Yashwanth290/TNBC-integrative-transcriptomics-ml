############################################
# 🧬 1️⃣ SETUP
############################################
dir.create("Cancer_Project", showWarnings = FALSE)
setwd("Cancer_Project")

dirs <- c("data","figures","results",
          "results/DEG","results/Enrichment",
          "results/GSEA","results/ML")

for(d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

############################################
# 🧬 2️⃣ LIBRARIES
############################################
packages <- c(
  "GEOquery","limma","ggplot2","pheatmap","dplyr",
  "clusterProfiler","org.Hs.eg.db","enrichplot",
  "pROC","randomForest","glmnet"
)

for(p in packages){
  if(!require(p, character.only=TRUE)){
    install.packages(p, dependencies=TRUE)
    library(p, character.only=TRUE)
  }
}

############################################
# 🧬 3️⃣ DOWNLOAD GEO DATA
############################################
gset <- getGEO("GSE38959", GSEMatrix=TRUE)
exprSet <- exprs(gset[[1]])
pheno <- pData(gset[[1]])
feature <- fData(gset[[1]])

############################################
# 🧬 4️⃣ SAMPLE FILTER
############################################
keep <- pheno$source_name_ch1 %in% c(
  "triple negative breast cancer cells",
  "normal mammary gland ductal cells"
)

expr <- exprSet[,keep]
pheno <- pheno[keep,]

group <- factor(pheno$source_name_ch1)
levels(group) <- c("Normal","TNBC")

colnames(expr) <- paste0(group,"_",1:length(group))

############################################
# 🧬 5️⃣ ANNOTATION (FINAL FIX — NO FAIL)
############################################
match_idx <- match(rownames(expr), feature$ID)
valid <- !is.na(match_idx)

expr <- expr[valid, ]
feature <- feature[match_idx[valid], ]

genes <- sapply(feature$GENE_SYMBOL, function(x){
  strsplit(as.character(x)," /// ")[[1]][1]
})

# 🔥 REMOVE BAD GENES (INCLUDING NA)
valid_gene <- !is.na(genes) & genes != "" & genes != "---"

expr <- expr[valid_gene, ]
genes <- genes[valid_gene]

# REMOVE DUPLICATES
dup <- duplicated(genes)
expr <- expr[!dup, ]
genes <- genes[!dup]

# ASSIGN GENE NAMES
rownames(expr) <- genes

# 🔥 FINAL CRITICAL FIX (prevents DEG error)
bad <- is.na(rownames(expr)) | rownames(expr) == ""
expr <- expr[!bad, ]
rownames(expr) <- make.names(rownames(expr), unique=TRUE)

cat("✅ FINAL GENES:", nrow(expr), "\n")

############################################
# 🧬 6️⃣ DEG ANALYSIS (WORKS 100%)
############################################
design <- model.matrix(~0+group)
colnames(design) <- c("Normal","TNBC")

fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, makeContrasts(TNBC-Normal, levels=design))
fit <- eBayes(fit)

deg <- topTable(fit, number=Inf)

write.csv(deg, "results/DEG/DEG_full.csv")

############################################
# 🧬 7️⃣ VOLCANO (UP/DOWN/NS)
############################################
deg$category <- "NS"
deg$category[deg$logFC > 1 & deg$adj.P.Val < 0.05] <- "UP"
deg$category[deg$logFC < -1 & deg$adj.P.Val < 0.05] <- "DOWN"

write.csv(deg, "results/DEG/DEG_annotated.csv")

png("figures/Volcano.png",1200,1000,res=150)
ggplot(deg, aes(logFC, -log10(adj.P.Val), color=category)) +
  geom_point(size=2) +
  scale_color_manual(values=c("blue","grey","red")) +
  theme_minimal() +
  ggtitle("Volcano Plot")
dev.off()

############################################
# 🧬 8️⃣ PCA
############################################
pca <- prcomp(t(expr))
pca_df <- data.frame(pca$x[,1:2], Group=group)

png("figures/PCA.png",1200,1000,res=150)
ggplot(pca_df, aes(PC1, PC2, color=Group)) +
  geom_point(size=4) +
  theme_minimal()
dev.off()

############################################
# 🧬 9️⃣ HEATMAP (TOP UP GENES)
############################################
top_up <- deg %>%
  filter(logFC > 1, adj.P.Val < 0.05) %>%
  arrange(desc(logFC)) %>%
  head(50)

heat <- expr[rownames(top_up), ]
heat <- t(scale(t(heat)))

write.csv(heat, "results/DEG/top50_expression.csv")

png("figures/Heatmap.png",1200,1000,res=150)
pheatmap(heat, show_colnames=TRUE, show_rownames=TRUE)
dev.off()

############################################
# 🧬 🔟 GO ENRICHMENT
############################################
sig <- deg %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)

gene_df <- bitr(rownames(sig),
                fromType="SYMBOL",
                toType="ENTREZID",
                OrgDb=org.Hs.eg.db)

ego <- enrichGO(gene = gene_df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP")

write.csv(as.data.frame(ego), "results/Enrichment/GO.csv")

############################################
# 🧬 1️⃣1️⃣ GO PLOTS
############################################
png("figures/GO_dot.png",1200,1000,res=150)
dotplot(ego, showCategory=20)
dev.off()

png("figures/GO_bar.png",1200,1000,res=150)
barplot(ego, showCategory=20)
dev.off()

############################################
# 🧬 1️⃣2️⃣ KEGG
############################################
ekegg <- enrichKEGG(gene = gene_df$ENTREZID, organism="hsa")

write.csv(as.data.frame(ekegg), "results/Enrichment/KEGG.csv")

png("figures/KEGG_dot.png",1200,1000,res=150)
dotplot(ekegg)
dev.off()
png("figures/KEGG_bar.png",1200,1000,res=150)
barplot(ekegg)
dev.off()

#############################################
# 🧬 1️⃣3️⃣ GSEA (FINAL FIXED — WORKING)
############################################

library(clusterProfiler)
library(org.Hs.eg.db)

# STEP 1: Prepare gene list
gene_df <- bitr(
  rownames(deg),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# 🔥 IMPORTANT: remove NA
gene_df <- gene_df[!is.na(gene_df$ENTREZID), ]

# STEP 2: match logFC
geneList <- deg$logFC
names(geneList) <- rownames(deg)

# keep only mapped genes
geneList <- geneList[gene_df$SYMBOL]

# assign ENTREZ IDs
names(geneList) <- gene_df$ENTREZID

# sort (VERY IMPORTANT)
geneList <- sort(geneList, decreasing = TRUE)

# 🚨 DEBUG CHECK
cat("Genes for GSEA:", length(geneList), "\n")
if(length(geneList) < 10){
  stop("❌ ERROR: Too few genes for GSEA")
}

# STEP 3: RUN GSEA
gsea <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# STEP 4: SAVE
dir.create("results/GSEA", showWarnings = FALSE)
write.csv(as.data.frame(gsea), "results/GSEA/GSEA_results.csv")

# STEP 5: PLOTS
dir.create("figures", showWarnings = FALSE)

png("figures/GSEA_dotplot.png",1200,1000)
dotplot(gsea, showCategory=20)
dev.off()

png("figures/GSEA_ridgeplot.png",1200,1000)
ridgeplot(gsea)
dev.off()

cat("✅ GSEA COMPLETED SUCCESSFULLY\n")

############################################
# 🧬 1️⃣4️⃣ ROC
############################################
top_genes <- rownames(top_up)[1:10]

roc_res <- data.frame()

for(g in top_genes){
  r <- roc(group, as.numeric(expr[g,]))
  roc_res <- rbind(roc_res,
                   data.frame(Gene=g, AUC=auc(r)))
}

write.csv(roc_res, "results/ML/ROC.csv")

############################################
# 🧬 1️⃣5️⃣ RANDOM FOREST
############################################
ml <- data.frame(t(expr[top_genes,]))
ml$Group <- group

rf <- randomForest(Group~., data=ml)

png("figures/RF.png",1200,1000,res=150)
varImpPlot(rf)
dev.off()

############################################
# 🧬 1️⃣6️⃣ LASSO
############################################
x <- as.matrix(t(expr[top_genes,]))
y <- ifelse(group=="TNBC",1,0)

lasso <- cv.glmnet(x,y,family="binomial")

png("figures/LASSO.png",1200,1000,res=150)
plot(lasso)
dev.off()

############################################
# ✅ DONE
############################################
cat("🎉 PIPELINE COMPLETED SUCCESSFULLY\n")


############################################
# 🧬 GENERATE ROC RESULTS (FIXED)
############################################

library(pROC)

# Ensure group
group <- factor(group, levels = c("Normal","TNBC"))

roc_res <- data.frame()

# Use top DEGs (not hub_genes — safer)
top_genes <- rownames(deg)[1:50]

for(g in top_genes){
  
  if(g %in% rownames(expr)){
    
    gene_exp <- as.numeric(expr[g, ])
    
    roc_obj <- roc(
      response = group,
      predictor = gene_exp,
      levels = c("Normal","TNBC"),
      direction = ">"
    )
    
    roc_res <- rbind(
      roc_res,
      data.frame(
        Gene = g,
        AUC = as.numeric(auc(roc_obj))
      )
    )
  }
}

# Sort by AUC
roc_res <- roc_res[order(-roc_res$AUC), ]

# SAVE (THIS WAS MISSING)
dir.create("results/ROC", recursive = TRUE, showWarnings = FALSE)
write.csv(roc_res, "results/ROC/ROC_results.csv", row.names = FALSE)

cat("✅ ROC results generated\n")


############################################
# 🧬 BEST ROC (TOP GENE)
############################################

library(pROC)
library(ggplot2)

# Ensure group is correct
group <- factor(group, levels = c("Normal", "TNBC"))

# Load ROC results
roc_res <- read.csv("results/ROC/ROC_results.csv")

# Select best gene
best_gene <- roc_res$Gene[1]

# ROC calculation
roc_obj <- roc(group, as.numeric(expr[best_gene, ]), direction = ">")

# Create dataframe for ggplot
roc_df <- data.frame(
  tpr = roc_obj$sensitivities,
  fpr = 1 - roc_obj$specificities
)

auc_val <- round(auc(roc_obj), 3)

# Plot (Publication quality)
p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(size = 1.5, color = "#D62728") +
  geom_abline(linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 16) +
  labs(
    title = paste("ROC Curve -", best_gene),
    subtitle = paste("AUC =", auc_val),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

# Save
ggsave("figures/Best_ROC.png", plot = p, width = 8, height = 6, dpi = 300)

cat("✅ Best ROC saved\n")


############################################
# 🧬 OVERALL ROC (TOP 5 GENES)
############################################

library(pROC)
library(ggplot2)

# Select top 5 genes
top_genes <- head(roc_res$Gene, 5)

roc_list <- list()

for(g in top_genes){
  
  roc_obj <- roc(group, as.numeric(expr[g, ]), direction=">")
  
  roc_df <- data.frame(
    tpr = roc_obj$sensitivities,
    fpr = 1 - roc_obj$specificities,
    gene = g,
    auc = round(auc(roc_obj), 3)
  )
  
  roc_list[[g]] <- roc_df
}

# Combine
roc_all <- do.call(rbind, roc_list)

# Plot
p_all <- ggplot(roc_all, aes(x = fpr, y = tpr, color = gene)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed", color = "black") +
  theme_classic(base_size = 16) +
  labs(
    title = "ROC Curves (Top 5 Biomarker Genes)",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Gene"
  )

# Save
ggsave("figures/Overall_ROC.png", plot = p_all, width = 8, height = 6, dpi = 300)

cat("✅ Overall ROC saved\n")


############################################
# 🧬 OVERALL ROC WITH AUC LABELS
############################################

roc_labels <- unique(roc_all[,c("gene","auc")])
roc_labels$label <- paste0(roc_labels$gene, " (AUC=", roc_labels$auc, ")")

roc_all$gene_label <- roc_labels$label[match(roc_all$gene, roc_labels$gene)]

p_all <- ggplot(roc_all, aes(x = fpr, y = tpr, color = gene_label)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed") +
  theme_classic(base_size = 16) +
  labs(
    title = "ROC Curves with AUC",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Genes"
  )

ggsave("figures/Overall_ROC_AUC.png", p_all, width=8, height=6, dpi=300)


############################################
# 🧬 PPI NETWORK (STRINGdb)
############################################

library(STRINGdb)
library(igraph)

# Initialize STRING
string_db <- STRINGdb$new(
  version="11.5",
  species=9606,
  score_threshold=400
)

# Use significant genes
sig_genes <- rownames(deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ])

# Map to STRING IDs
mapped <- string_db$map(
  data.frame(gene=sig_genes),
  "gene",
  removeUnmappedRows = TRUE
)

# Get interactions
ppi <- string_db$get_interactions(mapped$STRING_id)

# Convert to graph
graph <- graph_from_data_frame(
  ppi[, c("from","to")],
  directed = FALSE
)
dir.create("results/PPI", recursive = TRUE, showWarnings = FALSE)
# Save edge list
write.csv(ppi, "results/PPI/PPI_edges.csv", row.names=FALSE)

cat("✅ PPI network created\n")

dir.create("results/Drug", recursive = TRUE, showWarnings = FALSE)
dir.create("results/ROC", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

############################################
# 🧬 HUB GENES
############################################

# Calculate degree
deg_vals <- degree(graph)

# Top 20 hub genes (STRING IDs)
hub_ids <- names(sort(deg_vals, decreasing=TRUE))[1:20]

# Map back to gene symbols
hub_genes <- mapped$gene[match(hub_ids, mapped$STRING_id)]

# Remove NA
hub_genes <- na.omit(hub_genes)

# Save
write.csv(data.frame(Hub_Genes=hub_genes),
          "results/PPI/Hub_genes.csv",
          row.names=FALSE)

cat("🔥 Hub genes:", paste(hub_genes, collapse=", "), "\n")

############################################
# 🧬 PPI PLOT
############################################

V(graph)$color <- ifelse(
  names(V(graph)) %in% hub_ids,
  "red",
  "skyblue"
)

png("figures/PPI_network.png", width=1200, height=1000, res=150)

plot(
  graph,
  vertex.size=5,
  vertex.label=NA,
  edge.color="grey70"
)

legend(
  "topright",
  legend=c("Hub genes","Other genes"),
  col=c("red","skyblue"),
  pch=16
)

dev.off()

cat("✅ PPI plot saved\n")


############################################
# 🧬 HUB GENE HEATMAP
############################################

library(pheatmap)

hub_genes <- hub_genes[hub_genes %in% rownames(expr)]

heat_data <- expr[hub_genes, ]
heat_data <- t(scale(t(heat_data)))

png("figures/Hub_genes_heatmap.png",1200,1000,res=150)

pheatmap(
  heat_data,
  show_colnames = TRUE,
  show_rownames = TRUE,
  main = "Hub Genes Expression"
)

dev.off()


############################################
# 💊 DRUG–GENE INTERACTION
############################################

library(httr)
library(jsonlite)

dir.create("results/Drug", showWarnings = FALSE)

# Use top 5 hub genes (safe)
genes_use <- hub_genes[1:min(5,length(hub_genes))]

url <- paste0(
  "https://dgidb.org/api/v2/interactions.json?genes=",
  paste(genes_use, collapse=",")
)

res <- GET(url)
txt <- content(res, "text", encoding="UTF-8")

parsed <- tryCatch(fromJSON(txt), error=function(e) NULL)

if(!is.null(parsed) && length(parsed$matchedTerms) > 0){
  
  drug_df <- parsed$matchedTerms
  
} else {
  
  drug_df <- data.frame(
    geneName = genes_use,
    drugName = "Check manually",
    interactionTypes = "NA"
  )
}

# Save
write.csv(drug_df, "results/Drug/Drug_gene.csv", row.names=FALSE)

cat("✅ Drug-gene saved\n")


############################################
# 💊 DRUG–GENE PLOT
############################################

png("figures/Drug_gene_barplot.png",1200,1000)

barplot(
  table(drug_df$geneName),
  las=2,
  col="steelblue",
  main="Drug–Gene Interaction Count"
)

dev.off()

cat("✅ Drug plot saved\n")


if(!require(rms)) install.packages("rms")
if(!require(survival)) install.packages("survival")

library(rms)
library(survival)

if(!require(iml)) install.packages("iml")
library(iml)

# Use your ML dataset
ml_data <- data.frame(t(expr[top_genes, ]))
ml_data$Group <- group

# Train model
library(randomForest)
rf_model <- randomForest(Group ~ ., data = ml_data)

predictor <- Predictor$new(
  rf_model,
  data = ml_data[ , -ncol(ml_data)],
  y = ml_data$Group
)


shap <- Shapley$new(
  predictor,
  x.interest = ml_data[1, -ncol(ml_data)]
)


png("figures/SHAP_single.png", width=1800, height=1400, res=300)

plot(shap)

dev.off()

cat("✅ SHAP single sample saved\n")


imp <- FeatureImp$new(
  predictor,
  loss = "ce"
)

png("figures/SHAP_importance.png", width=1800, height=1400, res=300)

plot(imp)

dev.off()

cat("✅ SHAP importance saved\n")



