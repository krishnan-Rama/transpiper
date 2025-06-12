#!/usr/bin/env Rscript

# =========================
# GO Enrichment Full Pipeline
# Usage:
# Rscript go_enrichment_full.R <combined_annotation_csv> <edger_results_file> <output_dir>
# =========================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript go_enrichment_full.R <combined_annotation_csv> <edger_results_file> <output_dir>")
}

annotation_file <- args[1]
edgeR_file <- args[2]
outdir <- args[3]

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# --- Load packages ---
suppressMessages({
  library(data.table)
  library(clusterProfiler)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(GO.db)
  library(AnnotationDbi)
})

# === STEP 1: Create GO_annotations_expanded_full.csv ===
cat("📥 Reading annotation file...\n")
annot <- fread(annotation_file)

cat("🧬 Extracting GO terms...\n")
go_long <- annot[, .(GeneID, GO = unlist(strsplit(`Gene.Ontology..GO.`, ";\\s*"))), by=1:nrow(annot)]
go_long[, GO := sub(".*\\[(GO:\\d+)\\].*", "\\1", GO)]
go_long <- go_long[GO != ""]
go_df <- go_long[, .(GeneID, GO_Term = GO)]

go_annot_path <- file.path(outdir, "GO_annotations_expanded_full.csv")
fwrite(go_df, go_annot_path)
cat("✅ GO annotations saved to:", go_annot_path, "\n")

# === STEP 2: Filter DEGs from edgeR results ===
cat("📊 Reading edgeR results...\n")
deg <- fread(edgeR_file)
deg_filtered <- deg[abs(logFC) > 1 & FDR < 0.05]
deg_ids <- deg_filtered[[1]]  # First column is assumed to be GeneID

deg_path <- file.path(outdir, "OT8fold.csv")
writeLines(deg_ids, deg_path)
cat("✅ DE gene list saved to:", deg_path, "\n")

# === STEP 3: GO Enrichment ===
gene2go_df <- go_df
gene_universe <- unique(gene2go_df$GeneID)
deg_ids <- intersect(deg_ids, gene_universe)

TERM2GENE <- gene2go_df[, .(GO_Term, GeneID)]

# Get GO term names from GO.db
go_terms <- unique(gene2go_df$GO_Term)
go_names <- AnnotationDbi::select(GO.db, keys = go_terms, columns = "TERM", keytype = "GOID")
TERM2NAME <- data.frame(go_term = go_names$GOID, name = go_names$TERM)

cat("🚀 Running GO enrichment...\n")
ego_result <- enricher(
  gene = deg_ids,
  universe = gene_universe,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)

# === STEP 4: Save results ===
result_file <- file.path(outdir, "go_enrichment_results.csv")
fwrite(as.data.frame(ego_result), result_file)
cat("✅ GO enrichment results saved to:", result_file, "\n")

# === STEP 5: Save plots ===
png(file.path(outdir, "dotplot.png"), width = 1200, height = 900)
print(dotplot(ego_result, showCategory = 20) + ggtitle("GO Term Enrichment"))
dev.off()

png(file.path(outdir, "barplot.png"), width = 1200, height = 900)
print(barplot(ego_result, showCategory = 20) + ggtitle("GO Term Enrichment"))
dev.off()

png(file.path(outdir, "emapplot.png"), width = 1200, height = 900)
print(emapplot(pairwise_termsim(ego_result)))
dev.off()

cat("🎉 All done! Output written to:", outdir, "\n")
