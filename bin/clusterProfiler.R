#!/usr/bin/env Rscript

entrez_ids_to_symbol <- function(enrichment_table, entrez_ids) {

    entrez_map <- setNames(entrez_ids$SYMBOL, entrez_ids$ENTREZID)

    for (i in 1:nrow(enrichment_table)) {
        ids <- enrichment_table[i, 'geneID']
        ids <- unlist(strsplit(ids, "/"))

        converted_ids <- c()

        converted_ids <- sapply(ids, function(id) entrez_map[[id]])
        enrichment_table[i, 'geneID'] <- paste(converted_ids, collapse="/")

    }

    return(enrichment_table)
}

symbols_list_to_mirna_list <- function(enrichment_table, miRNA_gene_map) {

    enrichment_table$miRNA_list <- c()

    for (i in 1:nrow(enrichment_table)) {
        genes <- enrichment_table[i, 'geneID']
        genes <- unlist(strsplit(genes, "/"))

        # Define a set for the miRNAs
        mirnas <- c()

        for (gene in genes) {
            mirna_list <- miRNA_gene_map$miRNA_list[miRNA_gene_map$Target_Gene_Symbol == gene]
            mirna_list <- unlist(strsplit(mirna_list, "/"))
            mirnas <- union(mirnas, mirna_list)
        }

        enrichment_table[i, 'miRNA_list'] <- paste(mirnas, collapse="/")

    }

    return (enrichment_table)
}

######################################################################################
# Get inputs
######################################################################################
library(argparse)
parser <- ArgumentParser()

parser$add_argument("--up_miRNA_targets", help="Path to upregulated miRNA Targets", type="character")
parser$add_argument("--down_miRNA_targets", help="Path to downregulated miRNA Targets", type="character")

args <- parser$parse_args()

up_targets <- args$up_miRNA_targets
down_targets <- args$down_miRNA_targets

######################################################################################
# Load libraries
######################################################################################

library(logger)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

######################################################################################
# Start Processing Data
######################################################################################

# Set up log file
log_threshold(INFO)
log_formatter(formatter_paste)
log_appender(appender_file("clusterProfiler.log"))

log_info("Starting clusterProfiler analysis...")

# Process UPregulated miRNA targets
log_info("Processing upregulated miRNA targets from: ", up_targets)

# Make a directory for upregulated miRNA target enrichment results
dir.create("all_up_regulated_miRNA_enrichment")

up_mirna <- read.table(up_targets, header=TRUE, sep='\t')

# Only when the dataframe has genes, perfrom enrichment analysis
if (nrow(up_mirna) != 0) {
    entrez_ids <- bitr(up_mirna$Target_Gene_Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    up_entrez_ids_no_na <- na.omit(entrez_ids)
    write.table(up_entrez_ids_no_na, "all_up_regulated_miRNA_enrichment/up_symbol_to_entrezid.tsv", sep='\t', quote=FALSE, row.names = FALSE)

    na_entrez <- entrez_ids[is.na(entrez_ids$ENTREZID), ]
    log_info("There are ", nrow(na_entrez), " genes with no Entrez ID in the upregulated miRNA targets.")
    write.table(na_entrez, "all_up_regulated_miRNA_enrichment/up_no_entrezid.tsv", sep='\t', quote=FALSE, row.names = FALSE)

    # GO BP
    up_ego <- enrichGO(gene=up_entrez_ids_no_na$SYMBOL,
                    OrgDb= org.Hs.eg.db,
                    keyType= "SYMBOL",
                    ont= "BP",  # Biological Process
                    pAdjustMethod= "BH",
                    pvalueCutoff= 0.05)

    BP_go_results_up <- as.data.frame(up_ego)

    # GO MF
    up_go_mf <- enrichGO(gene=up_entrez_ids_no_na$SYMBOL,
                    OrgDb= org.Hs.eg.db,
                    keyType= "SYMBOL",
                    ont= "MF",
                    pAdjustMethod= "BH",
                    pvalueCutoff= 0.05)
    MF_go_results_up <- as.data.frame(up_go_mf)

    # GO CC
    up_go_cc <- enrichGO(gene=up_entrez_ids_no_na$SYMBOL,
                    OrgDb= org.Hs.eg.db,
                    keyType= "SYMBOL",
                    ont= "CC", # Cellular Component
                    pAdjustMethod= "BH",
                    pvalueCutoff= 0.05)
    CC_go_results_up <- as.data.frame(up_go_cc)

    # KEGG requires entrez ids
    up_kegg <- enrichKEGG(gene=up_entrez_ids_no_na$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    KEGG_results_up <- as.data.frame(up_kegg)

    # Reactome
    up_reactome <- enrichPathway(gene=up_entrez_ids_no_na$ENTREZID,
                                organism = "human",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")
    reactome_results_up <- as.data.frame(up_reactome)

    # WikiPathways
    up_wiki <- enrichWP(gene=up_entrez_ids_no_na$ENTREZID,
                        organism = "Homo sapiens",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    Wiki_results_up<- as.data.frame(up_wiki)

    BP_go_results_up <- symbols_list_to_mirna_list(BP_go_results_up, up_mirna)
    write.table(BP_go_results_up, "all_up_regulated_miRNA_enrichment/up_regulated_miRNA_GO_BP_results.tsv", sep='\t', quote=FALSE, row.names = FALSE)
    MF_go_results_up <- symbols_list_to_mirna_list(MF_go_results_up, up_mirna)
    write.table(MF_go_results_up, "all_up_regulated_miRNA_enrichment/up_regulated_miRNA_GO_MF_results.tsv", sep='\t', quote = FALSE, row.names = FALSE)
    CC_go_results_up <- symbols_list_to_mirna_list(CC_go_results_up, up_mirna)
    write.table(CC_go_results_up, "all_up_regulated_miRNA_enrichment/up_regulated_miRNA_GO_CC_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    if(nrow(KEGG_results_up) > 0){
        KEGG_results_up <- entrez_ids_to_symbol(KEGG_results_up, up_entrez_ids_no_na)
        KEGG_results_up <- symbols_list_to_mirna_list(KEGG_results_up, up_mirna)
    }
    write.table(KEGG_results_up, "all_up_regulated_miRNA_enrichment/up_regulated_miRNA_KEGG_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    if(nrow(Wiki_results_up) > 0){
        Wiki_results_up <- entrez_ids_to_symbol(Wiki_results_up, up_entrez_ids_no_na)
        Wiki_results_up <- symbols_list_to_mirna_list(Wiki_results_up, up_mirna)
    }
    write.table(Wiki_results_up, "all_up_regulated_miRNA_enrichment/up_regulated_miRNA_WikiPathways_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    if(nrow(reactome_results_up) > 0){
        reactome_results_up <- entrez_ids_to_symbol(reactome_results_up, up_entrez_ids_no_na)
        reactome_results_up <- symbols_list_to_mirna_list(reactome_results_up, up_mirna)
    }
    write.table(reactome_results_up, "all_up_regulated_miRNA_enrichment/up_regulated_miRNA_Reactome_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    log_info("Completed processing upregulated miRNA targets")
} else {
    log_info("No upregulated miRNA targets found. Skipping enrichment analysis.")
}

# Process DOWNregulated miRNA targets
log_info("Processing downregulated miRNA targets from: ", down_targets)

# Make a directory for downregulated miRNA target enrichment results
dir.create("all_down_regulated_miRNA_enrichment")

down_mirna <- read.table(down_targets, header=TRUE, sep='\t')

# Only when the dataframe has genes, perfrom enrichment analysis
if (nrow(down_mirna) != 0) {

    down_entrez_ids <- bitr(down_mirna$Target_Gene_Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    down_entrez_ids_no_na <- na.omit(down_entrez_ids)
    write.table(down_entrez_ids_no_na, "all_down_regulated_miRNA_enrichment/down_symbol_to_entrezid.tsv", sep='\t', quote=FALSE, row.names = FALSE)

    down_na_entrez <- down_entrez_ids[is.na(down_entrez_ids$ENTREZID), ]
    log_info("There are ", nrow(down_na_entrez), " genes with no Entrez ID in the downregulated miRNA targets.")

    write.table(down_na_entrez, "all_down_regulated_miRNA_enrichment/down_na_entrez.tsv", sep='\t', quote=FALSE, row.names = FALSE)

    # GO BP
    down_ego <- enrichGO(gene=down_entrez_ids_no_na$SYMBOL,
                    OrgDb= org.Hs.eg.db,
                    keyType= "SYMBOL",
                    ont= "BP",  # Biological Process
                    pAdjustMethod= "BH",
                    pvalueCutoff= 0.05)

    BP_go_results_down <- as.data.frame(down_ego)

    # GO MF
    down_go_mf <- enrichGO(gene=down_entrez_ids_no_na$SYMBOL,
                    OrgDb= org.Hs.eg.db,
                    keyType= "SYMBOL",
                    ont= "MF",
                    pAdjustMethod= "BH",
                    pvalueCutoff= 0.05)
    MF_go_results_down <- as.data.frame(down_go_mf)

    # GO CC
    down_go_cc <- enrichGO(gene=down_entrez_ids_no_na$SYMBOL,
                    OrgDb= org.Hs.eg.db,
                    keyType= "SYMBOL",
                    ont= "CC", # Cellular Component
                    pAdjustMethod= "BH",
                    pvalueCutoff= 0.05)
    CC_go_results_down <- as.data.frame(down_go_cc)

    # KEGG requires entrez ids
    down_kegg <- enrichKEGG(gene=down_entrez_ids_no_na$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    KEGG_results_down <- as.data.frame(down_kegg)

    # Reactome
    down_reactome <- enrichPathway(gene=down_entrez_ids_no_na$ENTREZID,
                                organism = "human",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")
    reactome_results_down <- as.data.frame(down_reactome)

    # WikiPathways
    down_wiki <- enrichWP(gene=down_entrez_ids_no_na$ENTREZID,
                        organism = "Homo sapiens",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    Wiki_results_down<- as.data.frame(down_wiki)

    BP_go_results_down <- symbols_list_to_mirna_list(BP_go_results_down, down_mirna)
    write.table(BP_go_results_down, "all_down_regulated_miRNA_enrichment/down_regulated_miRNA_GO_BP_results.tsv", sep='\t', quote=FALSE, row.names = FALSE)
    MF_go_results_down <- symbols_list_to_mirna_list(MF_go_results_down, down_mirna)
    write.table(MF_go_results_down, "all_down_regulated_miRNA_enrichment/down_regulated_miRNA_GO_MF_results.tsv", sep='\t', quote = FALSE, row.names = FALSE)
    CC_go_results_down <- symbols_list_to_mirna_list(CC_go_results_down, down_mirna)
    write.table(CC_go_results_down, "all_down_regulated_miRNA_enrichment/down_regulated_miRNA_GO_CC_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    if(nrow(KEGG_results_down) > 0){
        KEGG_results_down <- entrez_ids_to_symbol(KEGG_results_down, down_entrez_ids_no_na)
        KEGG_results_down <- symbols_list_to_mirna_list(KEGG_results_down, down_mirna)
    }
    write.table(KEGG_results_down, "all_down_regulated_miRNA_enrichment/down_regulated_miRNA_KEGG_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    if(nrow(Wiki_results_down) > 0){
        Wiki_results_down <- entrez_ids_to_symbol(Wiki_results_down, down_entrez_ids_no_na)
        Wiki_results_down <- symbols_list_to_mirna_list(Wiki_results_down, down_mirna)
    }
    write.table(Wiki_results_down, "all_down_regulated_miRNA_enrichment/down_regulated_miRNA_WikiPathways_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    if(nrow(reactome_results_down) > 0){
        reactome_results_down <- entrez_ids_to_symbol(reactome_results_down, down_entrez_ids_no_na)
        reactome_results_down <- symbols_list_to_mirna_list(reactome_results_down, down_mirna)
    }
    write.table(reactome_results_down, "all_down_regulated_miRNA_enrichment/down_regulated_miRNA_Reactome_results.tsv", sep = '\t', quote=FALSE, row.names = FALSE)

    log_info("Completed processing downregulated miRNA targets")
} else {
    log_info("No downregulated miRNA targets found. Skipping enrichment analysis.")
}

log_info("ClusterProfiler analysis completed successfully.")