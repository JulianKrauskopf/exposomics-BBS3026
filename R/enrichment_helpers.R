# Session 5 — Convert expression-result gene identifiers to Entrez IDs.
# Many enrichment packages require Entrez IDs rather than gene symbols.
map_results_to_entrez <- function(results) {
  required <- c("feature_id", "t", "P.Value", "adj.P.Val")
  if (!all(required %in% names(results))) {
    stop("Expected Session 4 full limma results with columns: ",
         paste(required, collapse = ", "))
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Run source('R/setup_packages.R') to install annotation packages.")
  }
  keys <- as.character(results$feature_id)
  # Remove an optional Ensembl version suffix.
  clean_keys <- sub("\\.[0-9]+$", "", keys)
  db <- org.Hs.eg.db::org.Hs.eg.db
  # Determine whether the expression row names are gene symbols or Ensembl IDs.
  symbol_ids <- AnnotationDbi::mapIds(
    db, keys = unique(keys), column = "ENTREZID",
    keytype = "SYMBOL", multiVals = "first"
  )
  ensembl_ids <- AnnotationDbi::mapIds(
    db, keys = unique(clean_keys), column = "ENTREZID",
    keytype = "ENSEMBL", multiVals = "first"
  )
  if (sum(!is.na(ensembl_ids)) > sum(!is.na(symbol_ids))) {
    mapped <- unname(ensembl_ids[clean_keys])
    key_type <- "ENSEMBL"
  } else {
    mapped <- unname(symbol_ids[keys])
    key_type <- "SYMBOL"
  }
  out <- results
  out$ENTREZID <- as.character(mapped)
  out <- out[!is.na(out$ENTREZID) & nzchar(out$ENTREZID) &
               is.finite(out$t) & is.finite(out$P.Value), , drop = FALSE]
  # A gene can have multiple expression features. Retain one ranking value
  # per Entrez gene, chosen by largest absolute t statistic.
  out <- out[order(-abs(out$t), out$feature_id), , drop = FALSE]
  out <- out[!duplicated(out$ENTREZID), , drop = FALSE]
  message("Mapped ", nrow(out), " unique genes using ", key_type,
          " identifiers for enrichment analysis.")
  if (!nrow(out)) stop("No gene identifiers could be mapped.")
  out
}
