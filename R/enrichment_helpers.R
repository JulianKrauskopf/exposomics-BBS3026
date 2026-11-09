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
  # Map gene symbols and Ensembl IDs independently. Query only IDs that
  # actually exist in each annotation key space: mapIds() errors if *none*
  # of its input keys are valid for that key type.
  valid_symbols <- intersect(unique(keys[!is.na(keys)]),
                             AnnotationDbi::keys(db, keytype = "SYMBOL"))
  valid_ensembl <- intersect(unique(clean_keys[!is.na(clean_keys)]),
                             AnnotationDbi::keys(db, keytype = "ENSEMBL"))

  symbol_ids <- if (length(valid_symbols)) {
    AnnotationDbi::mapIds(db, keys = valid_symbols, column = "ENTREZID",
                          keytype = "SYMBOL", multiVals = "first")
  } else {
    setNames(character(), character())
  }
  ensembl_ids <- if (length(valid_ensembl)) {
    AnnotationDbi::mapIds(db, keys = valid_ensembl, column = "ENTREZID",
                          keytype = "ENSEMBL", multiVals = "first")
  } else {
    setNames(character(), character())
  }

  # Prefer the appropriate identifier format; allow mixed tables.
  mapped_symbol <- unname(symbol_ids[keys])
  mapped_ensembl <- unname(ensembl_ids[clean_keys])
  mapped <- ifelse(!is.na(mapped_symbol), mapped_symbol, mapped_ensembl)
  key_type <- if (sum(!is.na(mapped_symbol)) >= sum(!is.na(mapped_ensembl))) {
    "SYMBOL (with Ensembl fallback)"
  } else {
    "ENSEMBL (with symbol fallback)"
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
