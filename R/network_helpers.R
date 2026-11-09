# Functions for Session 6: pathways and protein interaction networks.
# Run the notebook chunks; these functions perform the longer R steps.
#
# RCy3 controls Cytoscape from RStudio. Cytoscape Desktop must be open.
# Data are gene-expression regression estimates, NOT case/control log2 fold changes.

check_cytoscape <- function(apps=c("WikiPathways","stringApp")) {
  check_packages("RCy3")
  tryCatch(RCy3::cytoscapePing(), error=function(e) {
    stop("Open Cytoscape Desktop and try again. Details: ",conditionMessage(e))
  })
  for (app in apps) {
    current <- tolower(as.character(unlist(RCy3::getInstalledApps())))
    if (!any(grepl(tolower(app), current, fixed=TRUE))) {
      message("Installing the Cytoscape app: ",app)
      RCy3::installApp(app)
    }
  }
  invisible(TRUE)
}

# Match the expression features to human gene names and identifiers.
prepare_network_annotations <- function(results) {
  check_packages(c("AnnotationDbi","org.Hs.eg.db"))
  must <- c("feature_id","logFC","P.Value","adj.P.Val")
  if (!all(must %in% names(results)))
    stop("Use the limma results from Practical 3 (feature_id, logFC and p-values).")

  id <- as.character(results$feature_id)
  id_clean <- sub("\\.[0-9]+$", "", id)
  ens <- grepl("^ENSG[0-9]+$", id_clean)
  symbol <- rep(NA_character_,length(id))
  ensembl <- rep(NA_character_,length(id))
  db <- org.Hs.eg.db::org.Hs.eg.db

  if (any(ens)) {
    ensembl[ens] <- id_clean[ens]
    valid <- intersect(unique(ensembl[ens]), AnnotationDbi::keys(db,keytype="ENSEMBL"))
    if (length(valid)) {
      lookup <- AnnotationDbi::mapIds(db, keys=valid, keytype="ENSEMBL",
                                     column="SYMBOL", multiVals="first")
      symbol[ens] <- unname(lookup[match(ensembl[ens], names(lookup))])
    }
  }
  if (any(!ens)) {
    symbol[!ens] <- id[!ens]
    valid <- intersect(unique(symbol[!ens]), AnnotationDbi::keys(db,keytype="SYMBOL"))
    if (length(valid)) {
      lookup <- AnnotationDbi::mapIds(db, keys=valid, keytype="SYMBOL",
                                     column="ENSEMBL", multiVals="first")
      ensembl[!ens] <- unname(lookup[match(symbol[!ens], names(lookup))])
    }
  }

  # WikiPathways nodes may use Entrez/NCBI Gene IDs even when the
  # expression matrix contains gene symbols or Ensembl IDs.
  entrez <- rep(NA_character_,length(id))
  available_symbols <- unique(symbol[!is.na(symbol) & nzchar(symbol)])
  valid_symbols <- intersect(
    available_symbols,AnnotationDbi::keys(db,keytype="SYMBOL")
  )
  if (length(valid_symbols)) {
    entrez_lookup <- AnnotationDbi::mapIds(
      db,keys=valid_symbols,keytype="SYMBOL",
      column="ENTREZID",multiVals="first"
    )
    entrez <- unname(entrez_lookup[match(symbol,names(entrez_lookup))])
  }

  # Nominal p < 0.05 is the course-wide teaching/visualization criterion;
  # disease-set and transcriptome-wide BH FDR are retained for comparison.
  focused <- "FDR_disease_set" %in% names(results)
  analysis_fdr <- if (focused) results$FDR_disease_set else results$adj.P.Val
  tab <- data.frame(
    feature_id=as.character(results$feature_id),
    GeneName=symbol, Ensembl=ensembl, Entrez=entrez,
    logFC=as.numeric(results$logFC),
    moderated_t=if ("t" %in% names(results))
      as.numeric(results$t) else rep(NA_real_,nrow(results)),
    P.Value=as.numeric(results$P.Value),
    adj.P.Val=as.numeric(results$adj.P.Val),
    FDR_analysis=as.numeric(analysis_fdr),
    StudentNominal=as.integer(!is.na(results$P.Value) &
                               results$P.Value<0.05),
    StudentFDR=as.integer(!is.na(analysis_fdr) & analysis_fdr<0.05),
    stringsAsFactors=FALSE
  )
  message("Course network highlights use nominal p < 0.05. ",
          "Adjusted values remain visible for multiple-testing discussion.")
  tab <- tab[order(tab$P.Value, na.last=TRUE),,drop=FALSE]
  message("Annotated gene symbols: ",sum(!is.na(tab$GeneName)),
          " / ",nrow(tab), "; Ensembl IDs: ",sum(!is.na(tab$Ensembl)))
  tab
}

# Select nominally associated genes for the STRING network.
select_network_genes <- function(annotations, mode=c("nominal","exploratory_top"),
                                 n_top=30L, max_genes=80L) {
  mode <- match.arg(mode)
  tab <- annotations[!is.na(annotations$GeneName) &
                     nzchar(annotations$GeneName) &
                     is.finite(annotations$P.Value),,drop=FALSE]
  if (mode=="nominal")
    tab <- tab[!is.na(tab$P.Value) & tab$P.Value<0.05,,drop=FALSE]
  tab <- tab[order(tab$P.Value),,drop=FALSE]
  tab <- tab[!duplicated(tab$GeneName),,drop=FALSE]
  tab <- head(tab,if (mode=="nominal") max_genes else min(n_top,max_genes))
  message("Selected ",nrow(tab)," genes (mode: ",mode,").")
  if (mode=="nominal")
    message("Nominal p < 0.05 is an uncorrected EXPLORATORY network filter, ",
            "not FDR-significant evidence.")
  if (mode=="exploratory_top")
    message("These are exploratory top-ranked genes, NOT FDR-significant discoveries.")
  tab
}

# Match gene names to the identifiers in the Cytoscape network.
match_node_table <- function(data, network=NULL,
                             candidates=list(
                               c("Ensembl","Ensembl"),
                               c("XrefId","Ensembl"),
                               c("query term","GeneName"),
                               c("display name","GeneName"),
                               c("name","GeneName"))) {
  cols <- as.character(unlist(RCy3::getTableColumnNames("node",network=network)))
  best <- list(matches=0L)
  for (candidate in candidates) {
    col <- candidate[1]; input <- candidate[2]
    if (!col %in% cols) next
    vals <- RCy3::getTableColumns(table="node",columns=col,network=network)[[col]]
    n <- sum(!is.na(vals) & as.character(vals) %in% data[[input]])
    if (n>best$matches)
      best <- list(table_key=col,data_key=input,matches=n)
  }
  if (best$matches<1L)
    stop("No matching gene identifiers in this Cytoscape network. Check symbol/Ensembl IDs and the pathway annotations.")
  message("Matched ",best$matches," nodes with ",best$data_key,
          " ↔ ",best$table_key)
  best
}

# WikiPathways can identify genes by symbols, Ensembl or Entrez IDs.
match_pathway_nodes <- function(annotations, nodes) {
  required <- c("GeneName","Ensembl","Entrez","logFC","P.Value",
                "StudentNominal")
  if (!all(required %in% names(annotations)))
    stop("Missing analysis columns. Rerun prepare_network_annotations().")

  if (!"SUID" %in% names(nodes))
    stop("Cytoscape did not return node SUIDs. Check the imported pathway.")

  clean <- function(x) {
    x <- as.character(unlist(x,use.names=FALSE))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character())
    for (mark in c("[","]",intToUtf8(34L),"'"))
      x <- gsub(mark,"",x,fixed=TRUE)
    values <- trimws(unlist(strsplit(x,"[,;|]")))
    values[nzchar(values)]
  }
  lookup <- function(values, reference) {
    if (!length(values)) return(NA_integer_)
    values <- values[!is.na(values) & nzchar(values)]
    idx <- match(values,reference)
    idx <- idx[!is.na(idx)]
    if (length(idx)) idx[[1]] else NA_integer_
  }
  matching <- rep(NA_integer_,nrow(nodes))
  for (i in seq_len(nrow(nodes))) {
    cell <- function(col) {
      if (col %in% names(nodes)) clean(nodes[[col]][i]) else character()
    }
    ensembl <- cell("Ensembl")
    xref <- cell("XrefId")
    source <- tolower(paste(cell("XrefDatasource"),collapse=" "))
    ensembl <- c(ensembl, xref[grepl("^ENSG[0-9]+",xref)])
    if (grepl("ensembl",source,fixed=TRUE)) ensembl <- c(ensembl,xref)
    ensembl <- sub("\\.[0-9]+$","",ensembl)
    matching[i] <- lookup(ensembl,annotations$Ensembl)
    if (!is.na(matching[i])) next

    entrez <- if (grepl("entrez|ncbi gene|geneid",source)) xref else character()
    matching[i] <- lookup(entrez,annotations$Entrez)
    if (!is.na(matching[i])) next

    symbols <- c(cell("GeneName"),cell("display name"),cell("label"),
                 cell("name"))
    if (grepl("hgnc|gene symbol",source)) symbols <- c(symbols,xref)
    matching[i] <- lookup(symbols,annotations$GeneName)
  }
  keep <- which(!is.na(matching))
  if (!length(keep)) {
    message("No overlapping gene identifiers were found in this pathway. ",
            "It may not contain any of the modelled genes; try a more ",
            "biologically relevant human WikiPathways ID.")
    return(data.frame())
  }
  selected <- annotations[matching[keep],,drop=FALSE]
  result <- data.frame(
    OverlayNodeSUID=as.numeric(nodes$SUID[keep]),
    CourseGeneSymbol=selected$GeneName,
    logFC=selected$logFC,
    P.Value=selected$P.Value,
    adj.P.Val=selected$adj.P.Val,
    StudentNominal=selected$StudentNominal,
    stringsAsFactors=FALSE
  )
  if (anyDuplicated(result$OverlayNodeSUID))
    stop("Duplicate pathway node SUIDs; cannot safely add annotations.")
  message("Matched ",nrow(result)," of ",nrow(nodes),
          " pathway nodes using gene symbols/Ensembl/Entrez identifiers.")
  result
}

# Color the pathway using the full-transcriptome expression results.
# The individual disease-gene analysis is in notebook 04.
overlay_pathway <- function(annotations, pathway_id) {
  if (!grepl("^WP[0-9]+$",pathway_id))
    stop("Use a WikiPathways identifier such as WP179.")
  RCy3::commandsRun(paste0("wikipathways import-as-pathway id=",pathway_id))
  net <- RCy3::getNetworkSuid()
  if (RCy3::getNodeCount(network=net)<1L)
    stop("Imported pathway contains no nodes.")

  nodes <- RCy3::getTableColumns(table="node",network=net)
  mapping <- match_pathway_nodes(annotations,nodes)
  matched <- nrow(mapping)
  if (matched>0L) {
    RCy3::loadTableData(
      mapping,data.key.column="OverlayNodeSUID",
      table.key.column="SUID",network=net
    )
    lim <- max(abs(mapping$logFC),na.rm=TRUE)
    if (!is.finite(lim)||lim==0) lim <- 1
    RCy3::setNodeColorMapping("logFC",c(-lim,0,lim),
      c("#3B70B5","#FFFFFF","#D46B39"),
      default.color="#E6E6E6",style.name="WikiPathways",network=net)
    RCy3::setNodeBorderColorMapping("StudentNominal",c(0L,1L),
      c("#999999","#009900"),mapping.type="d",
      style.name="WikiPathways",network=net)
    RCy3::setNodeBorderWidthMapping("StudentNominal",c(0L,1L),
      c(1,6),mapping.type="d",style.name="WikiPathways",network=net)
    RCy3::setVisualStyle("WikiPathways",network=net)
  } else {
    warning("The pathway was imported but no genes could be annotated. ",
            "An uncolored pathway is not an analysis result. Try another ",
            "human pathway and check the available node identifiers.",
            call.=FALSE)
  }
  RCy3::fitContent(network=net)
  safe <- gsub("[^A-Za-z0-9_-]","_",pathway_id)
  RCy3::exportImage(file.path("outputs",paste0("pathway_",safe)),
                   type="PNG",zoom=250,network=net)
  message("Saved outputs/pathway_",safe,".png (",matched,
          " matched nodes).")
  list(network=net,matched_nodes=matched)
}

# Read protein interactions for calculating degree and betweenness.
string_graph_from_cytoscape <- function(network) {
  check_packages("igraph")
  node_table <- RCy3::getTableColumns(
    table="node",columns="name",network=network
  )
  edge_table <- RCy3::getTableColumns(
    table="edge",columns="name",network=network
  )
  if (!"name" %in% names(node_table) ||
      !"name" %in% names(edge_table))
    stop("STRING network node/edge tables do not have name columns.")

  nodes <- unique(as.character(node_table[["name"]]))
  nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
  edges <- as.character(edge_table[["name"]])
  edges <- edges[!is.na(edges) & nzchar(edges)]
  if (!length(nodes) || !length(edges))
    stop("STRING returned no usable nodes or interaction edges.")

  parsed <- regmatches(edges,regexec(
    "^(.+) \\([^)]*\\) (.+)$",edges
  ))
  if (any(lengths(parsed)!=3L)) {
    stop("Some STRING edge names do not use the usual Cytoscape ",
         "'source (interaction) target' format. Check the edge table; ",
         "these names cannot be converted safely.")
  }
  ends <- data.frame(
    from=vapply(parsed,function(x)x[[2L]],character(1)),
    to=vapply(parsed,function(x)x[[3L]],character(1)),
    stringsAsFactors=FALSE
  )
  if (!all(ends$from %in% nodes & ends$to %in% nodes)) {
    stop("STRING edge endpoints do not match the node name column. ",
         "Inspect Cytoscape's node and edge tables for their identifiers.")
  }
  igraph::graph_from_data_frame(
    ends,directed=FALSE,vertices=data.frame(name=nodes)
  )
}

# STRING network: import human proteins with no added interactors, calculate
# graph degree/betweenness for the largest component, then style in Cytoscape.
# The gene filter and green borders use NOMINAL p < 0.05 for teaching.
build_ppi_network <- function(input_genes,cutoff=0.9) {
  check_packages(c("RCy3","igraph"))
  stage <- "validating the input genes"
  tryCatch({
    if (!cutoff %in% c(0.4,0.7,0.9))
      stop("Choose STRING confidence 0.4, 0.7 or 0.9.")
    if (!is.data.frame(input_genes) ||
        !all(c("GeneName","logFC","P.Value","StudentNominal") %in%
             names(input_genes)))
      stop("Run the nominal-p network gene-selection chunk first.")
    if (nrow(input_genes)<3L) {
      message("Fewer than three genes; no STRING network was made. ",
              "Try selection_mode='exploratory_top' and label it exploratory.")
      return(NULL)
    }
    syms <- unique(as.character(input_genes$GeneName))
    syms <- syms[!is.na(syms) & grepl("^[A-Za-z0-9._-]+$",syms)]
    if (length(syms)<3L) {
      message("Fewer than three valid gene symbols for STRING.")
      return(NULL)
    }
    dir.create("outputs",showWarnings=FALSE)
    stage <- "sending the STRING query to Cytoscape"
    cmd <- paste0('string protein query query="',
                  paste(syms,collapse=","),'"',
                  ' cutoff=',cutoff,
                  ' limit=0 species="Homo sapiens"')
    imported <- tryCatch({
      RCy3::commandsPOST(cmd)
      TRUE
    },error=function(e) {
      message("STRING POST request failed: ",conditionMessage(e))
      FALSE
    })
    if (!imported) {
      if (nchar(cmd,type="bytes")>1800L)
        stop("STRING POST failed and the gene list is too long for ",
             "GET fallback. Try fewer genes (e.g., exploratory top 30).")
      message("Retrying the STRING import with commandsRun()...")
      RCy3::commandsRun(cmd)
    }

    stage <- "reading the imported STRING network"
    net <- RCy3::getNetworkSuid()
    if (RCy3::getNodeCount(network=net)<3L) {
      message("STRING returned fewer than three mapped proteins.")
      return(NULL)
    }
    if (RCy3::getEdgeCount(network=net)<1L) {
      message("STRING returned no interaction edges at cutoff ",cutoff,
              ". Try the scientifically justified sensitivity cutoff 0.7 or 0.4.")
      return(NULL)
    }

    stage <- "constructing the interaction graph from node and edge names"
    ig <- string_graph_from_cytoscape(net)
    if (igraph::ecount(ig)<1L) {
      message("The STRING network contains no interactions.")
      return(NULL)
    }

    stage <- "selecting the main connected component"
    parts <- igraph::components(ig)
    largest <- which.max(parts$csize)
    keep <- igraph::V(ig)$name[parts$membership==largest]
    if (length(keep)<3L) {
      message("Largest connected component has fewer than three proteins.")
      return(NULL)
    }
    gr <- igraph::induced_subgraph(ig,vids=keep)
    if (length(keep)<igraph::vcount(ig)) {
      net <- RCy3::createSubnetwork(
        nodes=as.character(keep),
        nodes.by.col="name",
        subnetwork.name="Exposomics PPI - main component",
        network=net
      )
    }
    stats <- data.frame(
      name=as.character(igraph::V(gr)$name),
      StudentDegree=as.integer(igraph::degree(gr)),
      StudentBetweenness=as.numeric(igraph::betweenness(
        gr,directed=FALSE,normalized=TRUE
      )),
      stringsAsFactors=FALSE
    )

    stage <- "retrieving protein labels"
    columns <- as.character(unlist(
      RCy3::getTableColumnNames("node",network=net)
    ))
    protein_label <- if ("display name" %in% columns)
      "display name" else "name"
    node_names <- RCy3::getTableColumns(
      table="node",columns=unique(c("name",protein_label)),network=net
    )
    stats$Protein <- as.character(
      node_names[[protein_label]][
        match(stats$name,as.character(node_names[["name"]]))
      ]
    )

    stage <- "loading degree, betweenness and exposure results"
    RCy3::loadTableData(
      stats,data.key.column="name",table.key.column="name",network=net
    )
    keys <- match_node_table(
      input_genes,network=net,
      candidates=list(c("query term","GeneName"),
                      c("display name","GeneName"),
                      c("name","GeneName"))
    )
    RCy3::loadTableData(
      input_genes,data.key.column=keys$data_key,
      table.key.column=keys$table_key,network=net
    )

    stage <- "styling the network"
    style <- "Exposomics PPI"
    if (!style %in% RCy3::getVisualStyleNames()) {
      RCy3::createVisualStyle(
        style,
        defaults=list(NODE_SIZE=50,NODE_FILL_COLOR="#E6E6E6",
                      EDGE_TRANSPARENCY=130),
        mappings=list()
      )
    }
    lim <- max(abs(input_genes$logFC),na.rm=TRUE)
    if (!is.finite(lim) || lim==0) lim <- 1
    RCy3::setNodeColorMapping(
      "logFC",c(-lim,0,lim),c("#3B70B5","#FFFFFF","#D46B39"),
      default.color="#E6E6E6",style.name=style,network=net
    )
    RCy3::setNodeSizeMapping(
      "StudentDegree",c(0,max(1L,stats$StudentDegree)),
      sizes=c(35,110),style.name=style,network=net
    )
    RCy3::setNodeBorderColorMapping(
      "StudentNominal",c(0L,1L),c("#999999","#009900"),
      mapping.type="d",style.name=style,network=net
    )
    labels <- as.character(unlist(
      RCy3::getTableColumnNames("node",network=net)
    ))
    RCy3::setNodeLabelMapping(
      if ("display name" %in% labels) "display name" else "name",
      style.name=style,network=net
    )
    RCy3::setVisualStyle(style,network=net)
    RCy3::layoutNetwork("force-directed",network=net)
    RCy3::fitContent(network=net)

    stage <- "exporting the PPI results"
    RCy3::exportImage(
      "outputs/PPI_network",type="PNG",zoom=250,network=net
    )
    write.csv(
      stats[order(-stats$StudentDegree),],
      "outputs/PPI_hub_statistics.csv",row.names=FALSE
    )
    message("Saved the STRING network and hub table. ",
            "Green borders indicate nominal p < 0.05 only.")
    list(network=net,stats=stats,matched_nodes=keys$matches,
         requested_genes=length(syms),
         mode="largest connected component")
  },error=function(e) {
    stop("STRING PPI failed while ",stage,": ",
         conditionMessage(e),
         ". If this continues, run traceback() immediately after ",
         "the error and check the active Cytoscape network.",
         call.=FALSE)
  })
}
