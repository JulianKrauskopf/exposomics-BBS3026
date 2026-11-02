# Shared functions for the Exposomics course notebooks.
# These functions load your data, prepare exposures, fit models,
# create plots, and help identify disease-associated genes.
# Usually you only need to edit the settings inside the notebooks.
check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Install packages: ", paste(missing, collapse=", "))
}
make_folders <- function() {
  dir.create("outputs", showWarnings=FALSE)
  dir.create("data/gene_sets", recursive=TRUE, showWarnings=FALSE)
}
# Generate a reproducible, unitless AIR POLLUTION teaching index.
# These values are simulated from participant age, cohort and random noise.
# They are NOT measured air pollution and have no known biological effects.
# Keep the measured chemicals in course$expo separate from this exercise.
simulate_air_pollution <- function(covar, seed=20261022L) {
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if (had_seed) previous_seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed",previous_seed,envir=.GlobalEnv)
    } else if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) {
      rm(".Random.seed",envir=.GlobalEnv)
    }
  },add=TRUE)
  set.seed(seed)
  age <- as.numeric(covar$age)
  age_center <- stats::median(age,na.rm=TRUE)
  age[!is.finite(age)] <- age_center
  groups <- sort(unique(as.character(stats::na.omit(covar$cohort))))
  cohort_term <- if (length(groups)>1L) {
    as.numeric(as.character(covar$cohort)==groups[2L])
  } else rep(0,nrow(covar))
  cohort_term[is.na(cohort_term)] <- 0
  # Designed to show a skewed distribution and illustrative associations
  # with age and cohort, not a measured exposure or health relationship.
  simulated_index <- round(2^(5 + 0.025*(age-age_center) +
                             0.65*cohort_term +
                             stats::rnorm(nrow(covar),sd=0.55)),3)
  result <- data.frame("Air pollution"=simulated_index,check.names=FALSE)
  rownames(result) <- rownames(covar)
  result
}

load_course_data <- function() {
  files <- c("data/covar.csv", "data/expo.csv")
  if (!all(file.exists(files))) stop("The course includes data/covar.csv and data/expo.csv. Make sure you opened the downloaded project folder in RStudio; see SETUP.md.")
  covar <- read.csv(files[1], row.names=1, check.names=FALSE)
  expo <- read.csv(files[2], row.names=1, check.names=FALSE)
  if (anyDuplicated(rownames(covar)) || anyDuplicated(rownames(expo)))
    stop("Duplicate sample IDs.")
  ids <- rownames(covar)[rownames(covar) %in% rownames(expo)]
  if (length(ids) < 10L) stop("Fewer than 10 shared participant IDs.")
  covar <- covar[ids,,drop=FALSE]; expo <- expo[ids,,drop=FALSE]
  for (v in intersect(c("sex","cohort","smoking_status","education","activity"),names(covar)))
    covar[[v]] <- factor(covar[[v]])
  list(covar=covar, expo=expo,
       simulated_expo=simulate_air_pollution(covar))
}
# Project-level participant selection. The study stores sex as F/M and
# cohort as E/N; course labels cohort 1=E and cohort 2=N.
# Keep exactly the same stratum throughout the individual project.
select_project_stratum <- function(covar, stratum="all") {
  valid <- c("all","male","female","cohort_1","cohort_2")
  if (!is.character(stratum) || length(stratum)!=1L ||
      is.na(stratum) || !stratum %in% valid)
    stop("Choose my_stratum from: ",paste(valid,collapse=", "))
  if (!all(c("sex","cohort") %in% names(covar)))
    stop("Participant data must include sex and cohort.")
  sex <- as.character(covar$sex)
  cohort <- as.character(covar$cohort)
  keep <- switch(stratum,
    all=rep(TRUE,nrow(covar)),
    male=!is.na(sex) & sex=="M",
    female=!is.na(sex) & sex=="F",
    cohort_1=!is.na(cohort) & cohort=="E",
    cohort_2=!is.na(cohort) & cohort=="N"
  )
  d <- droplevels(covar[keep,,drop=FALSE])
  if (nrow(d)<10L)
    stop("The selected stratum has fewer than 10 participants. ",
         "Choose a larger group before fitting models.")
  label <- switch(stratum,all="All participants",male="Males only (M)",
                  female="Females only (F)",cohort_1="Cohort 1 (E)",
                  cohort_2="Cohort 2 (N)")
  message("Study stratum: ",label,"; ",nrow(d)," of ",
          nrow(covar)," participants before model-specific missingness.")
  d
}

# A factor with only one observed category cannot be included in a
# regression design. Drop constant adjustments in sex/cohort strata;
# retain all remaining scientifically justified varying covariates.
project_adjustments <- function(data,
                                candidates=c("age","sex","cohort"),
                                exclude=character()) {
  candidates <- setdiff(unique(candidates),exclude)
  if (!all(candidates %in% names(data)))
    stop("Unknown adjustment variable(s): ",
         paste(setdiff(candidates,names(data)),collapse=", "))
  varying <- vapply(candidates,function(v) {
    z <- data[[v]]
    if (is.factor(z)) z <- as.character(z)
    length(unique(z[!is.na(z)]))>1L
  },logical(1))
  dropped <- candidates[!varying]
  if (length(dropped))
    message("Not adjusting for constant variable(s) in this stratum: ",
            paste(dropped,collapse=", "))
  candidates[varying]
}

prepare_exposure <- function(covar, expo, exposure) {
  if (!exposure %in% colnames(expo)) stop("Unknown exposure: ",exposure)
  if (!is.numeric(expo[[exposure]])) stop("Exposure is not numeric.")
  ids <- rownames(covar)[rownames(covar) %in% rownames(expo)]
  x <- covar[ids,,drop=FALSE]; values <- expo[ids,exposure]
  x$exposure_raw <- values
  x$log2_exposure <- ifelse(is.finite(values) & values>0, log2(pmax(values,.Machine$double.xmin)), NA_real_)
  if (any(!is.finite(values) | values<=0, na.rm=TRUE))
    warning("Nonpositive/nonfinite exposures excluded from log2 analysis; review detection limits.")
  x
}
run_exposure_models <- function(data, predictor, adjust_for=c("age","sex","cohort")) {
  adjust_for <- setdiff(unique(adjust_for),predictor)
  vars <- unique(c("log2_exposure",predictor,adjust_for))
  if (!all(vars %in% names(data))) stop("Unknown model column(s): ",paste(setdiff(vars,names(data)),collapse=", "))
  d <- droplevels(data[complete.cases(data[,vars,drop=FALSE]),,drop=FALSE])
  if (nrow(d)<10) stop("Not enough complete cases.")
  for (v in vars) if (length(unique(d[[v]]))<2) stop("Constant variable: ",v)
  m1 <- lm(reformulate(predictor,response="log2_exposure"),data=d)
  m2 <- lm(reformulate(c(predictor,adjust_for),response="log2_exposure"),data=d)
  if (anyNA(coef(m2))) warning("Rank-deficient model; review covariates.")
  list(crude=m1,adjusted=m2,data=d)
}
load_expression <- function(file="data/EGM_Transcriptomics_GeneSymbol_Full.Rdata",object_name="omics") {
  if (!file.exists(file)) stop("Gene-expression file missing: ",file, ". Download EGM_Transcriptomics_GeneSymbol_Full.Rdata from the Dropbox link in SETUP.md and save it inside data/.")
  e <- new.env(parent=baseenv()); load(file,envir=e)
  if (!exists(object_name,envir=e,inherits=FALSE)) stop("No object named ",object_name," found.")
  x <- as.matrix(get(object_name,envir=e,inherits=FALSE))
  if (!is.numeric(x) || is.null(rownames(x)) || is.null(colnames(x))) stop("Expected numeric genes-rows, participants-columns matrix with dimnames.")
  if (anyDuplicated(colnames(x))) stop("Duplicate expression sample IDs.")
  x
}
run_limma_exposure <- function(expression, exposure_data, adjust_for=c("age","sex","cohort"), normalize=FALSE) {
  check_packages("limma")
  vars <- c("log2_exposure",adjust_for)
  if (!all(vars %in% names(exposure_data))) stop("Missing model covariate.")
  ids <- rownames(exposure_data)[rownames(exposure_data) %in% colnames(expression)]
  if (length(ids)<10) stop("Insufficient matched expression samples.")
  d <- exposure_data[ids,,drop=FALSE]
  keep <- complete.cases(d[,vars,drop=FALSE])
  d <- droplevels(d[keep,,drop=FALSE])
  x <- expression[,rownames(d),drop=FALSE]
  stopifnot(identical(rownames(d),colnames(x)))
  if (any(!is.finite(x))) stop("Expression matrix has missing or infinite values; review preprocessing.")
  if (normalize) x <- limma::normalizeBetweenArrays(x)
  design <- model.matrix(reformulate(vars),data=d)
  if (qr(design)$rank<ncol(design) || nrow(design)<=ncol(design)+2)
    stop("Model not estimable: covariates are collinear or sample size is insufficient.")
  fit <- limma::eBayes(limma::lmFit(x,design))
  result <- limma::topTable(fit,coef="log2_exposure",number=Inf,adjust.method="BH",sort.by="P")
  result$feature_id <- rownames(result)
  list(results=result,fit=fit,design=design,metadata=d,expression=x)
}
# Match the selected disease-associated genes to the measured gene features.
# Choose the disease and retrieve its gene list before fitting models.
define_disease_gene_set <- function(expression, genes) {
  needed <- c("ensembl_id","gene_symbol","association_score")
  if (!all(needed %in% names(genes)))
    stop("Disease list is missing expected Open Targets gene columns.")
  features <- rownames(expression)
  if (is.null(features) || anyDuplicated(features))
    stop("Expression features must have unique gene identifiers.")
  if (is.null(nrow(genes)) || !nrow(genes))
    stop("Choose a nonempty Open Targets disease gene list before modelling.")

  gene_symbols <- as.character(genes$gene_symbol)
  ensembl_ids <- sub("\\.[0-9]+$", "", as.character(genes$ensembl_id))
  clean_features <- sub("\\.[0-9]+$", "", features)
  symbol_index <- match(features,gene_symbols)
  ensembl_index <- match(clean_features,ensembl_ids)
  match_index <- ifelse(!is.na(symbol_index),symbol_index,ensembl_index)
  selected <- which(!is.na(match_index))
  if (!length(selected))
    stop("None of the disease-associated genes match the expression matrix. ",
         "Check whether it uses gene symbols or Ensembl identifiers.")
  selected_gene <- match_index[selected]
  predefined <- data.frame(
    feature_id=features[selected],
    disease_gene_symbol=gene_symbols[selected_gene],
    ensembl_id=as.character(genes$ensembl_id[selected_gene]),
    association_score=as.numeric(genes$association_score[selected_gene]),
    stringsAsFactors=FALSE
  )
  if (anyDuplicated(predefined$feature_id))
    stop("The predefined gene set includes duplicate expression features.")
  message("Predefined hypothesis family: ",nrow(predefined),
          " measured expression features matching ",
          length(unique(selected_gene))," disease-associated genes.")
  predefined
}

# Interpret nominal p < 0.05 and retain adjusted p-values for comparison.
disease_focused_results <- function(all_results, predefined_genes) {
  needed <- c("feature_id","P.Value","adj.P.Val","logFC")
  if (!all(needed %in% names(all_results)))
    stop("Expected full-transcriptome limma results with p-values.")
  if (!all(c("feature_id","disease_gene_symbol","association_score") %in%
           names(predefined_genes)) || !nrow(predefined_genes))
    stop("First define your disease gene set before running the model.")
  if (anyDuplicated(predefined_genes$feature_id))
    stop("The selected hypothesis family contains duplicate features.")
  index <- match(predefined_genes$feature_id,all_results$feature_id)
  if (anyNA(index))
    stop("Some preselected features are absent from the fitted model results.")
  out <- all_results[index,,drop=FALSE]
  out$disease_gene_symbol <- predefined_genes$disease_gene_symbol
  out$association_score <- predefined_genes$association_score
  out$FDR_disease_set <- stats::p.adjust(out$P.Value,method="BH")
  out$nominal_p_0_05 <- !is.na(out$P.Value) & out$P.Value<0.05
  out$gene_set_tests <- nrow(predefined_genes)
  out <- out[order(out$P.Value,na.last=TRUE),,drop=FALSE]
  message("Disease-focused analysis: nominal p < 0.05 is the teaching ",
          "threshold across ",nrow(predefined_genes)," preselected features. ",
          "Disease-set and genome-wide BH FDR are retained for comparison.")
  out
}

# The volcano plot highlights genes with nominal p < 0.05.
plot_volcano <- function(results,title="Volcano plot",
                         highlight_column="P.Value") {
  check_packages("ggplot2")
  if (!all(c("logFC","P.Value",highlight_column) %in% names(results)))
    stop("Missing plot columns; check the selected p-value column.")
  d <- results
  d$minus_log10_p <- -log10(pmax(d$P.Value,.Machine$double.xmin))
  d$highlighted <- !is.na(d[[highlight_column]]) &
    d[[highlight_column]]<0.05
  description <- if (highlight_column=="P.Value") {
    "Nominal p < 0.05 (uncorrected)"
  } else if (highlight_column=="FDR_disease_set") {
    "Disease-set BH FDR < 0.05"
  } else if (highlight_column=="adj.P.Val") {
    "Genome-wide BH FDR < 0.05"
  } else {
    paste(highlight_column,"< 0.05")
  }
  ggplot2::ggplot(d,ggplot2::aes(x=logFC,y=minus_log10_p,
                                  color=highlighted))+
    ggplot2::geom_point(alpha=.65)+
    ggplot2::scale_color_manual(values=c("FALSE"="grey65","TRUE"="firebrick"))+
    ggplot2::labs(title=title,
                  x="Expression-scale coefficient per log2 exposure unit",
                  y="-log10(nominal p-value)",color=description)+
    ggplot2::theme_minimal()
}
ot_request <- function(query,variables=list()) {
  check_packages(c("httr","jsonlite"))
  response <- httr::RETRY("POST",
    "https://api.platform.opentargets.org/api/v4/graphql",
    body=list(query=query,variables=variables),encode="json",
    httr::timeout(45),times=3,pause_min=1,pause_cap=10,quiet=TRUE)
  httr::stop_for_status(response)
  p <- jsonlite::fromJSON(httr::content(response,as="text",encoding="UTF-8"),simplifyVector=FALSE)
  if (length(p$errors)) stop("Open Targets API error: ",paste(vapply(p$errors,function(e)e$message,character(1)),collapse="; "))
  p$data
}
ot_search_diseases <- function(search_term) {
  query <- 'query SearchDiseases($term: String!) { search(queryString:$term, entityNames:["disease"]) { hits { id name } } }'
  rows <- ot_request(query,list(term=search_term))$search$hits
  if (is.null(rows) || !length(rows)) return(data.frame(disease_id=character(),disease_name=character()))
  data.frame(disease_id=vapply(rows,function(x)x$id,character(1)),
             disease_name=vapply(rows,function(x)x$name,character(1)))
}
# Define a disease-associated gene set using external Open Targets evidence.
# Score cutoffs are teaching choices, not biological significance thresholds.
# Record the choice before looking at exposure-expression p-values.
ot_cache_path <- function(disease_id,min_score=0.2) {
  make_folders()
  disease <- gsub("[^A-Za-z0-9_-]","_",disease_id)
  if (min_score==0) {
    suffix <- "all"
  } else {
    suffix <- paste0("score_",gsub(".", "p",
                       format(min_score,scientific=FALSE,trim=TRUE),
                       fixed=TRUE))
  }
  file.path("data","gene_sets",paste0("ot_",disease,"_",suffix,".csv"))
}

get_disease_genes <- function(disease_id,allow_online=TRUE,
                              page_size=1000L,min_score=0.2) {
  if (!is.character(disease_id) || length(disease_id)!=1L ||
      is.na(disease_id) || !nzchar(trimws(disease_id)))
    stop("Enter a valid Open Targets disease ID.",call.=FALSE)
  if (length(min_score)!=1L || !is.numeric(min_score) ||
      is.na(min_score) || !is.finite(min_score) ||
      min_score<0 || min_score>1)
    stop("Choose an Open Targets score cutoff between 0 and 1.",
         call.=FALSE)
  if (length(page_size)!=1L || !is.numeric(page_size) ||
      is.na(page_size) || page_size<1 ||
      page_size!=as.integer(page_size) || page_size>3000)
    stop("page_size must be a whole number between 1 and 3000.",
         call.=FALSE)

  file <- ot_cache_path(disease_id,min_score)
  required <- c("disease_id","disease_name","ensembl_id",
                "gene_symbol","association_score")
  validate_genes <- function(genes,source) {
    if (!all(required %in% names(genes)) || !nrow(genes) ||
        anyNA(genes$ensembl_id) ||
        anyDuplicated(genes$ensembl_id) ||
        anyNA(genes$disease_id) ||
        any(genes$disease_id!=disease_id))
      stop("The saved gene list is invalid: ",source,
           ". Check the file and rerun the disease-gene chunk.",
           call.=FALSE)
    if (any(!is.finite(genes$association_score)) ||
        any(genes$association_score<min_score))
      stop("The saved list does not match the selected score cutoff: ",
           source,call.=FALSE)
    genes
  }
  if (file.exists(file)) {
    genes <- validate_genes(read.csv(file,stringsAsFactors=FALSE),file)
    message("Using ",nrow(genes)," saved disease-associated genes ",
            "(score >= ",min_score,"): ",file)
    return(genes)
  }

  # If a complete list is already cached, reuse it without a new download.
  all_file <- ot_cache_path(disease_id,0)
  if (min_score>0 && file.exists(all_file)) {
    genes <- read.csv(all_file,stringsAsFactors=FALSE)
    genes <- genes[is.finite(genes$association_score) &
                     genes$association_score>=min_score,,drop=FALSE]
    if (!nrow(genes))
      stop("No genes in the saved complete list meet score >= ",
           min_score,". Consider a lower score cutoff chosen ",
           "before examining exposure-expression associations.",
           call.=FALSE)
    genes <- validate_genes(genes,all_file)
    write.csv(genes,file,row.names=FALSE)
    message("Saved ",nrow(genes)," genes meeting score >= ",
            min_score," from the existing local gene list.")
    return(genes)
  }
  if (!allow_online)
    stop("No saved gene list for disease ",disease_id,
         " at score >= ",min_score,". Run the gene-list chunk in ",
         "notebook 04 with an internet connection first.",
         call.=FALSE)

  # Open Targets ranks results by association score (highest first).
  # Download only the pages needed to reach the chosen score threshold.
  query <- paste0(
    'query DiseaseTargets($id:String!,$index:Int!,$size:Int!) { ',
    'disease(efoId:$id) { id name ',
    'associatedTargets(page:{index:$index,size:$size}, ',
    'orderByScore:"score") { count rows { ',
    'score target { id approvedSymbol } } } } }'
  )
  all_rows <- list()
  count <- NULL
  disease_name <- NULL
  page_index <- 0L
  pages_read <- 0L
  fetched_rows <- 0L
  repeat {
    result <- tryCatch(
      ot_request(query,list(id=disease_id,index=page_index,
                            size=as.integer(page_size)))$disease,
      error=function(e)
        stop("Open Targets retrieval failed. Check your internet ",
             "connection and try again. Details: ",
             conditionMessage(e),call.=FALSE)
    )
    if (is.null(result))
      stop("Open Targets did not find disease ",disease_id,
           ". Check the selected disease ID.",call.=FALSE)
    response <- result$associatedTargets
    if (is.null(response) || is.null(response$count))
      stop("Open Targets returned no total association count.",
           call.=FALSE)
    reported_total <- as.integer(response$count)
    if (is.na(reported_total) || reported_total<0L)
      stop("Open Targets returned an invalid count.",call.=FALSE)
    if (is.null(count)) {
      count <- reported_total
      disease_name <- as.character(result$name)
      message("Open Targets reports ",count," association rows for ",
              disease_name,". Selecting score >= ",min_score,".")
    } else if (count!=reported_total) {
      stop("Open Targets changed its association count during this ",
           "download. Retry the gene-list chunk.",call.=FALSE)
    }
    if (count==0L)
      stop("Open Targets found no disease-associated genes.",call.=FALSE)

    rows <- response$rows
    if (is.null(rows) || !length(rows))
      stop("Open Targets returned an empty page before the expected ",
           "results were reached.",call.=FALSE)
    fetched_rows <- fetched_rows+length(rows)
    pages_read <- pages_read+1L
    if (fetched_rows>count)
      stop("Open Targets returned more association rows than its ",
           "reported count.",call.=FALSE)
    scores <- vapply(rows,function(x) {
      value <- x$score
      if (is.null(value) || length(value)!=1L || is.na(value))
        return(NA_real_)
      as.numeric(value)
    },numeric(1))
    chosen <- which(is.finite(scores) & scores>=min_score)
    if (length(chosen)) all_rows <- c(all_rows,rows[chosen])
    message("Read ",fetched_rows," of ",count,
            " association rows; ",length(all_rows),
            " rows meet the evidence-score cutoff.")
    if (fetched_rows==count) break
    # A score below the threshold means following lower-ranked pages
    # do not qualify; this avoids downloading thousands of weak hits.
    if (min_score>0 && any(is.finite(scores) & scores<min_score))
      break
    page_index <- page_index+1L
  }

  if (!length(all_rows))
    stop("No genes met the Open Targets score cutoff of ",min_score,
         ". Choose a lower cutoff based on disease evidence BEFORE ",
         "running the exposure-expression models.",call.=FALSE)
  valid <- Filter(function(x)
    !is.null(x$target) && !is.null(x$target$id) &&
    length(x$target$id)==1L && !is.na(x$target$id) &&
    nzchar(x$target$id),all_rows)
  if (!length(valid))
    stop("No retrieved associations had usable gene identifiers.",
         call.=FALSE)

  genes <- data.frame(
    disease_id=rep(disease_id,length(valid)),
    disease_name=rep(disease_name,length(valid)),
    ensembl_id=vapply(valid,function(x)
      as.character(x$target$id),character(1)),
    gene_symbol=vapply(valid,function(x) {
      value <- x$target$approvedSymbol
      if (is.null(value) || length(value)!=1L || is.na(value) ||
          !nzchar(value)) return(NA_character_)
      as.character(value)
    },character(1)),
    association_score=vapply(valid,function(x)
      as.numeric(x$score),numeric(1)),
    retrieved_on=rep(as.character(Sys.Date()),length(valid)),
    stringsAsFactors=FALSE
  )
  # Pagination can repeat an identifier at a page boundary. Keep the
  # strongest association score and report that duplication occurred.
  genes <- genes[order(-genes$association_score,
                       genes$ensembl_id,na.last=TRUE),,drop=FALSE]
  n_duplicated <- sum(duplicated(genes$ensembl_id))
  genes <- genes[!duplicated(genes$ensembl_id),,drop=FALSE]
  if (n_duplicated)
    warning(n_duplicated," repeated gene identifiers were returned by ",
            "Open Targets pages and removed. Page-boundary duplicates ",
            "can mean other genes were skipped by the API's ordering. ",
            "This selected list may not contain every eligible gene; ",
            "interpret it as a teaching gene set, not a complete ",
            "reference catalogue.",call.=FALSE)
  if (!nrow(genes))
    stop("No valid disease-associated genes to save.",call.=FALSE)

  write.csv(genes,file,row.names=FALSE)
  message("Saved ",nrow(genes)," distinct disease-associated genes ",
          "with score >= ",min_score," in ",file,".")
  message("Your expression analysis will use the genes measured ",
          "in the transcriptomics data.")
  genes
}
