# app.R
# WGCNA interactive Shiny app with DESeq2, fast QC plots, WGCNA, and GO/KEGG enrichment
# Author: LT (corrected & hardened)

# =============================
# Packages
# =============================
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinycssloaders)
  library(DT)
  library(tidyverse)
  library(magrittr)
  library(SummarizedExperiment)
  library(DESeq2)
  library(edgeR)
  library(ComplexHeatmap)
  library(DEGreport)
  library(IHW)
  library(WGCNA)
  library(BioNERO)
  library(matrixStats)
  library(grid)
  library(clusterProfiler)
  library(enrichplot)
})

options(shiny.maxRequestSize = 500*1024^2) # allow big uploads
allowWGCNAThreads()

# =============================
# Helpers
# =============================
read_delim_guess <- function(path){
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("rds")) return(readRDS(path))
  if (tolower(ext) %in% c("rdata","rda")) {
    e <- new.env(); load(path, envir=e); return(as.list(e))
  }
  if (tolower(ext) %in% c("csv")) return(read.csv(path, check.names = FALSE, row.names = 1))
  if (tolower(ext) %in% c("tsv","txt")) return(read.delim(path, check.names = FALSE, row.names = 1))
  validate(need(FALSE, paste0("Unsupported file type: ", ext)))
}

make_se <- function(counts, coldata){
  SummarizedExperiment(assays = list(counts = as.matrix(counts)), colData = coldata)
}

safe_plot <- function(expr){
  tryCatch({ expr }, error=function(e){ plot.new(); title(main=paste("Plot error:", e$message)) })
}

pick_power_auto <- function(sft, r2_cut=0.9){
  fi <- sft$fitIndices
  target <- -sign(fi[,3]) * fi[,2]
  ok <- which(target >= r2_cut)
  if (length(ok)>0) return(fi[min(ok),1])
  return(fi[which.max(target),1])
}

sanitize_expr <- function(mat){
  rgood <- apply(mat, 1, function(x) all(is.finite(x)))
  cgood <- apply(mat, 2, function(x) all(is.finite(x)))
  mat[rgood, cgood, drop=FALSE]
}

strip_ensembl_version <- function(x){ sub("\\.\\d+$", "", x) }

make_group_colors <- function(x){
  x <- as.factor(x); levs <- levels(x)
  cols <- grDevices::hcl(seq(15, 375, length.out = length(levs) + 1)[1:length(levs)], 100, 65)
  stats::setNames(cols, levs)
}

guess_species_from_ids <- function(ids){
  ids <- unique(ids); ids <- ids[!is.na(ids)]; n <- length(ids)
  if (n == 0) return(NA_character_)
  sample_ids <- ids[seq_len(min(n, 1000))]
  frac_hsa_ens <- mean(grepl("^ENSG[0-9]+", sample_ids))
  frac_mmu_ens <- mean(grepl("^ENSMUSG[0-9]+", sample_ids))
  if (isTRUE(frac_hsa_ens > 0.5)) return("hsa")
  if (isTRUE(frac_mmu_ens > 0.5)) return("mmu")
  NA_character_
}

org_map <- function(org_choice){
  # Returns list(OrgDb_string, kegg_code)
  switch(org_choice,
         "Human (Homo sapiens)" = list(OrgDb="org.Hs.eg.db", code="hsa"),
         "Mouse (Mus musculus)" = list(OrgDb="org.Mm.eg.db", code="mmu"),
         list(OrgDb="org.Hs.eg.db", code="hsa")
  )
}

# --- New: ID helpers for robust mapping ---
guess_id_type <- function(ids){
  ids <- unique(ids); ids <- ids[!is.na(ids)]; n <- length(ids)
  if (n == 0) return("SYMBOL")
  sample_ids <- ids[seq_len(min(n, 1000))]
  frac_ensembl <- mean(grepl("^ENS(MUS)?G[0-9]+", sample_ids))
  frac_entrez  <- mean(grepl("^[0-9]+$", sample_ids))
  if (isTRUE(frac_ensembl > 0.5)) return("ENSEMBL")
  if (isTRUE(frac_entrez  > 0.5)) return("ENTREZID")
  return("SYMBOL")
}

robust_map_to_entrez <- function(ids, id_type, OrgDb){
  ids <- unique(ids); ids <- ids[!is.na(ids)]
  strip_ver <- function(x) sub("\\.\\d+$", "", x)
  candidates <- unique(c(id_type, guess_id_type(ids), if(id_type!="SYMBOL") "SYMBOL", "ALIAS", if(id_type!="ENSEMBL") "ENSEMBL", "ENTREZID"))
  best <- list(df=NULL, coverage=0, fromType=NULL)
  for (ft in candidates){
    ids_ft <- if (ft == "ENSEMBL") strip_ver(ids) else ids
    m <- try(clusterProfiler::bitr(ids_ft, fromType = ft, toType = "ENTREZID", OrgDb = OrgDb), silent = TRUE)
    if (!inherits(m, "try-error") && !is.null(m) && nrow(m) > 0) {
      if (!ft %in% colnames(m)) next
      m <- m[!duplicated(m[[ft]]), , drop = FALSE]
      cov <- mean(ids_ft %in% m[[ft]])
      if (is.finite(cov) && cov > best$coverage) best <- list(df=m, coverage=cov, fromType=ft)
    }
  }
  best
}

# =============================
# UI
# =============================
ui <- navbarPage(
  title = "WGCNA Network Explorer",
  theme = shinytheme("flatly"),
  
  tabPanel("Data",
           sidebarLayout(
             sidebarPanel(width=4,
                          h4("Upload data"),
                          fileInput("counts_file", "Counts matrix (genes x samples)", accept=c('.csv','.tsv','.txt','.rds','.rdata','.rda')),
                          fileInput("meta_file", "Metadata (samples x variables)", accept=c('.csv','.tsv','.txt','.rds','.rdata','.rda')),
                          textInput("group_col", "Group column in metadata", value = "group"),
                          checkboxInput("has_gene_col", "First column are gene IDs (if not set as rownames)", value = FALSE),
                          hr(),
                          actionButton("load_example", "Load example (GSE54456.RData)", class = "btn-primary"),
                          hr(),
                          helpText("Counts must be raw integer counts. Sample names must match metadata rownames.")
             ),
             mainPanel(
               h4("Preview"),
               fluidRow(
                 column(6, h5("Counts"), withSpinner(DTOutput("counts_head"), type=6)),
                 column(6, h5("Metadata"), withSpinner(DTOutput("meta_head"), type=6))
               )
             )
           )
  ),
  
  tabPanel("QC",
           fluidRow(
             column(6, withSpinner(plotOutput("libsize_bar"), type=6)),
             column(6, withSpinner(plotOutput("heatmap_top500"), type=6))
           ),
           hr(),
           fluidRow(
             column(12, withSpinner(plotOutput("pca_plot"), type=6))
           )
  ),
  
  tabPanel("DESeq2",
           sidebarLayout(
             sidebarPanel(
               uiOutput("contrast_ui"),
               sliderInput("alpha", "FDR (alpha)", min=0.001, max=0.2, value=0.05, step=0.001),
               sliderInput("lfc", "Log2 FC threshold", min=0, max=3, value=1),
               checkboxInput("use_ihw", "Use IHW if feasible", value = TRUE),
               checkboxInput("use_shrink", "Shrink LFC (apeglm)", value = TRUE),
               actionButton("run_deseq", "Run DESeq2", class="btn-primary")
             ),
             mainPanel(
               withSpinner(DTOutput("res_table"), type=6)
             )
           )
  ),
  
  tabPanel("Variance filter",
           sidebarLayout(
             sidebarPanel(
               sliderInput("var_q", "Row variance quantile", min=0.5, max=0.99, value=0.95, step=0.01),
               radioButtons("plot_mode", "Plot mode", choices = c("Auto (fast)"="auto", "Fast box (precomputed)"="box", "Violin (subsample)"="violin"), selected = "auto"),
               numericInput("violin_n", "Violin gene subsample (per app)", value = 3000, min = 200, max = 20000, step = 100),
               actionButton("apply_var_filter", "Apply filter & VST", class="btn-primary")
             ),
             mainPanel(
               withSpinner(plotOutput("violin_norm"), type=6),
               hr(),
               withSpinner(DTOutput("norm_dims"), type=6)
             )
           )
  ),
  
  tabPanel("WGCNA",
           sidebarLayout(
             sidebarPanel(
               selectInput("net_type","Network type", choices = c("signed","signed hybrid","unsigned"), selected = "signed"),
               selectInput("cor_type","Correlation", choices = c("pearson","bicor"), selected = "pearson"),
               selectInput("tom_type","TOM type", choices = c("signed","unsigned"), selected = "signed"),
               checkboxInput("auto_power", "Auto-pick power to reach R^2 target", value = TRUE),
               sliderInput("r2_target", "Scale-free R^2 target", min=0.70, max=0.95, value=0.90, step=0.01),
               numericInput("power_max", "Max power to test", value = 20, min=6, max=50),
               numericInput("minModuleSize", "minModuleSize", 20, min=10),
               sliderInput("deepSplit", "deepSplit", min=0, max=4, value=2, step=1),
               checkboxInput("pamStage", "Use PAM stage", value = TRUE),
               numericInput("mergeCutHeight", "mergeCutHeight", 0.25, min=0, max=1, step=0.01),
               numericInput("maxBlockSize", "maxBlockSize", 4000, min=1000, max=50000, step=500),
               actionButton("run_sft", "Pick soft threshold", class="btn-secondary"),
               uiOutput("power_pick_ui"),
               actionButton("run_wgcna", "Run blockwiseModules", class="btn-primary")
             ),
             mainPanel(
               fluidRow(
                 column(6, withSpinner(plotOutput("sft_scale"), type=6)),
                 column(6, withSpinner(plotOutput("sft_conn"), type=6))
               ),
               hr(),
               withSpinner(plotOutput("dendro_colors"), type=6)
             )
           )
  ),
  
  tabPanel("Modules",
           sidebarLayout(
             sidebarPanel(
               uiOutput("module_select_ui")
             ),
             mainPanel(
               withSpinner(plotOutput("module_trait_heat"), type=6),
               hr(),
               withSpinner(plotOutput("module_lines"), type=6)
             )
           )
  ),
  
  tabPanel("Enrichment (GO/KEGG)",
           sidebarLayout(
             sidebarPanel(width = 4,
                          selectInput("enr_org", "Organism", choices = c("Human (Homo sapiens)", "Mouse (Mus musculus)"), selected = "Human (Homo sapiens)"),
                          selectInput("id_type", "Gene ID type", choices = c("SYMBOL","ENSEMBL","ENTREZID"), selected = "ENSEMBL"),
                          selectInput("go_ont", "GO Ontology", choices = c("BP","CC","MF"), selected = "BP"),
                          selectInput("enr_modules", "Modules", choices = NULL, multiple = TRUE),
                          radioButtons("universe_mode", "Background (universe)", choices = c("All network genes"="network", "All filtered genes"="filtered"), selected = "network"),
                          sliderInput("p_cut", "p-value cutoff", min=0.001, max=0.2, value=0.01, step=0.001),
                          sliderInput("q_cut", "q-value cutoff (FDR)", min=0.001, max=0.25, value=0.05, step=0.001),
                          numericInput("minGS", "Min gene set size", value = 10, min = 5, step = 1),
                          numericInput("maxGS", "Max gene set size", value = 500, min = 50, step = 50),
                          actionButton("run_enrich", "Run enrichment", class = "btn-primary"),
                          hr(),
                          actionButton("audit_ids", "Audit IDs"),
                          verbatimTextOutput("id_audit")
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("GO results", withSpinner(DTOutput("go_table"), type=6), hr(), plotOutput("go_dot")),
                 tabPanel("KEGG results", withSpinner(DTOutput("kegg_table"), type=6), hr(), plotOutput("kegg_dot"))
               ),
               hr(),
               downloadButton("dl_go", "Download GO CSV"),
               downloadButton("dl_kegg", "Download KEGG CSV")
             )
           )
  ),
  
  tabPanel("Network / TOM",
           sidebarLayout(
             sidebarPanel(
               sliderInput("tom_thresh", "TOM edge threshold", min=0, max=1, value=0.2, step=0.01),
               actionButton("build_tom", "Build TOM for selected modules", class="btn-primary"),
               downloadButton("dl_edges", "Download edge list")
             ),
             mainPanel(
               withSpinner(DTOutput("edges_tbl"), type=6)
             )
           )
  ),
  
  tabPanel("Downloads",
           downloadButton("dl_se", "SummarizedExperiment (RDS)"),
           downloadButton("dl_modules", "Module assignment (CSV)"),
           downloadButton("dl_me", "Module eigengenes (CSV)")
  ),
  
  tabPanel("About",
           fluidRow(
             column(8, offset=2,
                    h3("About this app"),
                    p("Interactive WGCNA app built from your R Markdown workflow. Steps: data upload, QC, DESeq2, variance filtering, soft-threshold selection, module detection, trait association, TOM edge export, and GO/KEGG enrichment per module."),
                    tags$ul(
                      tags$li("Inputs: raw counts matrix and metadata with a 'group' column (configurable)."),
                      tags$li("Uses DESeq2 VST for WGCNA input; variance filter by row variance quantile."),
                      tags$li("Soft-threshold picked via pickSoftThreshold; auto-picks power if enabled."),
                      tags$li("Modules via blockwiseModules with robust pre-checks (goodSamplesGenes)."),
                      tags$li("GO/KEGG enrichment with clusterProfiler, flexible IDs and organism settings.")
                    ),
                    p("Packages: WGCNA, DESeq2, BioNERO, IHW, ComplexHeatmap, DEGreport, edgeR, tidyverse, clusterProfiler, enrichplot.")
             )
           )
  )
)

# =============================
# Server
# =============================
server <- function(input, output, session){
  
  # -----------------
  # Data loading
  # -----------------
  rv <- reactiveValues(counts=NULL, meta=NULL, dds=NULL, vsd=NULL, expr_norm=NULL,
                       sft=NULL, power_vec=NULL, picked_power=NULL, net=NULL,
                       module_df=NULL, MEs0=NULL, module_order=NULL,
                       go_df=NULL, kegg_df=NULL)
  
  observeEvent(input$load_example, {
    if (!file.exists("GSE54456.RData")) {
      showNotification("GSE54456.RData not found in the working directory.", type = "error");
      return(NULL)
    }
    e <- new.env(); load("GSE54456.RData", envir = e)
    counts <- e$counts; metadat <- e$metadat
    if (is.null(counts) || is.null(metadat)) {
      showNotification("RData must contain objects named 'counts' (genes x samples) and 'metadat' (sample metadata)", type = "error");
      return(NULL)
    }
    counts <- as.matrix(counts)
    metadat <- as.data.frame(metadat)
    
    # Infer sample identifiers robustly
    sample_id <- NULL
    if ("name" %in% names(metadat))        sample_id <- metadat$name
    else if ("ind" %in% names(metadat))    sample_id <- metadat$ind
    else if (!is.null(rownames(metadat)))   sample_id <- rownames(metadat)
    else if ("SampleName" %in% names(metadat)) sample_id <- metadat$SampleName
    
    if (is.null(sample_id)) {
      showNotification("Couldn't infer sample IDs in 'metadat'. Please add a column 'name' or 'ind' with sample IDs matching counts column names.", type = "error");
      return(NULL)
    }
    
    # Harmonize and align to counts
    metadat$sample_id <- make.names(sample_id, unique = TRUE)
    colnames(counts) <- make.names(colnames(counts), unique = TRUE)
    common <- intersect(colnames(counts), metadat$sample_id)
    
    if (length(common) < 4) {
      showNotification(paste0("Found only ", length(common), " overlapping samples between counts and metadata. Check IDs."), type = "error");
      return(NULL)
    }
    
    counts <- counts[, common, drop = FALSE]
    metadat1 <- metadat[match(common, metadat$sample_id), , drop = FALSE]
    rownames(metadat1) <- metadat1$sample_id
    
    # Create a 'group' column if missing (fallback)
    if (!"group" %in% names(metadat1)) {
      if ("tissue_type" %in% names(metadat1)) {
        metadat1$group <- gsub("-.*", "", metadat1$tissue_type) %>% gsub("[.].*", "", .)
      } else {
        metadat1$group <- factor("group1")
      }
    }
    
    rv$counts <- counts
    rv$meta   <- metadat1
    showNotification(paste0("Loaded example with ", nrow(counts), " genes and ", ncol(counts), " samples."), type = "message")
  })
  
  observeEvent(input$counts_file, {
    obj <- read_delim_guess(input$counts_file$datapath)
    if (is.list(obj) && !is.data.frame(obj)) {
      if (!is.null(obj$counts)) obj <- obj$counts else validate(need(FALSE, "RData must contain object 'counts'"))
    }
    rv$counts <- as.matrix(obj)
  })
  
  observeEvent(input$meta_file, {
    obj <- read_delim_guess(input$meta_file$datapath)
    if (is.list(obj) && !is.data.frame(obj)) {
      if (!is.null(obj$metadat)) obj <- obj$metadat else validate(need(FALSE, "RData must contain object 'metadat' (metadata)"))
    }
    df <- as.data.frame(obj)
    if (!"name" %in% colnames(df) && !is.null(rownames(df))) df$name <- rownames(df)
    rownames(df) <- df$name
    rv$meta <- df
  })
  
  output$counts_head <- renderDT({ req(rv$counts); datatable(head(rv$counts), options=list(scrollX=TRUE)) })
  output$meta_head    <- renderDT({ req(rv$meta); datatable(head(rv$meta),    options=list(scrollX=TRUE)) })
  
  # -----------------
  # QC + DESeq2
  # -----------------
  dds_ready <- reactive({
    req(rv$counts, rv$meta)
    meta <- rv$meta
    grpcol <- req(input$group_col)
    validate(need(grpcol %in% colnames(meta), paste0("Column '", grpcol, "' not found in metadata")))
    common <- intersect(colnames(rv$counts), rownames(meta))
    validate(need(length(common)>=4, "Need >=4 overlapping samples between counts and metadata"))
    counts <- round(rv$counts[, common])
    coldata <- meta[common, , drop=FALSE]
    DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = as.formula(paste0("~", grpcol)))
  })
  
  observeEvent(input$run_deseq, {
    dds <- dds_ready()
    keep <- rowSums(edgeR::cpm(counts(dds)) > 0.5) >= 2
    dds <- dds[keep,]
    rv$dds <- DESeq(dds)
  })
  
  res_reactive <- eventReactive(input$run_deseq, {
    dds <- req(rv$dds)
    validate(need(input$grpA != input$grpB, "Choose two different groups."))
    
    contrast <- c(input$group_col, input$grpA, input$grpB)
    use_threshold <- isTRUE(input$lfc > 0)
    
    base_args <- list(object = dds, contrast = contrast, alpha = input$alpha, independentFiltering = TRUE)
    if (use_threshold) { base_args$lfcThreshold <- input$lfc; base_args$altHypothesis <- "greaterAbs" }
    
    res_base <- do.call(DESeq2::results, base_args)
    
    do_ihw <- isTRUE(input$use_ihw)
    can_ihw <- do_ihw && requireNamespace("IHW", quietly = TRUE) &&
      sum(!is.na(res_base$pvalue)) >= 1000 && sum(!is.na(res_base$baseMean)) >= 1000
    
    res_final <- res_base
    if (can_ihw) {
      ihw_args <- base_args; ihw_args$filterFun <- IHW::ihw
      res_try <- try(do.call(DESeq2::results, ihw_args), silent = TRUE)
      if (!inherits(res_try, "try-error")) res_final <- res_try else showNotification("IHW failed — using standard DESeq2 filtering.", type="warning")
    } else if (do_ihw) {
      msg <- if (!requireNamespace("IHW", quietly = TRUE)) "IHW not available" else "Dataset small or many NAs — skipping IHW"
      showNotification(paste0(msg, "."), type = "message")
    }
    
    if (isTRUE(input$use_shrink)) {
      if (requireNamespace("apeglm", quietly = TRUE)) {
        res_shrunk <- try(lfcShrink(dds, contrast = contrast, res = res_final, type = "apeglm"), silent = TRUE)
        if (!inherits(res_shrunk, "try-error")) return(res_shrunk) else showNotification("apeglm shrinkage failed — showing unshrunk results.", type="warning")
      } else showNotification("Package 'apeglm' not installed — showing unshrunk results.", type="message")
    }
    
    res_final
  })
  
  output$libsize_bar <- renderPlot({
    dds <- dds_ready(); sizes <- colSums(counts(dds))
    safe_plot({ barplot(sizes, las=2, main="Library sizes", ylab="Reads", col="gray") })
  })
  
  output$heatmap_top500 <- renderPlot({
    dds <- dds_ready(); mat <- log2(counts(dds) + 0.5); if (nrow(mat)>500) mat <- mat[1:500,]
    d <- as.data.frame(colData(dds))
    ann_df <- NULL; ann_col <- NULL
    if (!is.null(input$group_col) && input$group_col %in% names(d)){
      ann_df <- data.frame(group = d[[input$group_col]])
      ann_col <- list(group = make_group_colors(ann_df$group))
    }
    th <- if (is.null(ann_df)) NULL else HeatmapAnnotation(df = ann_df, col = ann_col)
    safe_plot({ Heatmap(mat, show_row_names = FALSE, top_annotation = th) })
  })
  
  output$pca_plot <- renderPlot({ dds <- dds_ready(); vsd <- vst(dds, blind=TRUE); safe_plot({ plotPCA(vsd, intgroup = input$group_col, ntop = 1000) }) })
  
  output$contrast_ui <- renderUI({
    meta <- req(rv$meta); grp <- req(input$group_col); validate(need(grp %in% colnames(meta), "Group column missing"))
    levs <- sort(unique(meta[[grp]]))
    fluidRow(column(12, selectInput("grpA", "Group A (treatment)", choices = levs), selectInput("grpB", "Group B (reference)", choices = levs, selected = if (length(levs)>=2) levs[2] else NULL)))
  })
  
  output$res_table <- renderDT({ res <- req(res_reactive()); df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id"); datatable(df, options=list(scrollX=TRUE), filter = "top") })
  
  # -----------------
  # Variance filtering & VST for WGCNA
  # -----------------
  observeEvent(input$apply_var_filter, {
    dds <- dds_ready(); vsd <- vst(dds, blind=TRUE); mat <- assay(vsd); rv$vsd <- vsd
    rv$expr_norm <- { rvv <- matrixStats::rowVars(mat); thr <- as.numeric(stats::quantile(rvv, input$var_q)); mat[rvv > thr, ] }
  })
  
  output$violin_norm <- renderPlot({
    mat <- req(rv$expr_norm); grpcol <- req(input$group_col); meta <- req(rv$meta)
    validate(need(grpcol %in% colnames(meta), paste0("Column '", grpcol, "' not found in metadata")))
    sample_names <- colnames(mat); group_vec <- as.character(meta[sample_names, grpcol]); names(group_vec) <- sample_names
    mode <- input$plot_mode; if (identical(mode, "auto")) { big <- nrow(mat)*ncol(mat) >= 2e6 || nrow(mat) >= 10000; mode <- if (big) "box" else "violin" }
    safe_plot({
      if (mode=="box"){
        q1 <- matrixStats::colQuantiles(mat, probs=0.25, na.rm=TRUE); med <- matrixStats::colMedians(mat, na.rm=TRUE); q3 <- matrixStats::colQuantiles(mat, probs=0.75, na.rm=TRUE)
        iqr <- q3-q1; lower <- pmax(matrixStats::colMins(mat, na.rm=TRUE), q1-1.5*iqr); upper <- pmin(matrixStats::colMaxs(mat, na.rm=TRUE), q3+1.5*iqr)
        sumDT <- tibble::tibble(name=sample_names, group=group_vec, ymin=lower, lower=q1, middle=med, upper=q3, ymax=upper)
        ggplot(sumDT, aes(x=name, ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax)) + geom_boxplot(stat="identity") +
          facet_grid(cols = vars(group), drop = TRUE, scales = "free_x") + coord_cartesian(ylim=c(0,NA)) + theme_bw() +
          theme(axis.text.x = element_text(angle=90), panel.spacing.x = unit(6, "pt")) + labs(title="Fast sample distributions (box)", x="Samples", y="VST expression")
      } else {
        set.seed(1); k <- min(input$violin_n, nrow(mat)); idx <- sample.int(nrow(mat), size = k)
        df <- as.data.frame(mat[idx,,drop=FALSE]) %>% tibble::rownames_to_column("Gene_id") %>% tidyr::pivot_longer(-Gene_id, names_to="name", values_to="value")
        df$group <- group_vec[df$name]
        ggplot(df, aes(x=name, y=value)) + geom_violin(trim=TRUE, adjust=1.5, scale="width") + stat_summary(fun=median, geom="point", size=0.6, alpha=0.8) +
          facet_grid(cols = vars(group), drop = TRUE, scales = "free_x") + coord_cartesian(ylim=c(0,NA)) + theme_bw() + theme(axis.text.x = element_text(angle=90)) +
          labs(title=paste0("Normalized expression (violin, n=", k, " genes)"), x="Samples", y="VST expression")
      }
    })
  })
  
  output$norm_dims <- renderDT({ mat <- req(rv$expr_norm); mode <- input$plot_mode; if (identical(mode, "auto")) { big <- nrow(mat)*ncol(mat) >= 2e6 || nrow(mat) >= 10000; mode <- if (big) "box" else "violin" }
  datatable(data.frame(n_genes=nrow(mat), n_samples=ncol(mat), plot_mode_used=mode, violin_subsample=ifelse(mode=="violin", input$violin_n, NA)), options=list(dom='t')) })
  
  # -----------------
  # WGCNA: soft-threshold and modules
  # -----------------
  output$power_pick_ui <- renderUI({ numericInput("picked_power", "Picked power", value = ifelse(is.null(rv$picked_power), 9, rv$picked_power), min=1, max=50) })
  
  observeEvent(input$run_sft, {
    mat <- req(rv$expr_norm); input_mat <- t(sanitize_expr(mat))
    gsg <- goodSamplesGenes(input_mat, verbose = 3); if (!gsg$allOK) { input_mat <- input_mat[gsg$goodSamples, gsg$goodGenes]; showNotification(paste0("Removed ", sum(!gsg$goodSamples), " samples and ", sum(!gsg$goodGenes), " genes after QC."), type = "message") }
    powers <- c(1:10, seq(12, input$power_max, by=2))
    corFnc <- if (input$cor_type == "bicor") "bicor" else "cor"
    corOpt <- if (input$cor_type == "bicor") list(maxPOutliers = 0.1) else list(use = 'p')
    sft <- pickSoftThreshold(input_mat, powerVector = powers, networkType = input$net_type, corFnc = corFnc, corOptions = corOpt, verbose = 5)
    rv$sft <- sft; rv$power_vec <- powers
    if (isTRUE(input$auto_power)) rv$picked_power <- pick_power_auto(sft, input$r2_target)
  })
  
  output$sft_scale <- renderPlot({ req(rv$sft); fi <- rv$sft$fitIndices; powers <- rv$power_vec; cex1 <- .9
  safe_plot({ plot(fi[,1], -sign(fi[,3]) * fi[,2], xlab="Soft Threshold (power)", ylab="Scale-free topology fit (signed R^2)", main="Scale independence"); text(fi[,1], -sign(fi[,3]) * fi[,2], labels=powers, cex=cex1, col="red"); abline(h = input$r2_target, col="red") }) })
  
  output$sft_conn <- renderPlot({ req(rv$sft); fi <- rv$sft$fitIndices; powers <- rv$power_vec; cex1 <- .9
  safe_plot({ plot(fi[,1], fi[,5], xlab="Soft Threshold (power)", ylab="Mean connectivity", type="n", main="Mean connectivity"); text(fi[,1], fi[,5], labels=powers, cex=cex1, col="red") }) })
  
  observeEvent(input$run_wgcna, {
    mat <- req(rv$expr_norm); input_mat <- t(sanitize_expr(mat))
    if (nrow(input_mat) < 10 || ncol(input_mat) < 100) showNotification("Few samples/genes for stable WGCNA. Consider lowering minModuleSize or variance filter.", type = "warning")
    gsg <- goodSamplesGenes(input_mat, verbose = 3); if (!gsg$allOK) { input_mat <- input_mat[gsg$goodSamples, gsg$goodGenes]; showNotification(paste0("Removed ", sum(!gsg$goodSamples), " samples and ", sum(!gsg$goodGenes), " genes after QC."), type = "message") }
    corType <- if (input$cor_type == "bicor") "bicor" else "pearson"
    corOptions <- if (input$cor_type == "bicor") list(maxPOutliers = 0.1) else list(use = 'p', method = 'pearson')
    params <- list(datExpr = input_mat, power = req(input$picked_power), networkType = input$net_type, TOMType = input$tom_type, corType = corType, corOptions = corOptions, deepSplit = input$deepSplit, pamStage = isTRUE(input$pamStage), pamRespectsDendro = FALSE, minModuleSize = input$minModuleSize, maxBlockSize = input$maxBlockSize, reassignThreshold = 0, mergeCutHeight = input$mergeCutHeight, saveTOMs = FALSE, numericLabels = TRUE, verbose = 3)
    attempt <- function(p){ try(do.call(WGCNA::blockwiseModules, p), silent = TRUE) }
    net <- attempt(params)
    no_modules <- function(net){ if (inherits(net, "try-error")) return(TRUE); cols <- net$colors; return(sum(cols > 0) == 0) }
    if (no_modules(net)) {
      params_fallback <- params; params_fallback$deepSplit <- 4; params_fallback$mergeCutHeight <- 0.1; params_fallback$minModuleSize <- max(10, floor(ncol(input_mat)/5))
      if (input$cor_type == "bicor") { params_fallback$corType <- "pearson"; params_fallback$corOptions <- list(use='p', method='pearson') }
      net2 <- attempt(params_fallback)
      if (!no_modules(net2)) { net <- net2; showNotification("Fallback parameters found modules (deepSplit=4, mergeCutHeight=0.1).", type = "warning", duration = 7) } else { showNotification("WGCNA could not detect modules with given settings.", type = "error", duration = 10); return(NULL) }
    }
    rv$net <- net
    mergedColors <- labels2colors(net$colors)
    rv$module_df <- data.frame(gene_id = names(net$colors), colors = mergedColors)
    rownames(rv$module_df) <- rv$module_df$gene_id
    MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes %>% orderMEs(); rv$MEs0 <- MEs0; rv$module_order <- gsub("ME","", colnames(MEs0))
    updateSelectInput(session, "enr_modules", choices = sort(unique(mergedColors)), selected = sort(unique(mergedColors))[1:min(3, length(unique(mergedColors)))])
  })
  
  output$dendro_colors <- renderPlot({ net <- req(rv$net); mergedColors <- labels2colors(net$colors)
  safe_plot({ plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) }) })
  
  # -----------------
  # Module-trait and module visuals
  # -----------------
  output$module_select_ui <- renderUI({ df <- req(rv$module_df); mods <- sort(unique(df$colors)); selectizeInput("modules_of_interest", "Modules of interest", choices = mods, multiple = TRUE, selected = head(mods,3)) })
  
  output$module_trait_heat <- renderPlot({ req(rv$MEs0, rv$module_order); MEs0 <- rv$MEs0; MEs0$sample <- rownames(MEs0)
  mME <- MEs0 %>% tibble::as_tibble() %>% tidyr::pivot_longer(-sample) %>% mutate(name = gsub("ME", "", name), name = factor(name, levels = rv$module_order))
  safe_plot({ ggplot(mME, aes(x=sample, y=name, fill=value)) + geom_tile() + theme_bw() + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limits=c(-1,1)) + theme(axis.text.x = element_text(angle=90)) + labs(title="Module–sample eigengene values", y="Modules", fill="ME") }) })
  
  output$module_lines <- renderPlot({ mat <- req(rv$expr_norm); dfm <- req(rv$module_df); sel <- req(input$modules_of_interest); genes <- dfm %>% dplyr::filter(colors %in% sel) %>% pull(gene_id); validate(need(length(genes)>0, "No genes in selected modules"))
  subexpr <- mat[intersect(rownames(mat), genes), ]; submod_df <- as.data.frame(subexpr) %>% tibble::rownames_to_column("gene_id") %>% tidyr::pivot_longer(-gene_id) %>% mutate(module = dfm[gene_id,]$colors)
  safe_plot({ ggplot(submod_df, aes(x=name, y=value, group=gene_id)) + geom_line(aes(color=module), alpha=0.2) + theme_bw() + theme(axis.text.x = element_text(angle=90)) + facet_grid(rows = vars(module)) + labs(x="Sample", y="VST expression") }) })
  
  # -----------------
  # Enrichment (GO/KEGG)
  # -----------------
  observe({
    df <- rv$module_df
    if (!is.null(df)) updateSelectInput(session, "enr_modules", choices = sort(unique(df$colors)))
  })
  
  output$id_audit <- renderText({
    dfm <- rv$module_df; mat <- rv$expr_norm
    validate(need(!is.null(dfm) && !is.null(mat), "Run WGCNA first."))
    universe_ids <- if (identical(input$universe_mode, "network")) dfm$gene_id else rownames(mat)
    id_type_guess <- guess_id_type(universe_ids)
    species_guess <- guess_species_from_ids(universe_ids)
    org <- org_map(input$enr_org)
    cov <- "(OrgDb not installed)"
    if (requireNamespace(org$OrgDb, quietly = TRUE)) {
      suppressPackageStartupMessages(require(org$OrgDb, character.only = TRUE))
      OrgDb_obj <- get(org$OrgDb)
      test_map <- robust_map_to_entrez(universe_ids, input$id_type, OrgDb_obj)
      cov <- if (!is.null(test_map$df)) sprintf("%.1f%%", 100*test_map$coverage) else "0%"
    }
    sample_show <- paste(utils::head(unique(universe_ids), 5), collapse=", ")
    paste0("ID audit:\n  • Sample IDs: ", sample_show,
           "\n  • Guessed ID type: ", id_type_guess,
           "\n  • Guessed species: ", ifelse(is.na(species_guess), "unknown", species_guess),
           "\n  • Mapping coverage (current settings): ", cov)
  })
  
  observeEvent(input$audit_ids, {
    dfm <- rv$module_df; mat <- rv$expr_norm; req(dfm, mat)
    ids <- if (identical(input$universe_mode, "network")) dfm$gene_id else rownames(mat)
    sp <- guess_species_from_ids(ids)
    if (!is.na(sp)) {
      if (sp == "hsa" && input$enr_org != "Human (Homo sapiens)") {
        updateSelectInput(session, "enr_org", selected = "Human (Homo sapiens)")
        showNotification("Switched organism to Human based on Ensembl prefixes.", type = "message")
      }
      if (sp == "mmu" && input$enr_org != "Mouse (Mus musculus)") {
        updateSelectInput(session, "enr_org", selected = "Mouse (Mus musculus)")
        showNotification("Switched organism to Mouse based on Ensembl prefixes.", type = "message")
      }
    } else {
      showNotification("Could not infer species from IDs (numeric or mixed). Ensure organism selection is correct.", type = "warning")
    }
  })
  
  enrich_run <- eventReactive(input$run_enrich, {
    dfm <- req(rv$module_df); mat <- req(rv$expr_norm)
    mods <- req(input$enr_modules)
    if (length(mods) == 0) { showNotification("Select at least one module.", type="warning"); return(NULL) }
    
    # Universe
    universe_ids <- if (input$universe_mode == "network") dfm$gene_id else rownames(mat)
    if (input$id_type == "ENSEMBL") universe_ids <- strip_ensembl_version(universe_ids)
    
    # Species sanity check from IDs (only if Ensembl-like)
    sp_guess <- guess_species_from_ids(universe_ids)
    
    org <- org_map(input$enr_org)
    if (!requireNamespace(org$OrgDb, quietly = TRUE)) {
      showNotification(paste0("Please install Bioconductor package '", org$OrgDb, "' for ID mapping."), type = "error", duration = 8)
      return(NULL)
    }
    suppressPackageStartupMessages(require(org$OrgDb, character.only = TRUE))
    OrgDb_obj <- get(org$OrgDb)
    
    if (!is.na(sp_guess) && sp_guess != org$code) {
      showNotification(sprintf("Your IDs look %s (based on Ensembl prefixes), but organism is set to %s. Please switch.",
                               ifelse(sp_guess=="hsa","Human","Mouse"), input$enr_org), type = "error", duration = 8)
      return(NULL)
    }
    
    # Map universe to ENTREZ robustly
    uni_map <- robust_map_to_entrez(universe_ids, input$id_type, OrgDb_obj)
    if (is.null(uni_map$df) || length(unique(uni_map$df$ENTREZID)) == 0) {
      showNotification("Failed to map universe IDs to ENTREZ. Check organism and ID type.", type = "error", duration = 8)
      return(NULL)
    }
    universe_entrez <- unique(uni_map$df$ENTREZID)
    if (length(universe_entrez) < input$minGS) {
      showNotification(sprintf("Universe is small after mapping (n=%d). Results may be unstable.", length(universe_entrez)), type="warning", duration=6)
    }
    showNotification(sprintf("Universe mapped using %s → ENTREZ (%.1f%% of %d).",
                             uni_map$fromType, 100*uni_map$coverage, length(universe_ids)),
                     type="message", duration=6)
    
    go_list <- list(); kegg_list <- list()
    
    for (m in mods) {
      genes_m <- dfm$gene_id[dfm$colors == m]
      if (input$id_type == "ENSEMBL") genes_m <- strip_ensembl_version(genes_m)
      mm <- robust_map_to_entrez(genes_m, input$id_type, OrgDb_obj)
      if (is.null(mm$df)) next
      genes_entrez <- unique(mm$df$ENTREZID)
      
      # GO
      eg <- try(enrichGO(gene = genes_entrez, OrgDb = OrgDb_obj, keyType = "ENTREZID", ont = input$go_ont,
                         universe = universe_entrez, pvalueCutoff = input$p_cut, qvalueCutoff = input$q_cut,
                         minGSSize = input$minGS, maxGSSize = input$maxGS, readable = TRUE), silent = TRUE)
      if (!inherits(eg, "try-error") && !is.null(eg) && nrow(as.data.frame(eg))>0) {
        df_go <- as.data.frame(eg); df_go$Module <- m; go_list[[m]] <- df_go
      }
      
      # KEGG (Entrez IDs)
      ek <- try(enrichKEGG(gene = genes_entrez, organism = org$code, keyType = "ncbi-geneid",
                           universe = universe_entrez, pvalueCutoff = input$p_cut, qvalueCutoff = input$q_cut,
                           minGSSize = input$minGS, maxGSSize = input$maxGS), silent = TRUE)
      if (!inherits(ek, "try-error") && !is.null(ek) && nrow(as.data.frame(ek))>0) {
        df_kegg <- as.data.frame(ek); df_kegg$Module <- m; kegg_list[[m]] <- df_kegg
      }
    }
    
    rv$go_df   <- if (length(go_list))   dplyr::bind_rows(go_list)   else NULL
    rv$kegg_df <- if (length(kegg_list)) dplyr::bind_rows(kegg_list) else NULL
    
    list(go=rv$go_df, kegg=rv$kegg_df)
  })
  
  output$go_table <- renderDT({ x <- enrich_run(); validate(need(!is.null(x$go), "No GO terms found.")); datatable(x$go, options=list(scrollX=TRUE), filter="top") })
  output$kegg_table <- renderDT({ x <- enrich_run(); validate(need(!is.null(x$kegg), "No KEGG pathways found.")); datatable(x$kegg, options=list(scrollX=TRUE), filter="top") })
  
  output$go_dot <- renderPlot({
    x <- enrich_run(); validate(need(!is.null(x$go), "Run enrichment first"))
    df <- x$go %>% dplyr::group_by(Module) %>% dplyr::slice_min(order_by = p.adjust, n = 10) %>% dplyr::ungroup()
    safe_plot({ ggplot(df, aes(x=Count, y=reorder(Description, p.adjust))) +
        geom_point(aes(size=Count, color=-log10(p.adjust))) +
        facet_wrap(~Module, scales="free_y") + theme_bw() + labs(x="Gene count", y=NULL, color="-log10(FDR)") })
  })
  
  output$kegg_dot <- renderPlot({
    x <- enrich_run(); validate(need(!is.null(x$kegg), "Run enrichment first"))
    df <- x$kegg %>% dplyr::group_by(Module) %>% dplyr::slice_min(order_by = p.adjust, n = 10) %>% dplyr::ungroup()
    safe_plot({ ggplot(df, aes(x=Count, y=reorder(Description, p.adjust))) +
        geom_point(aes(size=Count, color=-log10(p.adjust))) +
        facet_wrap(~Module, scales="free_y") + theme_bw() + labs(x="Gene count", y=NULL, color="-log10(FDR)") })
  })
  
  output$dl_go <- downloadHandler(
    filename = function(){ "GO_enrichment.csv" },
    content = function(file){ x <- enrich_run(); validate(need(!is.null(x$go), "Run enrichment first")); readr::write_csv(x$go, file) }
  )
  output$dl_kegg <- downloadHandler(
    filename = function(){ "KEGG_enrichment.csv" },
    content = function(file){ x <- enrich_run(); validate(need(!is.null(x$kegg), "Run enrichment first")); readr::write_csv(x$kegg, file) }
  )
  
  # -----------------
  # Network / TOM
  # -----------------
  edges_reactive <- eventReactive(input$build_tom, {
    mat <- req(rv$expr_norm); dfm <- req(rv$module_df); sel <- req(input$modules_of_interest)
    genes <- dfm %>% dplyr::filter(colors %in% sel) %>% pull(gene_id)
    expr_of_interest <- mat[intersect(rownames(mat), genes), ]
    validate(need(nrow(expr_of_interest) >= 2, "Not enough genes to build TOM."))
    picked <- req(input$picked_power)
    TOM <- TOMsimilarityFromExpr(t(expr_of_interest), power = picked, TOMType = input$tom_type)
    rownames(TOM) <- colnames(TOM) <- rownames(expr_of_interest)
    edge_list <- as.data.frame(TOM) %>% tibble::rownames_to_column("gene1") %>% tidyr::pivot_longer(-gene1, names_to = "gene2", values_to = "correlation") %>% dplyr::filter(gene1 != gene2, correlation >= input$tom_thresh) %>% distinct()
    edge_list$module1 <- dfm[edge_list$gene1,]$colors
    edge_list$module2 <- dfm[edge_list$gene2,]$colors
    edge_list
  })
  
  output$edges_tbl <- renderDT({ datatable(req(edges_reactive()), options=list(scrollX=TRUE), filter = "top") })
  output$dl_edges <- downloadHandler(
    filename = function(){ paste0("edges_thr", input$tom_thresh, ".csv") },
    content = function(file){ readr::write_csv(req(edges_reactive()), file) }
  )
  
  # -----------------
  # Downloads
  # -----------------
  output$dl_se <- downloadHandler(
    filename = function(){ "wgcna_input_se.rds" },
    content = function(file){ se <- make_se(round(counts(dds_ready())), as.data.frame(colData(dds_ready()))); saveRDS(se, file) }
  )
  
  output$dl_modules <- downloadHandler(
    filename = function(){ "module_assignment.csv" },
    content = function(file){ readr::write_csv(req(rv$module_df), file) }
  )
  
  output$dl_me <- downloadHandler(
    filename = function(){ "module_eigengenes.csv" },
    content = function(file){ MEs0 <- req(rv$MEs0); readr::write_csv(tibble::rownames_to_column(as.data.frame(MEs0), "sample"), file) }
  )
}

# For app.R entrypoint
shinyApp(ui, server)
