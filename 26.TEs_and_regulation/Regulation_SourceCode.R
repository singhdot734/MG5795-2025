suppressPackageStartupMessages({
  library(WGCNA)
  library(pheatmap)
  library(dplyr)
})

## Top-level WGCNA heatmap 
make_wgcna_heatmap <- function(MEs, traits, moduleColors, sampN,
                               pdf_file = "ME_Trait_heatmap.pdf",
                               cor_csv = "MEvsTraitCor.csv",
                               fdr_csv = "MEvsTraitFDR.csv",
                               trait_column = "tissue") {
  
  # ---- 1) Build a trait matrix ----
  traits2 <- data.frame(
    lapply(traits, function(x) if (is.character(x)) factor(x) else x),
    check.names = FALSE, row.names = rownames(traits)
  )
  traits2 <- traits2[, intersect(trait_column, names(traits2)), drop = FALSE]
  
  X <- model.matrix(~ 0 + ., data = traits2)
  colnames(X) <- make.names(colnames(X))
  
  is_const <- apply(X, 2, function(v) length(na.omit(unique(v))) == 1)
  traitMat <- X[, !is_const, drop = FALSE]
  
  stopifnot(identical(rownames(MEs), rownames(traitMat)))
  
  # ---- 2) (paste in your module ordering + renaming code here) ----
  moduleColors <- trimws(as.character(moduleColors))
  me_colors    <- trimws(sub("^ME", "", colnames(MEs)))
  color_sizes  <- sort(table(moduleColors), decreasing = TRUE)
  base_colors  <- names(color_sizes)
  present_colors <- intersect(base_colors, unique(me_colors))
  color_to_id <- setNames(seq_along(present_colors), present_colors)
  missing <- setdiff(unique(me_colors), present_colors)
  if (length(missing)) {
    warning("These ME colors did not match the size table: ",
            paste(missing, collapse = ", "),
            ". Assigning new IDs at the end.")
    extra_ids <- setNames(seq_along(missing) + length(present_colors), missing)
    color_to_id <- c(color_to_id, extra_ids)
  }
  me_ids <- paste0("M", unname(color_to_id[me_colors]))
  colnames(MEs) <- me_ids
  
  module_key <- data.frame(
    module_id    = paste0("M", seq_along(color_to_id)),
    module_color = names(color_to_id),
    size         = as.integer(color_sizes[names(color_to_id)]),
    row.names    = NULL
  )
  
  # ---- 3) Correlations + FDR (your existing code) ----
  MEvsTraitCor <- matrix(NA_real_, nrow = ncol(MEs), ncol = ncol(traitMat),
                         dimnames = list(colnames(MEs), colnames(traitMat)))
  
  for (j in seq_len(ncol(traitMat))) {
    y <- traitMat[, j]
    if (mad(y, na.rm = TRUE) == 0) {
      MEvsTraitCor[, j] <- cor(MEs, y, use = "pairwise.complete.obs", method = "pearson")
    } else {
      MEvsTraitCor[, j] <- WGCNA::bicor(MEs, y, maxPOutliers = 0.1,
                                        use = "pairwise.complete.obs")
    }
  }
  
  MEvsTraitP   <- WGCNA::corPvalueStudent(MEvsTraitCor, sampN)
  MEvsTraitFDR <- matrix(
    p.adjust(as.vector(MEvsTraitP), method = "BH"),
    nrow = nrow(MEvsTraitP),
    dimnames = dimnames(MEvsTraitP)
  )
  
  # ---- 4) Heatmap plotting (your pheatmap code, unchanged) ----
  fdr_to_stars <- function(q) ifelse(q < 0.001, "***",
                                     ifelse(q < 0.01,  "**",
                                            ifelse(q < 0.05,  "*", "")))
  stars <- matrix(fdr_to_stars(MEvsTraitFDR), nrow = nrow(MEvsTraitFDR),
                  dimnames = dimnames(MEvsTraitFDR))
  
  num_txt <- matrix(sprintf("%.2f", MEvsTraitCor), nrow = nrow(MEvsTraitCor),
                    dimnames = dimnames(MEvsTraitCor))
  
  disp <- matrix(paste0(num_txt, stars), nrow = nrow(num_txt),
                 dimnames = dimnames(num_txt))
  
  breaks <- seq(-1, 1, length.out = 101)
  
  p <- pheatmap(
    MEvsTraitCor,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(length(breaks)-1),
    breaks = breaks,
    cluster_rows = TRUE, cluster_cols = TRUE,
    display_numbers = disp, number_color = "black",
    fontsize_number = 9,
    main = "Module–Trait correlations (numbers=r, stars=FDR)",
    border_color = NA
  )
  
  pdf(pdf_file, width = 8.5, height = 6.5)
  print(p)
  dev.off()
  
  utils::write.csv(MEvsTraitCor, file = cor_csv)
  utils::write.csv(MEvsTraitFDR, file = fdr_csv)
  
  invisible(list(
    module_key   = module_key,
    MEvsTraitCor = MEvsTraitCor,
    MEvsTraitFDR = MEvsTraitFDR,
    MEs_renamed  = MEs
  ))
}

## Cytoscape and network analysis 
run_cytoscape_section <- function(datExpr, traits, MEs_input, moduleColors,
                                  modules            = "module",
                                  trait_column       = "tissue",
                                  tissue_label       = NULL,
                                  tissue_colors      = NULL,
                                  tom_top_quantile   = 0.9,
                                  tom_fixed_cutoff   = NULL,
                                  connect_cytoscape  = TRUE,
                                  out_dir_prefix     = "cyto",
                                  softPower          = 22,
                                  ltr_id             = NULL,
                                  top_k_partners     = 20,
                                  edge_types         = c("LTR-Gene"),
                                  skip_heatmap       = FALSE,
                                  export_sif         = TRUE,
                                  zip_bundle         = TRUE    
) {
  
  # ---- deps
  stopifnot(requireNamespace("WGCNA", quietly = TRUE))
  if (connect_cytoscape) stopifnot(requireNamespace("RCy3", quietly = TRUE))
  
  # ---- basic checks & alignment
  datExpr   <- as.data.frame(datExpr,   check.names = FALSE)
  traits    <- as.data.frame(traits,    check.names = FALSE)
  MEs_input <- as.data.frame(MEs_input, check.names = FALSE)
  stopifnot(!is.null(rownames(datExpr)), !is.null(rownames(traits)), !is.null(rownames(MEs_input)))
  samps <- Reduce(intersect, list(rownames(datExpr), rownames(traits), rownames(MEs_input)))
  stopifnot(length(samps) >= 3)
  datExpr   <- datExpr[samps, , drop = FALSE]
  traits    <- traits[samps, , drop = FALSE]
  MEs_input <- MEs_input[samps, , drop = FALSE]
  
  # ---- helper: cytoscape ping
  .wait_for_cy <- function(retries = 10, wait = 1) {
    for (i in seq_len(retries)) {
      ok <- tryCatch({ suppressMessages(RCy3::cytoscapePing()); TRUE }, error = function(e) FALSE)
      if (ok) return(TRUE)
      Sys.sleep(wait)
    }
    FALSE
  }
  
  if (is.null(tissue_colors)) {
    tissue_colors <- c(
      "Anther"    = "#7b3294",
      "Ear"       = "#c2a5cf",
      "Embryo"    = "#008837",
      "Endosperm" = "#a6dba0",
      "Leaf"      = "#fdae61",
      "Pollen"    = "#d01c8b",
      "Root"      = "#1f78b4",
      "Shoot"     = "#e6ab02",
      "Tassel"    = "#a6761d"
    )
  }
  
  stopifnot(is.character(names(tissue_colors)))
  
  # ---- normalize module selector
  modules    <- unique(trimws(as.character(modules)))
  want_is_me  <- grepl("^(ME?\\d+)$", modules, ignore.case = TRUE)
  want_meids  <- modules[want_is_me]
  want_cols   <- unique(tolower(modules[!want_is_me]))
  
  # ---- map ME ids <-> colors
  MEs_color <- WGCNA::moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
  stopifnot(identical(rownames(MEs_color), rownames(MEs_input)))
  C <- suppressWarnings(cor(MEs_input, MEs_color, use = "pairwise.complete.obs"))
  best_col_idx <- apply(abs(C), 1, which.max)
  inputME_to_color <- setNames(sub("^ME", "", colnames(MEs_color)[best_col_idx]), colnames(MEs_input))
  color_to_inputME <- vapply(unique(sub("^ME","",colnames(MEs_color))), function(colr) {
    hits <- names(inputME_to_color)[tolower(inputME_to_color) == tolower(colr)]
    if (length(hits)) hits[1] else NA_character_
  }, character(1))
  
  # ---- resolve modules
  resolved_colors <- character(0)
  resolved_mecols <- character(0)
  
  if (length(want_cols)) {
    want_cols <- unique(tolower(want_cols))
    resolved_colors <- c(resolved_colors, want_cols)
    if (exists("color_to_inputME")) {
      resolved_mecols <- c(resolved_mecols, na.omit(unname(color_to_inputME[want_cols])))
    }
  }
  
  if (length(want_meids)) {
    want_meids <- intersect(want_meids, colnames(MEs_input))
    if (!length(want_meids)) {
      stop("None of the requested ME ids exist in MEs_input.")
    }
    resolved_mecols <- c(resolved_mecols, want_meids)
    if (exists("inputME_to_color")) {
      resolved_colors <- c(resolved_colors,
                           tolower(na.omit(unname(inputME_to_color[want_meids]))))
    }
  }
  
  resolved_mecols <- unique(na.omit(resolved_mecols))
  resolved_colors <- unique(na.omit(resolved_colors))
  
  if (!length(resolved_mecols))
    stop("Could not resolve any ME columns from 'modules'.")
  if (!length(resolved_colors))
    stop("Could not map the requested ME ids to module colors.")
  
  # ---- select genes by color
  inMod <- tolower(moduleColors) %in% resolved_colors
  stopifnot(any(inMod))
  genesInMod <- colnames(datExpr)[inMod]
  
  # ---- TOM
  TOM <- WGCNA::TOMsimilarityFromExpr(
    datExpr[, inMod, drop = FALSE],
    power       = softPower,
    corType     = "bicor",
    maxPOutliers= 0.1,
    networkType = "signed",
    TOMType     = "signed"
  )
  diag(TOM) <- 0
  
  # ---- edge threshold
  ut_mask <- upper.tri(TOM)
  vals    <- as.numeric(TOM[ut_mask])
  stopifnot(length(vals) == sum(ut_mask))
  
  thr <- stats::quantile(vals, probs = tom_top_quantile, na.rm = TRUE, names = FALSE)
  idx <- which(ut_mask & (TOM >= thr), arr.ind = TRUE)
  cat("thr:", signif(thr,4), " | edges_pre_keep:", nrow(idx), "\n")
  
  # ---- classify edge types
  edge_types <- unique(match.arg(edge_types, choices = c("LTR-Gene","Gene-Gene","LTR-LTR"), several.ok = TRUE))
  
  isLTR <- grepl("LTR", genesInMod, ignore.case = TRUE)
  
  build_edges_df <- function(keep_mask) {
    idxK <- idx[keep_mask, , drop = FALSE]
    L    <- nrow(idxK)
    if (L == 0) {
      data.frame(source=character(0), target=character(0),
                 interaction=character(0), weight=numeric(0),
                 check.names=FALSE, stringsAsFactors=FALSE)
    } else {
      w <- as.numeric(TOM[cbind(idxK[,1], idxK[,2])])
      data.frame(
        source      = genesInMod[idxK[,1]],
        target      = genesInMod[idxK[,2]],
        interaction = rep("coexp", L),
        weight      = w,
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }
  }
  
  edge_buckets <- list()
  for (etype in edge_types) {
    if (etype == "LTR-Gene") {
      keep <- (isLTR[idx[,1]] & !isLTR[idx[,2]]) | (!isLTR[idx[,1]] & isLTR[idx[,2]])
    } else if (etype == "Gene-Gene") {
      keep <- (!isLTR[idx[,1]] & !isLTR[idx[,2]])
    } else if (etype == "LTR-LTR") {
      keep <- (isLTR[idx[,1]] & isLTR[idx[,2]])
    }
    edf <- build_edges_df(keep)
    cat(sprintf("edges_post_keep (%s): %d\n", etype, nrow(edf)))
    edge_buckets[[etype]] <- edf
  }
  
  edges_df <- if ("LTR-Gene" %in% names(edge_buckets)) edge_buckets[["LTR-Gene"]] else {
    edge_buckets[[ edge_types[1] ]]
  }
  
  # ---- kME
  kME_all <- WGCNA::bicor(datExpr, MEs_input, maxPOutliers = 0.1)
  stopifnot(all(resolved_mecols %in% colnames(kME_all)))
  kME_mod <- rowMeans(kME_all[, resolved_mecols, drop = FALSE], na.rm = TRUE)
  
  # ---- optional GS
  gs_name <- "GS"
  GS_vec  <- rep(NA_real_, ncol(datExpr)); names(GS_vec) <- colnames(datExpr)
  if (!is.null(tissue_label)) {
    tcol_idx <- which(tolower(names(traits)) == tolower(trait_column))
    stopifnot(length(tcol_idx) > 0)
    tcol <- names(traits)[tcol_idx[1]]
    labs <- tolower(traits[[tcol]])
    stopifnot(any(labs == tolower(tissue_label), na.rm = TRUE))
    bin_vec <- as.numeric(labs == tolower(tissue_label))
    safe_cor1 <- function(x, y) {
      ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
      if (length(unique(y)) < 2 || length(x) < 3) return(NA_real_)
      if (stats::mad(y, na.rm = TRUE) == 0 || stats::mad(x, na.rm = TRUE) == 0) stats::cor(x, y)
      else WGCNA::bicor(x, y, maxPOutliers = 0.1)
    }
    GS_vec <- apply(datExpr, 2, safe_cor1, y = bin_vec)
    gs_name <- paste0("GS_", gsub("\\s+", "", tissue_label))
    
    # --- NEW PATCH: GS significance for each element ---
    N_eff <- sum(is.finite(bin_vec))
    GS_p  <- WGCNA::corPvalueStudent(GS_vec, N_eff)
    GS_fdr <- p.adjust(GS_p, method = "BH")
    
    # helpful column names derived from the tissue label you already use
    gs_p_name   <- paste0(gs_name, "_p")
    gs_fdr_name <- paste0(gs_name, "_FDR")
    
    # Optional significance flag (tweak thresholds to taste)
    # Here: FDR < 0.05 and |GS| >= 0.3
    gs_sig_name <- paste0(gs_name, "_sig")
    GS_sig <- (GS_fdr < 0.05) & (abs(GS_vec) >= 0.3)
  }
  
  # Build base node table
  nodes_df <- data.frame(
    id       = as.character(genesInMod),
    name     = as.character(genesInMod),
    modColor = moduleColors[inMod],
    type     = ifelse(isLTR[inMod], "LTR", "Gene"),
    kME      = as.numeric(kME_mod[inMod]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Module label column (unchanged)
  if (!"Module" %in% names(nodes_df)) {
    if (exists("inputME_to_color")) {
      color_to_me <- setNames(names(inputME_to_color), tolower(inputME_to_color))
      nodes_df$Module <- unname(color_to_me[tolower(nodes_df$modColor)])
      nodes_df$Module[is.na(nodes_df$Module)] <- nodes_df$modColor
    } else {
      nodes_df$Module <- nodes_df$modColor
    }
  }
  
  # ---- ATTACH GS + p + FDR (only if GS was computed) ----
  if (!is.null(tissue_label)) {
    nodes_df[[gs_name]]     <- as.numeric(GS_vec[inMod])
    nodes_df[[gs_p_name]]   <- as.numeric(GS_p[inMod])
    nodes_df[[gs_fdr_name]] <- as.numeric(GS_fdr[inMod])
    nodes_df[[gs_sig_name]] <- as.logical(GS_sig[inMod])
  }
  
  N <- nrow(datExpr)  # number of samples (rows)
  
  kME_p <- WGCNA::corPvalueStudent(kME_mod, N)
  kME_q <- p.adjust(kME_p, method = "BH")
  
  # Optionally keep a stats table for debugging/exports
  kME_stats <- data.frame(
    id    = names(kME_mod),
    kME   = kME_mod,
    p_kME = kME_p,
    q_kME = kME_q,
    stringsAsFactors = FALSE
  )
  
  # Attach to nodes — only for genes in this module
  nodes_df$kME_p <- as.numeric(kME_p[inMod])
  nodes_df$kME_q <- as.numeric(kME_q[inMod])
  
  # IMPORTANT: copy AFTER adding the new columns
  nodes_df_full   <- nodes_df
  nodes_df_export <- nodes_df_full
  
  # ---- thresholds
  if (nrow(nodes_df) > 0) {
    kme_thr <- as.numeric(stats::quantile(nodes_df$kME, 0.9, na.rm = TRUE))
    cat("Top 10% kME threshold:", signif(kme_thr, 4), "\n")
  } else kme_thr <- NA
  
  if (nrow(edges_df) > 0) {
    w_thr <- as.numeric(stats::quantile(edges_df$weight, 0.9, na.rm = TRUE))
    cat("Top 10% TOM edge threshold:", signif(w_thr, 4), "\n")
  } else w_thr <- NA
  
  # ---- filtering
  prefilter_for_export <- TRUE
  kME_top_quantile     <- 0.9
  edge_top_quantile    <- 0.9
  min_degree           <- 1
  max_nodes            <- Inf
  
  if (prefilter_for_export && "Gene-Gene" %in% names(edge_buckets)) {
    gg_edges <- edge_buckets[["Gene-Gene"]]
    
    if (nrow(nodes_df) > 0) {
      kme_thr <- as.numeric(stats::quantile(nodes_df$kME, kME_top_quantile, na.rm = TRUE))
    } else kme_thr <- NA_real_
    
    if (is.finite(kme_thr)) {
      nodes_high <- subset(nodes_df, is.finite(kME) & kME >= kme_thr)
    } else {
      nodes_high <- nodes_df[0, , drop = FALSE]
    }
    keep_ids <- nodes_high$id
    
    if (nrow(gg_edges) > 0 && length(keep_ids) > 0) {
      gg_edges <- subset(gg_edges, source %in% keep_ids & target %in% keep_ids)
      
      if (nrow(gg_edges) > 0) {
        w_thr <- as.numeric(stats::quantile(gg_edges$weight, edge_top_quantile, na.rm = TRUE))
        gg_edges <- subset(gg_edges, is.finite(weight) & weight >= w_thr)
      } else {
        w_thr <- NA_real_
      }
      
      if (nrow(gg_edges) > 0 && min_degree > 0) {
        deg <- table(c(gg_edges$source, gg_edges$target))
        good <- names(deg)[deg >= min_degree]
        nodes_high <- subset(nodes_high, id %in% good)
        gg_edges   <- subset(gg_edges, source %in% nodes_high$id & target %in% nodes_high$id)
      }
      
      if (nrow(nodes_high) > max_nodes) {
        ord <- order(nodes_high$kME, decreasing = TRUE)
        nodes_high <- nodes_high[ord[seq_len(max_nodes)], , drop = FALSE]
        gg_edges   <- subset(gg_edges, source %in% nodes_high$id & target %in% nodes_high$id)
      }
      
    } else {
      w_thr <- NA_real_
      gg_edges <- gg_edges[0, , drop = FALSE]
      nodes_high <- nodes_df[0, , drop = FALSE]
    }
    
    edge_buckets[["Gene-Gene"]] <- gg_edges
    
    all_ids <- unique(unlist(lapply(edge_buckets, function(ed)
      if (is.data.frame(ed) && nrow(ed) > 0) c(ed$source, ed$target) else character(0)
    )))
    if (length(all_ids) > 0) {
      nodes_df_export <- unique(subset(nodes_df_full, id %in% all_ids))
      if (nrow(nodes_df_export) == 0) nodes_df_export <- nodes_df_full
    } else {
      nodes_df_export <- nodes_df_full[0, , drop = FALSE]
    }
    
    cat("Export (Gene-Gene) — kME_thr:", signif(kme_thr,4),
        " | w_thr:", ifelse(is.finite(w_thr), signif(w_thr,4), NA),
        " | nodes:", nrow(nodes_df_export),
        " | edges:", nrow(edge_buckets[["Gene-Gene"]]), "\n")
  }
  
  if (is.finite(kme_thr) && is.finite(w_thr)) {
    nodes_high_all <- subset(nodes_df_full, is.finite(kME) & kME >= kme_thr)
    edges_high <- subset(edges_df,
                         is.finite(weight) & weight >= w_thr &
                           source %in% nodes_high_all$id &
                           target %in% nodes_high_all$id)
    cat("High-confidence nodes:", nrow(nodes_high_all), " | edges:", nrow(edges_high), "\n")
  }
  
  # ---- export
  mod_label <- paste(
    c(if (length(want_cols))  paste0("C(",  paste(sort(unique(want_cols)),  collapse = "_"), ")"),
      if (length(want_meids)) paste0("ME(", paste(sort(unique(want_meids)), collapse = "_"), ")")),
    collapse = "_"
  )
  if (mod_label == "") mod_label <- "modules"
  
  outdir <- file.path(sprintf("%s_%s", out_dir_prefix, mod_label))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # 1) Nodes + edges as CSV (what you already had)
  nodes_csv <- file.path(outdir, sprintf("%s_nodes.csv", mod_label))
  utils::write.csv(nodes_df_export, nodes_csv, row.names = FALSE)
  
  edges_csv_files <- character(0)
  for (nm in names(edge_buckets)) {
    f <- file.path(outdir,
                   sprintf("%s_edges_%s.csv", mod_label, gsub("-", "_", nm)))
    utils::write.csv(edge_buckets[[nm]], f, row.names = FALSE)
    edges_csv_files <- c(edges_csv_files, f)
  }
  
  # 2) OPTIONAL: one SIF file combining all edges (source interaction target)
  sif_file <- NULL
  if (isTRUE(export_sif)) {
    all_edges <- do.call(rbind, lapply(edge_buckets, function(ed) {
      if (is.data.frame(ed) && nrow(ed) > 0) {
        data.frame(
          source      = ed$source,
          interaction = if ("interaction" %in% names(ed)) ed$interaction else "coexp",
          target      = ed$target,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))
    if (!is.null(all_edges) && nrow(all_edges) > 0) {
      sif_file <- file.path(outdir, sprintf("%s_network.sif", mod_label))
      # Cytoscape SIF is tab-delimited, 3 columns, no header
      utils::write.table(all_edges[, c("source", "interaction", "target")],
                         file = sif_file, quote = FALSE, sep = "\t",
                         row.names = FALSE, col.names = FALSE)
    }
  }
  
  # 3) OPTIONAL: zip everything into a bundle
  zip_file <- NULL
  if (isTRUE(zip_bundle)) {
    zip_file <- file.path(outdir, sprintf("%s_cytoscape_bundle.zip", mod_label))
    files_to_zip <- c(nodes_csv, edges_csv_files)
    if (!is.null(sif_file)) files_to_zip <- c(files_to_zip, sif_file)
    utils::zip(zipfile = zip_file, files = files_to_zip)
  }
  
  # ---- Cytoscape
  if (isTRUE(connect_cytoscape)) {
    if (!.wait_for_cy()) stop("Cytoscape (CyREST) not reachable.")
    for (nm in names(edge_buckets)) {
      title_txt <- sprintf("WGCNA_%s_%s_%s", mod_label, if (is.null(tissue_label)) "tissueNA" else tissue_label, gsub("-","",nm))
      edf <- edge_buckets[[nm]]
      if (nrow(edf) > 0) {
        RCy3::createNetworkFromDataFrames(nodes = nodes_df_export, edges = edf,
                                          title = title_txt, collection = "WGCNA")
      } else {
        warning(sprintf("No edges for '%s'; creating node-only network.", nm))
        RCy3::createNetworkFromDataFrames(nodes = nodes_df_export,
                                          title = title_txt, collection = "WGCNA")
      }
    }
  }
  
  # ============================================================================
  # HEATMAP SECTION - NOW WITH PROPER ERROR HANDLING
  # ============================================================================
  if (!skip_heatmap) {
    tryCatch({
      norm <- function(x) sub("^ME", "", toupper(as.character(x)))
      mods_norm <- norm(modules)
      
      mod_vec <- if ("Module" %in% names(nodes_df_full)) nodes_df_full$Module else nodes_df_full$modColor
      
      keep_mask <- norm(mod_vec) %in% mods_norm |
        tolower(nodes_df_full$modColor) %in% resolved_colors
      
      nodes_mod <- nodes_df_full[keep_mask, , drop = FALSE]
      
      if (nrow(nodes_mod) == 0) {
        warning(sprintf("No nodes found for modules: %s. Skipping heatmap.",
                        paste(modules, collapse = ", ")))
        return(invisible(list(nodes = nodes_df_export, edges = edges_df, edge_buckets = edge_buckets,
                              tom_threshold = thr, outdir = outdir)))
      }
      
      # Choose focal LTR
      if (is.null(ltr_id)) {
        ltr_interest <- nodes_mod %>%
          dplyr::filter(tolower(type) == "ltr") %>%
          dplyr::arrange(dplyr::desc(kME)) %>%
          dplyr::slice(1) %>%
          dplyr::pull(id)
        if (length(ltr_interest) == 0 || is.na(ltr_interest)) {
          warning("No LTR found in selected module(s). Skipping heatmap.")
          return(invisible(list(nodes = nodes_df_export, edges = edges_df, edge_buckets = edge_buckets,
                                tom_threshold = thr, outdir = outdir)))
        }
      } else {
        hit <- which(nodes_mod$id == ltr_id | nodes_mod$name == ltr_id)
        if (length(hit) == 0) {
          warning(sprintf("Requested LTR '%s' not found. Skipping heatmap.", ltr_id))
          return(invisible(list(nodes = nodes_df_export, edges = edges_df, edge_buckets = edge_buckets,
                                tom_threshold = thr, outdir = outdir)))
        }
        ltr_interest <- nodes_mod$id[hit[1]]
      }
      
      # Get edges within module - USE ALL EDGE TYPES, NOT JUST LTR-Gene
      all_edges_mod <- do.call(rbind, lapply(edge_buckets, function(eb) {
        if (is.data.frame(eb) && nrow(eb) > 0) {
          subset(eb, source %in% nodes_mod$id & target %in% nodes_mod$id)
        } else {
          data.frame(source=character(0), target=character(0), weight=numeric(0))
        }
      }))
      
      if (nrow(all_edges_mod) == 0) {
        warning("No edges found within the requested module(s). Skipping heatmap.")
        return(invisible(list(nodes = nodes_df_export, edges = edges_df, edge_buckets = edge_buckets,
                              tom_threshold = thr, outdir = outdir)))
      }
      
      # Get top-K neighbors
      top_neighbors <- all_edges_mod %>%
        dplyr::filter(source == ltr_interest | target == ltr_interest) %>%
        dplyr::transmute(
          partner = ifelse(source == ltr_interest, target, source),
          weight  = weight
        ) %>%
        dplyr::group_by(partner) %>%
        dplyr::summarise(weight = max(weight), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(weight)) %>%
        dplyr::slice_head(n = top_k_partners) %>%
        dplyr::pull(partner)
      
      gene_set <- unique(c(ltr_interest, top_neighbors))
      if (length(gene_set) < 2) {
        warning("LTR has no connected partners in this module. Skipping heatmap.")
        return(invisible(list(nodes = nodes_df_export, edges = edges_df, edge_buckets = edge_buckets,
                              tom_threshold = thr, outdir = outdir)))
      }
      
      # Build expression matrix
      expr_sub <- datExpr[, gene_set, drop = FALSE]
      
      trait_col <- traits[, trait_column, drop = FALSE]
      trait_labels <- trait_col[rownames(expr_sub), 1]
      
      if (length(trait_labels) != nrow(expr_sub)) {
        stop(sprintf("Trait label length mismatch: %d labels for %d samples.",
                     length(trait_labels), nrow(expr_sub)))
      }
      
      annotation_col <- data.frame(Tissue = trait_labels, check.names = FALSE)
      rownames(annotation_col) <- rownames(expr_sub)
      
      ann_levels <- unique(na.omit(as.character(annotation_col$Tissue)))
      
      missing <- setdiff(ann_levels, names(tissue_colors))
      if (length(missing)) {
        extra <- grDevices::hcl.colors(length(missing), "Dark3")
        names(extra) <- missing
        tissue_colors <- c(tissue_colors, extra)
        message("Added colors for missing tissues: ", paste(missing, collapse = ", "))
      }
      
      annotation_colors <- list(Tissue = tissue_colors[ann_levels])
      
      # Define a fixed color palette and breaks
      breaks_fixed <- seq(-2, 2, length.out = 101)
      color_fixed  <- colorRampPalette(c("cornflowerblue", "palegoldenrod", "brown3"))(100)
      
      pheatmap::pheatmap(
        t(expr_sub),
        scale = "row",
        color = color_fixed,
        breaks = breaks_fixed,
        clustering_method = "ward.D2",
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        show_rownames = TRUE,
        show_colnames = TRUE,
        main = sprintf("Expression: %s + top %d partners",
                       ltr_interest, length(gene_set) - 1),
        angle_col = 45
      )
      
      
    }, error = function(e) {
      warning(sprintf("Heatmap generation failed: %s", e$message))
    })
  }
  
  invisible(list(
    nodes       = nodes_df_export,
    edges       = edges_df,
    edge_buckets= edge_buckets,
    tom_threshold = thr,
    outdir        = outdir,
    sif_file      = sif_file,
    zip_file      = zip_file
  ))
}
