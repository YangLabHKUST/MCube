filter_genes <- function(
    counts, library_size = NULL,
    threshold = 5e-5, batch_size = 1000) {
  message(
    "filter_genes: filter genes based on expression threshold = ",
    threshold, "."
  )

  if (is.null(library_size)) {
    library_size <- rowSums(counts)
  }

  n_genes <- ncol(counts)
  gene_means <- numeric(n_genes)
  names(gene_means) <- colnames(counts)
  n_batches <- ceiling(n_genes / batch_size)
  for (j in 1:n_batches) {
    if (j < n_batches) {
      index_range <- (1:batch_size) + (j - 1) * batch_size
    } else {
      index_range <- (1 + (n_batches - 1) * batch_size):n_genes
    }
    norm_counts <- as.matrix(counts[, index_range, drop = FALSE]) / library_size
    gene_means[index_range] <- colMeans(norm_counts)
  }
  gene_list_total <- names(which(gene_means > threshold))
  mito_genes_idx <- c(
    grep("^MT-", gene_list_total),
    grep("^mt-", gene_list_total)
  )
  if (length(mito_genes_idx) > 0) {
    gene_list_total <- gene_list_total[-mito_genes_idx]
  }
  return(gene_list_total)
}

# filter_spots <- function(proportion, celltype_test, threshold = 0.9999) {
#   spot_names <- rownames(proportion)
#   spot_names <- spot_names[
#     rowSums(proportion[, celltype_test, drop = FALSE]) >= threshold
#   ]
#   return(spot_names)
# }

#' Aggregated Cauchy Assocaition Test
#'
#' A p-value combination method using the Cauchy distribution.
#'
#' @param weights a numeric vector/matrix of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
#' @param Pvals a numeric vector/matrix of p-values. When it is a matrix, each column of p-values is combined by ACAT.
#' @return The p-value(s) of ACAT.
#' @author Yaowu Liu
#' @examples p.values <- c(2e-02, 4e-04, 0.2, 0.1, 0.8)
#' ACAT(Pvals = p.values)
#' @examples ACAT(matrix(runif(1000), ncol = 10))
#' @references Liu, Y., & Xie, J. (2019). Cauchy combination test: a powerful test with analytic p-value calculation
#' under arbitrary dependency structures. \emph{Journal of American Statistical Association},115(529), 393-402. (\href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2018.1554485}{pub})
#' @export
ACAT <- function(Pvals, Weights = NULL, threshold = 5.55e-17) {
  #### check if there is NA
  if (sum(is.na(Pvals)) > 0) {
    stop("ACAT: Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals < 0) + sum(Pvals > 1)) > 0) {
    stop("ACAT: P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero <- (sum(Pvals == 0) >= 1)
  is.one <- (sum(Pvals == 1) >= 1)
  if (is.zero && is.one) {
    # stop("ACAT: Cannot have both 0 and 1 p-values!")
    return(NA)
  }
  if (is.zero) {
    Pvals[which(Pvals == 0)] <- 5.55e-17
  }
  if (is.one) {
    Pvals[which((1 - Pvals) < 1e-3)] <- 0.99
  }

  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)) {
    Weights <- rep(1 / length(Pvals), length(Pvals))
  } else if (length(Weights) != length(Pvals)) {
    stop("ACAT: The length of weights should be the same as that of the p-values")
  } else if (sum(Weights < 0) > 0) {
    stop("ACAT: All the weights must be positive!")
  } else {
    Weights <- Weights / sum(Weights)
  }

  #### check if there are very small non-zero p values
  is.small <- (Pvals < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(Weights * tan((0.5 - Pvals) * pi))
  } else {
    cct.stat <- sum((Weights[is.small] / Pvals[is.small]) / pi)
    cct.stat <- cct.stat + sum(Weights[!is.small] * tan((0.5 - Pvals[!is.small]) * pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat > 1e+15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - stats::pcauchy(cct.stat)
  }
  pval[pval < threshold] <- min(c(Pvals, threshold))

  return(pval)
}

get_gene_list_celltype <- function(
    celltype, celltype_all, gene_test,
    library_size, proportion, reference,
    reference_threshold = 0.5, platform_effects) {
  C <- 15
  N_cells <- colSums(proportion)[celltype]
  library_size_list <- library_size[
    which(proportion[, celltype] >= .99)
  ]
  if (length(library_size_list) < 10) {
    library_size_list <- library_size[
      which(proportion[, celltype] >= .80)
    ]
  }
  if (length(library_size_list) < 10) {
    library_size_list <- library_size[
      which(proportion[, celltype] >= .50)
    ]
  }
  if (length(library_size_list) < 10) {
    library_size_list <- library_size[
      which(proportion[, celltype] >= .01)
    ]
  }
  library_size_median <- median(library_size_list)
  expr_thresh <- C / (N_cells * library_size_median)
  reference_renorm <- exp(platform_effects[gene_test]) *
    reference[celltype_all, gene_test, drop = FALSE]
  gene_list_type <- setdiff(
    gene_test,
    gene_test[which(reference_renorm[celltype, ] < expr_thresh)]
  )
  # Compare with the reference of other cell-types
  celltype_means <- reference_renorm[celltype_all, gene_list_type, drop = FALSE]
  if (ncol(proportion) > 1) {
    celltype_mean_ratio <- celltype_means[celltype, ] / apply(celltype_means, 2, max)
    gene_list_type <- gene_list_type[which(celltype_mean_ratio >= reference_threshold)]
  }

  return(gene_list_type)
}

mcubePlotPvals <- function(
    pvalues_list, which_pvalue = "combined_pvalue",
    type = "qqplot", minus_log10p_max = NULL,
    under_null = FALSE, ci = 0.95,
    nrow = 1) {
  if (is.null(minus_log10p_max)) {
    # minus_log10p_max <- ceiling(-log10(0.05 / max(sapply(pvalues_list, length))))
    minus_log10p_max <- 5
  }
  pvalues_long <- do.call(
    rbind,
    lapply(
      names(pvalues_list),
      FUN = function(x) {
        order_x <- order(pvalues_list[[x]][, which_pvalue])
        pvalues_x <- pvalues_list[[x]][, which_pvalue][order_x]
        n_palues <- length(pvalues_x)
        pvalues_theoretical_x <- ppoints(n_palues)
        minus_log10p_x <- -log10(pvalues_x)
        minus_log10p_x[minus_log10p_x > minus_log10p_max] <- minus_log10p_max
        minus_log10p_theoretical_x <- -log10(pvalues_theoretical_x)
        lower_bound <- -log10(qbeta(
          p = (1 - ci) / 2, shape1 = 1:n_palues, shape2 = n_palues:1
        ))
        upper_bound <- -log10(qbeta(
          p = (1 + ci) / 2, shape1 = 1:n_palues, shape2 = n_palues:1
        ))

        data.frame(
          celltype = x,
          gene = rownames(pvalues_list[[x]])[order_x],
          pvalue = pvalues_x,
          pvalues_theoretical = pvalues_theoretical_x,
          minus_log10p = minus_log10p_x,
          minus_log10p_theoretical = minus_log10p_theoretical_x,
          lower_bound = lower_bound, upper_bound = upper_bound
        )
      }
    )
  )
  if (type == "qqplot") {
    expected_minus_log10p_lab <- expression(paste("Expected -log"[10], plain(P)))
    observed_minus_log10p_lab <- expression(paste("Observed -log"[10], plain(P)))
    p <- ggplot2::ggplot(
      data = pvalues_long,
      ggplot2::aes(x = minus_log10p_theoretical, y = minus_log10p)
    )

    if (under_null) {
      p <- p + ggplot2::geom_ribbon(
        aes(
          x = minus_log10p_theoretical,
          ymin = lower_bound, ymax = upper_bound
        ),
        fill = "darkgrey", alpha = 0.5
      ) +
        ggplot2::geom_abline(
          intercept = 0, slope = 1,
          linewidth = 2, col = "white"
        )
    } else {
      p <- p + ggplot2::geom_abline(
        intercept = 0, slope = 1,
        linetype = "longdash", linewidth = 2, color = "black"
      )
    }

    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = celltype),
      size = 5, alpha = 0.5
    ) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::facet_wrap(celltype ~ ., nrow = nrow) +
      ggplot2::labs(
        title = "P-value QQ plots",
        x = expected_minus_log10p_lab, y = observed_minus_log10p_lab
      ) +
      ggplot2::scale_color_discrete(name = "Cell-type") +
      ggplot2::theme_classic() +
      ggplot2::theme(strip.text = ggplot2::element_blank())
  } else if (type == "histogram") {
    p <- ggplot2::ggplot(
      data = pvalues_long,
      ggplot2::aes(x = pvalue, fill = celltype)
    ) +
      ggplot2::geom_histogram() +
      ggplot2::facet_wrap(celltype ~ ., scales = "free", nrow = nrow) +
      ggplot2::labs(
        title = "P-value histograms",
        x = "P values", y = "Frequency"
      ) +
      ggplot2::scale_fill_discrete(name = "Cell-type") +
      ggplot2::theme_classic() +
      ggplot2::theme(strip.text = ggplot2::element_blank())
  } else {
    stop("mcubePlot: type must be either 'qqplot' or 'histogram'!")
  }

  return(p)
}

get_platform_effects <- function(counts, library_size, proportion, reference, spot_effects) {
  bulk_vec <- colSums(counts)
  weight_celltype <- colSums(
    library_size * exp(spot_effects) * proportion
  ) / sum(library_size)
  weight_avg <- colSums(
    reference * weight_celltype / sum(weight_celltype),
  )
  target_means <- bulk_vec / sum(library_size)
  platform_effects <- log(target_means / weight_avg)
  names(platform_effects) <- colnames(reference)
  return(platform_effects)
}

choose_celltypes <- function(
    proportion, celltype_test = NULL,
    proportion_threshold = 0.05, celltype_threshold = 100) {
  proportion[proportion < proportion_threshold] <- 0
  celltype_default <- names(which(colSums(proportion) >= celltype_threshold))

  if (!is.null(celltype_test)) {
    diff_celltypes <- setdiff(celltype_test, celltype_default)
    if (length(diff_celltypes) > 0) {
      message(
        "choose_celltypes: cell-types ",
        paste0(diff_celltypes, collapse = ", "),
        " have less than the minimum celltype_threshold = ", celltype_threshold,
        ". To include these cell-types, please reduce the celltype_threshold."
      )
    }
    celltype_test <- intersect(celltype_test, celltype_default)
  } else {
    celltype_test <- celltype_default
  }
  if (length(celltype_test) == 0) {
    stop("choose_celltypes: no cell-types occure greater than celltype_threshold!")
  }

  message(
    "choose_celltypes: cell-types ",
    paste(celltype_default, collapse = ", "),
    " pass the celltype_threshold = ", celltype_threshold, "."
  )

  return(celltype_test)
}

filter_spots_celltypes <- function(proportion, celltype_test, threshold = 0.75) {
  spot_names <- rownames(proportion)
  spot_names <- spot_names[
    rowSums(proportion[, celltype_test, drop = FALSE]) >= threshold
  ]
  if (length(spot_names) == 0) {
    stop("filter_spots_celltypes: No spots remain after filtering based on cell type proportions!")
  }
  return(spot_names)
}

# get_doublet_proportion <- function(proportion) {
#   doublet_idx <- apply(
#     proportion,
#     MARGIN = 1,
#     FUN = function(x) {
#       order(x, decreasing = TRUE)[1:2]
#     }
#   )
#   doublet_idx <- cbind(
#     rep(1:nrow(proportion), each = 2),
#     as.vector(doublet_idx)
#   )
#   proportion_doublet <- Matrix::sparseMatrix(
#     i = doublet_idx[, 1],
#     j = doublet_idx[, 2],
#     x = proportion[doublet_idx],
#     dims = dim(proportion)
#   )
#   rownames(proportion_doublet) <- rownames(proportion)
#   colnames(proportion_doublet) <- colnames(proportion)
#   return(proportion_doublet)
# }

rotation_matrix_2d <- function(theta) {
  theta_rad <- theta * (pi / 180) # 将角度转换为弧度
  matrix(
    c(
      cos(theta_rad), -sin(theta_rad),
      sin(theta_rad), cos(theta_rad)
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
}

mcubePlot_converge <- function(
    pvalues_list, converge,
    which_pvalue = "combined_pvalue", type = "qqplot") {
  pvalues_long <- do.call(
    rbind,
    lapply(
      names(pvalues_list),
      FUN = function(x) {
        pvalues_x <- pvalues_list[[x]][, which_pvalue]
        minus_log10p_x <- -log10(pvalues_x)
        data.frame(
          celltype = x,
          pvalue = pvalues_x,
          minus_log10p = minus_log10p_x
        )
      }
    )
  )[converge, ]
  if (type == "qqplot") {
    expected_minus_log10p_lab <- expression(paste("Expected -log"[10], plain(P)))
    observed_minus_log10p_lab <- expression(paste("Observed -log"[10], plain(P)))
    p <- ggplot2::ggplot(
      data = pvalues_long,
      ggplot2::aes(sample = -log10(pvalue), color = celltype)
    ) +
      ggplot2::stat_qq(distribution = qexp, dparams = list(rate = log(10))) +
      ggplot2::geom_abline(
        intercept = 0, slope = 1,
        linetype = "dashed", color = "black", linewidth = 2
      ) +
      ggplot2::facet_wrap(celltype ~ ., scales = "free") +
      ggplot2::labs(
        title = "P-value QQ plots",
        x = expected_minus_log10p_lab, y = observed_minus_log10p_lab
      ) +
      ggplot2::scale_color_discrete(name = "Cell-type") +
      ggplot2::theme_bw()
  } else if (type == "histogram") {
    p <- ggplot2::ggplot(
      data = pvalues_long,
      ggplot2::aes(x = pvalue, fill = celltype)
    ) +
      ggplot2::geom_histogram() +
      ggplot2::facet_wrap(celltype ~ ., scales = "free") +
      ggplot2::labs(
        title = "P-value histograms",
        x = "P values", y = "Frequency"
      ) +
      ggplot2::scale_fill_discrete(name = "Cell-type") +
      ggplot2::theme_bw()
  } else {
    stop("mcubePlot: type must be either 'qqplot' or 'histogram'!")
  }

  return(p)
}

mcubeGetSigGenes <- function(
    pvalues_list, which_pvalue = "combined_pvalue",
    adjust_method = "BH", alpha = 0.05) {
  all_methods <- c(
    "holm", "hochberg", "hommel", "bonferroni",
    "BH", "BY", "fdr", "none"
  )
  if (!(adjust_method %in% all_methods)) {
    stop(
      "mcubeGetSigGenes: adjust_method must be one of ",
      paste(all_methods, collapse = ", "), "!"
    ) # End
  }
  message(
    "mcubeGetSigGenes: Set adjust_method as ", adjust_method,
    " and alpha as ", alpha, "."
  )

  sig_genes_list <- lapply(
    pvalues_list,
    FUN = function(x) {
      x <- x[order(x[, which_pvalue]), ]
      adjusted_pvalues <- p.adjust(x[, which_pvalue], method = adjust_method)
      sig_genes_idx <- which(adjusted_pvalues <= alpha)
      data.frame(
        pvalue = x[sig_genes_idx, which_pvalue],
        adjusted_pvalue = adjusted_pvalues[sig_genes_idx],
        row.names = rownames(x)[sig_genes_idx]
      )
    }
  )
  names(sig_genes_list) <- names(pvalues_list)

  return(sig_genes_list)
}

mcubePlotExprCellType <- function(
    object, celltype, gene, normalize = TRUE,
    he_image = NULL, size = 1, palettes = rev(pals::brewer.spectral(10)),
    xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
  pair_name <- paste(celltype, gene, sep = "_")
  null_model_results <- object@null_models[[pair_name]]
  spot_names <- null_model_results$spots
  expr_level <- null_model_results$u[, celltype]
  if (normalize) {
    expr_level <- (expr_level - min(expr_level)) /
      (max(expr_level) - min(expr_level))
  }
  plot_df <- data.frame(
    x = object@coordinates[spot_names, 1],
    y = object@coordinates[spot_names, 2],
    expr_level = expr_level
  )

  p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = x, y = y))

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = 0, xmax = max(object@coordinates[, 1]),
      ymin = 0, ymax = max(object@coordinates[, 2])
    )
  }

  p <- p + ggplot2::geom_point(ggplot2::aes(color = expr_level), size = size) +
    ggplot2::scale_colour_gradientn(
      name = "Level", colors = palettes,
      values = scales::rescale(
        c(min(expr_level), median(expr_level), max(expr_level))
      )
    )

  if (is.null(xlim)) {
    xlim <- c(0, max(object@coordinates[, 1]) + 1)
  }
  if (is.null(ylim)) {
    ylim <- c(0, max(object@coordinates[, 2]) + 1)
  }
  p <- p + ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::coord_fixed(ratio = ratio)

  title <- ifelse(is.null(title),
    ifelse(normalize,
      paste("Normalized expression of", gene, "in", celltype),
      paste("Expression of", gene, "in", celltype)
    ),
    title
  )
  p <- p + ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_classic()

  return(p)
}

# mcubePlotExprCellTypeBinary <- function(
#     object, celltype, gene, spots = NULL,
#     upper_percentile = 15, lower_percentile = 15,
#     he_image = NULL, size = 1, palettes = c("#32CD32", "#c5c5c5", "#ff69b4"),
#     xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
#   pair_name <- paste(celltype, gene, sep = "_")
#   null_model_results <- object@null_models[[pair_name]]
#   if (is.null(spots)) {
#     spots <- null_model_results$spots
#   } else {
#     spots <- intersect(spots, null_model_results$spots)
#   }
#   u <- null_model_results$u[spots, celltype]
#   upper_threshold <- quantile(u, 1 - upper_percentile / 100)
#   lower_threshold <- quantile(u, lower_percentile / 100)
#   rel_expr_level <- factor(
#     ifelse(
#       u >= upper_threshold, 1,
#       ifelse(u <= lower_threshold, -1, 0)
#     )
#   )
#   plot_df <- data.frame(
#     x = object@coordinates[spots, 1],
#     y = object@coordinates[spots, 2],
#     rel_expr_level = rel_expr_level
#   )

#   p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = x, y = y))

#   # H&E image background
#   if (!is.null(he_image)) {
#     p <- p + ggplot2::annotation_raster(he_image,
#       xmin = 0, xmax = max(object@coordinates[, 1]),
#       ymin = -max(object@coordinates[, 2]), ymax = 0
#     )
#   }

#   p <- p + ggplot2::geom_point(ggplot2::aes(color = rel_expr_level), size = size) +
#     ggplot2::scale_color_manual(
#       name = "Level",
#       values = c("-1" = palettes[1], "0" = palettes[2], "1" = palettes[3])
#     )

#   if (is.null(xlim)) {
#     xlim <- c(0, max(object@coordinates[, 1]) + 1)
#   }
#   if (is.null(ylim)) {
#     ylim <- c(0, max(object@coordinates[, 2]) + 1)
#   }
#   p <- p + ggplot2::scale_x_continuous(limits = xlim) +
#     ggplot2::scale_y_continuous(trans = "reverse", limits = rev(ylim)) +
#     ggplot2::coord_fixed(ratio = ratio)

#   title <- ifelse(
#     is.null(title),
#     paste("Relative expression of", gene, "in", celltype),
#     title
#   )
#   p <- p + ggplot2::labs(title = title, x = NULL, y = NULL) +
#     ggplot2::theme_classic()

#   return(p)
# }

mcubePlotExprCellTypeBinary <- function(
    object, celltype, gene, spots = NULL,
    he_image = NULL, size = 1, palettes = c("#32CD32", "#ff69b4"),
    xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
  pair_name <- paste(celltype, gene, sep = "_")
  null_model_results <- object@null_models[[pair_name]]
  if (is.null(spots)) {
    spots <- null_model_results$spots
  } else {
    spots <- intersect(spots, null_model_results$spots)
  }
  u <- null_model_results$u[spots, celltype]
  rel_expr_level <- factor(ifelse(u >= 0, 1, 0))
  plot_df <- data.frame(
    x = object@coordinates[spots, 1],
    y = object@coordinates[spots, 2],
    rel_expr_level = rel_expr_level
  )

  p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = x, y = y))

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = 0, xmax = max(object@coordinates[, 1]),
      ymin = 0, ymax = max(object@coordinates[, 2])
    )
  }

  p <- p + ggplot2::geom_point(ggplot2::aes(color = rel_expr_level), size = size) +
    ggplot2::scale_color_manual(
      name = "Level",
      values = c("0" = palettes[1], "1" = palettes[2])
    )

  if (is.null(xlim)) {
    xlim <- c(0, max(object@coordinates[, 1]) + 1)
  }
  if (is.null(ylim)) {
    ylim <- c(0, max(object@coordinates[, 2]) + 1)
  }
  p <- p + ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::coord_fixed(ratio = ratio)

  title <- ifelse(
    is.null(title),
    paste("Relative expression of", gene, "in", celltype),
    title
  )
  p <- p + ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_classic()

  return(p)
}

mcubePlotExprCellType3D <- function(
    object, celltype, gene,
    spots = NULL, proportion_threshold = 0.5,
    plot_method = "plotly", # "plotly" or "rgl"
    palettes = c("#32CD32", "#FF69B4"),
    opacity_target = 0.9, opacity_background = 0.1,
    axis_rescale = c(1, 1, 1), spot_size = 4, spot_radius = 0.5,
    plotly_center = list(x = 0, y = 0, z = 0),
    plotly_eye = list(x = 1.25, y = 1.25, z = 1.25),
    rgl_um = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)) {
  pair_name <- paste(celltype, gene, sep = "_")
  null_model_results <- object@null_models[[pair_name]]

  spots_all <- rownames(object@proportion)
  if (is.null(spots)) {
    spots_target <- null_model_results$spots
  } else {
    spots_target <- intersect(spots, null_model_results$spots)
  }
  spots_target <- intersect(
    spots_target,
    spots_all[object@proportion[, celltype] >= proportion_threshold]
  )
  spots_background <- setdiff(spots_all, spots_target)

  u <- null_model_results$u[spots_target, celltype]
  rel_expr_level <- ifelse(u >= 0, 1, -1)
  rel_expr_level_all <- c(-1, 1)

  target_df <- data.frame(
    x = object@coordinates[spots_target, 1],
    y = object@coordinates[spots_target, 2],
    z = object@coordinates[spots_target, 3],
    rel_expr_level = rel_expr_level
  )
  background_df <- data.frame(
    x = object@coordinates[spots_background, 1],
    y = object@coordinates[spots_background, 2],
    z = object@coordinates[spots_background, 3]
  )

  # 3D plot
  if (plot_method == "plotly") {
    plotly_df <- data.frame(
      x = c(target_df$x, background_df$x) * axis_rescale[1],
      y = c(target_df$y, background_df$y) * axis_rescale[2],
      z = c(target_df$z, background_df$z) * axis_rescale[3],
      col_group = factor(
        c(target_df$rel_expr_level, rep(0, nrow(background_df))),
        levels = c(-1, 1, 0)
      ),
      opacity = c(
        rep(opacity_target, nrow(target_df)),
        rep(opacity_background, nrow(background_df))
      )
    )
    plotly::plot_ly(
      plotly_df,
      type = "scatter3d", mode = "markers",
      x = ~x, y = ~y, z = ~z,
      color = ~col_group, colors = c(palettes, "gray"),
      opacity = ~opacity, marker = list(size = spot_size)
    ) %>% plotly::layout(
      scene = list(
        aspectmode = "data",
        camera = list(center = plotly_center, eye = plotly_eye),
        xaxis = list(
          showline = TRUE, zeroline = FALSE, showgrid = FALSE,
          showaxeslabels = FALSE, showticklabels = FALSE
        ),
        yaxis = list(
          showline = TRUE, zeroline = FALSE, showgrid = FALSE,
          showaxeslabels = FALSE, showticklabels = FALSE
        ),
        zaxis = list(
          showline = TRUE, zeroline = FALSE, showgrid = FALSE,
          showaxeslabels = FALSE, showticklabels = FALSE
        )
      ),
      showlegend = FALSE
    )
  } else if (plot_method == "rgl") {
    rgl::open3d(windowRect = c(0, 0, 720, 720))
    rgl::par3d(persp)
    rgl::view3d(userMatrix = matrix(rgl_um, byrow = TRUE, nrow = 4))
    # Target spots
    for (l in 1:length(rel_expr_level_all)) {
      rgl::spheres3d(
        target_df[target_df$rel_expr_level == rel_expr_level_all[l], ]$x * axis_rescale[1],
        target_df[target_df$rel_expr_level == rel_expr_level_all[l], ]$y * axis_rescale[2],
        target_df[target_df$rel_expr_level == rel_expr_level_all[l], ]$z * axis_rescale[3],
        col = palettes[l], radius = spot_radius, alpha = opacity_target
      )
    }
    # Background in grey
    rgl::spheres3d(
      background_df$x * axis_rescale[1],
      background_df$y * axis_rescale[2],
      background_df$z * axis_rescale[3],
      col = "gray", radius = spot_radius, alpha = opacity_background
    )
    rgl::decorate3d()
  }
}

mcubePlotExpr <- function(
    counts, coordinates, gene, spots = NULL, normalize = TRUE,
    he_image = NULL, size = 1, palettes = pals::brewer.blues(20),
    xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
  # Check sample names
  if (!all.equal(rownames(counts), rownames(coordinates))) {
    stop("mcubePlotExpr: rownames of counts and coordinates do not match!")
  }
  expr <- counts[, gene]
  if (is.null(spots)) {
    spots <- rownames(counts)[expr > 0]
  } else {
    spots <- intersect(spots, rownames(counts)[expr > 0])
  }
  if (normalize) {
    expr <- log(expr + 1)
    expr <- (expr - min(expr)) / (max(expr) - min(expr))
    value_lim <- c(0, 1)
  } else {
    value_lim <- c(0, max(expr))
  }
  expr <- pmax(pmin(expr, value_lim[2] - 1e-8), value_lim[1] + 1e-8)
  plot_df <- data.frame(
    x = coordinates[spots, 1],
    y = coordinates[spots, 2],
    expr = expr[spots]
  )

  p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = x, y = y))

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = 0, xmax = max(coordinates[, 1]),
      ymin = 0, ymax = max(coordinates[, 2])
    )
  }

  p <- p + ggplot2::geom_point(ggplot2::aes(color = expr), size = size) +
    ggplot2::scale_colour_gradientn(
      name = "Level", colors = palettes,
      # values = scales::rescale(c(min(expr), median(expr[expr > 0]), max(expr)))
    )

  if (is.null(xlim)) {
    xlim <- c(0, max(coordinates[, 1]) + 1)
  }
  if (is.null(ylim)) {
    ylim <- c(0, max(coordinates[, 2]) + 1)
  }
  p <- p + ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::coord_fixed(ratio = ratio)

  if (is.null(title)) {
    title <- ifelse(
      normalize,
      paste("Normalized expression of", gene),
      paste("Expression of", gene)
    )
  }
  p <- p + ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_classic()

  return(p)
}

mcubePlotPropCellType <- function(
    proportion, coordinates, celltype, spots = NULL,
    he_image = NULL, size = 1, palettes = pals::brewer.orrd(22)[3:22],
    xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
  # Check sample names
  if (!all.equal(rownames(proportion), rownames(coordinates))) {
    stop("mcubePlotPropCellType: rownames of proportion and coordinates do not match!")
  }
  if (is.null(spots)) {
    spots <- rownames(proportion)
  } else {
    spots <- intersect(spots, rownames(proportion))
  }
  plot_df <- data.frame(
    x = coordinates[spots, 1],
    y = coordinates[spots, 2],
    prop = proportion[spots, celltype]
  )
  plot_df <- plot_df[plot_df$prop > 0, , drop = FALSE]

  p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = x, y = y))

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = 0, xmax = max(coordinates[, 1]),
      ymin = 0, ymax = max(coordinates[, 2])
    )
  }

  p <- p + ggplot2::geom_point(ggplot2::aes(color = prop), size = size) +
    ggplot2::scale_colour_gradientn(name = NULL, colors = palettes)

  if (is.null(xlim)) {
    xlim <- c(-1, max(coordinates[, 1] + 1))
  }
  if (is.null(ylim)) {
    ylim <- c(-1, max(coordinates[, 2] + 1))
  }
  p <- p + ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::coord_fixed(ratio = ratio)

  if (is.null(title)) {
    title <- paste("Proportion of", celltype)
  }
  p <- p + ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_classic()

  return(p)
}

mcubePlotPropHeatmap <- function(
    proportion, spots = NULL, celltypes = NULL,
    palettes = pals::brewer.purd(20), title = NULL, filename = NA) {
  # Check sample names
  if (!is.null(spots)) {
    spots <- intersect(spots, rownames(proportion))
    if (length(spots) == 0) {
      stop("mcubePlotProp: No samples shared between proportion and spots!")
    }
    proportion <- proportion[spots, , drop = FALSE]
  }

  # Check cell-type names
  if (!is.null(celltypes)) {
    celltypes <- intersect(celltypes, colnames(proportion))
    if (length(celltypes) == 0) {
      stop("mcubePlotProp: No cell-types shared between proportion and celltypes!")
    }
    proportion <- proportion[, celltypes, drop = FALSE]
  }

  # Heatmap of the proportion matrix
  p <- pheatmap::pheatmap(t(proportion),
    scale = "none",
    color = palettes, border_color = NA,
    clustering_method = "complete", cluster_row = TRUE, cluster_col = TRUE,
    treeheight_row = 3, treeheight_col = 10,
    show_rownames = TRUE, show_colnames = FALSE,
    cellwidth = 0.05, cellheight = 16, fontsize = 12,
    filename = filename
  )

  return(p)
}
