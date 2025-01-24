#' Filter out lowly expressed genes in spatial transcriptomics data
#'
#' @description
#' Fuction for data preprocessing.
#' Credit goes to the R package `spacexr` (\url{https://github.com/dmcable/spacexr}).
#'
#' @param counts A matrix. Each row represents a spot and each column represents a gene.
#' @param library_size A numeric vector containing library sizes of all spots.
#' @param threshold A numeric value. Genes with average relative expression below this threshold are filtered out.
#' @param batch_size An integer.
#' The number of genes to process in each batch.
#'
#' @return A character vector containing the names of genes that pass the threshold.
#'
#' @export
mcubeFilterGenes <- function(
    counts, library_size = NULL,
    threshold = 5e-5, batch_size = 1000) {
  message(
    "mcubeFilterGenes: Filter genes based on relative expression",
    " with threshold = ", threshold, "."
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

#' Filter out lowly expressed genes for a specific cell type
#'
#' @description
#' Fuction for data preprocessing.
#' Credit goes to the R package `spacexr` (\url{https://github.com/dmcable/spacexr}).
#'
#' @importFrom stats median
#'
#' @param celltype A character specifying the cell type.
#' @param celltype_all A character vector containing all cell types.
#' @param gene_test A character vector specifying the genes to test.
#' @param library_size A numeric vector containing library sizes of all spots.
#' @param proportion A numeric matrix containing the proportion of each cell type at each spot.
#' @param reference A numeric matrix containing the reference expression of genes for each cell type.
#' @param reference_threshold A numeric value. Genes with relative expression below this threshold will be filtered out.
#' @param platform_effects A numeric vector containing the platform effects for each gene.
#'
#' @return A character vector containing the names of genes that pass the threshold.
#'
#' @export
mcubeFilterGenesCellType <- function(
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
  reference_renorm <- sweep(
    reference[celltype_all, gene_test, drop = FALSE],
    MARGIN = 2,
    STATS = exp(platform_effects[gene_test]),
    FUN = "*"
  )
  gene_list_type <- setdiff(
    gene_test,
    gene_test[which(reference_renorm[celltype, ] < expr_thresh)]
  )
  # Compare with the reference of other cell types
  celltype_means <- reference_renorm[celltype_all, gene_list_type, drop = FALSE]
  if (ncol(proportion) > 1) {
    celltype_mean_ratio <- celltype_means[celltype, ] / apply(celltype_means, 2, max)
    gene_list_type <- gene_list_type[which(celltype_mean_ratio >= reference_threshold)]
  }

  return(gene_list_type)
}

#' Get the platform effects for each gene
#'
#' @description
#' Fuction for data preprocessing.
#' Credit goes to the R package `spacexr` (\url{https://github.com/dmcable/spacexr}).
#'
#' @param counts A matrix. Each row represents a spot and each column represents a gene.
#' @param library_size A numeric vector containing library sizes of all spots.
#' @param proportion A numeric matrix containing the proportion of each cell type at each spot.
#' @param reference A numeric matrix containing the reference expression of genes for each cell type.
#' @param spot_effects A numeric vector containing the spot effects for each spot.
#'
#' @return A numeric vector containing the platform effects for each gene.
#'
#' @export
mcubeGetPlatformEffects <- function(
    counts, library_size, proportion,
    reference, spot_effects) {
  bulk_vec <- colSums(counts)
  weight_celltype <- colSums(
    library_size * exp(spot_effects) * proportion
  ) / sum(library_size)
  weight_avg <- colSums(
    reference * weight_celltype / sum(weight_celltype),
  )
  target_means <- bulk_vec / sum(library_size)
  platform_effects <- log(target_means / weight_avg)
  platform_effects[!is.finite(platform_effects)] <- 0
  names(platform_effects) <- colnames(reference)
  return(platform_effects)
}

#' Get main cell types
#'
#' @description
#' Fuction for data preprocessing.
#' Credit goes to the R package `spacexr` (\url{https://github.com/dmcable/spacexr}).
#'
#' @param proportion A numeric matrix containing the proportion of each cell type at each spot.
#' @param celltype_test A character vector specifying the cell types to test.
#' @param proportion_threshold A numeric value. The proportions below this threshold will be ignored.
#' @param celltype_threshold A numeric value. Cell types with column sums less than this number will be filtered out.
#'
#' @return A character vector containing the names of cell types that pass the threshold.
#'
#' @export
mcubeFilterCellTypes <- function(
    proportion, celltype_test = NULL,
    proportion_threshold = 0.1, celltype_threshold = 100) {
  proportion[proportion < proportion_threshold] <- 0
  celltype_default <- names(which(colSums(proportion) >= celltype_threshold))

  if (!is.null(celltype_test)) {
    diff_celltypes <- setdiff(celltype_test, celltype_default)
    if (length(diff_celltypes) > 0) {
      message(
        "mcubeFilterCellTypes: Cell types ",
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
    stop(
      "mcubeFilterCellTypes: No cell types occure greater than",
      " celltype_threshold ", celltype_threshold, "!"
    )
  }

  message(
    "mcubeFilterCellTypes: Cell types ",
    paste(celltype_default, collapse = ", "),
    " pass the celltype_threshold = ", celltype_threshold, "."
  )

  return(celltype_test)
}

#' Get cell-type-specific spatially variable genes
#'
#' @description Get spatially variable genes specific to a cell type based on the adjusted p-values.
#'
#' @param pvalues_list A list of data frames. Each data frame contains p-values of genes for a certain cell type.
#' @param which_pvalue A character specifying the column name of p-values in each data frame.
#' @param adjust_method A character specifying the method for multiple testing correction.
#' @param alpha A numeric value. The significance level.
#'
#' @return A list of data frames. Each data frame contains the significant genes for a certain cell type.
#'
#' @export
mcubeGetSigGenes <- function(
    pvalues_list, which_pvalue = "combined_pvalue",
    adjust_method = "BH", alpha = 0.05) {
  all_methods <- c(
    "holm", "hochberg", "hommel", "bonferroni",
    "BH", "BY", "fdr", "none"
  )
  if (!(adjust_method %in% all_methods)) {
    stop(
      "mcubeGetSigGenes: Adjust method must be one of ",
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

#' Aggregated Cauchy Assocaition Test
#'
#' @description
#' A p-value combination method using the Cauchy distribution.
#' Credit goes to the R package `ACAT` (\url{https://github.com/yaowuliu/ACAT}).
#'
#' @param weights a numeric vector/matrix of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
#' @param Pvals a numeric vector/matrix of p-values. When it is a matrix, each column of p-values is combined by ACAT.
#'
#' @return The p-value(s) of ACAT.
#'
#' @examples p.values <- c(2e-02, 4e-04, 0.2, 0.1, 0.8)
#' ACAT(Pvals = p.values)
#' @examples ACAT(matrix(runif(1000), ncol = 10))
#'
#' @references Liu, Y., & Xie, J. (2019). Cauchy combination test: a powerful test with analytic p-value calculation
#' under arbitrary dependency structures. \emph{Journal of American Statistical Association},115(529), 393-402. (\href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2018.1554485}{pub})
#'
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
