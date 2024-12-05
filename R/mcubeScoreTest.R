mcubeTest <- function(
    object,
    kernels_list = NULL, standardize = TRUE, max_cores = 1) {
  if (is.null(kernels_list)) {
    num_kernels <- 2
    length_scale_seq <- c(sqrt(3), 3) *
      mcubeLengthScale(object@coordinates, standardize = standardize)
    object@kernels <- mapply(
      FUN = function(kernel_type, length_scale) {
        mcubeKernel(
          coordinates = object@coordinates[object@spots, , drop = FALSE],
          standardize = standardize,
          kernel_type = kernel_type,
          length_scale = length_scale
        )
      },
      kernel_type = c(
        "linear",
        rep("Gaussian", num_kernels), rep("Gaussian_transformed", num_kernels)
      ),
      length_scale = c(0, length_scale_seq, length_scale_seq),
      SIMPLIFY = FALSE
    )
    names(object@kernels) <- c(
      "linear",
      paste0("Gaussian_", 1:num_kernels),
      paste0("Gaussian_transformed_", 1:num_kernels)
    )
  } else if (!is.list(kernels_list)) {
    stop("mcubeTest: Please store the kernel matrices in a list!") # End
  } else if (!all(
    sapply(kernels_list,
      FUN = function(kernel_mat) {
        identical(rownames(object@counts), rownames(kernel_mat)) &&
          identical(rownames(object@counts), colnames(kernel_mat))
      }
    )
  )) {
    stop("mcubeTest: The rownames and colnames of the kernel matrices must match the sample names of the counts!") # End
  } else {
    if (is.null(kernels_list)) {
      names(kernels_list) <- paste0("kernel_", length(object@kernels))
    }
    object@kernels <- kernels_list
  }

  if (max_cores == 1) {
    pvalues_df <- data.frame()
    for (i in 1:nrow(object@celltype_gene_test_pairs)) {
      pvalues_i <- mcubeTestSinglePairMultiKernels(
        null_model_results = object@null_models[[i]],
        X = object@covariates[object@spots, , drop = FALSE],
        kernels_list = object@kernels,
        celltype = object@celltype_gene_test_pairs$celltype[i]
      )
      if (!is.null(pvalues_i)) {
        pvalues_i <- data.frame(
          celltype = object@celltype_gene_test_pairs$celltype[i],
          gene = object@celltype_gene_test_pairs$gene[i],
          t(pvalues_i)
        )
        pvalues_df <- rbind(pvalues_df, pvalues_i)
      }
    }
  } else if (max_cores > 1) {
    num_cores <- parallel::detectCores()
    num_cores <- ifelse(num_cores <= max_cores, num_cores - 1, max_cores)
    message("Number of cores used: ", num_cores, ".")
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    pvalues_df <- foreach::foreach(
      null_model_results_i = object@null_models,
      celltype_i = iterators::iter(object@celltype_gene_test_pairs$celltype),
      gene_i = iterators::iter(object@celltype_gene_test_pairs$gene),
      .packages = "Matrix",
      .export = c(
        "mcubeTestSinglePairSingleKernel",
        "mcubeTestSinglePairMultiKernels",
        "ACAT"
      ),
      .combine = "rbind"
    ) %dopar% {
      pvalues_i <- mcubeTestSinglePairMultiKernels(
        null_model_results = null_model_results_i,
        X = object@covariates[object@spots, , drop = FALSE],
        kernels_list = object@kernels, celltype = celltype_i
      )
      if (!is.null(pvalues_i)) {
        pvalues_i <- data.frame(
          celltype = celltype_i,
          gene = gene_i,
          t(pvalues_i)
        )
      }
      pvalues_i
    }

    parallel::stopCluster(cl)
  } else {
    stop("mcubeTest: max_cores must be a positive integer!") # End
  }
  pvalues_df <- split(pvalues_df, pvalues_df$celltype)
  pvalues_df <- lapply(pvalues_df,
    FUN = function(pvalues_celltype) {
      pvalues_celltype$celltype <- NULL
      rownames(pvalues_celltype) <- pvalues_celltype$gene
      pvalues_celltype$gene <- NULL
      return(pvalues_celltype)
    }
  )
  object@pvalues <- pvalues_df

  return(object)
}

mcubeTestSinglePairMultiKernels <- function(
    null_model_results, X,
    kernels_list, celltype) {
  if (!(celltype %in% null_model_results$celltype)) {
    return(NULL)
  }

  kernels_list <- lapply(
    kernels_list,
    FUN = function(kernel_mat) {
      kernel_mat[null_model_results$spots, null_model_results$spots, drop = FALSE]
    }
  )
  X <- X[null_model_results$spots, , drop = FALSE]

  # Compute the projection matrix
  Sigma_inv_X_mat <- null_model_results$Sigma_inv * X
  X_Sigma_inv_X_mat <- crossprod(X, Sigma_inv_X_mat)
  P_mat <- diag(null_model_results$Sigma_inv) -
    tcrossprod(
      Sigma_inv_X_mat %*% chol2inv(chol(X_Sigma_inv_X_mat)),
      Sigma_inv_X_mat
    )
  P_Y_tilde <- as.vector(P_mat %*% null_model_results$Y_tilde)

  pvalues <- sapply(
    kernels_list,
    FUN = function(kernel_mat) {
      mcubeTestSinglePairSingleKernel(
        null_model_results, P_mat, P_Y_tilde,
        kernel_mat, celltype
      )
    }
  )

  if (length(kernels_list) > 1) {
    pvalues <- c(pvalues, combined_pvalue = ACAT(pvalues))
  }
  names(pvalues)[1:length(kernels_list)] <- names(kernels_list)

  return(pvalues)
}

mcubeTestSinglePairSingleKernel <- function(
    null_model_results, P_mat, P_Y_tilde,
    kernel_mat, celltype) {
  if (!(celltype %in% null_model_results$celltype)) {
    return(NULL)
  }

  membership <- null_model_results$membership[, celltype]
  kernel_mat <- kernel_mat * (membership %o% membership)
  P_kernel_mat <- P_mat %*% kernel_mat

  teststat_var <- 0.5 * sum(P_kernel_mat^2) # variance of test statistic
  teststat_mean <- 0.5 * sum(diag(P_kernel_mat)) # mean of test statistic
  scaled_para <- 0.5 * teststat_var / teststat_mean # scaled parameter
  df_para <- 2.0 * teststat_mean^2 / teststat_var # degrees of freedom

  P_kernel_P_Y_tilde <- as.vector(P_kernel_mat %*% P_Y_tilde)
  teststat <- 0.5 *
    sum(null_model_results$Y_tilde * P_kernel_P_Y_tilde) # test statistic

  pvalue <- pchisq(q = teststat / scaled_para, df = df_para, lower.tail = FALSE)

  return(pvalue)
}
