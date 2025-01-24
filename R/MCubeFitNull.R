#' Fit the null models
#'
#' @description Fit the null MMM models for all celltype-gene pairs using a PQL-based approach.
#'
#' @param object An \code{\linkS4class{MCUBE}} object.
#' @param reference_threshold A numeric value between 0 and 1.
#' The minimum relative expression level of a gene to be considered for a cell type when fitting the null model.
#' @param safeguard A numeric value.
#' A small positive number to avoid numerical issues when computing the inverse of a matrix.
#' @param iter_max A positive integer.
#' The maximum number of iterations for the optimization algorithm.
#' @param tol A numeric value.
#' The convergence criteria for the optimization algorithm.
#' @param verbose A logical value.
#' Whether to print the progress of the optimization algorithm. Default is FALSE.
#' @param max_cores A positive integer.
#' The maximum number of cores to use for parallel computing. Default is 1.
#'
#' @return An \code{\linkS4class{MCUBE}} object with fitted null models for all celltype-gene pairs.
#'
#' @export
mcubeFitNull <- function(
    object, reference_threshold = 0.25,
    safeguard = 1e-6, iter_max = 100, tol = 1e-6,
    verbose = FALSE, max_cores = 1) {
  object@config$reference_threshold_fit <- reference_threshold
  object@config$safeguard <- safeguard
  object@config$iter_max <- iter_max
  object@config$tol <- tol

  if (max_cores == 1) {
    object@null_models <- list()
    for (i in 1:nrow(object@celltype_gene_test_pairs)) {
      null_model_results_i <- try(mcubeFitNullSinglePair(
        Y = object@counts[object@spots, object@celltype_gene_test_pairs[i, "gene"]],
        library_size = object@library_size[object@spots],
        X = object@covariates[object@spots, , drop = FALSE],
        batch_id = object@batch_id[object@spots],
        proportion = object@proportion[object@spots, , drop = FALSE],
        reference = object@reference[, object@celltype_gene_test_pairs[i, "gene"]],
        used_for_deconvolution = object@used_for_deconvolution[object@celltype_gene_test_pairs[i, "gene"]],
        spot_effects = object@spot_effects[object@spots],
        platform_effect = object@platform_effects[, object@celltype_gene_test_pairs[i, "gene"]],
        celltype_test = object@celltype_gene_test_pairs[i, "celltype"],
        proportion_threshold = object@config$proportion_threshold,
        reference_threshold = object@config$reference_threshold_fit,
        safeguard = object@config$safeguard,
        iter_max = object@config$iter_max,
        tol = object@config$tol,
        verbose = verbose
      ))
      object@null_models <- append(object@null_models, list(null_model_results_i))
    }
  } else if (max_cores > 1) {
    num_cores <- parallel::detectCores()
    num_cores <- ifelse(num_cores <= max_cores, num_cores - 1, max_cores)
    message("Number of cores used: ", num_cores, ".")
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    object@null_models <- foreach::foreach(
      Y_i = iterators::iter(
        as.matrix(object@counts[object@spots, object@celltype_gene_test_pairs$gene, drop = FALSE]),
        by = "column"
      ),
      reference_i = iterators::iter(
        object@reference[, object@celltype_gene_test_pairs$gene, drop = FALSE],
        by = "column"
      ),
      used_for_deconvolution_i = iterators::iter(
        object@used_for_deconvolution[object@celltype_gene_test_pairs$gene]
      ),
      platform_effect_i = iterators::iter(
        object@platform_effects[, object@celltype_gene_test_pairs$gene, drop = FALSE],
        by = "column"
      ),
      celltype_test_i = iterators::iter(
        object@celltype_gene_test_pairs$celltype
      ),
      .export = "mcubeFitNullSinglePair",
      .combine = "c"
    ) %dopar% {
      null_model_results_i <- try(mcubeFitNullSinglePair(
        Y = as.vector(Y_i),
        library_size = object@library_size[object@spots],
        X = object@covariates[object@spots, , drop = FALSE],
        batch_id = object@batch_id[object@spots],
        proportion = object@proportion[object@spots, , drop = FALSE],
        reference = as.vector(reference_i),
        used_for_deconvolution = used_for_deconvolution_i,
        spot_effects = object@spot_effects[object@spots],
        platform_effect = as.vector(platform_effect_i),
        celltype_test = celltype_test_i,
        proportion_threshold = object@config$proportion_threshold,
        reference_threshold = object@config$reference_threshold_fit,
        safeguard = object@config$safeguard,
        iter_max = object@config$iter_max,
        tol = object@config$tol,
        verbose = verbose
      ))
      list(null_model_results_i)
    }

    parallel::stopCluster(cl)
  } else {
    stop("max_cores must be a positive integer!") # End
  }

  names(object@null_models) <- paste(
    object@celltype_gene_test_pairs$celltype,
    object@celltype_gene_test_pairs$gene,
    sep = "_"
  )
  error <- sapply(
    object@null_models,
    function(x) {
      # inherits(x, "try-error") || !(x$converge)
      inherits(x, "try-error")
    }
  )
  if (any(error)) {
    object@null_models <- object@null_models[!error]
    object@celltype_gene_test_pairs <- object@celltype_gene_test_pairs[!error, , drop = FALSE]
  }

  return(object)
}

#' Fit the null MMM model
#'
#' @description Fit the null MMM model for a single celltype-gene pair using a PQL-based approach.
#'
#' @param Y A numeric vector containing gene expression counts of all spots.
#' @param library_size A numeric vector containing library sizes of all spots.
#' @param X A numeric matrix containing covariates of all spots.
#' Each row represents a spot and each column represents a covariate.
#' If `NULL`, a matrix with one column of all 1s will be used as the covariate matrix. Default is `NULL`.
#' @param proportion A numeric matrix containing cell type proportions of all spots.
#' Each row represents a spot and each column represents a cell type.
#' @param batch_id A character/factor vector indicating which batch each spot comes from.
#' It's applicable to the case of multiple samples/replicates/slices and specific gene platform effects required.
#' If `NULL`, all spots will be assumed to come from the same batch and share the same gene platform effects. Default is `NULL`.
#' @param platform_effect A numeric value or vector.
#' In the single batch case, a numeric value is expected.
#' When in the case of multiple batches and specific platform effects required, a vector is expected with each element corresponding to a batch.
#' If `NULL`, the platform effect will be estimated from data with zero initialization. Default is `NULL`.
#' @param celltype_test A character specifying the cell type to test after fitting the null model.
#' @param proportion_threshold A numeric value between 0 and 1.
#' The minimum proportion of a cell type at a spot to be considered.
#' @param reference_threshold A numeric value between 0 and 1.
#' The minimum relative gene expression level in a cell type to be considered for fitting the null model.
#' @param safeguard A numeric value.
#' A small positive number to avoid numerical issues when computing the inverse of a matrix.
#' @param iter_max A positive integer.
#' The maximum number of iterations for the optimization algorithm.
#' @param tol A numeric value.
#' The convergence criteria for the optimization algorithm.
#' @param verbose A logical value.
#' Whether to print the progress of the optimization algorithm. Default is TRUE.
#'
#' @return A list containing the fitted null model results.
#'
#' @export
mcubeFitNullSinglePair <- function(
    Y, library_size, X = NULL, proportion, batch_id = NULL,
    reference, used_for_deconvolution = TRUE,
    spot_effects = NULL, platform_effect = NULL,
    celltype_test, proportion_threshold = 0.1, reference_threshold = 0.25,
    safeguard = 1e-6, iter_max = 100, tol = 1e-6, verbose = TRUE) {
  spot_names <- rownames(proportion)
  celltype_names <- colnames(proportion)

  spots_filter <- which(proportion[, celltype_test] >= proportion_threshold)
  if (length(spots_filter) == 0) {
    stop(
      "mcubeFitNullSinglePair: No spot has proportion >= ",
      proportion_threshold, " for the cell type to test!"
    ) # End
  } else if (length(spots_filter) < length(spot_names)) {
    Y <- Y[spots_filter]
    library_size <- library_size[spots_filter]
    X <- X[spots_filter, , drop = FALSE]
    proportion <- proportion[spots_filter, , drop = FALSE]
    if (!is.null(batch_id)) {
      batch_id <- batch_id[spots_filter]
    }
    if (!is.null(spot_effects)) {
      spot_effects <- spot_effects[spots_filter]
    }
    spot_names <- spot_names[spots_filter]
  }
  num_spots <- length(spot_names)
  num_covariates <- ncol(X) # Inclding intercept

  celltype_minor <- which(reference / max(reference) < reference_threshold)
  if (length(celltype_minor) > 0) {
    if (celltype_test %in% celltype_names[celltype_minor]) {
      stop(
        "mcubeFitNullSinglePair: The gene to test is lowly expressed in the cell type to test!"
      ) # End
    }
    reference_minor <- reference[celltype_minor]
    proportion_minor <- proportion[, celltype_minor, drop = FALSE]
    Y_mean_minor_vec <- library_size *
      as.vector(proportion_minor %*% reference_minor)

    celltype_names <- celltype_names[-celltype_minor]
    reference <- reference[-celltype_minor]
    proportion <- proportion[, -celltype_minor, drop = FALSE]
  } else {
    Y_mean_minor_vec <- rep(0, num_spots)
  }
  num_celltypes <- length(celltype_names)

  membership_mat <- sweep(
    proportion,
    MARGIN = 2,
    STATS = reference,
    FUN = "*"
  )
  membership_mat <- membership_mat / rowSums(membership_mat)
  MMT_vec <- rowSums(membership_mat^2) # M * M^T, diagnoal martix

  if ((!used_for_deconvolution) || is.null(spot_effects)) {
    spot_effects <- rep(0, num_spots)
  }
  if (is.null(platform_effect)) {
    spot_platform_effects <- spot_effects
  } else if (length(levels(batch_id)) > 1) {
    spot_platform_effects <- spot_effects +
      as.vector(model.matrix(~ batch_id - 1) %*% platform_effect)
  } else {
    spot_platform_effects <- spot_effects + platform_effect
  }
  Y_mean_minor_vec <- exp(spot_platform_effects) * Y_mean_minor_vec

  # Cell type level intercept: log reference
  log_reference <- matrix(
    log(reference),
    nrow = num_spots, ncol = num_celltypes,
    byrow = TRUE
  )

  # Initialization
  tau <- 1e-5 # all cell types share same tau
  xi <- rep(0, num_covariates) # P * K
  u_mat <- matrix(0, nrow = num_spots, ncol = num_celltypes) # I * K

  for (step in 1:iter_max) {
    ### Step 1: Compute Y_tilde with current estimate of eta(xi, u)

    # Matrix of spot * cell type
    eta_mat <- log_reference + u_mat
    Y_mean_derivative_mat <- library_size *
      exp(as.vector(X %*% xi) + spot_platform_effects) *
      (proportion * exp(eta_mat))
    Y_mean_vec <- rowSums(Y_mean_derivative_mat) +
      exp(as.vector(X %*% xi)) * Y_mean_minor_vec

    # membership_mat <- Y_mean_derivative_mat / Y_mean_vec
    score_vec <- Y - Y_mean_vec
    W_vec <- Y_mean_vec
    W_inv_vec <- 1 / W_vec

    ## Compute Y_tilde
    Y_tilde <- W_inv_vec * score_vec +
      as.vector(X %*% xi) + rowSums(membership_mat * u_mat)

    ### Step 2: Variance components estimation, update tau using AI algorithm

    # Compute noise matrix of IK * IK
    # all cell types share the same tau
    # MMT_vec <- rowSums(membership_mat^2)
    Sigma_vec <- W_inv_vec + tau * MMT_vec
    Sigma_inv_vec <- 1 / (Sigma_vec + safeguard)
    Sigma_inv_X_mat <- Sigma_inv_vec * X
    X_Sigma_inv_X_mat <- crossprod(X, Sigma_inv_X_mat)

    P_mat <- diag(Sigma_inv_vec) -
      tcrossprod(
        Sigma_inv_X_mat %*% chol2inv(chol(X_Sigma_inv_X_mat)),
        Sigma_inv_X_mat
      )
    P_Y_tilde <- as.vector(P_mat %*% Y_tilde)

    # partial Sigma / partial tau = M * M^T, diagnoal matrix
    deriv_first_vec <- sum(MMT_vec * P_Y_tilde^2) / 2 -
      sum(diag(P_mat * MMT_vec)) / 2

    # Compute second order derivative
    par_Sigma_par_tau_P_Y_tilde <- MMT_vec * P_Y_tilde
    deriv_sec <- as.vector(t(par_Sigma_par_tau_P_Y_tilde) %*% P_mat %*%
                             par_Sigma_par_tau_P_Y_tilde) / 2

    # Update tau
    tau_new <- tau + deriv_first_vec / deriv_sec
    step_size <- 1.0
    while (tau_new <= 0.0 || tau_new >= 5) {
      step_size <- step_size * 0.5
      tau_new <- tau +
        step_size * deriv_first_vec / deriv_sec
    }

    ### Step 3: Update xi and u using current Y_tilde and new tau

    Sigma_vec <- W_inv_vec + tau_new * MMT_vec
    Sigma_inv_vec <- 1 / (Sigma_vec + safeguard)

    # Boost computation through matrix multiplication tricks
    Sigma_inv_X_mat <- Sigma_inv_vec * X
    X_Sigma_inv_X_mat <- crossprod(X, Sigma_inv_X_mat)
    Sigma_inv_Y_vec <- Sigma_inv_vec * Y_tilde
    X_Sigma_inv_Y_mat <- crossprod(X, Sigma_inv_Y_vec)

    xi_new <- as.vector(
      chol2inv(chol(X_Sigma_inv_X_mat)) %*% X_Sigma_inv_Y_mat
    )
    u_mat_new <- tau_new * Sigma_inv_vec * membership_mat *
      as.vector(Y_tilde - X %*% xi_new)

    gap <- max(
      c(
        abs(as.vector(xi) - xi_new) /
          (abs(as.vector(xi)) + abs(xi_new)),
        abs(as.vector(tau - tau_new)) /
          (abs(as.vector(tau)) + abs(tau_new))
      )
    )

    if (verbose) {
      message(
        "Iteration: ", step,
        ", intercept: ", xi_new[1], ", tau: ", tau_new[1],
        ", gap: ", gap, "."
      )
    }
    if (gap < tol) {
      break
    }

    tau <- tau_new
    xi <- xi_new
    u_mat <- u_mat_new
  }

  null_model_results <- list(
    tau = tau_new, xi = xi_new, u = u_mat_new,
    Y_tilde = Y_tilde, membership = membership_mat,
    W = W_vec, Sigma_inv = Sigma_inv_vec,
    spots = spot_names, celltypes = celltype_names
  )
  names(null_model_results$Y_tilde) <-
    rownames(null_model_results$membership) <-
    names(null_model_results$W) <-
    names(null_model_results$Sigma_inv) <-
    rownames(null_model_results$u) <- spot_names
  colnames(null_model_results$membership) <-
    colnames(null_model_results$u) <- celltype_names
  null_model_results$converge <- ifelse(gap < tol, TRUE, FALSE)
  return(null_model_results)
}

# # A sparse version of the `mcubeFitNullSinglePair` function.
# # Suitable for the case of a sparse cell type proportion matrix.
# mcubeFitNullSinglePair <- function(
    #     Y, library_size, X = NULL, proportion, batch_id = NULL,
#     reference, used_for_deconvolution = TRUE,
#     spot_effects = NULL, platform_effect = NULL,
#     celltype_test, proportion_threshold = 0.1, reference_threshold = 0.25,
#     safeguard = 1e-6, iter_max = 100, tol = 1e-6, verbose = TRUE) {
#   spot_names <- rownames(proportion)
#   celltype_names <- colnames(proportion)

#   spots_filter <- which(proportion[, celltype_test] >= proportion_threshold)
#   if (length(spots_filter) == 0) {
#     stop(
#       "mcubeFitNullSinglePair: No spot has proportion >= ",
#       proportion_threshold, " for the cell type to test!"
#     ) # End
#   } else if (length(spots_filter) < length(spot_names)) {
#     Y <- Y[spots_filter]
#     library_size <- library_size[spots_filter]
#     X <- X[spots_filter, , drop = FALSE]
#     proportion <- proportion[spots_filter, , drop = FALSE]
#     if (!is.null(batch_id)) {
#       batch_id <- batch_id[spots_filter]
#     }
#     if (!is.null(spot_effects)) {
#       spot_effects <- spot_effects[spots_filter]
#     }
#     spot_names <- spot_names[spots_filter]
#   }
#   num_spots <- length(spot_names)
#   num_covariates <- ncol(X) # Inclding intercept

#   celltype_minor <- which(reference / max(reference) < reference_threshold)
#   if (length(celltype_minor) > 0) {
#     if (celltype_test %in% celltype_names[celltype_minor]) {
#       stop(
#         "mcubeFitNullSinglePair: The gene to test is lowly expressed in the cell type to test!"
#       ) # End
#     }
#     reference_minor <- reference[celltype_minor]
#     proportion_minor <- proportion[, celltype_minor, drop = FALSE]
#     Y_mean_minor_vec <- library_size *
#       as.vector(proportion_minor %*% reference_minor)

#     celltype_names <- celltype_names[-celltype_minor]
#     reference <- reference[-celltype_minor]
#     proportion <- proportion[, -celltype_minor, drop = FALSE]
#   } else {
#     Y_mean_minor_vec <- rep(0, num_spots)
#   }
#   num_celltypes <- length(celltype_names)

#   # "row" means the spot, "col" means the cell-type
#   not_zero_proportion <- which(proportion != 0, arr.ind = TRUE)

#   membership_mat <- sweep(
#     proportion,
#     MARGIN = 2,
#     STATS = reference,
#     FUN = "*"
#   )
#   membership_mat <- membership_mat / rowSums(membership_mat)
#   membership_mat <- Matrix::sparseMatrix(
#     i = not_zero_proportion[, "row"],
#     j = not_zero_proportion[, "col"],
#     x = membership_mat[not_zero_proportion],
#     dims = c(num_spots, num_celltypes)
#   )
#   MMT_vec <- rowSums(membership_mat^2) # M * M^T, diagnoal martix

#   if ((!used_for_deconvolution) || is.null(spot_effects)) {
#     spot_effects <- rep(0, num_spots)
#   }
#   if (is.null(platform_effect)) {
#     spot_platform_effects <- spot_effects
#   } else if (length(levels(batch_id)) > 1) {
#     spot_platform_effects <- spot_effects +
#       as.vector(model.matrix(~ batch_id - 1) %*% platform_effect)
#   } else {
#     spot_platform_effects <- spot_effects + platform_effect
#   }
#   Y_mean_minor_vec <- exp(spot_platform_effects) * Y_mean_minor_vec

#   # Cell type level intercept: log reference
#   log_reference <- Matrix::sparseMatrix(
#     i = not_zero_proportion[, "row"],
#     j = not_zero_proportion[, "col"],
#     x = (matrix(
#       log(reference),
#       nrow = num_spots, ncol = num_celltypes,
#       byrow = TRUE
#     ))[not_zero_proportion],
#     dims = c(num_spots, num_celltypes)
#   )

#   # Initialization
#   tau <- 1e-5 # all cell types share same tau
#   xi <- rep(0, num_covariates) # P * K
#   u <- Matrix::sparseMatrix(
#     i = not_zero_proportion[, "row"],
#     j = not_zero_proportion[, "col"],
#     x = 0,
#     dims = c(num_spots, num_celltypes)
#   ) # I * K

#   for (step in 1:iter_max) {
#     ### Step 1: Compute Y_tilde with current estimate of eta(xi, u)

#     # Matrix of spot * cell type
#     eta_mat <- log_reference + u
#     Y_mean_derivative_mat <- library_size *
#       exp(as.vector(X %*% xi) + spot_platform_effects) *
#       (proportion * exp(eta_mat))
#     Y_mean_vec <- rowSums(Y_mean_derivative_mat) +
#       exp(as.vector(X %*% xi)) * Y_mean_minor_vec

#     # membership_mat <- Y_mean_derivative_mat / Y_mean_vec
#     score_vec <- Y - Y_mean_vec
#     W_vec <- Y_mean_vec
#     W_inv_vec <- 1 / W_vec

#     ## Compute Y_tilde
#     Y_tilde <- W_inv_vec * score_vec +
#       as.vector(X %*% xi) + rowSums(membership_mat * u)

#     ### Step 2: Variance components estimation, update tau using AI algorithm

#     # Compute noise matrix of IK * IK
#     # all cell types share the same tau
#     # MMT_vec <- rowSums(membership_mat^2)
#     Sigma_vec <- W_inv_vec + tau * MMT_vec
#     Sigma_inv_vec <- 1 / (Sigma_vec + safeguard)
#     Sigma_inv_X_mat <- Sigma_inv_vec * X
#     X_Sigma_inv_X_mat <- crossprod(X, Sigma_inv_X_mat)

#     P_mat <- diag(Sigma_inv_vec) -
#       tcrossprod(
#         Sigma_inv_X_mat %*% chol2inv(chol(X_Sigma_inv_X_mat)),
#         Sigma_inv_X_mat
#       )
#     P_Y_tilde <- as.vector(P_mat %*% Y_tilde)

#     # partial Sigma / partial tau = M * M^T, diagnoal matrix
#     deriv_first_vec <- sum(MMT_vec * P_Y_tilde^2) / 2 -
#       sum(diag(P_mat * MMT_vec)) / 2

#     # Compute second order derivative
#     par_Sigma_par_tau_P_Y_tilde <- MMT_vec * P_Y_tilde
#     deriv_sec <- as.vector(t(par_Sigma_par_tau_P_Y_tilde) %*% P_mat %*%
#       par_Sigma_par_tau_P_Y_tilde) / 2

#     # Update tau
#     tau_new <- tau + deriv_first_vec / deriv_sec
#     step_size <- 1.0
#     while (tau_new <= 0.0 || tau_new >= 5) {
#       step_size <- step_size * 0.5
#       tau_new <- tau +
#         step_size * deriv_first_vec / deriv_sec
#     }

#     ### Step 3: Update xi and u using current Y_tilde and new tau

#     Sigma_vec <- W_inv_vec + tau_new * MMT_vec
#     Sigma_inv_vec <- 1 / (Sigma_vec + safeguard)

#     # Modifiy to boost computation
#     Sigma_inv_X_mat <- Sigma_inv_vec * X
#     X_Sigma_inv_X_mat <- crossprod(X, Sigma_inv_X_mat)
#     Sigma_inv_Y_vec <- Sigma_inv_vec * Y_tilde
#     X_Sigma_inv_Y_mat <- crossprod(X, Sigma_inv_Y_vec)

#     xi_new <- as.vector(
#       chol2inv(chol(X_Sigma_inv_X_mat)) %*% X_Sigma_inv_Y_mat
#     )
#     u_new <- tau_new * Sigma_inv_vec * membership_mat *
#       as.vector(Y_tilde - X %*% xi_new)
#     u_new <- Matrix::sparseMatrix(
#       i = not_zero_proportion[, "row"],
#       j = not_zero_proportion[, "col"],
#       x = u_new[not_zero_proportion],
#       dims = c(num_spots, num_celltypes)
#     )

#     gap <- max(
#       c(
#         abs(as.vector(xi) - xi_new) /
#           (abs(as.vector(xi)) + abs(xi_new)),
#         abs(as.vector(tau - tau_new)) /
#           (abs(as.vector(tau)) + abs(tau_new))
#       )
#     )

#     if (verbose) {
#       message(
#         "Iteration: ", step,
#         ", intercept: ", xi_new[1], ", tau: ", tau_new[1],
#         ", gap: ", gap, "."
#       )
#     }
#     if (gap < tol) {
#       break
#     }

#     tau <- tau_new
#     xi <- xi_new
#     u <- u_new
#   }

#   null_model_results <- list(
#     tau = tau_new, xi = xi_new, u = u_new,
#     Y_tilde = Y_tilde, membership = membership_mat,
#     W = W_vec, Sigma_inv = Sigma_inv_vec,
#     spots = spot_names, celltypes = celltype_names
#   )
#   names(null_model_results$Y_tilde) <-
#     rownames(null_model_results$membership) <-
#     names(null_model_results$W) <-
#     names(null_model_results$Sigma_inv) <-
#     rownames(null_model_results$u) <- spot_names
#   colnames(null_model_results$membership) <-
#     colnames(null_model_results$u) <- celltype_names
#   null_model_results$converge <- ifelse(gap < tol, TRUE, FALSE)
#   return(null_model_results)
# }
