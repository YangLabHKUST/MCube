#' Construct the kernel matrix
#'
#' @importFrom stats median dist
#'
#' @param coordinates A matrix of spatial coordinates.
#' Each row represents a spot and each column represents a spatial dimension.
#' @param spots A character vector specifying the spots to include in the kernel matrix.
#' If `NULL`, all spots in the coordinates will be used. Default is `NULL`.
#' @param standardize A logical value.
#' If `TRUE`, the coordinates will be standardized. Default is `TRUE`.
#' @param kernel_type A character string specifying the kernel type.
#' It must be one of "linear", "Gaussian", "Cauchy", "periodic", or "Gaussian_transformed". Default is "Gaussian".
#' @param length_scale A numeric value specifying the length scale of the kernel.
#' If `NULL`, the length scale will be estimated from the coordinates. Default is `NULL`.
#'
#' @return A kernel matrix.
#'
#' @export
mcubeKernel <- function(
    coordinates, spots = NULL, standardize = TRUE,
    kernel_type = "Gaussian", length_scale = NULL) {
    if (is.null(spots)) {
        spots <- rownames(coordinates)
    } else {
        spots <- intersect(rownames(coordinates), spots)
        if (length(spots) == 0) {
            stop("mcubeKernel: No common spots between coordinates and spots!")
        }
    }

    if (kernel_type == "linear") {
        if (standardize) {
            coordinates <- scale(coordinates, center = TRUE, scale = TRUE)
        }
        kernel_mat <- tcrossprod(coordinates[spots, , drop = FALSE])
    } else {
        if (standardize) {
            coordinates <- scale(coordinates, center = TRUE, scale = FALSE)
            coordinates <- coordinates / sqrt(stats::median(rowSums(coordinates^2)))
        }

        # Set the length scale
        if (is.null(length_scale)) {
            length_scale <- mcubeLengthScale(coordinates)
        } else if (length_scale <= 0) {
            stop("mcubeKernel: length_scale must be positive!") # End
        }
        message(
            "mcubeKernel: length scale is set as ", length_scale,
            " for the ", kernel_type, " kernel."
        )

        coordinates <- coordinates[spots, , drop = FALSE]

        if (kernel_type == "Gaussian") {
            kernel_mat <- exp(
                -(as.matrix(stats::dist(coordinates)))^2 /
                    (2 * length_scale^2)
            )
        } else if (kernel_type == "Cauchy") {
            kernel_mat <- 1 /
                (1 + (as.matrix(stats::dist(coordinates)))^2 /
                    (2 * length_scale^2))
        } else if (kernel_type == "periodic") {
            kernel_mat <- exp(
                -sin(pi * (as.matrix(stats::dist(coordinates))))^2 /
                    (2 * length_scale^2)
            )
        } else if (kernel_type == "Gaussian_transformed") {
            coordinates <- exp(-coordinates^2 / (2 * length_scale^2))
            coordinates <- scale(coordinates, center = TRUE, scale = FALSE)
            kernel_mat <- coordinates %*% solve(crossprod(coordinates)) %*% t(coordinates)
        }
    }

    # if (!is.null(sparse)) {
    #     kernel_mat[kernel_mat < sparse] <- 0
    #     kernel_mat <- as(kernel_mat, "sparseMatrix")
    # }

    rownames(kernel_mat) <- colnames(kernel_mat) <- rownames(coordinates)

    return(kernel_mat)
}

#' Estimate the length scale of the kernel
#'
#' @importFrom stats dist median
#'
#' @param coordinates A matrix of spatial coordinates.
#' Each row represents a spot and each column represents a spatial dimension.
#' @param standardize A logical value.
#' If `TRUE`, the coordinates will be standardized. Default is `TRUE`.
#'
#' @return A numeric value specifying the length scale of the kernel.
#'
#' @export
mcubeLengthScale <- function(coordinates, standardize = TRUE) {
    if (standardize) {
        coordinates <- scale(coordinates, center = TRUE, scale = FALSE)
        coordinates <- coordinates / sqrt(stats::median(rowSums(coordinates^2)))
    }
    dist_mat <- as.matrix(stats::dist(coordinates))
    dist_min <- apply(
        dist_mat,
        MARGIN = 2,
        FUN = function(x) {
            sort(x, decreasing = FALSE)[2]
        }
    )
    length_scale <- stats::median(dist_min) * sqrt(6)

    if (is.na(length_scale)) {
        warning("mcubeLengthScale: length_scale estimatet from the data is NA!")
        length_scale <- 0.1
        warning("mcubeLengthScale: length_scale will be set as ", length_scale, ".\n")
    }

    return(length_scale)
}
