mcubeKernel <- function(
    coordinates, standardize = TRUE,
    kernel_type = "Gaussian", length_scale = NULL, sparse = NULL) {
    if (kernel_type == "linear") {
        if (standardize) {
            coordinates <- scale(coordinates, center = TRUE, scale = TRUE)
        }
        kernel_mat <- tcrossprod(coordinates)
    } else {
        if (standardize) {
            coordinates <- scale(coordinates, center = TRUE, scale = FALSE)
            # coordinates <- coordinates / sqrt(median(rowSums(coordinates^2)))
            coordinates <- coordinates / sqrt(median(rowSums(coordinates^2)))
        }

        # Set the length scale
        if (is.null(coordinates) && is.null(length_scale)) {
            stop("mcubeKernel: length scale is either estimated from coordinates or provided by the user!") # End
        } else if (is.null(length_scale)) {
            length_scale <- mcubeLengthScale(coordinates)
        } else if (length_scale <= 0) {
            stop("mcubeKernel: length_scale must be positive!") # End
        }
        message("mcubeKernel: length scale is set as ", length_scale, " for the ", kernel_type, " kernel.")

        if (kernel_type == "Gaussian") {
            kernel_mat <- exp(
                -(as.matrix(dist(coordinates)))^2 / (2 * length_scale^2)
            )
        } else if (kernel_type == "Cauchy") {
            kernel_mat <- 1 /
                (1 + (as.matrix(dist(coordinates)))^2 / (2 * length_scale^2))
        } else if (kernel_type == "periodic") {
            kernel_mat <- exp(
                -sin(pi * (as.matrix(dist(coordinates))))^2 / (2 * length_scale^2)
            )
        } else if (kernel_type == "Gaussian_transformed") {
            coordinates <- exp(-coordinates^2 / (2 * length_scale^2))
            coordinates <- scale(coordinates, center = TRUE, scale = FALSE)
            kernel_mat <- coordinates %*% solve(crossprod(coordinates)) %*% t(coordinates)
        }
    }

    if (!is.null(sparse)) {
        kernel_mat[kernel_mat < sparse] <- 0
        kernel_mat <- as(kernel_mat, "sparseMatrix")
    }

    rownames(kernel_mat) <- colnames(kernel_mat) <- rownames(coordinates)

    return(kernel_mat)
}

mcubeLengthScale <- function(coordinates, standardize = TRUE) {
    if (standardize) {
        coordinates <- scale(coordinates, center = TRUE, scale = FALSE)
        coordinates <- coordinates / sqrt(median(rowSums(coordinates^2)))
    }
    dist <- as.matrix(dist(coordinates))
    dist_min <- apply(
        dist,
        MARGIN = 2,
        FUN = function(x) {
            sort(x, decreasing = FALSE)[2]
        }
    )
    length_scale <- median(dist_min)

    if (is.na(length_scale)) {
        warning("mcubeLengthScale: length_scale estimatet from the data is NA!")
        length_scale <- 0.1
        warning("mcubeLengthScale: length_scale will be set as ", length_scale, ".\n")
    }

    return(length_scale)
}

# mcubeKernel <- function(
#     counts = NULL, coordinates, standardize = TRUE,
#     kernel_type = "Gaussian", length_scale = NULL,
#     sparse = NULL) {
#     if (kernel_type == "linear") {
#         if (standardize) {
#             coordinates <- scale(coordinates, center = TRUE, scale = TRUE)
#         }
#         kernel_mat <- tcrossprod(coordinates)
#     } else {
#         # Set the length scale
#         if (is.null(counts) && is.null(length_scale)) {
#             stop("mcubeKernel: length scale is either estimated from counts or provided by the user!") # End
#         } else if (is.null(length_scale)) {
#             length_scale <- mcubeLengthScale(counts)
#         } else if (length_scale <= 0) {
#             stop("mcubeKernel: length_scale must be positive!") # End
#         }
#         message("mcubeKernel: length scale is set as ", length_scale, " for the ", kernel_type, " kernel.")

#         if (standardize) {
#             coordinates <- scale(coordinates, center = TRUE, scale = FALSE)
#             coordinates <- coordinates / sqrt(median(rowSums(coordinates^2)))
#         }

#         if (kernel_type == "Gaussian") {
#             kernel_mat <- exp(
#                 -as.matrix(dist(coordinates)^2) / length_scale
#             )
#         } else if (kernel_type == "Cauchy") {
#             kernel_mat <- 1 /
#                 (1 + as.matrix(dist(coordinates)^2) / length_scale)
#         } else if (kernel_type == "periodic") {
#             kernel_mat <- exp(
#                 -sin(pi * as.matrix(dist(coordinates)))^2 / length_scale
#             )
#         }
#     }

#     if (!is.null(sparse)) {
#         kernel_mat[kernel_mat < sparse] <- 0
#     }

#     rownames(kernel_mat) <- colnames(kernel_mat) <- rownames(coordinates)

#     return(kernel_mat)
# }

# mcubeLengthScale <- function(counts, method = NULL) {
#     n_spots <- nrow(counts)
#     # log_counts <- log(1e+3 * counts / rowSums(counts) + 1)
#     log_counts <- log(counts + 1)
#     if (is.null(method)) {
#         method <- ifelse(n_spots > 10000, "Silverman", "SJ")
#     }

#     if (method == "SJ") {
#         length_scale_SJ <- apply(
#             log_counts,
#             MARGIN = 2,
#             FUN = function(x) {
#                 tryCatch(
#                     {
#                         bw.SJ(x, method = "dpi")
#                     },
#                     error = function(e) {
#                         return(NA)
#                     }
#                 )
#             }
#         )
#         length_scale <- median(length_scale_SJ, na.rm = TRUE)
#     } else if (method == "Silverman") {
#         length_scale_Silverman <- apply(
#             log_counts,
#             MARGIN = 2,
#             FUN = function(x) {
#                 tryCatch(
#                     {
#                         bw.nrd0(x)
#                     },
#                     error = function(e) {
#                         return(NA)
#                     }
#                 )
#             }
#         )
#         length_scale <- median(length_scale_Silverman, na.rm = TRUE)
#     }

#     if (is.na(length_scale)) {
#         warning("mcubeLengthScale: length_scale estimatet from the data is NA!\n")
#         warning("mcubeLengthScale: length_scale will be set as 0.1!\n")
#         length_scale <- 0.1
#     }

#     return(length_scale)
# }
