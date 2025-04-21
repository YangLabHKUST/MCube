#' An S4 class to represent spatial transcriptomic and single-cell RNA-sequencing reference data
#'
#' @description Each \code{\linkS4class{mcube}} object has a number of slots which store information. Key slots to access are listed below.
#'
#' @importFrom methods setClass
#'
#' @slot counts A matrix, Matrix::dgCMatrix, or Matrix::dgTMatrix.
#' Each row represents a spot and each column represents a gene.
#' @slot coordinates A matrix.
#' Each row represents a spot and each column represents a spatial dimension.
#' @slot proportions A cell type proportion matrix.
#' Each row represents a spot and each column represents a cell type.
#' @slot library_sizes A numeric vector containing the library sizes of all spots.
#' @slot covariates A matrix.
#' Each row represents a spot and each column represents a covariate.
#' @slot batch_id A factor indicating which batch each spot comes from.
#' @slot spots A character vector containing the names of spots to analyze.
#' @slot reference A matrix of average gene expression calculated from scRNA-seq reference data.
#' Each row represents a cell type and each column represents a gene.
#' @slot used_for_deconvolution A logical vector indicating whether each gene has been used for cell type deconvolution.
#' @slot spot_effects A numeric vector containing the spot effects of all spots.
#' @slot platform_effects A matrix containing the gene-specific platform effects.
#' Each row represents a batch and each column represents a gene.
#' @slot config A list recording the configuration of the \code{\linkS4class{mcube}} object.
#' @slot null_models A listing recording the fitted results of the null models with each element corresponding to a celltype-gene pair.
#' @slot kernels A list recording the kernel matrices used for score testing.
#' @slot celltype_test A character vector specifying the cell types to test.
#' @slot gene_test A character vector specifying the genes to test.
#' @slot celltype_gene_test_pairs A data.frame specifying the celltype-gene pairs to test.
#' @slot pvalues A list recording the p-values from the SVG test.
#' Each element is a data.frame corresponding to a cell type.
#' @slot project A character recording the name of the project.
#'
#' @export
setClass(
    # Set the name for the class
    "mcube",

    # Define the slots
    slots = c(
        counts = "ANY",
        coordinates = "matrix",
        proportions = "matrix",
        library_sizes = "numeric",
        covariates = "matrix",
        batch_id = "factor",
        spots = "character",
        reference = "matrix",
        used_for_deconvolution = "logical",
        spot_effects = "numeric",
        platform_effects = "matrix",
        config = "list",
        null_models = "list",
        kernels = "list",
        celltype_test = "character",
        gene_test = "character",
        celltype_gene_test_pairs = "data.frame",
        pvalues = "list",
        project = "character"
    ),

    # Assign the default prototypes
    prototype = list(
        project = "MCube"
    )
)

#' Create the \code{\linkS4class{mcube}} object
#'
#' @importFrom methods new
#'
#' @param counts A matrix, Matrix::dgCMatrix, or Matrix::dgTMatrix.
#' Each Row represents a spot and each column represents a gene.
#' @param coordinates A matrix.
#' Each row represents a spot and each column represents a spatial dimension.
#' @param proportions A cell type proportion matrix.
#' Each row represents a spot and each column represents a cell type.
#' @param library_sizes A numeric vector containing the library sizes of all spots.
#' @param covariates A matrix.
#' Each row represents a spot and each column represents a covariate.
#' @param batch_id A character/factor vector indicating which batch each spot comes from.
#' It's applicable to the case of multiple samples/replicates/slices and specific gene platform effects required.
#' If `NULL`, all spots will be assumed to come from the same batch and share the same gene platform effects. Default is `NULL`.
#' @param reference A matrix of average gene expression calculated from scRNA-seq reference data.
#' Each row represents a cell type and each column represents a gene.
#' @param used_for_deconvolution A logical vector indicating whether each gene has been used for cell type deconvolution.
#' If `NULL`, all genes will be assumed have been used for cell type deconvolution. Default is `NULL`.
#' @param spots A character vector specifying the names of spots to analyze.
#' @param library_size_min A numeric value.
#' The minimum library size to filter out spots.
#' @param spot_effects A numeric vector containing the spot effects of all spots.
#' @param platform_effects A numeric vector or matrix containing the gene-specific platform effects.
#' In the single batch case, a numeric vector is expected with each element corresponding to a gene.
#' When in the case of multiple batches and specific platform effects required, a matrix is expected with rows corresponding to batches and columns corresponding to genes.
#' If `NULL`, the platform effects will be estimated from the data. Default is `NULL`.
#' @param celltype_test A character vector specifying the cell types to test.
#' If `NULL`, all cell types that pass the filtering criteria will be included. Default is `NULL`.
#' @param gene_test A character vector specifying the genes to test.
#' If `NULL`, all genes that pass the filtering criteria will be included. Default is `NULL`.
#' @param celltype_threshold A numeric value.
#' The minimum proportion sum across spots of a cell type to be considered.
#' @param gene_threshold A numeric value.
#' The minimum average expression level across spots of a gene to be considered.
#' @param proportion_threshold A numeric value between 0 and 1.
#' The minimum proportion of a cell type to be considered at a spot.
#' @param reference_threshold A numeric value between 0 and 1.
#' The minimum relative expression level of a gene to be considered for a cell type in testing.
#' @param project A character recording the name of the project.
#'
#' @return An \code{\linkS4class{mcube}} object.
#'
#' @export
createMCube <- function(
    counts, coordinates, proportions,
    library_sizes = NULL, covariates = NULL, batch_id = NULL,
    reference, used_for_deconvolution = NULL,
    spots = NULL, library_size_min = 10,
    spot_effects = NULL, platform_effects = NULL,
    celltype_test = NULL, gene_test = NULL,
    celltype_threshold = 100, gene_threshold = 5e-5,
    proportion_threshold = 0.1, reference_threshold = 0.5,
    project = "MCube") {
    # Check spot names
    if (!identical(rownames(counts), rownames(coordinates))) {
        stop("Spot names of counts and coordinates do not match!") # End
    }
    if (!identical(rownames(counts), rownames(proportions))) {
        stop("Spot names of counts and proportions do not match!") # End
    }
    if (!is.null(library_sizes) &&
        !identical(rownames(counts), names(library_sizes))) {
        stop("Spot names of counts and library_sizes do not match!") # End
    } else if (is.null(library_sizes)) {
        library_sizes <- rowSums(as.matrix(counts))
        names(library_sizes) <- rownames(counts)
    }
    if (!is.null(covariates) &&
        !identical(rownames(counts), rownames(covariates))) {
        stop("Spot names of counts and covariates do not match!") # End
    }
    if (!is.null(batch_id) &&
        !identical(rownames(counts), names(batch_id))) {
        stop("Spot names of counts and batch_id do not match!") # End
    } else if (!is.null(batch_id)) {
        batch_id <- as.factor(batch_id)
    } else {
        message(
            "The batch_id is not provided!\n",
            "All spots are assumed to be from the same batch ",
            "and share the same gene platform effects."
        )
        batch_id <- factor(rep("batch_1", nrow(counts)))
        names(batch_id) <- rownames(counts)
    }
    if (!is.null(spot_effects) &&
        !identical(rownames(counts), names(spot_effects))) {
        stop("Spot names of counts and spot_effects do not match!") # End
    } else if (is.null(spot_effects)) {
        spot_effects <- rep(0, nrow(counts))
        names(spot_effects) <- rownames(counts)
    }

    # Check covariates (including intercept term)
    if (is.null(covariates)) {
        covariates <- matrix(1, nrow = nrow(counts), ncol = 1)
        rownames(covariates) <- rownames(counts)
    } else {
        all_ones_indices <- apply(
            covariates,
            MARGIN = 2,
            function(col) {
                all(col == 1)
            }
        )
        # Add intercept term
        if (!any(all_ones_indices)) {
            covariates <- cbind(1, covariates)
        }
    }

    # Check spots to analyze
    if (is.null(spots)) {
        # Filter out spots with low library size
        spots <- rownames(counts)[library_sizes > library_size_min]
    } else {
        spots <- intersect(
            rownames(counts)[library_sizes > library_size_min], spots
        )
        if (length(spots) == 0) {
            stop("The spots to analyze do not match the input data!") # End
        }
    }
    counts <- counts[spots, , drop = FALSE]
    coordinates <- coordinates[spots, , drop = FALSE]
    proportions <- proportions[spots, , drop = FALSE]
    library_sizes <- library_sizes[spots]
    covariates <- covariates[spots, , drop = FALSE]
    batch_id <- batch_id[spots]
    spot_effects <- spot_effects[spots]

    # Check cell type names
    if (!identical(colnames(proportions), rownames(reference))) {
        stop("Cell type names of proportions and reference do not match!")
    } # End

    # Check/assign celltype_test
    if (is.null(celltype_test)) {
        celltype_test <- rownames(reference)
    } else {
        celltype_test <- intersect(rownames(reference), celltype_test)
        if (length(celltype_test) == 0) {
            stop("The cell types in celltype_test do not match the input data!") # End
        }
    }
    celltype_test <- mcubeFilterCellTypes(
        proportions = proportions,
        celltype_test = celltype_test,
        proportion_threshold = proportion_threshold,
        celltype_threshold = celltype_threshold
    )

    # Take the intersection of genes in counts, reference, and (set up) gene_test
    if (is.null(gene_test)) {
        gene_test <- intersect(colnames(counts), colnames(reference))
    } else {
        gene_test <- intersect(
            intersect(colnames(counts), colnames(reference)),
            gene_test
        )
    }
    if (length(gene_test) == 0) {
        stop("No common genes between counts, reference, and gene_test!") # End
    }
    counts <- counts[, gene_test, drop = FALSE]

    # Filter out lowly expressed genes
    gene_test <- mcubeFilterGenes(as.matrix(counts), library_sizes, gene_threshold)
    if (length(gene_test) == 0) {
        stop("No genes remain after filtering based on expression level!") # End
    }
    counts <- counts[, gene_test, drop = FALSE]
    reference <- reference[, gene_test, drop = FALSE]

    # Check/calculate platform effects
    if (is.null(platform_effects)) {
        message("The platform effects are not provided and need to be estimated from data!")
        # Platform effects matrix
        # Each row is a batch and each column is a gene
        platform_effects <- t(sapply(
            levels(batch_id),
            FUN = function(x) {
                spots_used <- intersect(spots, names(batch_id)[batch_id == x])
                if (length(spots_used) == 0) {
                    return(rep(0, length(gene_test)))
                } else {
                    mcubeGetPlatformEffects(
                        counts = as.matrix(counts[spots_used, , drop = FALSE]),
                        library_sizes = library_sizes[spots_used],
                        proportions = proportions[spots_used, , drop = FALSE],
                        reference = reference,
                        spot_effects = spot_effects[spots_used]
                    )
                }
            }
        ))
        rownames(platform_effects) <- levels(batch_id)
        colnames(platform_effects) <- gene_test
    } else {
        if (is.vector(platform_effects)) {
            if (length(levels(batch_id)) > 1) {
                stop(
                    "The platform effects are provided as a vector for one batch, ",
                    "but multiple batches are detected!"
                ) # End
            }
            platform_effects <- t(as.matrix(platform_effects))
            rownames(platform_effects) <- levels(batch_id)
        } else {
            if (!all(levels(batch_id) %in% rownames(platform_effects))) {
                stop("The row names of platform_effects do not match the batch_id!") # End
            }
            platform_effects <- platform_effects[levels(batch_id), , drop = FALSE]
        }

        if (!all(gene_test %in% colnames(platform_effects))) {
            stop("The genes contained in platform_effects do not match the counts/reference data!") # End
        }
        platform_effects <- platform_effects[, gene_test, drop = FALSE]
    }

    # Get the platform effects by considering all spots together
    platform_effects_all <- mcubeGetPlatformEffects(
        counts = as.matrix(counts),
        library_sizes = library_sizes,
        proportions = proportions,
        reference = reference,
        spot_effects = spot_effects
    )

    # Get the gene list to test for each cell type
    gene_test_each_celltype_list <- lapply(
        celltype_test,
        FUN = function(celltype) {
            mcubeFilterGenesCellType(
                celltype = celltype,
                celltype_all = celltype_test,
                gene_test = gene_test,
                library_sizes = library_sizes,
                proportions = proportions,
                reference = reference,
                reference_threshold = reference_threshold,
                platform_effects = platform_effects_all
            )
        }
    )
    # Delete the cell type with no genes can be tested after filtering
    empty_celltype <- sapply(gene_test_each_celltype_list, length) == 0
    if (all(empty_celltype)) {
        stop("No genes can be tested for any cell type in celltype_test!") # End
    }
    gene_test_each_celltype_list <-
        gene_test_each_celltype_list[!empty_celltype]
    names(gene_test_each_celltype_list) <-
        celltype_test <- celltype_test[!empty_celltype]

    # Record the spots that will actually be analyzed
    spots <- spots[
        apply(
            proportions[, celltype_test, drop = FALSE],
            MARGIN = 1,
            FUN = function(row) {
                any(row >= proportion_threshold)
            }
        )
    ]

    # Update gene_test and crop the data
    gene_test <- gene_test[
        match(unique(unlist(gene_test_each_celltype_list)), gene_test)
    ]
    if (!is.null(used_for_deconvolution)) {
        used_for_deconvolution <- gene_test %in% used_for_deconvolution
    } else {
        used_for_deconvolution <- rep(TRUE, length(gene_test))
    }
    names(used_for_deconvolution) <- gene_test
    counts <- counts[, gene_test, drop = FALSE]
    reference <- reference[, gene_test, drop = FALSE]
    platform_effects <- platform_effects[, gene_test, drop = FALSE]

    celltype_gene_test_pairs <- do.call(
        "rbind",
        lapply(
            celltype_test,
            function(celltype) {
                data.frame(
                    celltype = celltype,
                    gene = gene_test_each_celltype_list[[celltype]]
                )
            }
        )
    )

    message(
        "Preprocessed data description: ",
        nrow(counts), " spots and ",
        ncol(proportions), " cell types in total. ",
        length(spots), " spots, ",
        length(gene_test), " genes, and ",
        length(celltype_test), " cell types to analyze."
    )

    # Inheriting
    object <- methods::new(
        Class = "mcube",
        counts = counts,
        coordinates = coordinates,
        proportions = proportions,
        library_sizes = library_sizes,
        covariates = covariates,
        batch_id = batch_id,
        spots = spots,
        reference = reference,
        used_for_deconvolution = used_for_deconvolution,
        spot_effects = spot_effects,
        platform_effects = platform_effects,
        config = list(
            celltype_threshold = celltype_threshold,
            gene_threshold = gene_threshold,
            proportion_threshold = proportion_threshold,
            reference_threshold_test = reference_threshold
        ),
        celltype_test = celltype_test,
        gene_test = gene_test,
        celltype_gene_test_pairs = celltype_gene_test_pairs,
        project = project
    )

    return(object)
}
