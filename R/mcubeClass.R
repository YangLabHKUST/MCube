#' Define the matrixORdgCMatrix class as the union of matrix and Matrix::dgCMatrix.
#'
#' @import Matrix
#'
setClassUnion(
    name = "matrixORdgCMatrix",
    members = c("matrix", "dgCMatrix", "dgTMatrix")
)

#' Each MCUBE object has a number of slots which store information. Key slots to access are listed below.
#'
#' @import Matrix
#'
#' @slot counts a matrix or Matrix::dgCMatrix. The main data matrix of N samples and M features.
#'
#' @return MCUBE class.
#' @export
setClass(
    # Set the name for the class
    "MCUBE",

    # Define the slots
    slots = c(
        counts = "matrixORdgCMatrix",
        coordinates = "matrix",
        proportion = "matrixORdgCMatrix",
        library_size = "numeric",
        covariates = "matrix",
        spots = "character",
        reference = "matrix",
        used_for_deconvolution = "logical",
        spot_effects = "numeric",
        platform_effects = "numeric",
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
        project = "MCUBE"
    )
)

#' Create the MCUBE object.
#'
#' @import Matrix
#' @importFrom methods new
#'
#' @param counts a matrix or Matrix::dgCMatrix.
#'
#' @return Returns MCUBE object.
#' @export
createMCUBE <- function(
    counts, coordinates, proportion,
    library_size = NULL, covariates = NULL, spots = NULL,
    reference, used_for_deconvolution = NULL,
    spot_effects = NULL, platform_effects = NULL,
    celltype_test = NULL, gene_test = NULL,
    celltype_threshold = 100, gene_threshold = 5e-5,
    proportion_threshold = 0.1, reference_threshold = 0.5,
    project = "MCUBE") {
    # Check sample names
    if (!identical(rownames(counts), rownames(coordinates))) {
        stop("Sample names of counts and coordinates do not match!") # End
    }
    if (!identical(rownames(counts), rownames(proportion))) {
        stop("Sample names of counts and proportion do not match!") # End
    }
    if (!is.null(library_size) &&
        !identical(rownames(counts), names(library_size))) {
        stop("Sample names of counts and library_size do not match!") # End
    } else if (is.null(library_size)) {
        library_size <- rowSums(counts)
        names(library_size) <- rownames(counts)
    }
    if (!is.null(covariates) &&
        !identical(rownames(counts), rownames(covariates))) {
        stop("Sample names of counts and covariates do not match!") # End
    }
    if (!is.null(spot_effects) &&
        !identical(rownames(counts), names(spot_effects))) {
        stop("Sample names of counts and spot_effects do not match!") # End
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
            function(col) all(col == 1)
        )
        if (any(all_ones_indices)) {
            covariates <- covariates[, !all_ones_indices, drop = FALSE]
        }
        covariates <- cbind(1, covariates)
    }

    # Check spots to analyze
    if (is.null(spots)) {
        spots <- rownames(counts) # Analyze all spots
    } else {
        spots <- intersect(rownames(counts), spots)
        if (length(spots) == 0) {
            stop("The spots to analyze do not match the input data!") # End
        }
    }

    # Check cell type names
    if (!identical(colnames(proportion), rownames(reference))) {
        stop("Celltype names of proportion and reference do not match!")
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
    celltype_test <- choose_celltypes(
        proportion = proportion[spots, , drop = FALSE],
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
    gene_test <- filter_genes(
        counts[spots, , drop = FALSE], library_size[spots], gene_threshold
    )
    if (length(gene_test) == 0) {
        stop("No genes remain after filtering based on expression level!") # End
    }
    counts <- counts[, gene_test, drop = FALSE]
    reference <- reference[, gene_test, drop = FALSE]

    # Check/calculate platform effects
    if (!is.null(platform_effects)) {
        if (all(gene_test %in% names(platform_effects))) {
            platform_effects <- platform_effects[gene_test]
        } else {
            stop("The platform_effects does not match the counts/reference data!") # End
        }
    } else {
        platform_effects <- get_platform_effects(
            counts = counts[spots, , drop = FALSE],
            library_size = library_size[spots],
            proportion = proportion[spots, , drop = FALSE],
            reference = reference,
            spot_effects = spot_effects[spots]
        )
    }

    # Get the gene list to test for each cell type
    gene_test_each_celltype_list <- lapply(
        celltype_test,
        FUN = function(celltype) {
            get_gene_list_celltype(
                celltype = celltype, celltype_all = celltype_test,
                gene_test = gene_test,
                library_size = library_size[spots],
                proportion = proportion[spots, , drop = FALSE],
                reference = reference, reference_threshold = reference_threshold,
                platform_effects = platform_effects
            )
        }
    )
    # Delete the cell type with no genes can be tested after filtering
    empty_celltype <- sapply(gene_test_each_celltype_list, length) == 0
    if (all(empty_celltype)) {
        stop("No genes can be tested for any cell type in celltype_test!") # End
    } else {
        gene_test_each_celltype_list <-
            gene_test_each_celltype_list[!empty_celltype]
        names(gene_test_each_celltype_list) <-
            celltype_test <- celltype_test[!empty_celltype]
    }
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
    platform_effects <- platform_effects[gene_test]

    celltype_gene_test_pairs <- do.call(
        rbind,
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
        nrow(counts), " spots, ", ncol(counts), " genes, and ",
        ncol(proportion), " cell types in total. ",
        length(spots), " spots and ",
        length(celltype_test), " cell types to analyze."
    )

    # Inheriting
    object <- methods::new(
        Class = "MCUBE",
        counts = counts,
        coordinates = coordinates,
        proportion = proportion,
        library_size = library_size,
        covariates = covariates,
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
