#' Plot the spatial distribution of gene expression
#'
#' @param counts A matrix of gene expression counts.
#' Each row represents a spot and each column represents a gene.
#' @param coordinates A matrix of spatial coordinates.
#' Each row represents a spot and each column represents a spatial dimension.
#' @param gene A character specifying the gene to plot.
#' @param spots A character vector specifying the spots to plot.
#' If `NULL`, all spots will be plotted. Default is `NULL`.
#' @param normalize A logical value.
#' If `TRUE`, the gene expression will be log-transformed and then normalized to [0, 1]. Default is `TRUE`.
#' @param he_image H&E image background.
#' A H&E image in TIFF format can be read into a raster array using `tiff::readTIFF`. Default is `NULL`.
#' @param background A logical value.
#' If `TRUE`, the spots not pass the expression threshold will be plotted in grey as background. Default is `TRUE`.
#' @param expr_threshold A numeric value specifying the expression threshold.
#' Spots with expression values greater than or equal to the threshold will be plotted.
#' Other spots will be used as background if `background` is `TRUE`. Default is 1.
#' @param spot_size A numeric value specifying the size of spots. Default is 1.
#' @param palettes A vector of color palettes. Default is `pals::brewer.blues(20)`.
#' @param opacity_target A numeric value specifying the opacity of target spots. Default is 1.
#' @param opacity_background A numeric value specifying the opacity of background spots. Default is 0.2.
#' @param xlim A numeric vector specifying the x-axis limits. Default is `NULL`.
#' @param ylim A numeric vector specifying the y-axis limits. Default is `NULL`.
#' @param ratio A numeric value specifying the aspect ratio. Default is 1.
#' @param title A character specifying the plot title. Default is `NULL`.
#'
#' @return A ggplot object.
#'
#' @export
mcubePlotExpr <- function(
    counts, coordinates, gene, spots = NULL, normalize = TRUE,
    he_image = NULL, background = TRUE, expr_threshold = 1,
    spot_size = 1, palettes = pals::brewer.blues(20),
    opacity_target = 1, opacity_background = 0.2,
    xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
  # Check sample names
  if (!all.equal(rownames(counts), rownames(coordinates))) {
    stop("mcubePlotExpr: rownames of counts and coordinates do not match!")
  }
  expr <- counts[, gene]
  if (is.null(spots)) {
    spots_target <- rownames(counts)[expr >= expr_threshold]
  } else {
    spots_target <- intersect(spots, rownames(counts)[expr >= expr_threshold])
  }
  if (normalize) {
    expr <- log(expr + 1)
    expr <- (expr - min(expr)) / (max(expr) - min(expr))
    value_lim <- c(0, 1)
  } else {
    value_lim <- c(0, max(expr))
  }
  expr <- pmax(pmin(expr, value_lim[2] - 1e-8), value_lim[1] + 1e-8)
  target_df <- data.frame(
    x = coordinates[spots_target, 1],
    y = coordinates[spots_target, 2],
    expr = expr[spots_target]
  )

  p <- ggplot2::ggplot(data = target_df, ggplot2::aes(x = x, y = y))

  # Axis ranges
  if (is.null(xlim)) {
    xlim <- c(min(coordinates[, 1]) - 1, max(coordinates[, 1]) + 1)
  }
  if (is.null(ylim)) {
    ylim <- c(min(coordinates[, 2]) - 1, max(coordinates[, 2]) + 1)
  }

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = xlim[1], xmax = xlim[2],
      ymin = -ylim[2], ymax = -ylim[1]
    )
  }
  # Other spots as background
  if (background) {
    spots_background <- setdiff(rownames(coordinates), spots_target)
    if (length(spots_background) > 0) {
      background_df <- data.frame(
        x = coordinates[spots_background, 1],
        y = coordinates[spots_background, 2]
      )
      p <- p + ggplot2::geom_point(
        data = background_df, ggplot2::aes(x = x, y = y),
        color = "gray", size = spot_size, alpha = opacity_background
      )
    }
  }

  p <- p + ggplot2::geom_point(
    ggplot2::aes(color = expr),
    size = spot_size, alpha = opacity_target
  ) +
    ggplot2::scale_colour_gradientn(
      name = "Level", colors = palettes,
      # values = scales::rescale(c(min(expr), median(expr[expr > 0]), max(expr)))
    ) +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(trans = "reverse", limits = rev(ylim)) +
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

#' Plot the spatial distribution of cell type proportions
#'
#' @param proportion A matrix of cell type proportions.
#' Each row represents a spot and each column represents a cell type.
#' @param coordinates A matrix of spatial coordinates.
#' Each row represents a spot and each column represents a spatial dimension.
#' @param celltype A character specifying the cell type to plot.
#' @param spots A character vector specifying the spots to plot.
#' If `NULL`, all spots will be plotted. Default is `NULL`.
#' @param he_image H&E image background.
#' A H&E image in TIFF format can be read into a raster array using `tiff::readTIFF`. Default is `NULL`.
#' @param background A logical value.
#' If `TRUE`, the spots with low proportions will be removed. Default is `TRUE`.
#' @param proportion_threshold A numeric value between 0 and 1.
#' Spots with proportions greater than or equal to the threshold will be plotted. Default is 0.01.
#' @param spot_size A numeric value specifying the size of spots. Default is 1.
#' @param palettes A vector of color palettes. Default is `pals::brewer.orrd(22)[3:22]`.
#' @param xlim A numeric vector specifying the x-axis limits. Default is `NULL`.
#' @param ylim A numeric vector specifying the y-axis limits. Default is `NULL`.
#' @param ratio A numeric value specifying the aspect ratio. Default is 1.
#' @param title A character specifying the plot title. Default is `NULL`.
#'
#' @return A ggplot object.
#'
#' @export
mcubePlotPropCellType <- function(
    proportion, coordinates, celltype, spots = NULL,
    he_image = NULL, background = TRUE, proportion_threshold = 0.01,
    spot_size = 1, palettes = pals::brewer.orrd(22)[3:22],
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
  if (!background) {
    # Remove spots with low proportion
    plot_df <- plot_df[plot_df$prop >= proportion_threshold, , drop = FALSE]
  }

  p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = x, y = y))

  # Axis ranges
  if (is.null(xlim)) {
    xlim <- c(min(coordinates[, 1]) - 1, max(coordinates[, 1]) + 1)
  }
  if (is.null(ylim)) {
    ylim <- c(min(coordinates[, 2]) - 1, max(coordinates[, 2] + 1))
  }

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = xlim[1], xmax = xlim[2],
      ymin = -ylim[2], ymax = -ylim[1]
    )
  }

  p <- p + ggplot2::geom_point(ggplot2::aes(color = prop), size = spot_size) +
    ggplot2::scale_colour_gradientn(name = NULL, colors = palettes) +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(trans = "reverse", limits = rev(ylim)) +
    ggplot2::coord_fixed(ratio = ratio)

  if (is.null(title)) {
    title <- paste("Proportion of", celltype)
  }
  p <- p + ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_classic()

  return(p)
}

#' Plot the spatial distribution of cell type proportions in 3D
#'
#' @param proportion A matrix of cell type proportions.
#' Each row represents a spot and each column represents a cell type.
#' @param coordinates A matrix of spatial coordinates.
#' Each row represents a spot and each column represents a spatial dimension.
#' @param celltype A character specifying the cell type to plot.
#' @param proportion_threshold A numeric value between 0 and 1.
#' Spots with proportions greater than or equal to the threshold will be plotted. Default is 0.1.
#' @param spots A character vector specifying the spots to plot.
#' If `NULL`, all spots will be plotted. Default is `NULL`.
#' @param plot_method A character specifying the plotting method.
#' Either "plotly" or "rgl". Default is "plotly".
#' @param color_target A character specifying the color of target spots. Default is "blue".
#' @param color_background A character specifying the color of background spots. Default is "gray".
#' @param opacity_target A numeric value specifying the opacity of target spots. Default is 0.5.
#' @param opacity_background A numeric value specifying the opacity of background spots. Default is 0.05.
#' @param axis_rescale A numeric vector specifying the rescale factors for the spatial coordinates. Default is c(1, 1, 1).
#' @param spot_size A numeric value specifying the size of spots. Default is 3.
#' @param spot_radius A numeric value specifying the radius of spots. Default is 0.5.
#' @param plotly_center A list specifying the center of the plotly scene. Default is list(x = 0, y = 0, z = 0).
#' @param plotly_eye A list specifying the eye of the plotly scene. Default is list(x = 1.25, y = 1.25, z = 1.25).
#' @param rgl_um A numeric vector specifying the user matrix for rgl. Default is c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1).
#'
#' @return A plotly or rgl object.
#'
#' @export
mcubePlotPropCellType3D <- function(
    proportion, coordinates, celltype,
    proportion_threshold = 0.1, spots = NULL,
    plot_method = "plotly", # "plotly" or "rgl"
    color_target = "blue", color_background = "gray",
    opacity_target = 0.5, opacity_background = 0.05,
    axis_rescale = c(1, 1, 1), spot_size = 3, spot_radius = 0.5,
    plotly_center = list(x = 0, y = 0, z = 0),
    plotly_eye = list(x = 1.25, y = 1.25, z = 1.25),
    rgl_um = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)) {
  spots_all <- rownames(proportion)
  spots_target <- rownames(proportion)[
    proportion[, celltype] >= proportion_threshold
  ]
  if (!is.null(spots)) {
    spots_target <- intersect(spots_target, spots)
  }
  spots_background <- setdiff(spots_all, spots_target)

  target_df <- data.frame(
    x = coordinates[spots_target, 1] * axis_rescale[1],
    y = coordinates[spots_target, 2] * axis_rescale[2],
    z = coordinates[spots_target, 3] * axis_rescale[3],
    proportion = proportion[spots_target, celltype]
  )
  background_df <- data.frame(
    x = coordinates[spots_background, 1] * axis_rescale[1],
    y = coordinates[spots_background, 2] * axis_rescale[2],
    z = coordinates[spots_background, 3] * axis_rescale[3]
  )

  # 3D plot
  if (plot_method == "plotly") {
    target_color_rgb <- 0.3 +
      (target_df$proportion - min(target_df$proportion)) * (1 - 0.3) /
        (max(target_df$proportion) - min(target_df$proportion))
    plotly::plot_ly(type = "scatter3d", mode = "markers") |>
      plotly::add_markers(
        data = target_df,
        x = ~x, y = ~y, z = ~z,
        marker = list(
          color = rgb(1 - target_color_rgb, 1 - target_color_rgb, 1),
          size = spot_size,
          opacity = opacity_target
        )
      ) |>
      plotly::add_markers(
        data = background_df,
        x = ~x, y = ~y, z = ~z,
        marker = list(
          color = color_background,
          size = spot_size,
          opacity = opacity_background
        )
      ) |>
      plotly::layout(
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
    rgl::spheres3d(
      target_df$x, target_df$y, target_df$z,
      col = color_target, radius = spot_radius, alpha = target_df$proportion
    )
    # Background in grey
    rgl::spheres3d(
      background_df$x, background_df$y, background_df$z,
      col = color_background, radius = spot_radius, alpha = opacity_background
    )
    rgl::decorate3d()
  }
}

#' Plot the heatmap of cell type proportions
#'
#' @param proportion A matrix of cell type proportions.
#' Each row represents a spot and each column represents a cell type.
#' @param spots A character vector specifying the spots to plot.
#' If `NULL`, all spots will be plotted. Default is `NULL`.
#' @param celltypes A character vector specifying the cell types to plot.
#' If `NULL`, all cell types will be plotted. Default is `NULL`.
#' @param palettes A vector of color palettes. Default is `pals::brewer.purd(20)`.
#' @param title A character specifying the plot title. Default is `NULL`.
#' @param filename A character specifying the filename to save the plot.
#' If `NA`, the plot will not be saved. Default is `NA`.
#'
#' @return A pheatmap object.
#'
#' @export
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

#' Plot function for p-values
#'
#' @param pvalues_list A list of data frames containing p-values.
#' Each data frame corresponds to a cell type. Each row represents a gene and each column represents a kernel matrix/combined p-value.
#' @param which_pvalue A character specifying the p-value to plot. Default is "combined_pvalue".
#' @param type A character specifying the plot type.
#' Either "qqplot" or "histogram". Default is "qqplot".
#' @param minus_log10p_max A numeric value specifying the maximum value of -log10(p-value). Default is `NULL`.
#' @param under_null A logical value.
#' If `TRUE`, the confidence interval under the null hypothesis will be plotted. Default is `FALSE`.
#' @param ci A numeric value specifying the confidence interval. Default is 0.95.
#' @param nrow A positive integer specifying the number of rows of subplots.
#'
#' @return A ggplot object.
#'
#' @export
mcubePlotPvalues <- function(
    pvalues_list, which_pvalue = "combined_pvalue",
    type = "qqplot", minus_log10p_max = NULL,
    under_null = FALSE, ci = 0.95,
    nrow = 1) {
  if (is.null(minus_log10p_max)) {
    # minus_log10p_max <- ceiling(
    #   -log10(0.05 / max(sapply(pvalues_list, length)))
    # )
    minus_log10p_max <- 5
  }
  pvalues_long <- do.call(
    rbind,
    lapply(
      names(pvalues_list),
      FUN = function(x) {
        pvalues_order_x <- order(pvalues_list[[x]][, which_pvalue])
        pvalues_x <- pvalues_list[[x]][, which_pvalue][pvalues_order_x]
        n_pvalues <- length(pvalues_x)
        pvalues_theoretical_x <- ppoints(n_pvalues)
        minus_log10p_x <- -log10(pvalues_x)
        minus_log10p_x[minus_log10p_x > minus_log10p_max] <- minus_log10p_max
        minus_log10p_theoretical_x <- -log10(pvalues_theoretical_x)
        lower_bound <- -log10(qbeta(
          p = (1 - ci) / 2, shape1 = 1:n_pvalues, shape2 = n_pvalues:1
        ))
        upper_bound <- -log10(qbeta(
          p = (1 + ci) / 2, shape1 = 1:n_pvalues, shape2 = n_pvalues:1
        ))

        data.frame(
          celltype = x,
          gene = rownames(pvalues_list[[x]])[pvalues_order_x],
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
      ggplot2::scale_color_discrete(name = "Cell type") +
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
      ggplot2::scale_fill_discrete(name = "Cell type") +
      ggplot2::theme_classic() +
      ggplot2::theme(strip.text = ggplot2::element_blank())
  } else {
    stop("mcubePlot: type must be either 'qqplot' or 'histogram'!")
  }

  return(p)
}

#' Plot the spatial distribution of gene expression variations
#'
#' @param object An \code{\linkS4class{MCUBE}} object with fitted results of the null model.
#' @param celltype A character specifying the cell type to plot.
#' @param gene A character specifying the gene to plot.
#' @param spots A character vector specifying the spots to plot.
#' If `NULL`, all spots will be plotted. Default is `NULL`.
#' @param normalize A logical value.
#' If `TRUE`, the gene expression variations will be normalized to [0, 1]. Default is `TRUE`.
#' @param he_image H&E image background.
#' A H&E image in TIFF format can be read into a raster array using `tiff::readTIFF`. Default is `NULL`.
#' @param background A logical value.
#' If `TRUE`, the spots not pass the expression threshold will be plotted in grey as background. Default is `FALSE`.
#' @param spot_size A numeric value specifying the size of spots. Default is 1.
#' @param palettes A vector of color palettes. Default is `rev(pals::brewer.piyg(10))`.
#' @param opacity_target A numeric value specifying the opacity of target spots. Default is 1.
#' @param opacity_background A numeric value specifying the opacity of background spots. Default is 0.2.
#' @param xlim A numeric vector specifying the x-axis limits. Default is `NULL`.
#' @param ylim A numeric vector specifying the y-axis limits. Default is `NULL`.
#' @param ratio A numeric value specifying the aspect ratio. Default is 1.
#' @param title A character specifying the plot title. Default is `NULL`.
#'
#' @return A ggplot object.
#'
#' @export
mcubePlotExprCellType <- function(
    object, celltype, gene, spots = NULL, normalize = TRUE,
    he_image = NULL, background = FALSE,
    spot_size = 1, palettes = rev(pals::brewer.piyg(10)),
    opacity_target = 1, opacity_background = 0.2,
    xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
  pair_name <- paste(celltype, gene, sep = "_")
  null_model_results <- object@null_models[[pair_name]]
  if (is.null(spots)) {
    spots_target <- null_model_results$spots
  } else {
    spots_target <- intersect(spots, null_model_results$spots)
  }
  expr_level <- null_model_results$u[spots_target, celltype]
  if (normalize) {
    expr_level <- (expr_level - min(expr_level)) /
      (max(expr_level) - min(expr_level))
  }
  target_df <- data.frame(
    x = object@coordinates[spots_target, 1],
    y = object@coordinates[spots_target, 2],
    expr_level = expr_level
  )

  p <- ggplot2::ggplot(data = target_df, ggplot2::aes(x = x, y = y))

  # Axis ranges
  if (is.null(xlim)) {
    xlim <- c(
      min(object@coordinates[, 1]) - 1,
      max(object@coordinates[, 1]) + 1
    )
  }
  if (is.null(ylim)) {
    ylim <- c(
      min(object@coordinates[, 2]) - 1,
      max(object@coordinates[, 2]) + 1
    )
  }

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = xlim[1], xmax = xlim[2],
      ymin = -ylim[2], ymax = -ylim[1]
    )
  }
  # Other spots as background
  if (background) {
    spots_background <- setdiff(rownames(object@coordinates), spots_target)
    background_df <- data.frame(
      x = object@coordinates[spots_background, 1],
      y = object@coordinates[spots_background, 2]
    )
    p <- p + ggplot2::geom_point(
      data = background_df, ggplot2::aes(x = x, y = y),
      color = "gray", size = spot_size, alpha = opacity_background
    )
  }

  p <- p + ggplot2::geom_point(
    ggplot2::aes(color = expr_level),
    size = spot_size, alpha = opacity_target
  ) +
    ggplot2::scale_colour_gradientn(
      name = "Level", colors = palettes,
      values = scales::rescale(
        c(min(expr_level), median(expr_level), max(expr_level))
      )
    ) +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(trans = "reverse", limits = rev(ylim)) +
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

#' Plot the spatial distribution of gene expression variations in binary mode
#'
#' @param object An \code{\linkS4class{MCUBE}} object with fitted results of the null model.
#' @param celltype A character specifying the cell type to plot.
#' @param gene A character specifying the gene to plot.
#' @param spots A character vector specifying the spots to plot.
#' If `NULL`, all spots will be plotted. Default is `NULL`.
#' @param he_image H&E image background.
#' A H&E image in TIFF format can be read into a raster array using `tiff::readTIFF`. Default is `NULL`.
#' @param background A logical value.
#' If `TRUE`, the spots not pass the expression threshold will be plotted in grey as background. Default is `FALSE`.
#' @param spot_size A numeric value specifying the size of spots. Default is 1.
#' @param palettes A vector of color palettes. Default is `c("#32CD32", "#FF69B4")`.
#' @param opacity_target A numeric value specifying the opacity of target spots. Default is 1.
#' @param opacity_background A numeric value specifying the opacity of background spots. Default is 0.2.
#' @param xlim A numeric vector specifying the x-axis limits. Default is `NULL`.
#' @param ylim A numeric vector specifying the y-axis limits. Default is `NULL`.
#' @param ratio A numeric value specifying the aspect ratio. Default is 1.
#' @param title A character specifying the plot title. Default is `NULL`.
#'
#' @return A ggplot object.
#'
#' @export
mcubePlotExprCellTypeBinary <- function(
    object, celltype, gene, spots = NULL,
    he_image = NULL, background = FALSE,
    spot_size = 1, palettes = c("#32CD32", "#FF69B4"),
    opacity_target = 1, opacity_background = 0.3,
    xlim = NULL, ylim = NULL, ratio = 1, title = NULL) {
  pair_name <- paste(celltype, gene, sep = "_")
  null_model_results <- object@null_models[[pair_name]]
  if (is.null(spots)) {
    spots_target <- null_model_results$spots
  } else {
    spots_target <- intersect(spots, null_model_results$spots)
  }
  u <- null_model_results$u[spots_target, celltype]
  rel_expr_level <- factor(ifelse(u >= 0, "Higher", "Lower"))
  target_df <- data.frame(
    x = object@coordinates[spots_target, 1],
    y = object@coordinates[spots_target, 2],
    rel_expr_level = rel_expr_level
  )

  p <- ggplot2::ggplot(data = target_df, ggplot2::aes(x = x, y = y))

  # Axis ranges
  if (is.null(xlim)) {
    xlim <- c(
      min(object@coordinates[, 1]) - 1,
      max(object@coordinates[, 1]) + 1
    )
  }
  if (is.null(ylim)) {
    ylim <- c(
      min(object@coordinates[, 2]) - 1,
      max(object@coordinates[, 2]) + 1
    )
  }
  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = xlim[1], xmax = xlim[2],
      ymin = -ylim[2], ymax = -ylim[1]
    )
  }
  # Other spots as background
  if (background) {
    spots_background <- setdiff(rownames(object@coordinates), spots_target)
    background_df <- data.frame(
      x = object@coordinates[spots_background, 1],
      y = object@coordinates[spots_background, 2]
    )
    p <- p + ggplot2::geom_point(
      data = background_df, ggplot2::aes(x = x, y = y),
      color = "gray", size = spot_size, alpha = opacity_background
    )
  }

  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(color = rel_expr_level),
      size = spot_size, alpha = opacity_target
    ) +
    ggplot2::scale_color_manual(
      name = "Level",
      values = c("Lower" = palettes[1], "Higher" = palettes[2])
    ) +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(trans = "reverse", limits = rev(ylim)) +
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

#' Plot the spatial distribution of gene expression variations in 3D
#'
#' @param object An \code{\linkS4class{MCUBE}} object with fitted results of the null model.
#' @param celltype A character specifying the cell type to plot.
#' @param gene A character specifying the gene to plot.
#' @param spots A character vector specifying the spots to plot.
#' If `NULL`, all spots will be plotted. Default is `NULL`.
#' @param proportion_threshold A numeric value between 0 and 1.
#' Spots with proportions greater than or equal to the threshold will be plotted. Default is 0.1.
#' @param plot_method A character specifying the plotting method.
#' Either "plotly" or "rgl". Default is "plotly".
#' @param palettes A vector of color palettes. Default is `c("#32CD32", "#FF69B4")`.
#' @param opacity_target A numeric value specifying the opacity of target spots. Default is 0.5.
#' @param opacity_background A numeric value specifying the opacity of background spots. Default is 0.05.
#' @param axis_rescale A numeric vector specifying the rescale factors for the spatial coordinates. Default is c(1, 1, 1).
#' @param spot_size A numeric value specifying the size of spots. Default is 3.
#' @param spot_radius A numeric value specifying the radius of spots. Default is 0.5.
#' @param plotly_center A list specifying the center of the plotly scene. Default is list(x = 0, y = 0, z = 0).
#' @param plotly_eye A list specifying the eye of the plotly scene. Default is list(x = 1.25, y = 1.25, z = 1.25).
#' @param rgl_um A numeric vector specifying the user matrix for rgl. Default is c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1).
#'
#' @return A plotly or rgl object.
#'
#' @export
mcubePlotExprCellType3D <- function(
    object, celltype, gene,
    spots = NULL, proportion_threshold = 0.1,
    plot_method = "plotly", # "plotly" or "rgl"
    palettes = c("#32CD32", "#FF69B4"),
    opacity_target = 0.8, opacity_background = 0.05,
    axis_rescale = c(1, 1, 1), spot_size = 3, spot_radius = 0.5,
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
    x = object@coordinates[spots_target, 1] * axis_rescale[1],
    y = object@coordinates[spots_target, 2] * axis_rescale[2],
    z = object@coordinates[spots_target, 3] * axis_rescale[3],
    rel_expr_level = rel_expr_level
  )
  background_df <- data.frame(
    x = object@coordinates[spots_background, 1] * axis_rescale[1],
    y = object@coordinates[spots_background, 2] * axis_rescale[2],
    z = object@coordinates[spots_background, 3] * axis_rescale[3]
  )

  # 3D plot
  if (plot_method == "plotly") {
    plotly_df <- data.frame(
      x = c(target_df$x, background_df$x),
      y = c(target_df$y, background_df$y),
      z = c(target_df$z, background_df$z),
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
    ) |> plotly::layout(
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
        target_df[target_df$rel_expr_level == rel_expr_level_all[l], ]$x,
        target_df[target_df$rel_expr_level == rel_expr_level_all[l], ]$y,
        target_df[target_df$rel_expr_level == rel_expr_level_all[l], ]$z,
        col = palettes[l], radius = spot_radius, alpha = opacity_target
      )
    }
    # Background in grey
    rgl::spheres3d(
      background_df$x, background_df$y, background_df$z,
      col = "gray", radius = spot_radius, alpha = opacity_background
    )
    rgl::decorate3d()
  }
}
