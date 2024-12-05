mcubePlotPvals <- function(
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

mcubePlotExprCellType <- function(
    object, celltype, gene, spots = NULL, normalize = TRUE,
    he_image = NULL, background = FALSE,
    size = 1, palettes = rev(pals::brewer.piyg(10)),
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

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = 0, xmax = max(object@coordinates[, 1]),
      ymin = 0, ymax = max(object@coordinates[, 2])
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
      color = "gray", size = size, alpha = opacity_background
    )
  }

  p <- p + ggplot2::geom_point(
    ggplot2::aes(color = expr_level),
    size = size, alpha = opacity_target
  ) +
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
#     he_image = NULL, size = 1, palettes = c("#32CD32", "#C5C5C5", "#FF69B4"),
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
    he_image = NULL, background = FALSE,
    size = 1, palettes = c("#32CD32", "#FF69B4"),
    opacity_target = 1, opacity_background = 0.2,
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

  # H&E image background
  if (!is.null(he_image)) {
    p <- p + ggplot2::annotation_raster(
      he_image,
      xmin = 0, xmax = max(object@coordinates[, 1]),
      ymin = 0, ymax = max(object@coordinates[, 2])
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
      color = "gray", size = size, alpha = opacity_background
    )
  }

  p <- p + ggplot2::geom_point(
    ggplot2::aes(color = rel_expr_level),
    size = size, alpha = opacity_target
  ) + ggplot2::scale_color_manual(
    name = "Level",
    values = c("Lower" = palettes[1], "Higher" = palettes[2])
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
    spots = NULL, proportion_threshold = 0.1,
    plot_method = "plotly", # "plotly" or "rgl"
    palettes = c("#32CD32", "#FF69B4"),
    opacity_target = 0.8, opacity_background = 0.05,
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

mcubePlotPropCellType3D <- function(
    proportion, coordinates, celltype,
    proportion_threshold = 0.1, spots = NULL,
    plot_method = "plotly", # "plotly" or "rgl"
    color_target = "blue", color_background = "gray",
    opacity_target = 0.5, opacity_background = 0.05,
    axis_rescale = c(1, 1, 1), spot_size = 4, spot_radius = 0.5,
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
    plotly::plot_ly(type = "scatter3d", mode = "markers") %>%
      plotly::add_markers(
        data = target_df,
        x = ~x, y = ~y, z = ~z,
        marker = list(
          color = rgb(1 - target_color_rgb, 1 - target_color_rgb, 1),
          size = spot_size,
          opacity = opacity_target
        )
      ) %>%
      plotly::add_markers(
        data = background_df,
        x = ~x, y = ~y, z = ~z,
        marker = list(
          color = color_background,
          size = spot_size,
          opacity = opacity_background
        )
      ) %>%
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
