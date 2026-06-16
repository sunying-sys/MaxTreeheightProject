viviNetwork_new <- function(
    mat,
    intThreshold = NULL,
    intLims = NULL,
    impLims = NULL,
    intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n = 100)),
    impPal = rep("#D73027", 100),
    removeNode = FALSE,
    layout = igraph::layout_in_circle,
    cluster = NULL,
    nudge_x = 0.05,
    nudge_y = 0.03,
    edgeWidths = 1:4
) {
  nnodes <- nrow(mat)
  if (nnodes == 1) stop("Only one node provided, no graph drawn")
  if (is.numeric(cluster) && length(cluster) != nnodes) cluster <- NULL
  if (is.numeric(layout) && !identical(dim(layout), as.integer(c(nnodes, 2)))) layout <- igraph::layout_in_circle
  
  df <- vivid:::as.data.frame.vivid(mat)
  dfImp <- df[df$Measure == "Vimp", ]
  dfInt <- df[df$Measure == "Vint", ]
  dfInt <- dfInt[-which(dfInt$Row < dfInt$Col), ]
  dfInt <- dfInt[!is.na(dfInt$Value), ]
  dfInt <- dfInt[with(dfInt, order(Value)), ]
  
  if (is.null(impLims)) {
    impLimits <- range(dfImp$Value, na.rm = TRUE)
    if (impLimits[1] == impLimits[2]) {
      impLimits[1] <- impLimits[1] - impLimits[1] / 4
      impLimits[2] <- impLimits[2] + impLimits[2] / 4
    }
    impLimits <- range(labeling::rpretty(impLimits[1], impLimits[2]))
  } else {
    impLimits <- impLims
  }
  
  if (is.null(intLims)) {
    intLimits <- range(dfInt$Value)
    if (intLimits[1] == intLimits[2]) {
      intLimits[1] <- intLimits[1] - intLimits[1] / 4
      intLimits[2] <- intLimits[2] + intLimits[2] / 4
    }
    intLimits <- range(labeling::rpretty(intLimits[1], intLimits[2]))
  } else {
    intLimits <- intLims
  }
  
  dfInt1 <- dfInt
  
  if (!is.null(intThreshold)) {
    if (intThreshold > max(dfInt$Value) | intThreshold < min(dfInt$Value)) {
      warning("Interaction threshold value is outside range of interaction values and will be ignored")
      intThreshold <- NULL
    }
    if (!is.null(intThreshold)) dfInt1 <- dfInt[dfInt$Value > intThreshold, ]
  }
  
  g <- igraph::make_empty_graph(nnodes, directed = FALSE)
  g <- igraph::add_edges(
    graph = g,
    edges = as.vector(t(dfInt1[c("Row", "Col")]))
  )
  igraph::E(g)$weight <- dfInt1$Value
  
  if (removeNode) {
    rnode <- igraph::degree(g) == 0
    g <- igraph::delete_vertices(g, rnode)
    dfImp <- dfImp[!rnode, ]
    if (is.numeric(cluster)) cluster <- cluster[!rnode]
    if (is.numeric(layout)) layout <- layout[!rnode, , drop = FALSE]
  }
  
  if (is.function(layout)) {
    glayout <- layout(g)
    if (identical(layout, igraph::layout_in_circle)) glayout <- glayout[, 2:1]
  } else {
    glayout <- layout
  }
  
  m1 <- apply(glayout, 2, min)
  r <- apply(glayout, 2, max) - m1
  glayout <- -1 + 2 * scale(glayout, m1, r)
  attr(glayout, "scaled:scale") <- NULL
  attr(glayout, "scaled:center") <- NULL
  
  if (r[1] == 0) glayout[, 1] <- seq(-1, 1, length.out = nrow(glayout))
  if (r[2] == 0) glayout[, 2] <- seq(-1, 1, length.out = nrow(glayout))
  
  mapinto <- function(x, lims, v) {
    x[is.na(x)] <- lims[1]
    x <- pmin(pmax(x, lims[1]), lims[2])
    i <- cut(
      x,
      breaks = seq(lims[1], lims[2], length = length(v) + 1),
      include.lowest = TRUE
    )
    v[i]
  }
  
  edgeCols <- mapinto(dfInt1$Value, intLimits, intPal)
  edgeWidthScaled <- mapinto(dfInt1$Value, intLimits, sort(edgeWidths))
  
  glayout[abs(glayout) < 1e-04] <- 0
  nudged <- sign(glayout)
  nudged[nudged[, 2] == 0, 2] <- 1
  
  nodeSize <- mapinto(dfImp$Value, impLimits, seq(1, 2.4, length.out = 10))
  nudged[, 1] <- nudged[, 1] * nodeSize * nudge_x
  nudged[, 2] <- nudged[, 2] * nodeSize * nudge_y
  
  xlim <- c(-0.05 + min(nudged[, 1]), 1.05 + max(nudged[, 1]))
  ylim <- c(-0.05 + min(nudged[, 2]), 1.05 + max(nudged[, 2]))
  
  suppressMessages(
    p <- GGally::ggnet2(
      g,
      mode = glayout,
      size = 0,
      edge.label = NULL,
      edge.size = edgeWidthScaled,
      edge.color = edgeCols
    ) +
      ggplot2::xlim(xlim) +
      ggplot2::ylim(ylim) +
      ggplot2::geom_label(
        ggplot2::aes(label = dfImp$Variable_1),
        size = 4.5,
        nudge_x = nudged[, 1],
        nudge_y = nudged[, 2],
        hjust = "middle",
        vjust = "middle",
        linewidth = 0
      )
  )
  
  p$scales$scales <- Filter(
    function(s) !"size" %in% s$aesthetics,
    p$scales$scales
  )
  
  for (i in seq_along(p$layers)) {
    if ("size" %in% names(p$layers[[i]]$mapping)) {
      p$layers[[i]]$mapping$size <- NULL
      p$layers[[i]]$aes_params$size <- 0
    }
  }
  
  if (!is.null(cluster)) {
    if (!is.numeric(cluster)) cluster <- cluster(g)$membership
    colPal <- rainbow(length(unique(cluster)))
    colCluster <- colPal[cluster]
    p <- p +
      geom_encircle_vivi(
        ggplot2::aes(group = cluster),
        spread = 0.01,
        alpha = 0.2,
        expand = 0.03,
        fill = colCluster
      )
  }
  
  node_df <- data.frame(
    x = (glayout[, 1] + 1) / 2,
    y = (glayout[, 2] + 1) / 2,
    Vimp = as.numeric(dfImp$Value)
  )
  
  dummy_df <- data.frame(
    x = 0,
    y = 0,
    Vint = mean(intLimits)
  )
  
  p <- p +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = x, y = y, size = Vimp),
      inherit.aes = FALSE,
      fill = impPal[length(impPal)],
      colour = "transparent",
      shape = 21,
      show.legend = TRUE
    ) +
    ggplot2::scale_size_continuous(
      name = "Vimp",
      range = c(1, 8),
      limits = range(node_df$Vimp, na.rm = TRUE),
      breaks = pretty(node_df$Vimp, n = 4)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        title = "Vimp",
        order = 2,
        override.aes = list(
          fill = impPal[length(impPal)],
          colour = impPal[length(impPal)],
          shape = 21
        )
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      data = dummy_df,
      ggplot2::aes(x = x, y = y, fill = Vint),
      inherit.aes = FALSE,
      size = -1
    ) +
    ggplot2::scale_fill_gradientn(
      name = "Vint",
      colors = intPal,
      limits = intLimits,
      guide = ggplot2::guide_colorbar(
        order = 1,
        frame.colour = "black",
        ticks.colour = "black"
      ),
      oob = scales::squish
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1
    )
  
  p
}
