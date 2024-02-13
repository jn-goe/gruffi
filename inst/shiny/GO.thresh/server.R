server <- shiny::shinyServer(function(input, output, session) {
  count.not.null <- 0
  if (!is.null(stress.ident1)) count.not.null <- count.not.null + 1
  if (!is.null(stress.ident2)) count.not.null <- count.not.null + 1
  if (!is.null(notstress.ident3)) count.not.null <- count.not.null + 1
  if (!is.null(notstress.ident4)) count.not.null <- count.not.null + 1

  w <- 1000
  if (count.not.null == 3) w <- 750
  if (count.not.null == 2) w <- 500
  if (count.not.null == 1) w <- 250

  # dataframes for stress assignment per GO term
  if (!is.null(stress.ident1)) {
    df.data.stress.ident1 <- shiny::reactive({
      i <- input$"t.stress.ident1"
      df <- data.frame(av.stress.ident1 = average.vec$"av.stress.ident1", assignment = F)
      df$assignment[df$"av.stress.ident1" > i] <- T
      df
    })
  }
  if (!is.null(stress.ident2)) {
    df.data.stress.ident2 <- shiny::reactive({
      i <- input$"t.stress.ident2"
      df <- data.frame(av.stress.ident2 = average.vec$"av.stress.ident2", assignment = F)
      df$assignment[df$"av.stress.ident2" > i] <- T
      df
    })
  }
  if (!is.null(notstress.ident3)) {
    df.data.notstress.ident3 <- shiny::reactive({
      i <- input$"t.notstress.ident3"
      df <- data.frame(av.notstress.ident3 = average.vec$"av.notstress.ident3", assignment = F)
      df$assignment[df$"av.notstress.ident3" > i] <- T
      df
    })
  }
  if (!is.null(notstress.ident4)) {
    df.data.notstress.ident4 <- shiny::reactive({
      i <- input$"t.notstress.ident4"
      df <- data.frame(av.notstress.ident4 = average.vec$"av.notstress.ident4", assignment = F)
      df$assignment[df$"av.notstress.ident4" > i] <- T
      df
    })
  }

  # stress assignment
  c.data <- shiny::reactive({
    if (!is.null(stress.ident1)) i.stress.ident1 <- input$"t.stress.ident1"
    if (!is.null(stress.ident2)) i.stress.ident2 <- input$"t.stress.ident2"
    if (!is.null(notstress.ident3)) i.notstress.ident3 <- input$"t.notstress.ident3"
    if (!is.null(notstress.ident4)) i.notstress.ident4 <- input$"t.notstress.ident4"

    c <- obj

    if (!is.null(stress.ident1)) {
      i1.bool <- as.numeric(levels(c@meta.data[ , idents$"stress.ident1"]))[c@meta.data[ , idents$"stress.ident1"]] > i.stress.ident1
      if (!is.null(stress.ident2)) {
        i2.bool <- as.numeric(levels(c@meta.data[ , idents$"stress.ident2"]))[c@meta.data[ , idents$"stress.ident2"]] > i.stress.ident2
        stress.bool <- i1.bool | i2.bool
      } else {
        stress.bool <- i1.bool
      }
    } else {
      if (!is.null(stress.ident2)) {
        i2.bool <- as.numeric(levels(c@meta.data[ , idents$"stress.ident2"]))[c@meta.data[ , idents$"stress.ident2"]] > i.stress.ident2
        stress.bool <- i2.bool
      }
    }

    notstress.bool <- NULL
    if (!is.null(notstress.ident3)) {
      i3.bool <- as.numeric(levels(c@meta.data[ , idents$"notstress.ident3"]))[c@meta.data[ , idents$"notstress.ident3"]] > i.notstress.ident3
      if (!is.null(notstress.ident4)) {
        i4.bool <- as.numeric(levels(c@meta.data[ , idents$"notstress.ident4"]))[c@meta.data[ , idents$"notstress.ident4"]] > i.notstress.ident4
        notstress.bool <- i3.bool | i4.bool
      } else {
        notstress.bool <- i3.bool
      }
    } else {
      if (!is.null(notstress.ident4)) {
        i4.bool <- as.numeric(levels(c@meta.data[ , idents$"notstress.ident4"]))[c@meta.data[ , idents$"notstress.ident4"]] > i.notstress.ident4
        notstress.bool <- i4.bool
      }
    }

    if (!is.null(notstress.bool)) {
      c$is.Stressed <- stress.bool & !notstress.bool
    } else {
      c$is.Stressed <- stress.bool
    }
    # c$is.Stressed[c$is.Stressed == FALSE] <- F
    # c$is.Stressed[c$is.Stressed == TRUE] <- T
    c
  })
  c.data.part.GO <- shiny::reactive({
    if (!is.null(stress.ident1)) i.stress.ident1 <- input$"t.stress.ident1"
    if (!is.null(stress.ident2)) i.stress.ident2 <- input$"t.stress.ident2"
    if (!is.null(notstress.ident3)) i.notstress.ident3 <- input$"t.notstress.ident3"
    if (!is.null(notstress.ident4)) i.notstress.ident4 <- input$"t.notstress.ident4"

    c <- obj

    if (!is.null(stress.ident1)) {
      c$"stress.ident1.thresh_cluster" <- as.numeric(levels(c@meta.data[ , idents$"stress.ident1"]))[c@meta.data[ , idents$"stress.ident1"]] > i.stress.ident1
      c$"stress.ident1.thresh_cluster"[c$"stress.ident1.thresh_cluster" == FALSE] <- F
      c$"stress.ident1.thresh_cluster"[c$"stress.ident1.thresh_cluster" == TRUE] <- T
    }

    if (!is.null(stress.ident2)) {
      c$"stress.ident2.thresh_cluster" <- as.numeric(levels(c@meta.data[ , idents$"stress.ident2"]))[c@meta.data[ , idents$"stress.ident2"]] > i.stress.ident2
      c$"stress.ident2.thresh_cluster"[c$"stress.ident2.thresh_cluster" == FALSE] <- F
      c$"stress.ident2.thresh_cluster"[c$"stress.ident2.thresh_cluster" == TRUE] <- T
    }

    if (!is.null(notstress.ident3)) {
      c$"notstress.ident3.thresh_cluster" <- as.numeric(levels(c@meta.data[ , idents$"notstress.ident3"]))[c@meta.data[ , idents$"notstress.ident3"]] > i.notstress.ident3
      c$"notstress.ident3.thresh_cluster"[c$"notstress.ident3.thresh_cluster" == FALSE] <- F
      c$"notstress.ident3.thresh_cluster"[c$"notstress.ident3.thresh_cluster" == TRUE] <- T
    }

    if (!is.null(notstress.ident4)) {
      c$"notstress.ident4.thresh_cluster" <- as.numeric(levels(c@meta.data[ , idents$"notstress.ident4"]))[c@meta.data[ , idents$"notstress.ident4"]] > i.notstress.ident4
      c$"notstress.ident4.thresh_cluster"[c$"notstress.ident4.thresh_cluster" == FALSE] <- F
      c$"notstress.ident4.thresh_cluster"[c$"notstress.ident4.thresh_cluster" == TRUE] <- T
    }

    c
  })

  #------ Output Plots --------#
  # GO score histograms
  if (!is.null(stress.ident1)) {
    output$hist.stress.ident1 <- shiny::renderPlot({
      df <- df.data.stress.ident1()
      ggplot2::ggplot(df, ggplot2::aes(x = av.stress.ident1, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = c(Seurat.utils::gg_color_hue(2)[2], Seurat.utils::gg_color_hue(2)[1])) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(av.stress.ident1) / length(av.stress.ident1)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(av.stress.ident1, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(av.stress.ident1, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$min.x.stress.ident1, sliders$max.x.stress.ident1)) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }
  if (!is.null(stress.ident2)) {
    output$hist.stress.ident2 <- shiny::renderPlot({
      df <- df.data.stress.ident2()
      ggplot2::ggplot(df, ggplot2::aes(x = av.stress.ident2, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = c(Seurat.utils::gg_color_hue(2)[2], Seurat.utils::gg_color_hue(2)[1])) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(av.stress.ident2) / length(av.stress.ident2)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(av.stress.ident2, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(av.stress.ident2, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$min.x.stress.ident2, sliders$max.x.stress.ident2)) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }
  if (!is.null(notstress.ident3)) {
    output$hist.notstress.ident3 <- shiny::renderPlot({
      df <- df.data.notstress.ident3()
      ggplot2::ggplot(df, ggplot2::aes(x = av.notstress.ident3, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = c(Seurat.utils::gg_color_hue(2)[1], Seurat.utils::gg_color_hue(2)[2])) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(av.notstress.ident3) / length(av.notstress.ident3)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(av.notstress.ident3, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(av.notstress.ident3, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$min.x.notstress.ident3, sliders$max.x.notstress.ident3)) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }
  if (!is.null(notstress.ident4)) {
    output$hist.notstress.ident4 <- shiny::renderPlot({
      df <- df.data.notstress.ident4()
      ggplot2::ggplot(df, ggplot2::aes(x = av.notstress.ident4, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = c(Seurat.utils::gg_color_hue(2)[1], Seurat.utils::gg_color_hue(2)[2])) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(av.notstress.ident4) / length(av.notstress.ident4)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(av.notstress.ident4, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(av.notstress.ident4, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$min.x.notstress.ident4, sliders$max.x.notstress.ident4)) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }

  # stress assignment UMAP
  output$stress.umap <- shiny::renderPlot({
    c <- c.data()
    Seurat.utils::clUMAP(obj = c, ident = "is.Stressed", save.plot = F) + Seurat::NoAxes() +
      ggplot2::scale_color_manual(values = c(Seurat.utils::gg_color_hue(2)[2], Seurat.utils::gg_color_hue(2)[1]))
  })

  # stress assignment barplot cell numbers
  output$count.bar <- shiny::renderPlot({
    c <- c.data()
    ident.plot <- plot.cluster.shiny
    ggplot2::ggplot(c@meta.data, ggplot2::aes(x = ident.plot, fill = is.Stressed, stat = "count")) +
      ggplot2::scale_fill_manual(values = c(Seurat.utils::gg_color_hue(2)[2], Seurat.utils::gg_color_hue(2)[1])) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 1)) +
      ggplot2::geom_bar() +
      ggplot2::xlab("")
  })

  output$go.score <- shiny::renderPlot(
    {
      c <- c.data.part.GO()
      ptlist <- list()
      count <- 1
      if (!is.null(stress.ident1)) {
        ptlist[[count]] <- Seurat::FeaturePlot(obj = obj, features = ww.convert.GO_term.2.score(paste0(substr(x = strsplit(stress.ident1, "cl.av")[[1]][2], 2, nchar(strsplit(stress.ident1, "cl.av")[[1]][2])))), min.cutoff = "q01", max.cutoff = "q99") + Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(stress.ident2)) {
        ptlist[[count]] <- Seurat::FeaturePlot(obj = obj, features = ww.convert.GO_term.2.score(paste0(substr(x = strsplit(stress.ident2, "cl.av")[[1]][2], 2, nchar(strsplit(stress.ident2, "cl.av")[[1]][2])))), min.cutoff = "q01", max.cutoff = "q99") + Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident3)) {
        ptlist[[count]] <- Seurat::FeaturePlot(obj = obj, features = ww.convert.GO_term.2.score(paste0(substr(x = strsplit(notstress.ident3, "cl.av")[[1]][2], 2, nchar(strsplit(notstress.ident3, "cl.av")[[1]][2])))), min.cutoff = "q01", max.cutoff = "q99") + Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident4)) {
        ptlist[[count]] <- Seurat::FeaturePlot(obj = obj, features = ww.convert.GO_term.2.score(paste0(substr(x = strsplit(notstress.ident4, "cl.av")[[1]][2], 2, nchar(strsplit(notstress.ident4, "cl.av")[[1]][2])))), min.cutoff = "q01", max.cutoff = "q99") + Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(stress.ident1)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = c, ident = "stress.ident1.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = c(Seurat.utils::gg_color_hue(2)[2], Seurat.utils::gg_color_hue(2)[1])) + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(stress.ident2)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = c, ident = "stress.ident2.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = c(Seurat.utils::gg_color_hue(2)[2], Seurat.utils::gg_color_hue(2)[1])) + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident3)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = c, ident = "notstress.ident3.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = c(Seurat.utils::gg_color_hue(2)[1], Seurat.utils::gg_color_hue(2)[2])) + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident4)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = c, ident = "notstress.ident4.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = c(Seurat.utils::gg_color_hue(2)[1], Seurat.utils::gg_color_hue(2)[2])) + Seurat::NoAxes()
        count <- count + 1
      }
      gridExtra::grid.arrange(grobs = ptlist, ncol = ceiling(length(ptlist) / 2), title = "GO score")
    },
    height = 400,
    width = w
  )

  shiny::observeEvent(input$save_inputs, {
    obj <- c.data()
    if (!is.null(stress.ident1)) {
      i.stress.ident1 <- input$"t.stress.ident1"
      obj@misc$gruffi$thresh.stress.ident1 <- i.stress.ident1
    }
    if (!is.null(stress.ident2)) {
      i.stress.ident2 <- input$"t.stress.ident2"
      obj@misc$gruffi$thresh.stress.ident2 <- i.stress.ident2
    }
    if (!is.null(notstress.ident3)) {
      i.notstress.ident3 <- input$"t.notstress.ident3"
      obj@misc$gruffi$thresh.notstress.ident3 <- i.notstress.ident3
    }
    if (!is.null(notstress.ident4)) {
      i.notstress.ident4 <- input$"t.notstress.ident4"
      obj@misc$gruffi$thresh.notstress.ident4 <- i.notstress.ident4
    }
    shiny::stopApp(returnValue = obj)
  })
})
