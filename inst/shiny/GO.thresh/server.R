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
      df <- data.frame(gr.av.stress.scores1 = average.vec$"gr.av.stress.scores1", assignment = F)
      df$assignment[df$"gr.av.stress.scores1" > i] <- T
      df
    })
  }
  if (!is.null(stress.ident2)) {
    df.data.stress.ident2 <- shiny::reactive({
      i <- input$"t.stress.ident2"
      df <- data.frame(gr.av.stress.scores2 = average.vec$"gr.av.stress.scores2", assignment = F)
      df$assignment[df$"gr.av.stress.scores2" > i] <- T
      df
    })
  }
  if (!is.null(notstress.ident3)) {
    df.data.notstress.ident3 <- shiny::reactive({
      i <- input$"t.notstress.ident3"
      df <- data.frame(gr.av.notstress.scores3 = average.vec$"gr.av.notstress.scores3", assignment = F)
      df$assignment[df$"gr.av.notstress.scores3" > i] <- T
      df
    })
  }
  if (!is.null(notstress.ident4)) {
    df.data.notstress.ident4 <- shiny::reactive({
      i <- input$"t.notstress.ident4"
      df <- data.frame(gr.av.notstress.scores4 = average.vec$"gr.av.notstress.scores4", assignment = F)
      df$assignment[df$"gr.av.notstress.scores4" > i] <- T
      df
    })
  }

  # Stress Filtering & Assignment
  # Cells are considered stressed if their scores are above the threshold for any of the stress identifiers
  # and not above the threshold for any of the non-stress identifiers
  obj2.data <- shiny::reactive({
    if (!is.null(stress.ident1)) i.stress.ident1 <- input$"t.stress.ident1"
    if (!is.null(stress.ident2)) i.stress.ident2 <- input$"t.stress.ident2"
    if (!is.null(notstress.ident3)) i.notstress.ident3 <- input$"t.notstress.ident3"
    if (!is.null(notstress.ident4)) i.notstress.ident4 <- input$"t.notstress.ident4"

    obj2 <- obj
    meta2 <- obj2@meta.data

    gr.av.scores.1 <- meta2[, idents$"stress.ident1"]
    gr.av.scores.2 <- meta2[, idents$"stress.ident2"]
    if (!is.null(stress.ident1)) {
      i1.bool <- as.numeric(levels(gr.av.scores.1))[gr.av.scores.1] > i.stress.ident1
      if (!is.null(stress.ident2)) {
        i2.bool <- as.numeric(levels(gr.av.scores.2))[gr.av.scores.2] > i.stress.ident2
        stress.bool <- i1.bool | i2.bool # Combine both boolean vectors for stress determination
      } else {
        stress.bool <- i1.bool # Use only stress.ident1 for stress determination
      }
    } else {
      # Process only the second stress identifier if the first one is null
      if (!is.null(stress.ident2)) {
        i2.bool <- as.numeric(levels(gr.av.scores.2))[gr.av.scores.2] > i.stress.ident2
        stress.bool <- i2.bool # Use only stress.ident2 for stress determination
      }
    }

    # Process not stress identifiers similarly
    notstress.bool <- NULL
    gr.av.scores.3 <- meta2[, idents$"notstress.ident3"]
    gr.av.scores.4 <- meta2[, idents$"notstress.ident4"]

    # Determine not stressed cells based on the thresholds for notstress.ident3 and potentially notstress.ident4
    if (!is.null(notstress.ident3)) {
      i3.bool <- as.numeric(levels(gr.av.scores.3))[gr.av.scores.3] > i.notstress.ident3
      if (!is.null(notstress.ident4)) {
        i4.bool <- as.numeric(levels(gr.av.scores.4))[gr.av.scores.4] > i.notstress.ident4
        notstress.bool <- i3.bool | i4.bool # Combine boolean vectors for notstress.ident3 and notstress.ident4
      } else {
        notstress.bool <- i3.bool # Use only notstress.ident3 for not stress determination
      }
    } else {
      if (!is.null(notstress.ident4)) {
        i4.bool <- as.numeric(levels(gr.av.scores.4))[gr.av.scores.4] > i.notstress.ident4
        notstress.bool <- i4.bool
      }
    }

    # Final assignment of stressed cells
    if (!is.null(notstress.bool)) {
      obj2$"is.Stressed" <- stress.bool & !notstress.bool
    } else {
      obj2$"is.Stressed" <- stress.bool
    }

    obj2
  })

  obj2.data.part.GO <- shiny::reactive({
    if (!is.null(stress.ident1)) i.stress.ident1 <- input$"t.stress.ident1"
    if (!is.null(stress.ident2)) i.stress.ident2 <- input$"t.stress.ident2"
    if (!is.null(notstress.ident3)) i.notstress.ident3 <- input$"t.notstress.ident3"
    if (!is.null(notstress.ident4)) i.notstress.ident4 <- input$"t.notstress.ident4"

    obj2 <- obj
    meta2 <- obj2@meta.data

    if (!is.null(stress.ident1)) {
      gr.av.scores.1 <- meta2[, idents$"stress.ident1"]
      obj2$"stress.ident1.thresh_cluster" <- as.numeric(levels(gr.av.scores.1))[gr.av.scores.1] > i.stress.ident1
      obj2$"stress.ident1.thresh_cluster"[obj2$"stress.ident1.thresh_cluster" == FALSE] <- F
      obj2$"stress.ident1.thresh_cluster"[obj2$"stress.ident1.thresh_cluster" == TRUE] <- T
    }

    if (!is.null(stress.ident2)) {
      gr.av.scores.2 <- meta2[, idents$"stress.ident2"]
      obj2$"stress.ident2.thresh_cluster" <- as.numeric(levels(gr.av.scores.2))[gr.av.scores.2] > i.stress.ident2
      obj2$"stress.ident2.thresh_cluster"[obj2$"stress.ident2.thresh_cluster" == FALSE] <- F
      obj2$"stress.ident2.thresh_cluster"[obj2$"stress.ident2.thresh_cluster" == TRUE] <- T
    }

    if (!is.null(notstress.ident3)) {
      gr.av.scores.3 <- meta2[, idents$"notstress.ident3"]
      obj2$"notstress.ident3.thresh_cluster" <- as.numeric(levels(gr.av.scores.3))[gr.av.scores.3] > i.notstress.ident3
      obj2$"notstress.ident3.thresh_cluster"[obj2$"notstress.ident3.thresh_cluster" == FALSE] <- F
      obj2$"notstress.ident3.thresh_cluster"[obj2$"notstress.ident3.thresh_cluster" == TRUE] <- T
    }

    if (!is.null(notstress.ident4)) {
      gr.av.scores.4 <- meta2[, idents$"notstress.ident4"]
      obj2$"notstress.ident4.thresh_cluster" <- as.numeric(levels(gr.av.scores.4))[gr.av.scores.4] > i.notstress.ident4
      obj2$"notstress.ident4.thresh_cluster"[obj2$"notstress.ident4.thresh_cluster" == FALSE] <- F
      obj2$"notstress.ident4.thresh_cluster"[obj2$"notstress.ident4.thresh_cluster" == TRUE] <- T
    }

    obj2
  })

  #------ Output Plots --------#
  colorz <- Seurat.utils::gg_color_hue(2)[2:1]

  # GO score histograms
  if (!is.null(stress.ident1)) {
    output$"hist.stress.ident1" <- shiny::renderPlot({
      df <- df.data.stress.ident1()
      ggplot2::ggplot(df, ggplot2::aes(x = gr.av.stress.scores1, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = colorz) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(gr.av.stress.scores1) / length(gr.av.stress.scores1)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(gr.av.stress.scores1, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(gr.av.stress.scores1, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$"min.x.stress.ident1", sliders$"max.x.stress.ident1")) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }
  if (!is.null(stress.ident2)) {
    output$"hist.stress.ident2" <- shiny::renderPlot({
      df <- df.data.stress.ident2()
      ggplot2::ggplot(df, ggplot2::aes(x = gr.av.stress.scores2, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = colorz) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(gr.av.stress.scores2) / length(gr.av.stress.scores2)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(gr.av.stress.scores2, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(gr.av.stress.scores2, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$"min.x.stress.ident2", sliders$"max.x.stress.ident2")) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }
  if (!is.null(notstress.ident3)) {
    output$"hist.notstress.ident3" <- shiny::renderPlot({
      df <- df.data.notstress.ident3()
      ggplot2::ggplot(df, ggplot2::aes(x = gr.av.notstress.scores3, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = rev(colorz)) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(gr.av.notstress.scores3) / length(gr.av.notstress.scores3)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(gr.av.notstress.scores3, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(gr.av.notstress.scores3, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$"min.x.notstress.ident3", sliders$"max.x.notstress.ident3")) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }
  if (!is.null(notstress.ident4)) {
    output$"hist.notstress.ident4" <- shiny::renderPlot({
      df <- df.data.notstress.ident4()
      ggplot2::ggplot(df, ggplot2::aes(x = gr.av.notstress.scores4, fill = factor(assignment))) +
        ggplot2::scale_fill_manual(values = rev(colorz)) +
        ggplot2::geom_histogram(binwidth = (2 * IQR(gr.av.notstress.scores4) / length(gr.av.notstress.scores4)^(1 / 3))) +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(gr.av.notstress.scores4, 0.9)), colour = "black") +
        ggplot2::ylab("count") +
        ggplot2::geom_text(hjust = -.1, vjust = 10, mapping = ggplot2::aes(x = quantile(gr.av.notstress.scores4, 0.9), y = Inf, label = "90% quantile")) +
        ggplot2::coord_cartesian(xlim = c(sliders$"min.x.notstress.ident4", sliders$"max.x.notstress.ident4")) +
        ggplot2::xlab("GO Score") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "assignment"))
    })
  }

  # stress assignment UMAP
  output$stress.umap <- shiny::renderPlot({
    obj2 <- obj2.data()
    Seurat.utils::clUMAP(obj = obj2, ident = "is.Stressed", save.plot = F) + Seurat::NoAxes() +
      ggplot2::scale_color_manual(values = colorz)
  })

  # stress assignment barplot cell numbers
  output$count.bar <- shiny::renderPlot({
    obj2 <- obj2.data()
    meta2 <- obj2@meta.data

    ident.plot <- plot.cluster.shiny
    ggplot2::ggplot(meta2, ggplot2::aes(x = ident.plot, fill = is.Stressed, stat = "count")) +
      ggplot2::scale_fill_manual(values = colorz) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 1)) +
      ggplot2::geom_bar() +
      ggplot2::xlab("")
  })

  output$"go.score" <- shiny::renderPlot(
    {
      obj2 <- obj2.data.part.GO()
      ptlist <- list()
      count <- 1
      if (!is.null(stress.ident1)) {
        GOterm <- ww.parse.GO(stress.ident1)
        GOscore <- ww.convert.GO_term.2.score(GOterm)
        ptlist[[count]] <- Seurat::FeaturePlot(object = obj, features = GOscore, min.cutoff = "q01", max.cutoff = "q99") +
          Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(stress.ident2)) {
        GOterm <- ww.parse.GO(stress.ident2)
        GOscore <- ww.convert.GO_term.2.score(GOterm)
        ptlist[[count]] <- Seurat::FeaturePlot(object = obj, features = GOscore, min.cutoff = "q01", max.cutoff = "q99") +
          Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident3)) {
        GOterm <- ww.parse.GO(notstress.ident3)
        GOscore <- ww.convert.GO_term.2.score(GOterm)
        ptlist[[count]] <- Seurat::FeaturePlot(object = obj, features = GOscore, min.cutoff = "q01", max.cutoff = "q99") +
          Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident4)) {
        GOterm <- ww.parse.GO(notstress.ident4)
        GOscore <- ww.convert.GO_term.2.score(GOterm)
        ptlist[[count]] <- Seurat::FeaturePlot(object = obj, features = GOscore, min.cutoff = "q01", max.cutoff = "q99") +
          Seurat::NoLegend() + Seurat::NoAxes()
        count <- count + 1
      }


      if (!is.null(stress.ident1)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = obj2, ident = "stress.ident1.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = colorz) + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(stress.ident2)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = obj2, ident = "stress.ident2.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = colorz) + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident3)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = obj2, ident = "notstress.ident3.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = rev(colorz)) + Seurat::NoAxes()
        count <- count + 1
      }
      if (!is.null(notstress.ident4)) {
        ptlist[[count]] <- Seurat.utils::clUMAP(obj = obj2, ident = "notstress.ident4.thresh_cluster", save.plot = F) +
          ggplot2::ggtitle(ggplot2::element_blank()) +
          ggplot2::scale_color_manual(values = rev(colorz)) + Seurat::NoAxes()
        count <- count + 1
      }
      gridExtra::grid.arrange(grobs = ptlist, ncol = ceiling(length(ptlist) / 2), title = "GO score")
    },
    height = 400,
    width = w
  )

  shiny::observeEvent(input$save_inputs, {
    obj <- obj2.data()
    if (!is.null(stress.ident1)) {
      i.stress.ident1 <- input$"t.stress.ident1"
      obj@misc$gruffi$"thresh.stress.ident1" <- i.stress.ident1
    }
    if (!is.null(stress.ident2)) {
      i.stress.ident2 <- input$"t.stress.ident2"
      obj@misc$gruffi$"thresh.stress.ident2" <- i.stress.ident2
    }
    if (!is.null(notstress.ident3)) {
      i.notstress.ident3 <- input$"t.notstress.ident3"
      obj@misc$gruffi$"thresh.notstress.ident3" <- i.notstress.ident3
    }
    if (!is.null(notstress.ident4)) {
      i.notstress.ident4 <- input$"t.notstress.ident4"
      obj@misc$gruffi$"thresh.notstress.ident4" <- i.notstress.ident4
    }
    shiny::stopApp(returnValue = obj)
  })
})
