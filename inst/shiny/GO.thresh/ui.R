ui <- shiny::shinyUI(shiny::fluidPage(
  shiny::titlePanel(
    paste0("Threshold proposal:",
           if(!is.null(stress.ident1)) paste0(" GO",substr(x = strsplit(stress.ident1, "GO")[[1]][2], 1,8),": ", round(thresh.stress.ident1,2)),
           if(!is.null(stress.ident2)) paste0(", GO",substr(x = strsplit(stress.ident2, "GO")[[1]][2], 1,8),": ", round(thresh.stress.ident2,2)),
           if(!is.null(notstress.ident3)) paste0(", GO",substr(x = strsplit(notstress.ident3, "GO")[[1]][2], 1,8),": ", round(thresh.notstress.ident3,2)),
           if(!is.null(notstress.ident4)) paste0(", GO",substr(x = strsplit(notstress.ident4, "GO")[[1]][2], 1,8),": ", round(thresh.notstress.ident4,2)))),
  shiny::sidebarLayout(position = "left",
                       shiny::sidebarPanel(
                        if(!is.null(stress.ident1)) {
                          shiny::sliderInput("t.stress.ident1", paste0("Threshold for GO",substr(x = strsplit(stress.ident1, "GO")[[1]][2], 1,8),":"),
                                         min=sliders$min.x.stress.ident1, max=sliders$max.x.stress.ident1, value=thresh.stress.ident1, step=sliders$step.stress.ident1,
                                         animate=shiny::animationOptions(100))
                        },
                        if(!is.null(stress.ident2)) {
                          shiny::sliderInput("t.stress.ident2", paste0("Threshold for GO",substr(x = strsplit(stress.ident2, "GO")[[1]][2], 1,8),":"),
                                         min=sliders$min.x.stress.ident2, max=sliders$max.x.stress.ident2, value=thresh.stress.ident2, step=sliders$step.stress.ident2,
                                         animate=shiny::animationOptions(100))
                        },
                        if(!is.null(notstress.ident3)) {
                          shiny::sliderInput("t.notstress.ident3", paste0("Threshold for GO",substr(x = strsplit(notstress.ident3, "GO")[[1]][2], 1,8),":"),
                                                min=sliders$min.x.notstress.ident3, max=sliders$max.x.notstress.ident3, value=thresh.notstress.ident3, step=sliders$step.notstress.ident3,
                                                animate=shiny::animationOptions(100))
                        },
                        if(!is.null(notstress.ident4)) {
                          shiny::sliderInput("t.notstress.ident4", paste0("Threshold for GO",substr(x = strsplit(notstress.ident4, "GO")[[1]][2],1,8),":"),
                                                min=sliders$min.x.notstress.ident4, max=sliders$max.x.notstress.ident4, value=thresh.notstress.ident4, step=sliders$step.notstress.ident4,
                                                animate=shiny::animationOptions(100))
                        },
                        shiny::actionButton('save_inputs', 'Save asad thresholds')
                       ),
                shiny::mainPanel(shiny::tabsetPanel(
                  if(!is.null(stress.ident1)) shiny::tabPanel(paste0("Histogram GO",substr(x = strsplit(stress.ident1, "GO")[[1]][2],1,8)), shiny::plotOutput("hist.stress.ident1")),
                  if(!is.null(stress.ident2)) shiny::tabPanel(paste0("Histogram GO",substr(x = strsplit(stress.ident2, "GO")[[1]][2],1,8)), shiny::plotOutput("hist.stress.ident2")),
                  if(!is.null(notstress.ident3)) shiny::tabPanel(paste0("Histogram GO",substr(x = strsplit(notstress.ident3, "GO")[[1]][2],1,8)), shiny::plotOutput("hist.notstress.ident3")),
                  if(!is.null(notstress.ident4)) shiny::tabPanel(paste0("Histogram GO",substr(x = strsplit(notstress.ident4, "GO")[[1]][2],1,8)), shiny::plotOutput("hist.notstress.ident4")),
                  shiny::tabPanel("Stress Classification", shiny::plotOutput("stress.umap")),
                  shiny::tabPanel("Counts", shiny::plotOutput("count.bar"))
                ))
  ),
  shiny::mainPanel(shiny::plotOutput("go.score"))
))
