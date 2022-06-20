



contrastPlotUI <- function(id) {
  tagList(
    # shiny::fluidRow(
    #   column(3, tagList(
    #     shinyWidgets::pickerInput(NS(id, "trans_a"), label = "", choices = NA, inline = TRUE, width = "fit"),
    #     shinyWidgets::pickerInput(NS(id, "alpha_a"), label = "alpha", choices = NA, inline = TRUE, width = "fit")
    #   )),
    #   column(1, p(" vs ")),
    #   column(3,tagList(
    #     shinyWidgets::pickerInput(NS(id, "trans_b"), label = "", choices = NA, inline = TRUE, width = "fit"),
    #     shinyWidgets::pickerInput(NS(id, "alpha_b"), label = "alpha", choices = NA, inline = TRUE, width = "fit")
    #   ))
    # ),
    fluidRow(
      column(8, offset = 2, h4(
        div(shinyWidgets::pickerInput(NS(id, "trans_a"), label = "", choices = NA, inline = TRUE, 
                                  width = "fit"), class = "inline-transformation-picker"), 
        "(α=", 
        div(shinyWidgets::pickerInput(NS(id, "alpha_a"), label = "", choices = NA, inline = TRUE, width = "fit"), class = "inline-transformation-picker"),
        ")",
        p(" vs ", style = "display: inline;"),
        div(shinyWidgets::pickerInput(NS(id, "trans_b"), label = "", choices = NA, inline = TRUE, 
                                  width = "fit"), class = "inline-transformation-picker"),
        "(α=", 
        div(shinyWidgets::pickerInput(NS(id, "alpha_b"), label = "", choices = NA, inline = TRUE, width = "fit"), class = "inline-transformation-picker"),
        ")"
      ))
    ),
    div_with_floating_gear(
      shinycssloaders::withSpinner(plotOutput(NS(id, "contrast")), 4, hide.ui = FALSE),
      tagList(
        menu_content = shinyWidgets::awesomeCheckbox(NS(id, "zoom_in"), "zoom in"),
        shinyWidgets::awesomeCheckbox(NS(id, "relative_knn"), "Relative k-NN")
      )
    )
  )
}

contrastPlotServer <- function(id, data, metric = overlap, pcadim_sel = reactive(NULL), knn_sel = reactive(NULL), dataset_sel = reactive(NULL)) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("pca_dim", "knn", "dataset", "transformation", "replicate") %in% colnames(data)))
  stopifnot(! is.null(data %>% pull({{metric}})))
  stopifnot(is.reactive(pcadim_sel))
  stopifnot(is.reactive(knn_sel))
  stopifnot(is.reactive(dataset_sel))
  
  trans_choices <- trans_families %>%
    left_join(enframe(trans_labels_plain, name = "transformation")) %>%
    group_by(family) %>%
    group_map(~ set_names(list(deframe(.x[,c("value", "transformation")])), .y)) %>%
    flatten() 
  
  moduleServer(id, function(input, output, session) {
    
    shinyWidgets::updatePickerInput(session, "trans_a", choices = trans_choices, selected = "pearson_clip")
    shinyWidgets::updatePickerInput(session, "trans_b", choices = trans_choices, selected = "logp1")
    
    
    filtered_dat <- reactive({
      filtered_dat <- data
      filtered_dat <- filter_data_with_pca_sel(filtered_dat, pcadim_sel())
      filtered_dat <- filter_data_with_dataset_sel(filtered_dat, dataset_sel())
      if(! is.null(knn_sel())) filtered_dat <- filter(filtered_dat, knn == knn_sel())
      filtered_dat
    })
    observe({
      choices <- unique(filter(filtered_dat(), transformation == input$trans_a)$alpha)
      shinyWidgets::updatePickerInput(session, "alpha_a", choices = choices)
    })
    observe({
      choices <- unique(filter(filtered_dat(), transformation == input$trans_b)$alpha)
      shinyWidgets::updatePickerInput(session, "alpha_b", choices = choices)
    })
    
    output$contrast <- renderPlot({
      req(input$trans_a, input$trans_b, input$alpha_a, input$alpha_b)
      if(nrow(filtered_dat()) == 0){

      }else if(input$relative_knn){
        dat <- filtered_dat() %>%
          group_by(dataset, knn, pca_dim, replicate) %>%
          mutate(reference_performance = mean({{metric}})) %>%
          mutate({{metric}} := {{metric}} / reference_performance) %>%
          left_join(trans_families, by = "transformation") %>%
          filter((transformation == input$trans_a & alpha == input$alpha_a) | 
                   (transformation == input$trans_b & alpha == input$alpha_b))

        if(input$trans_a == input$trans_b && input$alpha_a == input$alpha_b){
          stop(shiny::safeError("Left and right must differ"))
        }else if(input$trans_a == input$trans_b){
          a <- paste0(input$trans_a, " (alpha=", input$alpha_a, ")")
          b <- paste0(input$trans_b, " (alpha=", input$alpha_b, ")")
          dat <- mutate(dat, transformation = ifelse(transformation == input$trans_a & alpha == input$alpha_a, a, b))
        }else{
          a <- input$trans_a
          b <- input$trans_b
        }
        
        dabestr::dabest(dat, x = transformation, y = {{metric}}, idx = c(a, b)) %>%
          dabestr::mean_diff(ci = 95, reps = 5000) %>%
            plot(rawplot.ylabel = "Relative k-NN overlap", color.column = dataset, palette = "Accent")
      }else{
        dat <- filtered_dat() %>%
          left_join(trans_families, by = "transformation") %>%
          filter((transformation == input$trans_a & alpha == input$alpha_a) | 
                   (transformation == input$trans_b & alpha == input$alpha_b)) |>
          suppressWarnings()

        if(input$trans_a == input$trans_b && input$alpha_a == input$alpha_b){
          stop(shiny::safeError("Left and right must differ"))
        }else if(input$trans_a == input$trans_b){
          a <- paste0(input$trans_a, " (alpha=", input$alpha_a, ")")
          b <- paste0(input$trans_b, " (alpha=", input$alpha_b, ")")
          dat <- mutate(dat, transformation = ifelse(as.character(transformation) == input$trans_a & alpha == input$alpha_a, a, b))
        }else{
          a <- input$trans_a
          b <- input$trans_b
        }
        
        dabestr::dabest(dat, x = transformation, y = {{metric}}, idx = c(a, b)) %>%
          dabestr::mean_diff(ci = 95, reps = 5000) %>%
          plot(rawplot.ylabel = "k-NN overlap", color.column = dataset, palette = "Accent") |>
          suppressWarnings()
      }
    }, res = 96)
  })
}

contrastPlotApp <- function() {
  
  ui <- fluidPage(
    includeCSS("www/main.css"),
    contrastPlotUI("contrast"),
    optionPaneUI("consistency_pca_sel", show_alpha_selector = FALSE),
  )
  server <- function(input, output, session) {
    op <- optionPaneServer("consistency_pca_sel", data = res, metric = overlap, knn_sel_init = 50, pca_sel_init = 50)
    contrastPlotServer("contrast", data = res, metric = overlap, pcadim_sel = op$pca_sel, 
                       knn_sel = op$knn_sel, dataset_sel = op$dataset_sel)
  }
  shinyApp(ui, server)  
}
# contrastPlotApp()
