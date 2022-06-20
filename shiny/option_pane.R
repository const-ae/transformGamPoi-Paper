
filter_data_with_dataset_sel <- function(data, dataset_sel){
  if(is.null(dataset_sel) || nrow(data) == 0){
    # Do nothing
    data
  }else if(str_starts(dataset_sel, "All ")){
    benchmark <- str_remove(dataset_sel, "All ")
    filter(data, benchmark == .env$benchmark) 
  }else{
    filter(data, dataset  == dataset_sel) 
  }
}

filter_data_with_pca_sel <- function(data, pcadim_sel, dataset_var = dataset){
  stopifnot("pca_dim" %in% colnames(data))
  datasets <- unique(pull(data, {{dataset_var}}))
  
  if(is.null(pcadim_sel) || nrow(data) == 0){
    # Do nothing
    data
  }else if(! is.null(names(pcadim_sel))){
    stopifnot(all(names(pcadim_sel) %in% datasets))
    inner_join(data, enframe(pcadim_sel, name = "dataset", value = "pca_dim"), by = set_names(c("pca_dim", "dataset"), nm = c("pca_dim", ensym(dataset_var))))
  }else if(length(pcadim_sel) == 1){
    filter(data, pca_dim == pcadim_sel)
  }
}


optionPaneUI <- function(id, show_pcadim_selector = TRUE, show_detailed_pcadim_selector = TRUE, 
                         show_knn_selector = TRUE, show_alpha_selector = TRUE, show_dataset_selector = TRUE) {
  res <- tagList(
    dataset_selector = if(show_dataset_selector) shinyWidgets::pickerInput(NS(id, "dataset_sel"), "Dataset", choices = NA),
    alpha_selector = if(show_alpha_selector) shinyWidgets::checkboxGroupButtons(NS(id, "alpha_sel"), "Overdispersion", choices = c("TRUE", "0", "0.05"), selected = c("TRUE", "0")),
    knn_selector = if(show_knn_selector) shinyWidgets::sliderTextInput(NS(id, "knn_sel"), "k-NN", choices = NA, grid = TRUE),
    pcadim_selector = if(show_pcadim_selector){
      if(show_detailed_pcadim_selector){
        random_id <- paste0(sample(letters[1:6], replace = TRUE, size = 8), collapse = "")
        collapse_button <- a(class = "", role="button", `data-toggle`="collapse",
                             href= paste0("#collapse_", random_id), `aria-expanded`="false", 
                             `aria-controls`=paste0("collapse_", random_id), "Details")
        shinyWidgets::sliderTextInput(NS(id, "global_pca_sel"), div("PCA (", collapse_button, ")"), choices = NA, grid = TRUE)
      }else{
        shinyWidgets::sliderTextInput(NS(id, "global_pca_sel"), "PCA", choices = NA, grid = TRUE)
      }
    }
  ) %>%
    purrr::discard(~ is.null(.x)) %>%
    purrr::map(~ column(3, .x))
  
  if(show_detailed_pcadim_selector){
    detail_view <- div(div_with_floating_gear(
      shinycssloaders::withSpinner(plotOutput(NS(id, "pcaplot"), click = NS(id, "pcaplot_click")), 4, hide.ui = FALSE),
      menu_content = shinyWidgets::awesomeCheckbox(NS(id, "zoom_in"), "zoom in")
    ), id = paste0("collapse_", random_id), class="collapse")
    tagList(fluidRow(res), detail_view)
  }else{
    fluidRow(res)
  }
}


optionPaneServer <- function(id, data, metric = overlap, pca_sel_init = 50, knn_sel_init = 50, 
                             alpha_sel_init = c("TRUE", "FALSE"), dataset_sel_init = "All consistency") {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("pca_dim", "knn", "dataset", "transformation") %in% colnames(data)))
  stopifnot(! is.null(data %>% pull({{metric}})))
  stopifnot(!is.reactive(pca_sel_init))
  stopifnot(!is.reactive(knn_sel_init))
  stopifnot(!is.reactive(alpha_sel_init))
  stopifnot(!is.reactive(dataset_sel_init))
  
  pca_steps <- sort(unique(data$pca_dim))
  datasets <- sort(unique(data$dataset))
  
  moduleServer(id, function(input, output, session) {
    # Initial widget update
    data %>%
      dplyr::group_by(benchmark) %>%
      group_map(~ set_names(list(c(paste0("All ", .y),unique(.x$dataset))), .y)) %>%
      purrr::flatten() %>%
      shinyWidgets::updatePickerInput(session, "dataset_sel", choices = ., selected = dataset_sel_init)
    shinyWidgets::updateSliderTextInput(session, "knn_sel", choices = c(0, knn_sel_init), selected = knn_sel_init)
    shinyWidgets::updateSliderTextInput(session, "global_pca_sel", choices = c(0, knn_sel_init), selected = pca_sel_init)
    isolate(shinyWidgets::updateCheckboxGroupButtons(session, "alpha_sel", choices = as.character(alpha_sel_init), selected = alpha_sel_init))
    # shinyWidgets::updateCheckboxGroupButtons(session, "alpha_sel", choices = c("TRUE", "FALSE")) # This is somehow broken...
    
    # Update clicking event
    pca_sel <- reactiveVal(deframe(data.frame(datasets, pca_sel_init)))
    observeEvent(input$pcaplot_click, {
      vals <- pca_sel()
      vals[input$pcaplot_click$panelvar1] <- pca_steps[which.min(abs(pca_steps - input$pcaplot_click$x))]
      pca_sel(vals)
    }, ignoreNULL = TRUE)
    
    observeEvent(input$global_pca_sel, {
      pca_sel(deframe(tibble(datasets, input$global_pca_sel)))
    })
    
    data_sel_data <- reactive({
      data_sel_data <- data
      if(! is.null(input$dataset_sel)) data_sel_data <- filter_data_with_dataset_sel(data, input$dataset_sel)
      data_sel_data
    })
    
    # Dynamic widget update
    observeEvent(data_sel_data(), {
      shinyWidgets::updateSliderTextInput(session, "knn_sel", choices = sort(unique(data_sel_data()$knn)))
      shinyWidgets::updateSliderTextInput(session, "global_pca_sel", choices = sort(unique(data_sel_data()$pca_dim)))
      shinyWidgets::updateCheckboxGroupButtons(session, "alpha_sel", choices = sort(unique(as.character(data_sel_data()$alpha))),
                                               selected = input$alpha_sel)
    })
    
    # Update detailed PCA selection plot
    filtered_dat <- reactive({
      filtered_dat <- data_sel_data()
      if(! is.null(input$alpha_sel)) filtered_dat <- filter(filtered_dat, alpha %in% input$alpha_sel)
      if(! is.null(input$knn_sel)) filtered_dat <- filter(filtered_dat, knn ==  input$knn_sel)
      filtered_dat
    })
      
    output$pcaplot <- renderPlot({
      if(nrow(filtered_dat()) > 0){
        filtered_dat() %>%
          group_by(dataset, pca_dim, knn, transformation, alpha) %>%
          summarize({{metric}} := mean({{metric}})) %>%
          ggplot(aes(x = pca_dim, y = {{metric}})) +
            geom_line(aes(group = paste0(transformation, "-", alpha))) +
            scale_x_log10() +
            ggh4x::facet_wrap2(vars(dataset), ncol = 5, drop = TRUE, scales = if(input$zoom_in) "free" else "fixed") +
            (if(! is.null(pca_sel())) geom_vline(data = enframe(pca_sel(), name = "dataset") %>% 
                                                   filter(dataset %in% unique(filtered_dat()$dataset)), 
                                                 aes(xintercept = value))) 
      }
    }, res = 96)
    
    # Return everything
    list(pca_sel = pca_sel, knn_sel = reactive(input$knn_sel), dataset_sel = reactive(input$dataset_sel), alpha_sel = reactive(input$alpha_sel))
  })
}








optionPaneApp <- function() {
  ui <- fluidPage(
    includeCSS("www/main.css"),
    optionPaneUI("dataset_1"),
  )
  server <- function(input, output, session) {
    optionPaneServer("dataset_1", data = res, pca_sel = 20)
  }
  shinyApp(ui, server)  
}
# optionPaneApp()
