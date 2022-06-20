

durationPlotUI <- function(id) {
  tagList(
    div_with_floating_gear(
      shinycssloaders::withSpinner(plotOutput(NS(id, "duration")), 4, hide.ui = FALSE),
      menu_content = tagList(
        shinyWidgets::switchInput(NS(id, "show_cputime"), label = "Show", onLabel = "CPU Time", offLabel = "Elapsed", value = TRUE),
        shinyWidgets::awesomeCheckbox(NS(id, "relative_duration"), "Relative Duration"),
      )
    )
  )
}

durationPlotServer <- function(id, data, pcadim_sel = reactive(NULL), knn_sel = reactive(NULL), alpha_sel = reactive(NULL), dataset_sel = reactive(NULL)) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("pca_dim", "knn", "dataset", "transformation", "replicate", "cputime_sec", "elapsed_sec") %in% colnames(data)))
  stopifnot(is.reactive(pcadim_sel))
  stopifnot(is.reactive(knn_sel))
  stopifnot(is.reactive(alpha_sel))
  stopifnot(is.reactive(dataset_sel))
  
  moduleServer(id, function(input, output, session) {
    
    filtered_dat <- reactive({
      filtered_dat <- data 
      filtered_dat <- filter_data_with_pca_sel(filtered_dat, pcadim_sel())
      filtered_dat <- filter_data_with_dataset_sel(filtered_dat, dataset_sel())
      if(! is.null(knn_sel())) filtered_dat <- filter(filtered_dat, knn == knn_sel())
      if(! is.null(alpha_sel())) filtered_dat <- filter(filtered_dat, alpha %in% alpha_sel())
      filtered_dat
    })
    
    output$duration <- renderPlot({
      if(is.null(input$show_cputime)){
        metric <- "cputime_sec"
        label <- "CPU Time"
      }else if(input$show_cputime){
        metric <- "cputime_sec"
        label <- "CPU Time"
      }else{
        metric <- "elapsed_sec"
        label <- "Elapsed Time"
      }
      
      if(nrow(filtered_dat()) == 0){
        
      }else if(input$relative_duration){
        dat <- filtered_dat() %>%
          group_by(dataset, knn, pca_dim, replicate) %>%
          mutate(reference_performance = mean(.data[[metric]][transformation == "logp1"])) %>%
          mutate(normed_dur = .data[[metric]] / reference_performance) %>%
          left_join(trans_families, by = "transformation")
        
        ggplot(dat, aes(x = normed_dur, y = transformation, color = family, shape = alpha)) +
          geom_vline(xintercept = 1, size = 0.3, linetype = 2) +
          ggbeeswarm::geom_quasirandom(color = "grey", size = 0.3, alpha = 0.7, groupOnX = FALSE) +
          stat_summary(geom = "point", position = position_dodge2(width = 0.3), fun.data = mean_cl_boot, size = 1.8) +
          scale_y_grouped_discrete(grouping = ~ trans_families_labels[deframe(trans_families)[.x]], gap_size = 1.7, limits = rev,
                                   labels = trans_labels, add_group_label = TRUE) +
          scale_color_manual(values = trans_families_colors, labels = trans_families_labels, guide = "none") +
          scale_x_log10() +
          labs(y = "", x = paste0("Relative ", label), shape = "Overdispersion")
      }else{
        time_minor_ticks <- tibble(ticks = c(seq(0, 60, by = 10), # seconds
                                             seq(0, 60 * 60, by = 10 * 60), # 10 minutes
                                             seq(0, 24 * 60 * 60, by = 6 * 60 * 60),  # 6 hours
                                             seq(0, 7 * 24 * 60 * 60, by = 24 * 60 * 60), # 1 day
                                             seq(0, 8* 7 * 24 * 60 * 60, by = 5 * 24 * 60 * 60))) # 4 weeks
        time_major_ticks <- tibble(ticks = c(1, 60, 60 * 60, 24 * 60 * 60, 7 * 24 * 60 * 60))
        
        
        
        dat <- filtered_dat()  %>%
          left_join(trans_families, by = "transformation")
        
        ggplot(dat, aes(x = .data[[metric]], y = transformation, color = family)) +
          ggbeeswarm::geom_quasirandom(aes(shape = alpha), color = "grey", size = 0.3, alpha = 0.7, groupOnX = FALSE) +
          stat_summary(aes(shape = alpha), geom = "point", position = position_dodge2(width = 0.3), fun.data = mean_cl_boot, size = 1.8) +
          scale_y_grouped_discrete(grouping = ~ trans_families_labels[deframe(trans_families)[.x]], gap_size = 1.7, limits = rev,
                                   labels = trans_labels, add_group_label = TRUE) +
          scale_color_manual(values = trans_families_colors, labels = trans_families_labels, guide = "none") +
          geom_rug(data = time_minor_ticks, aes(x = ticks, y = 0), sides = "b", color = "black") +
          geom_rug(data = time_major_ticks, aes(x = ticks, y = 0), sides = "b", color = "black") +
          scale_x_log10(breaks = c(0.001, 1, 60, 60 * 60, 24 * 60 * 60, 7 * 24 * 60 * 60),
                        labels = c("1ms", "1sec", "1min", "1hour", "1day", "1week"), limits = c(1, NA),
                        expand = expansion(mult = 0.01, 0.05)) +
          labs(y = "", x = label, shape = "Overdispersion")
      }
    }, res = 96)
  })
}




durationPlotApp <- function() {
  ui <- fluidPage(
    includeCSS("www/main.css"),
    durationPlotUI("duration"),
    optionPaneUI("duration_option_pane", show_detailed_pcadim_selector = FALSE),
  )
  server <- function(input, output, session) {
    op <- optionPaneServer("duration_option_pane", data = res, knn_sel = 50, pca_sel = 50)
    durationPlotServer("duration", data = res, 
                       pcadim_sel = op$pca_sel, knn_sel = op$knn_sel, alpha_sel = op$alpha_sel, dataset_sel = op$dataset_sel)
  }
  shinyApp(ui, server)  
}
# durationPlotApp()
