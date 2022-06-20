


reducedDimPlotUI <- function(id) {
  tagList(
    div_with_floating_gear(
      shinycssloaders::withSpinner(plotOutput(NS(id, "reduced_dim_plot")), 4, hide.ui = FALSE),
      menu_content = shinyWidgets::awesomeCheckbox(NS(id, "zoom_in"), "zoom in")
    ),
    shinyWidgets::pickerInput(inputId = NS(id, "dataset"), label = "Dataset", choices = ""),
    shinyWidgets::switchInput(inputId = NS(id, "show_tsne"), label = "Dimension Reduction", onLabel = "tSNE", offLabel = "PCA"),
    shinyWidgets::switchInput(inputId = NS(id, "color_cluster"), label = "Color by", onLabel = "cluster", offLabel = "Seq. Depth"),
  )
}


reducedDimPlotServer <- function(id, data) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("axis1", "axis2") %in% colnames(data)))
  stopifnot(all(c("dim_red_method", "dataset", "origin", "downsampling", "benchmark") %in% colnames(data)))
  stopifnot(all(c("cluster", "col_sums") %in% colnames(data)))
  
  
  moduleServer(id, function(input, output, session) {
    shinyWidgets::updatePickerInput(session, "dataset", choices = sort(unique(data$dataset)))
    
    filtered_dat <- reactive({
      req(input$dataset)
      data %>%
        filter(dataset == input$dataset, dim_red_method == (if(input$show_tsne) "tsne" else "pca"))
    })
    
    output$reduced_dim_plot <- renderPlot({
      req(input$dataset)
      if(nrow(filtered_dat()) > 0){
        benchmark <- unique(filtered_dat()$benchmark)[1]
        pl <- filtered_dat() %>%
          group_by(downsampling, origin) %>%
          # mutate(scaled_axis1 = axis1 / max(axis1, axis2),
          #        scaled_axis2 = axis2 / max(axis1, axis2)) %>%
          group_by(groups = if(.env$benchmark == "downsampling") downsampling
                   else if(.env$benchmark == "simulation") origin
                   else 1) %>%
          group_map(function(dat, key){
            key <- key[[1]][1]
            ggplot(dat, aes(x = axis1, y = axis2)) +
              geom_point(aes(color = if(input$color_cluster) as.factor(cluster) else col_sums)) +
              coord_fixed() +
              (if(input$color_cluster) scale_color_discrete(name = "Cluster")
               else scale_color_viridis_c(name = "Seq. Depth", trans = "log10")) +
              labs(subtitle = if(key == 1) "" else key,
                   x = if(input$show_tsne) "tSNE 1" else "PCA 1",
                   y = if(input$show_tsne) "tSNE 2" else "PCA 2")
          })
        cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_text(input$dataset, fontface = "bold"),
                           cowplot::plot_grid(plotlist = pl, nrow = 1, align = "vh"),
                           ncol = 1, rel_heights = c(0.2, 1)
                           )
      }
    }, res = 96)
    
  })
}






reducedDimPlotApp <- function() {
  
  tmp <- read_rds("../benchmark/output/benchmark_results/dataset_plot_data.RDS")
  source("../notebooks/annotation_helper.R")
  reduced_dim_data <- bind_rows(
    bind_rows(tmp$downsampling) %>%
      pivot_longer(-c(name, cluster, col_sums_full, col_sums_reduced), names_pattern = "(tsne|pca)_(ground_truth|log_counts)_(full|reduced)_(axis\\d)", names_to = c("dim_red_method", "origin", "downsampling", ".value")) %>%
      pivot_longer(c(col_sums_full, col_sums_reduced), names_pattern = "(.+)_(full|reduced)", names_to = c(".value", "downsampling2")) %>%
      filter(downsampling == downsampling2) %>% 
      dplyr::select(-downsampling2)%>%
      dplyr::rename(dataset = name),
    bind_rows(tmp$consistency) %>%
      pivot_longer(-c(name, cluster, col_sums), names_pattern = "(tsne|pca)_(ground_truth|log_counts)_(axis\\d)", names_to = c("dim_red_method", "origin", ".value")) %>%
      dplyr::rename(dataset = name),
    bind_rows(tmp$simulation) %>%
      pivot_longer(-c(simulator, cluster, col_sums), names_pattern = "(tsne|pca)_(ground_truth|log_counts)_(axis\\d)", names_to = c("dim_red_method", "origin", ".value")) %>%
      dplyr::rename(dataset = simulator)
  ) %>%
    left_join(enframe(dataset_benchmark, name = "dataset", value = "benchmark"))

  ui <- fluidPage(
    includeCSS("www/main.css"),
    reducedDimPlotUI("reduced_dim_plot"),
  )
  server <- function(input, output, session) {
    reducedDimPlotServer("reduced_dim_plot", data = reduced_dim_data)
  }
  shinyApp(ui, server)  
}
reducedDimPlotApp()
