
######### Custom ggplot2 theme #########

font_size <- 8
font_size_small <- 6
font_size_tiny <- 5
font_size_large <- 10
publ_theme <- cowplot::theme_cowplot(font_size = font_size, rel_small = font_size_small / font_size,
                                     rel_tiny = font_size_tiny / font_size, rel_large = font_size_large / font_size,
                                     line_size = 0.3) +
  theme(plot.title = element_text(size = font_size),
        axis.title = element_text(size = font_size_small))
theme_set(publ_theme)


######### Custom plotting functions #########

convert_dims <- function( width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1){
  units <- match.arg(units)
  if(units == "inches"){
    units <- "in"
  }
  to_inches <- function(x) x/c(`in` = 1, cm = 2.54, mm = 2.54 * 
                                 10, px = dpi)[units]
  to_inches(c(width, height)) * scale
}

my_pdf <- function(filename, width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1, ...){
  dim <- convert_dims(width, height, units, dpi, scale)
  grDevices::pdf(filename, width = dim[1], height = dim[2], useDingbats = FALSE, ...)
}


my_tikz <- function(filename, width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1, stand_alone = TRUE, ...){
  dim <- convert_dims(width, height, units, dpi, scale)
  tikzDevice::tikz(filename, width = dim[1], height = dim[2], standAlone = stand_alone, 
                   documentDeclaration = c(getOption("tikzDocumentDeclaration"), r"(\renewcommand{\familydefault}{\sfdefault})", r"(\usepackage{helvet})"),  # Use sans serif font Helvetica
                   packages = c(options("tikzLatexPackages")$tikzLatexPackages, "\\usepackage{amssymb}", "\\usepackage{amsmath}", "\\usepackage{bm}"),  ...)
}

save_plot <- function(filename, plot = ggplot2::last_plot(), width = 6.2328, height = 3.71, units = c("inches", "cm", "mm", "px"), dpi = 300, scale = 1, latex_support = FALSE, ...){
  
  old_dev <- grDevices::dev.cur()
  if(latex_support){
    filename <- if(stringr::str_ends(filename, "\\.pdf")){
      paste0(stringr::str_sub(filename, end  = -5L), ".tex")
    }
    my_tikz(filename, width = width, height = height, units = units, dpi = dpi, scale = scale, stand_alone = TRUE)
  }else{
    dim <- convert_dims(width, height, units, dpi, scale)
    dev <- ggplot2:::plot_dev(NULL, filename, dpi = dpi)
    dev(filename = filename, width = dim[1], height = dim[2], ...)
    on.exit(utils::capture.output({
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
    }))
  }
  
  grid::grid.draw(plot)
  
  if(latex_support){
    grDevices::dev.off()
    if (old_dev > 1) grDevices::dev.set(old_dev)
    
    withr::with_dir(dirname(filename), {
      tools::texi2pdf(basename(filename), clean = TRUE)
      # Remove .tex file
      file.remove(basename(filename)) 
      raw_file_name <- tools::file_path_sans_ext(basename(filename))
      # rastered images
      ras_files <- list.files(pattern = paste0("^", raw_file_name, "_ras\\d*.png"))
      if(length(ras_files) > 0){
        file.remove(ras_files)
      }
    })
  }
  
  invisible(filename)  
}

plot_assemble <- function(..., .plot_objs = NULL, width = 6.2328, height = 3.71, units = c("inches", "cm", "mm", "px"), latex_support = FALSE, show_grid_lines = TRUE, filename = NULL){
  units <- match.arg(units)
  
  plots <- if(is.null(.plot_objs)){
    list(...)
  }else{
    .plot_objs
  }
  
  if(show_grid_lines){
    x_breaks <- scales::breaks_pretty(n = 10)(seq(0, width, length.out = 100))
    y_breaks <- scales::breaks_pretty(n = 10)(seq(0, height, length.out = 100))
  }else{
    x_breaks <- c(0,Inf)
    y_breaks <- c(0,Inf)
  }
  
  if(! is.null(filename)){
    old_dev <- grDevices::dev.cur()
    if(latex_support){
      filename <- if(stringr::str_ends(filename, "\\.pdf")){
        paste0(stringr::str_sub(filename, end  = -5L), ".tex")
      }
      my_tikz(filename, width = width, height = height, units = units, stand_alone = TRUE)
    }else{
      my_pdf(filename, width = width, height = height, units = units)
      on.exit(utils::capture.output({
        grDevices::dev.off()
        if (old_dev > 1) grDevices::dev.set(old_dev)
      }))
    }
  }
  
  
  plotgardener::pageCreate(width = width, height = height, default.units = units, xgrid = diff(x_breaks)[1], ygrid = diff(y_breaks)[1], showGuides = show_grid_lines)
  for(obj in plots){
    if(is.ggplot(obj)){
      plotgardener::plotGG(obj, x = 0, y = 0, width = width, height = height, default.units = units)
    }else if(grid::is.grob(obj)){
      grid::grid.draw(obj)
    }else if(is.list(obj)){
      stopifnot(! is.null(names(obj)))
      stopifnot("plot" %in% names(obj))
      .x <- obj$x %||% 0
      .y <- obj$y %||% 0
      .width <- obj$width %||% width
      .height <- obj$height %||% height
      .units <- obj$units %||% units
      plotgardener::plotGG(obj$plot, x = .x, y = .y, width = .width, height = .height, default.units = .units)
    }
  }
  
  if(! is.null(filename)){
    if(latex_support){
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
      withr::with_dir(dirname(filename), {
        tools::texi2pdf(basename(filename), clean = TRUE)
        # Remove .tex file
        file.remove(basename(filename)) 
        raw_file_name <- tools::file_path_sans_ext(basename(filename))
        # rastered images
        ras_files <- list.files(pattern = paste0("^", raw_file_name, "_ras\\d*.png"))
        if(length(ras_files) > 0){
          file.remove(ras_files)
        }
      })
    }
  }
}

annotate_text <- function(label, x = 0, y = 0, fontsize = 12, hjust = 0, vjust = 0, ...){
  list(plot = cowplot::ggdraw() + cowplot::draw_label(label, size = fontsize, hjust = hjust, vjust = vjust, ...), x = x, y = y, width =0, height = 0)
}

# Note that x and y are from the lower left corner (instead of upper left :/)
annotate_graphic <- function(filename, x = 0, y = 0, units = c("inches", "cm", "mm", "px", "user")){
  stopifnot(file.exists(filename))
  units <- match.arg(units)
  abs_filepath <- tools::file_path_as_absolute(filename)
  tikzDevice::grid.tikzNode(x = x, y = y, units = units, opts = "draw=none,fill=none",  
                            content = paste0(r"(\includegraphics{")", abs_filepath, r"("})"),
                            draw = FALSE)
}

