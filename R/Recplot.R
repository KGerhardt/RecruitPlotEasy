#' @import reticulate
#' @import ggplot2
#' @import shiny
#' @import data.table
#' @import plotly
#' @import cowplot
#' @import shinyBS
#' @import hms
#' @import easycsv
#' @import shinyalert
#' @import htmlwidgets



#Currently removed.
#enveomics.R

# # #Dev
# library(reticulate)
# library(ggplot2)
# library(shiny)
# library(data.table)
# library(plotly)
# library(cowplot)
# library(shinyBS)
# library(hms)
# library(easycsv)
# library(shinyalert)
# library(htmlwidgets)
# #library(enveomics.R)

#Helper functions
{
  #This will download whatever the current python script is. You have to run it before landing page.

  get_python <- function(){

    recplot_py <- import("RecruitPlotEasy")

    return(recplot_py)

  }

  #Checks for necessary setup steps and makes the requisite installs as needed.

  prepare_environment <- function(){

    cat("Checking for Miniconda and installing if necessary...\n")
    try({
      install_miniconda()
    })

    #Checking for first-time use of recplots
    if(!"recruitment_plots" %in% conda_list()$name){
      cat("Creating Miniconda environment: 'recruitment_plots'\n")
      conda_create(envname = "recruitment_plots")
    }

    use_miniconda(condaenv = "recruitment_plots", required = T)

    if(!py_module_available("numpy")){
      cat("Attempting to install NumPy... ")
      try({

        py_install(packages = "numpy", envname = "recruitment_plots", pip = T)
        cat("Done!\n")
      })
    }

    if(py_module_available("RecruitPlotEasy")){
      cat("The python component of RecruitPlotEasy is already installed. You probably shouldn't be seeing this warning. Did you call prepare_environment() twice?\n")
    }else{
      cat("Attempting to install RecruitPlotEasy python... ")
      try({
        py_install(packages = "RecruitPlotEasy", envname = "recruitment_plots", pip = T)
        cat("Done!\n")
      })
    }

    recplot_py <- get_python()

    return(recplot_py)

  }

  #Prepares the background miniconda env if necessary; otherwise, sets the environment and loads the python script functions

  initiate <- function(){
    cat("Initiating recruitment plot environment. Please wait a moment. ")

    tryCatch({

      use_miniconda(condaenv = "recruitment_plots", required = T )

      recplot_py <- get_python()

    }, error = function(cond){

      cat("\nPerforming first-time setup. Wait a moment, please.\n")

      recplot_py <- prepare_environment()


    })

    return(recplot_py)

  }

  choose_directory = function(caption = 'Select data directory') {
    if (exists('choose.dir')) {
      choose.dir(caption = caption)
    } else {
      easycsv::choose_dir()
    }
  }

  #Wish I could add captions to the unix/osx versions
  choose_file = function(caption) {
    if (exists('choose.files')) {
      choose.files(caption = caption, multi = F)
    } else {
      file.choose(new = F)
    }
  }

  path_simplifier <- function(working_directory, file_path){

    file_path_t <- gsub("\\\\", "/", file_path)
    working_directory <- gsub("\\\\", "/", working_directory)
    working_directory <- paste0(working_directory, "/")

    if(substr(working_directory, 1, 12) == "Working in: "){
      working_directory <- substr(working_directory, 13, nchar(working_directory))
    }

    if(grepl(working_directory, file_path_t)){
      return(gsub(working_directory, "", file_path_t))
    }else{
      return(file_path)
    }

  }

}

warning_plot = function(){
  warning_plot <- ggplot(data = NULL, aes(x = 1, y = 1, label = "There is no plot ready to load.\n\nYou may not have selected a database or the database may be empty.\n\nPlease build or select a RecruitPlotEasy database and try again.\nIf you just switched from a still plot to an interactive one, just hit the 'load genome' button again."))+
    geom_text(size = 6) +
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  return(warning_plot)
}

static_plot = function(database_handle){

  x_ax = as.numeric(unlist(database_handle$plot_x_axis))
  y_ax = as.numeric(as.character(database_handle$describe_y))

  bin_sizes = as.numeric(as.character(database_handle$widths))
  norm_factor = min(bin_sizes)/ bin_sizes

  main = data.table(t(database_handle$plot_ready))
  colnames(main) = as.character(database_handle$describe_y)
  main$x_axis = x_ax
  main$bin_size = bin_sizes
  main$norm = norm_factor
  main = melt(main, id.vars = c("x_axis", "bin_size", "norm"))
  colnames(main)[4:5] = c("y_axis", "bp_count")
  main$y_axis = as.numeric(as.character(main$y_axis))
  main[, normalized_count := bp_count * norm, ]

  depth_in = data.table(depth = database_handle$depth_chart_data$in_group, label = "depth.in", pos = x_ax)
  depth_out = data.table(depth = database_handle$depth_chart_data$out_group, label = "depth.out", pos = x_ax)
  #zero depth vals
  depth_zin = depth_in[depth  == 0, ]
  depth_zout = depth_out[depth  == 0, ]
  depth_zin$label = "depth.zin"
  depth_zout$label = "depth.zout"
  depth_zin$depth = -Inf
  depth_zout$depth = -Inf

  depth_in$depth[depth_in$depth == 0] = NA
  depth_out$depth[depth_out$depth == 0] = NA

  total_depth = rbindlist(list(depth_in, depth_out))
  zero_depth = rbindlist(list(depth_zin, depth_zout))

  #order of magnitude
  oom = floor(log10(max(x_ax)))
  #Select unit
  unit = ceiling(oom/2)

  denom = c(1, 1e3, 1e6, 1e9)[unit]
  unit_name = c("(bp)", "(kbp)","(mbp)","(gbp)")[unit]

  #Find good break positions
  position_labels = pretty(c(0, max(x_ax)), 6)
  #Make sure they've got no decimals
  position_labels = as.integer(position_labels[1:(length(position_labels)-1)])

  placement = position_labels

  #format names nicely
  position_labels = position_labels / denom

  #Add readability as needed.
  position_labels = format(position_labels, scientific = F, nsmall = 1, big.mark = ",")

  main_plot = ggplot(main, aes(x = x_axis, y = y_axis, fill=log10(normalized_count)))+
    scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
    ylab("Percent Identity") +
    xlab(paste("Position in Genome", unit_name)) +
    scale_y_continuous(expand = c(0, 0), limits = c(NA, 100 + database_handle$height/2)) +
    scale_x_continuous(expand = c(0, 0), breaks = placement, labels = position_labels) +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    geom_tile()


  main_plot = main_plot + annotate("rect", xmin = 0, xmax = max(x_ax),
                    ymin = database_handle$in_grp,
                    ymax = 100+database_handle$height/2, fill = "darkblue", alpha = .15)


  contig_starts = lapply(database_handle$describe_x, function(x){
    return(x[[1]])
  })

  contig_ends = lapply(database_handle$describe_x, function(x){
    return(x[[2]])
  })

  #On a static plot, this is useless info. The horiz. res is too low
  #main_plot = main_plot + geom_vline(xintercept = ends$V1/bp_div[-nrow(ends)], col = "#AAAAAA40")

  group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.zin = "darkblue", depth.zout = "lightblue")

  depth_chart = ggplot(total_depth, aes(x = pos, y = depth, color = label)) +
    geom_step(alpha = 0.75) +
    scale_x_continuous(expand=c(0,0), breaks = placement, labels = position_labels)+
    scale_y_continuous(limits = c(min(total_depth$depth, na.rm=T)-0.75, NA), )+
    theme(legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = group.colors) +
    ylab("Avg. Depth")

  #If we have zero depth data, add it
  if(nrow(zero_depth)> 0){

    depth_chart = depth_chart +
      geom_segment(data = zero_depth, aes(x = pos, xend=pos, y = -Inf, yend = min(total_depth$depth ,na.rm=T)-0.5, color = label))

  }

  hist_plot = data.table(ct = database_handle$bases_by_pct_id, y = y_ax)
  buffer = hist_plot[1,]
  buffer$y = buffer$y + database_handle$height
  hist_plot = rbindlist(list(buffer, hist_plot))

  #breaks = scales::pretty_breaks(n = 3)
  bp_placement = pretty(c(0, max(hist_plot$ct)), n = 3)[1:3]
  bp_labs = format(as.integer(bp_placement), scientific = F, nsmall = 1, big.mark = ",")

  if(database_handle$plot_linear){

    bp_count_hist = ggplot(hist_plot, aes(x = y - database_handle$height/2, y = ct)) +
      geom_step()+
      geom_point(aes(x = y, y = ct), inherit.aes = F)+
      #geom_density()+
      #geom_line()+
      scale_y_continuous(expand = c(0,0), limits = c(0, 1.01*max(hist_plot$ct)), breaks = bp_placement, labels = bp_labs) +
      scale_x_continuous(expand = c(0, 0), limits = c(NA, 100+database_handle$height/2)) +
      theme(legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.background = element_rect(fill = "#EEF7FA"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14))+
      ylab("Base Pair Count by % ID")+
      annotate("rect", xmin = database_handle$in_grp-database_handle$height/2, xmax = 100+database_handle$height/2, ymin = 0, ymax = Inf, fill = "darkblue", alpha = .15)+
      coord_flip()
  }else{

    bp_count_hist = ggplot(hist_plot, aes(x = y - database_handle$height/2, y = ct)) +
      geom_step()+
      geom_point(aes(x = y, y = ct), inherit.aes = F)+
      scale_y_continuous(expand = c(0,0), limits = c(1, 1.01*max(hist_plot$ct)), trans='log10') +
      scale_x_continuous(expand = c(0, 0), limits = c(NA, 100+database_handle$height/2)) +
      theme(legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.background = element_rect(fill = "#EEF7FA"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14))+
      ylab("Log 10 Base Pair Count by % ID")+
      annotate("rect", xmin = database_handle$in_grp-database_handle$height/2, xmax = 100+database_handle$height/2, ymin = 0, ymax = Inf, fill = "darkblue", alpha = .15)+
      coord_flip()
  }


  total_depth = total_depth[!is.na(total_depth$depth),]

  #Reverse4
  #total_depth$label = as.factor(total_depth$label)
  #levels(total_depth$label) = c("depth.out", "depth.in", "depth.zout", "depth.zin")
  group.colors = c(depth.out = "lightblue", depth.in = "darkblue", depth.zout = "lightblue", depth.zin = "darkblue")
  #total_depth = total_depth[nrow(total_depth):1,]
  group_alphas = c(depth.in = 1, depth.out = 0.7)

  seqdepth.lim = range(total_depth$depth)
  hist_binwidth = (seqdepth.lim[2] - seqdepth.lim[1])/199

  depth_hist = ggplot(total_depth, aes(x = depth, fill = label, alpha = label)) +
    geom_histogram(binwidth = hist_binwidth, position = "identity") +
    scale_fill_manual(values = group.colors) +
    scale_alpha_manual(values = group_alphas)+
    scale_x_continuous(limits = c(min(total_depth$depth, na.rm=T)-0.75, NA))+
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "white"),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    xlab(label = element_blank()) +
    ylab(label = element_blank())

  #TODO
  #This code is temporarily disabled so that we can fix an enveomics.R dependency
  # if(database_handle$show_peaks){
  #
  #   depth_in = data.table(depth = database_handle$depth_chart_data$in_group, label = "depth.in", pos = x_ax)
  #
  #
  #   h.breaks <- seq(seqdepth.lim[1], seqdepth.lim[2], length.out = 200)
  #   h.mids <- (h.breaks[-1] + h.breaks[-length(h.breaks)])/2
  #
  #   min_info <- list()
  #
  #   #min_info$pos.breaks = c(0, base$contiguous_end[match(ddSave$seq_pos[ddSave$group_label == "depth.in"], base$seq_pos)])
  #   #min_info$pos.counts.in = ddSave$V1[ddSave$group_label == "depth.in"]
  #
  #   ends = lapply(database_handle$describe_x, function(x){
  #     return(x[[2]])
  #   })
  #
  #   if(length(ends) > 1){
  #
  #     for(i in 2:length(ends)){
  #
  #       ends[[i]] = ends[[i]] + max(ends[[i-1]])
  #
  #     }
  #   }
  #
  #   ends = unname(unlist(ends))
  #   ends = c(0, ends)
  #
  #   min_info$pos.breaks = ends
  #   min_info$pos.counts.in = as.integer(depth_in$depth * bin_sizes)
  #
  #   ### Internal ancilliary function (see `enve.RecPlot2.Peak`).
  #   enve.recplot2.__peakHist <- function(x, mids, counts=TRUE){
  #     d.o <- x$param.hat
  #     if(length(x$log)==0) x$log <- FALSE
  #     if(x$log){
  #       d.o$x <- log(mids)
  #     }else{
  #       d.o$x <- mids
  #     }
  #     prob  <- do.call(paste('d', x$dist, sep=''), d.o)
  #     if(!counts) return(prob)
  #     if(length(x$values)>0) return(prob*length(x$values)/sum(prob))
  #     return(prob*x$n.hat/sum(prob))
  #   }
  #
  #
  #   #This is finicky
  #
  #   try({
  #
  #     peaks <- enve.recplot2.findPeaks(min_info, method = "em")
  #
  #   })
  #
  #   if(length(peaks) > 0){
  #
  #   try({
  #     dpt <- signif(as.numeric(lapply(peaks, function(x) x$seq.depth)), 2)
  #     frx <- signif(100 * as.numeric(lapply(peaks,function(x) ifelse(length(x$values) == 0, x$n.hat, length(x$values))/x$n.total)), 2)
  #
  #     if (peaks[[1]]$err.res < 0) {
  #       err <- paste(", LL:", signif(peaks[[1]]$err.res,3))
  #     }  else {
  #       err <- paste(", err:", signif(as.numeric(lapply(peaks, function(x) x$err.res)), 2))
  #     }
  #     labels <- paste(letters[1:length(peaks)], ". ", dpt, "X (", frx, "%", err, ")", sep = "")
  #
  #     peak_counts <- lapply(peaks, enve.recplot2.__peakHist, h.mids)
  #
  #     plot_breaks = h.breaks[-length(h.breaks)]
  #
  #     gg_peak_info <- data.table(plot_breaks = rep(plot_breaks, length(peak_counts)), count = unlist(peak_counts), grp = rep(labels, each = length(plot_breaks)))
  #
  #     depth_hist = depth_hist + geom_line(data = gg_peak_info, aes(x = plot_breaks, y = count, color = grp, group = grp), inherit.aes = F, color = "red", lwd = 1.13)
  #
  #   })
  #
  # }
  # }

  depth_hist = depth_hist + coord_flip()

  #Same with this.
  # if(database_handle$show_peaks){
  #   if(length(peaks) > 0){
  #
  #   try({
  #     o_max <- max(table(findInterval(depth_in$depth, h.breaks)))*.68
  #     x_start = max(h.breaks)
  #
  #     if(length(labels) > 0){
  #     labels = paste(labels, collapse="\n")
  #     }
  #     depth_hist = depth_hist + annotate("text", label = labels, y = o_max, x = x_start, size = 5, hjust = 1, vjust = 1)
  #
  #   })
  #
  #   }
  # }


  total_plot = plot_grid(depth_chart, depth_hist, main_plot, bp_count_hist, align = "hv", ncol = 2, rel_widths = c(2.7, 1), rel_heights = c(1, 2.3))


  #Format data for printout
  total_depth = rbindlist(list(total_depth, zero_depth))
  total_depth[label == "depth.zout", label := "depth.out"]
  total_depth[label == "depth.zin", label := "depth.in"]
  total_depth[depth < 0, depth := 0, ]
  total_depth[is.na(depth), depth := 0, ]
  colnames(total_depth) = c("Avg_Depth_of_Coverage", "In_Group_or_Out_Group", "Position in Genome")

  colnames(main) = c("Genome_Pos_Midpt", "Bin_width", "Normalization_Multiplier", "Pct_ID_Bin", "Base_Pair_count", "Normalized_Base_Count")
  colnames(total_depth) = c("Avg_Depth_of_Coverage", "In_Group_or_Out_Group", "Genome_Pos_Midpt")
  colnames(hist_plot) = c("Count_of_bp_at_Pct_ID", "Pct_ID_Bin")
  #Remove buffer
  hist_plot = hist_plot[2:nrow(hist_plot),]

  query = database_handle$constructed_query

  mags = unlist(database_handle$contigs_in_mag)

  query = paste(query, "; Genome ID List:", paste(names(mags), ":", unname(mags)))

  writeable_data = list(main, total_depth, hist_plot, query)

  return(list(total_plot, main, total_depth, hist_plot, query))

}

interactive_plot = function(database_handle){
  x_ax = as.numeric(unlist(database_handle$plot_x_axis))
  y_ax = as.numeric(as.character(database_handle$describe_y))

  bin_sizes = as.numeric(as.character(database_handle$widths))
  norm_factor = min(bin_sizes)/ bin_sizes

  main = data.table(t(database_handle$plot_ready))
  colnames(main) = as.character(database_handle$describe_y)
  main$x_axis = x_ax
  main$bin_size = bin_sizes
  main$norm = norm_factor
  main = melt(main, id.vars = c("x_axis", "bin_size", "norm"))
  colnames(main)[4:5] = c("y_axis", "bp_count")
  main$y_axis = as.numeric(as.character(main$y_axis))
  main[, normalized_count := bp_count * norm, ]
  main$annot = unlist(database_handle$main_annotations)

  depth_in = data.table(depth = database_handle$depth_chart_data$in_group, label = "depth.in", pos = x_ax, annot = database_handle$depth_annotations$in_group)
  depth_out = data.table(depth = database_handle$depth_chart_data$out_group, label = "depth.out", pos = x_ax, annot = database_handle$depth_annotations$out_group)
  #zero depth vals
  depth_zin = depth_in[depth  == 0, ]
  depth_zout = depth_out[depth  == 0, ]
  depth_zin$label = "depth.zin"
  depth_zout$label = "depth.zout"
  depth_zin$depth = -Inf
  depth_zout$depth = -Inf

  depth_in$depth[depth_in$depth == 0] = NA
  depth_out$depth[depth_out$depth == 0] = NA

  total_depth = rbindlist(list(depth_in, depth_out))
  zero_depth = rbindlist(list(depth_zin, depth_zout))

  #order of magnitude
  oom = floor(log10(max(x_ax)))
  #Select unit
  unit = ceiling(oom/2)

  denom = c(1, 1e3, 1e6, 1e9)[unit]
  unit_name = c("(bp)", "(kbp)","(mbp)","(gbp)")[unit]

  #Find good break positions
  position_labels = pretty(c(0, max(x_ax)), 6)
  #Make sure they've got no decimals
  position_labels = as.integer(position_labels[1:(length(position_labels)-1)])

  placement = position_labels

  #format names nicely
  position_labels = position_labels / denom

  #Add readability as needed.
  position_labels = format(position_labels, scientific = F, nsmall = 1, big.mark = ",")

  main_plot = ggplot(main, aes(x = x_axis, y = y_axis, fill=log10(normalized_count),  text = annot))+
    scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
    ylab("Percent Identity") +
    xlab(paste("Position in Genome", unit_name)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), breaks = placement, labels = position_labels) +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    annotate("rect", xmin = 0, xmax = max(x_ax),
               ymin = database_handle$in_grp,
               ymax = 100+database_handle$height/2, fill = "darkblue", alpha = .15)+
    geom_tile()

  main_plot = ggplotly(main_plot, dynamicTicks = T, tooltip = c("text")) %>%
    layout(plot_bgcolor = "grey90") %>%
    style(hoverinfo = "none", traces = c(1))

  contig_starts = lapply(database_handle$describe_x, function(x){
    return(x[[1]])
  })

  contig_ends = lapply(database_handle$describe_x, function(x){
    return(x[[2]])
  })

  #On a static plot, this is useless info. The horiz. res is too low
  #main_plot = main_plot + geom_vline(xintercept = ends$V1/bp_div[-nrow(ends)], col = "#AAAAAA40")

  group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.zin = "darkblue", depth.zout = "lightblue")

  depth_chart = ggplot(total_depth, aes(x = pos, y = depth, color = label, group = label, text = annot)) +
    geom_step(alpha = 0.75) +
    scale_x_continuous(expand=c(0,0), breaks = placement, labels = position_labels)+
    scale_y_continuous(limits = c(min(total_depth$depth, na.rm=T)-0.75, NA), )+
    theme(legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = group.colors) +
    ylab("Avg. Depth")

  #If we have zero depth data, add it
  if(nrow(zero_depth)> 0){

    depth_chart = depth_chart +
      geom_segment(data = zero_depth, aes(x = pos, xend=pos, y = -Inf, yend = min(total_depth$depth ,na.rm=T)-0.5, color = label))

  }

  a <- list(
    range = c(min(total_depth$depth, na.rm = T), max(total_depth$depth, na.rm = T)),
    showticklabels = TRUE,
    exponentformat = "e"
  )

  depth_chart = ggplotly(depth_chart, dynamicTicks = T, tooltip = c("text")) %>%
    layout(plot_bgcolor = "grey90", yaxis = a)

  hist_plot = data.table(ct = database_handle$bases_by_pct_id, y = y_ax)
  buffer = hist_plot[1,]
  buffer$y = buffer$y + database_handle$height
  hist_plot = rbindlist(list(buffer, hist_plot))

  #COnvert to interactive.

  #Format data for printout
  total_depth = rbindlist(list(total_depth, zero_depth))
  total_depth[label == "depth.zout", label := "depth.out"]
  total_depth[label == "depth.zin", label := "depth.in"]
  total_depth[depth < 0, depth := 0, ]
  total_depth[is.na(depth), depth := 0, ]

  total_depth[, annot := NULL]
  main[, annot := NULL]
  colnames(main) = c("Genome_Pos_Midpt", "Bin_width", "Normalization_Multiplier", "Pct_ID_Bin", "Base_Pair_count", "Normalized_Base_Count")
  colnames(total_depth) = c("Avg_Depth_of_Coverage", "In_Group_or_Out_Group", "Genome_Pos_Midpt")
  colnames(hist_plot) = c("Count_of_bp_at_Pct_ID", "Pct_ID_Bin")
  #Rempve buffer
  hist_plot = hist_plot[2:nrow(hist_plot),]

  query = database_handle$constructed_query

  mags = unlist(database_handle$contigs_in_mag)

  query = paste(query, "; Genome ID List:", paste(names(mags), ":", unname(mags)))

  overplot <- subplot(list(depth_chart, main_plot), nrows = 2, shareX = T, heights = c(1/3, 2/3))

  return(list(overplot, main, total_depth, hist_plot, query))

}

additional_stats = function(database_handle){

  breadths = database_handle$breadth_for_this_genome
  depths = database_handle$depth_for_this_genome
  tad_n = database_handle$tad_middle
  anir = database_handle$ani_r

  basic_set = unlist(breadths)

  name_id = names(basic_set)

  contigs = unlist(lapply(name_id, function(x){

    return(strsplit(x, split = "[.]")[[1]][1])

  }))

  lens = nchar(contigs)

  pct_ids = as.numeric(substr(name_id, lens+2, nchar(name_id)))

  combined_data = data.table(contig = contigs, y = pct_ids, breadth = unname(unlist(breadths)), depth = unname(unlist(depths)), anir = unname(unlist(anir)))
  combined_data[, breadth := breadth * 100, ]

  combined_data = combined_data[y >= database_handle$in_grp ,]

  combined_data = combined_data[, .SD[which.min(y)], by = contig]
  combined_data[, breadth := round(breadth, 2),]
  combined_data[, depth := round(depth, 2),]
  combined_data[, anir := round(anir, 2),]

  combined_data[, y := NULL]

  b = combined_data$breadth
  d = combined_data$depth
  ani = combined_data$anir
  n = combined_data$contig

  b = paste("Breadths:", paste(b, collapse = ", "))
  d = paste(paste0("TAD-", database_handle$tad_middle), paste(d, collapse = ", "))
  ani = paste("ANIr:", paste(ani, collapse = ", "))
  n = paste("Genome/Contigs:", paste(n, collapse = ', '))

  combined_data = NULL

  pretty = paste0("The current within-population group (dark blue) has the following:\n\n", n, "\n\n", ani, "\n\n", b, "\n\n", d)

  return(pretty)

}


#This is the GUI function

recplot_UI <- function(){
  recplot_py <- get_python()

  #System/choices issue
  {

    platform = import("platform")
    system <- platform$system()

    format_choices <- c("Tabular BLAST" = "blast", "SAM" = "sam", "BAM" = "bam")

    gene_choices <- c("Prodigal GFF" = "prodigal")
  }

  #UI
  {
    ui <- fluidPage(

      tags$head(
        tags$style(
          HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
          )
        )
      ),

      useShinyalert(),

      tabsetPanel(id = "tabs",
                  tabPanel("Database Creation",
                           fluidRow(
                             #todo

                             column(1),
                             column(4,
                                    h3("Change your working directory"),
                                    actionButton('dir', 'Choose Directory', icon = icon("folder-open")),
                                    textInput("cur_dir",label = NULL, value = paste("Working in:", getwd()), width = '100%'),

                                    h3("Create a new database"),
                                    textInput("dbname",label = "Name a new database to create", value = "Enter name here.", width = '100%'),
                                    actionButton('db' , "Create the database", icon("coins")),

                                    h3("Or work with an existing database"),
                                    actionButton('exist_db', 'Select an existing DB', icon = icon("coins")),
                                    verbatimTextOutput("exist_dbname"),


                                    h3("Add genomes to your database"),
                                    actionButton('contigs', 'Select Reference Genomes', icon = icon("file-upload")),
                                    actionButton('add_contigs', 'Add to DB', icon = icon("plus-square")),
                                    actionButton("what_are_mags", "Info", icon = icon("question-circle")),
                                    textInput("contig_file",label = NULL, value = "No genomes selected.", width = '100%'),
                                    checkboxInput("one_mag", label = "This is a single binned genome.", value = F),

                                    h3("Add reads to your database"),
                                    actionButton('reads', 'Select Mapped Reads', icon = icon("file-upload")),
                                    actionButton('add_reads', 'Add to DB', icon = icon("plus-square")),
                                    actionButton("read_info", "Info", icon = icon("question-circle")),
                                    textInput("read_file",label = NULL, value = "No mapped read file selected.", width = '100%'),

                                    h3("Add genes to your database"),
                                    actionButton('genes', 'Select a Prodigal GFF to add genes', icon = icon("file-upload")),
                                    actionButton('add_genes', 'Add to DB', icon = icon("plus-square")),
                                    actionButton("gene_info", "Info", icon = icon("question-circle")),
                                    textInput("gene_file",label = NULL, value = "No gene file selected.", width = '100%'),

                                    #Tooltips

                                    bsTooltip("dir", "(Optional) Select a working directory. The database and any saved plots will be placed here.", placement = "right"),

                                    bsTooltip("dbname", "Name your database. A .db extension will be added to the end of the name you give it.", placement = "right"),
                                    bsTooltip("db", "Click me after selecting all input files to create your database.", placement = "right"),

                                    bsTooltip("exist_db", "Click me to load an existing RecruitPlotEasy database.", placement = "right"),


                                    bsTooltip("contigs", "Select a FASTA format file containing multiple genome sequences or multiple contigs belonging to a single MAG.", placement = "right"),
                                    bsTooltip("one_mag", "If you binned an assembly and have a single multi-FASTA containing the sequences for this MAG as your reference genomes, check this box, leave the association file blank, and they will all appear on a single plot.", placement = "right"),
                                    bsTooltip("what_are_mags", "Click me for a longer explanation on MAGs.", placement = "right"),
                                    bsTooltip("add_contigs", "Add genomes or a MAG to the DB. You must click this or the genomes will not be added.", placement = "right"),


                                    bsTooltip("reads", "Select a mapped read file. These reads should be mapped to the genomes in the reference genomes file selected above.", placement = "right"),
                                    bsTooltip("add_reads", "Add reads to the DB. You must click this or the reads will not be added.", placement = "right"),
                                    bsTooltip("read_info", "Click me for more info on reads", placement = "right"),

                                    bsTooltip("genes", "Select genes to add to the DB. These must be produced with Prodigal using the '-f gff' format option. Only works if you've also added the genomes the genes were predicted from.", placement = "right"),
                                    bsTooltip("add_genes", "Add genes to the DB. You must click this or the genes will not be added.", placement = "right"),
                                    bsTooltip("gene_info", "Click me for more info on genes.", placement = "right"),
                             ),
                             column(1),

                             column(5,
                                    h2("Build Report"),
                                    verbatimTextOutput("message")
                             ),
                             column(1)

                           )
                  ),

                  tabPanel("Database Management",
                           fluidRow(
                             column(1),
                             column(4,
                                    h2("Current Database:"),
                                    verbatimTextOutput("db_stat"),
                                    br(),
                                    actionButton('show_samples', 'Show samples in this database', icon = icon("coins"), width = '65%'),

                                    h2("Set Advanced Plot Settings"),
                                    br(),
                                    selectInput("plot_genes", "Plot genes?", choices = c("Genome" = F, "Genes" = T), selected = F, width = '65%'),
                                    br(),
                                    actionButton('current_args', 'Show me the current plotting parameters', icon = icon("coins"), width = '65%'),
                                    br(),
                                    br(),
                                    selectInput("local_id", "Global or local percent ID", choices = c("Local" = T, "Global" = F), selected = T, width = '65%'),

                                    #min_align_len = 50, min_pct_align = 90,

                                    selectInput("multiple_aligns", "Allow multiple alignments per read?", choices = c("No" = T, "Yes" = F), selected = T, width = '65%'),


                                    selectInput("bh_crit", "Select best hit criteria", choices = c("Local Percent ID" = "local_percent_ID", "Global Percent ID" = "global_percent_ID", "Alignment Length" = "alignment_length", "Percent Alignment" = "percent_alignment"), selected = "local_percent_ID", width = '65%'),

                                    #min_align_len = 50, min_pct_align = 90,
                                    numericInput("min_align_len", "Minimum aligned bases", value = 50, min = 0, max = 500, step = 5, width = '65%'),
                                    numericInput("min_pct_align", "Minimum percent alignment", value = 90, min = 0, max = 100, step = 5, width = '65%'),

                                    numericInput("tadn", "Set the value of TAD-N", value = 80, min = 0.5, max = 100, step = 5, width = '65%'),

                                    actionButton('export_reads', "Export current read cart", width = '65%'),
                                    bsTooltip('export_reads', "You can add reads to a cart on the plotting page and then save them to files with this button. Does nothing if the cart is empty.", placement = 'right'),

                                    bsTooltip("plot_genes", "Select between plotting your genome/MAG in evenly sized windows, or use gene starts and stops as the windows. You must have already added genes for the genome you wish to plot to see genes."),

                                    bsTooltip("current_args", "Report the current advanced plotting settings", placement = "right"),

                                    bsTooltip("local_id", "Local %ID = (# bases matching reference) / (# bases aligned) Global %ID = (# bases matching reference) / (length of read)", placement = "right"),
                                    bsTooltip("multiple_aligns", "Reads can align either to multiple locations in the genome or to multiple genomes. If multiple alignment is not allowed (No), then a particular read will only be allowed to align once across both of these possibilities.", placement = "right"),
                                    bsTooltip("bh_crit", "Select the criteria to determine the the best-matching alignment for each read. Has no effect if multiple alignments are allowed.", placement = 'right'),

                                    bsTooltip("min_align_len", "Minimum number of bases aligning to the reference genome to include a read.", placement = 'right'),
                                    bsTooltip("min_pct_align", "Minimum percent of the read which must align to the reference to be included.", placement = 'right'),
                                    bsTooltip("tadn", "TAD-N (Truncated Average Depth) the average depth of coverage over the middle N quantiles, excluding the lower and upper (100-N)/2 percentile depth loci, i.e. the lowest (incl. zero-depth) and highest depth positions in the genomes are excluded.", placement = "right"),

                             ),
                             column(1),
                             column(5,
                                    h2("Database Status"),
                                    verbatimTextOutput("message2")
                             ),
                             column(1)

                           )
                  ),

                  tabPanel("Recruitment Plot",
                           sidebarPanel(
                             width = 3,

                             h4("Choose a Genome"),

                             selectInput("samples", "(1) Select a sample in the database", selected = NULL, choices = NULL),

                             selectInput("mags_in_db", "(2) Select a genome in the sample", selected = NULL, choices = NULL),

                             h4("Select Bin Resolution"),

                             numericInput("height", "(3) Pct. ID Resolution", min = 0.1, max = 3, value = 0.5, step = 0.1),

                             numericInput("width", "(4) Genome Resolution", min = 75, max = 5000, value = 1000),

                             h4("Fine Tuning"),

                             numericInput("in_group_min_stat", "(5) In-Group Pct. ID", min = 50, max = 100, value = 95),

                             selectInput("linear_stat", "(6) BP Histogram Scale", choices = c("Linear" = T, "Logarithmic" = F), selected = T),

                             #TODO
                             #This is also disabled and numbers are shifted below to match.
                             #checkboxInput("show_peaks", "(7) Display Depth Peaks?"),

                             h4("Load Selected Genome"),

                             #actionButton('get_a_mag', '(8) View Selected Genome', icon = icon("jedi-order")),
                             actionButton('get_a_mag', '(7) View Selected Genome', icon = icon("jedi-order")),
                             actionButton('addtl_stats_static', 'Show additional stats', icon = icon("eye")),
                             bsTooltip('addtl_stats_static', 'Calculates the average percent identity, breadth of coverage, and truncated average depth for the current plot\'s within-population group', placement = 'right'),

                             h4("Output Data"),

                             textInput("pdf_name", "Name and save current plot."),
                             actionButton("print_stat", "Save this plot", icon = icon("save")),

                             actionButton("add_cart_stat", "Add reads to cart", icon = icon("save")),
                             bsTooltip("add_cart_stat", "Add reads from the current within-population (dark blue) to the cart.", placement = "right"),


                             bsTooltip("samples", "This menu contains a list of samples within the database. Select one, and the genome field will be populated with the genomes found in that sample.", placement = "right"),
                             bsTooltip("mags_in_db", "This menu contains the set of genomes in currently selected sample. Select one, then select resolution parameters.", placement = "right"),
                             bsTooltip("width", "Approximate number of base pairs in each genome window. The Recruitment Plot attempts to normalize bin width for each contig to this size. Lower values = higher resolution, but is slower. Higher values = lower resolution, but is faster.", placement = "right"),
                             bsTooltip("height", "Controls the resolution of percent identity to the reference. Lower values here will result in finer resolution, but will be slower. Hint: The default 0.5% window means a resolution of 1 base pair mismatch per 200 bases; finer resolution is probably uneccessary.", placement = "right"),
                             #bsTooltip("low_bound", "Reads mapping below this percent identity will not be included in the current recruitment plot.", placement = "right"),
                             bsTooltip("get_a_mag", "Click this to load the current Genome into the viewer and plot it. Please wait for the plot to appear after clicking this. This loads the data for all tabs.", placement = "right"),
                             bsTooltip("in_group_min_stat", "Controls the lower edge of the shaded region in the recruitment plot's lower panels. Reads mapping at or above this percent identity are regarded as the \"in-group\" for the Recruitment Plot, and are represented by the dark blue lines in the upper panels.", placement = "right"),
                             bsTooltip("linear_stat", "Causes the lower right panel to display base pair counts per percent identity bin in linear scale or log scale.", placement = "right"),
                             bsTooltip("print_stat", "After loading a plot (meaning you should be able to see it), add a name in the associated text box and then click this to print a PDF and save data for the current plot", placement = "right"),


                             #Same disabling as above.
                             #bsTooltip("show_peaks", "Calculate and overlay peaks for the depth of coverage histogram (top right panel)", placement = "right"),

                             #bsTooltip("recplot_main", "Bottom left panel: a 2-D histogram of the counts of base pairs falling into a bin defined by position in the genome (x-axis) and percent identity (y-axis) Bins are as wide as the genome resolution parameter                             if viewing contigs, and cover genes & intergenic regions in contiguous chunks if viewing genes. The shaded section is the current in-group, which is the dark blue line on the top two plots Top left panel: Average sequencing depth for each x-axis bin on in the bottom left panel. Dark blue corresponds to depth of coverage for bins in the in-group, and light blue to the out-group. Segments at the bottom of the plot have zero coverage. Top right panel: a histogram of the depths of coverage observed in the corresponding in/out group in the sequencing depths chart (top left). If peaks are selected, they correspond to the estimates of the genome's average sequencing depth. Bottom right panel: a histogram of the number of bases falling into each percent identity bin across the entire genome, displayed in linear or log scale depending on your selection.", trigger = "click", placement = "left")
                           ),
                           mainPanel(id = "recplot_main",
                                     #Spacing
                                     fluidRow(
                                       column(12,
                                              div(style = "height:60px;background-color: white;", "")
                                       )
                                     ),
                                     fluidRow(

                                       plotOutput("read_recruitment_plot", height = "850px")

                                     )

                           )
                  ),
                  tabPanel("Interactive Plot",
                           sidebarPanel(
                             width = 3,

                             h3("Choose a Genome"),

                             selectInput("samples_interact", "(1) Select a sample in the database", selected = NULL, choices = NULL),

                             selectInput("mags_in_db_interact", "(2) Select a Genome in the sample", selected = NULL, choices = NULL),

                             h3("Select Bin Resolution"),



                             numericInput("height_interact", "(3) Pct. ID Resolution", min = 0.05, max = 3, value = 0.5),


                             numericInput("width_interact", "(4) Genome Resolution", min = 75, max = 5000, value = 1000),


                             h3("Fine Tuning"),

                             numericInput("in_group_min_interact", "(5) In-Group Pct. ID", min = 50, max = 99.5, value = 95),

                             h3("Load Selected Genome"),

                             actionButton('get_a_mag_interact', '(6) View selected Genome', icon = icon("jedi-order")),
                             actionButton('addtl_stats_interact', 'Show additional stats', icon = icon("eye")),
                             bsTooltip('addtl_stats_interact', 'Calculates the average percent identity, breadth of coverage, and truncated average depth for the current plot\'s within-population group', placement = 'right'),

                             h4("Output Data"),

                             textInput("pdf_name_interact", "Name and save current plot."),
                             actionButton("print_interact", "Save interactive HTML", icon = icon("save")),

                             actionButton("add_cart_interact", "Add reads to cart", icon = icon("save")),
                             bsTooltip("add_cart_interact", "Add reads from the current within-population (dark blue) to the cart.", placement = "right"),

                             bsTooltip("pdf_name_interact", "Name the current interactive plot as an HTML page before saving it.", placement = "right"),
                             bsTooltip("print_interact", "After loading a plot (meaning you should be able to see it), add a name in the associated text box and then click this to save an HTML and data for the current plot", placement = "right"),

                             bsTooltip("samples_interact", "This menu contains a list of samples within the database. Select one, and the genome field will be populated with the genomes found in that sample.", placement = "right"),
                             bsTooltip("mags_in_db_interact", "This menu contains the set of genomes in currently selected sample. Select one, then select resolution parameters.", placement = "right"),
                             bsTooltip("width_interact", "Approximate number of base pairs in each genome window. The Recruitment Plot attempts to normalize bin width for each contig to this size. Lower values = higher resolution, but is slower. Higher values = lower resolution, but is faster.", placement = "right"),
                             bsTooltip("height_interact", "Controls the resolution of percent identity to the reference. Lower values here will result in finer resolution, but will be slower. Hint: The default 0.5% window means a resolution of 1 base pair mismatch per 200 bases; finer resolution is probably uneccessary.", placement = "right"),

                             bsTooltip("get_a_mag_interact", "Click this to load the current genome into the viewer and plot it. Please wait for the plot to appear after clicking this.", placement = "right"),
                             bsTooltip("in_group_min_interact", "Controls the lower edge of the shaded region in the recruitment plot's lower panels. Reads mapping at or above this percent identity are regarded as the \"in-group\" for the Recruitment Plot, and are represented by the dark blue lines in the upper panel.", placement = "right")
                           ),
                           mainPanel(
                             column(12,
                                    plotlyOutput("Plotly_interactive", height = "850px")
                             )
                           )

                  #),
                  # tabPanel("Additional Statistics", mainPanel(
                  #   column(12,
                  #     plotOutput("additional_stats", height = "850px")
                  #   )
                  #
                  # )

                  )
      )
      #End fluid page
    )

  }

  return(ui)
}

recplot_server <- function(input, output, session) {

  initial_message <- "Welcome to RecruitPlotEasy!\n\nThis page is meant to allow you to build a new database.\nYou must build or select a database before you can plot anything.\n\nTo build your database, first create the database using the create database inputs,\nthen add any genomes, MAGs, reads, and genes using the other buttons.\nYou can add as many files as you want, here."
  initial_message2 <- "This page is for taking a queek peek at your database and managing more complicated plotting options.\nIf you just created a database, it should already be loaded here."
  dbstat = "No database selected."

  output$message <- renderText(initial_message)
  output$message2 <- renderText(initial_message2)
  output$db_stat <- renderText(dbstat)
  output$exist_dbname <- renderText(dbstat)

  mod_handle = ""
  plot_handle = ""

  directory <- "No directory selected. Try again?"
  reads <- "No reads selected. Try again?"
  contigs <- "No genomes selected. Try again?"
  mags <- "No association file selected. Try again?"
  new_samp <- "No new sample. Try again?"
  db <- "No existing database selected. Try again?"
  new_genes <- "No genes selected. Try again?"

  stat_plot = NA
  interact_plot = NA
  main_dat = NA
  depth_dat = NA
  hist_dat = NA
  query = NA

  readied = FALSE

  exist_db <- "No existing database selected. Try again?"

  samples_in_db <- "No database selected or built yet."

  recplot_py <- get_python()

  #TODO - this is the contig end cutoff that removes the first & last 75 bp from consideration in avg. depth.
  #These will need to be later incorporated into the interactive inputs
  trunc_degree <- as.integer(75)
  wp = warning_plot()

  output$read_recruitment_plot <- renderPlot(wp)
  #output$additional_stats <- renderPlot(wp)

  wp2 = ggplotly(wp)
  output$Plotly_interactive <- renderPlotly(wp2)

  #Database building
  #Directories
  observeEvent(input$dir, {

    tryCatch({
      directory <<- choose_directory()
    },
    error = function(cond){
      directory <<- "No directory selected. Try again?"
      updateTextInput(session, "cur_dir", value = paste(directory))
    })

    tryCatch({
      setwd(directory)
    },
    error = function(cond){
      directory <<- "No directory selected. Try again?"
      updateTextInput(session, "cur_dir", value = paste(directory))
    })

    if(length(directory) == 0){
      directory <<- "No directory selected. Try again?"
    }


    if(directory == "No directory selected. Try again?"){
      updateTextInput(session, "cur_dir", value = paste(directory))
    }else{
      updateTextInput(session, "cur_dir", value = paste("Working in:", directory))
    }

  })

  #Create a database or select and/or modify it
  observeEvent(input$db, {

    ready_to_make = TRUE

      if(input$dbname == "Enter name here." || nchar(input$dbname) == 0){
        initial_message <<- paste0(initial_message, "\nThe database needs a name. Please give it one!")

        output$message <- renderText(initial_message)
        ready_to_make <- F

      }

      if(ready_to_make){

        name = paste0(gsub(" ", "_", input$dbname), ".db")

        if(file.exists(name)){
          mod_handle <<- recplot_py$reads_database(name)
          is_correct = mod_handle$check_valid()
          if(!is_correct){
            shinyalert("That file already exists, but doens't appear to be a RecruitPlotEasy database. I won't try to modify it.\nPlease select another file.")
            ready_to_make = F
            updateTextInput(session, "dbname", value = "Enter name here.")
            mod_handle <<- ""
          }
        }else{
          mod_handle <<- recplot_py$reads_database(name)
          mod_handle$initialize_db()
        }

        if(ready_to_make){
          plot_handle <<- recplot_py$recruitment_plot_data(name)

          updateTextInput(session, "dbname", value = name)

          initial_message <<- paste0(initial_message, "\nDatabase Built! Let's add some data! Use the reads, genomes, and genes buttons to add your data.")

          initial_message2 <<- paste0(initial_message2, "\n\nI have the database you just built loaded in.\nYou can go to view plots now, or check it out here.\n")

          output$message <- renderText(initial_message)
          output$message2 <- renderText(initial_message2)

          exist_db <<- name
          output$exist_dbname = renderText(name)
          output$db_stat = renderText(name)

        }
      }


  })

  observeEvent(input$exist_db,{

    tryCatch({
      db <- choose_file(caption = "Select an Existing Recruitment Plot Database")
    },
    error = function(cond){
      db <- "No existing database selected. Try again?"
      return(db)
    })

    if(length(db) == 0){
      db <- "No existing database selected. Try again?"
      return(db)
    }else{

      plot_handle <<- recplot_py$recruitment_plot_data(db)
      usable = plot_handle$check_valid()

        if(!usable){
          shinyalert("This doesn't look like a RecruitPlotEasy database to me!", "Did you make this file with the database creation tab?", type = "error")
          plot_handle <<- ""
          db <- "No existing database selected. Try again?"

        }else{
          mod_handle <<- recplot_py$reads_database(db)

          plot_handle$get_samples()

          samps = unlist(plot_handle$samples)

          samples_in_db <- unlist(samps)

          names(samples_in_db) = samps

          if(length(samps)==0){
            shinyalert("","I didn't find any samples in this database, but it looks like a RecruitPlotEasy database. It's possible that reads haven't yet been added to it. Reads must be added to the database before anything can be plotted.", type = "info")
          }else{
            updateSelectInput(session, "samples", choices = samples_in_db)
            updateSelectInput(session, "samples_interact", choices = samples_in_db)
          }

        }

    }

    exist_db <<- db
    output$exist_dbname = renderText(db)
    output$db_stat = renderText(db)
    #updateTextInput(session, "exist_dbname", value = db)


  })

  #Add genomes
  observeEvent(input$contigs, {

    tryCatch({
      contigs <- choose_file(caption = "Select Reference Genomes")
    },
    error = function(cond){
      contigs <- "No genomes selected. Try again?"
      updateTextInput(session, "contig_file", value = contigs)
    })

    if(length(contigs) == 0){
      contigs <- "No genomes selected. Try again?"
      updateTextInput(session, "contig_file", value = contigs)
    }else{
      if(!file.exists(contigs)){
        contigs <- "No genomes selected. Try again?"
      }

      fmt = recplot_py$detect_file_format(contigs)
      if(fmt == "fasta"){
        updateTextInput(session, "contig_file", value = contigs)
      }else{
        shinyalert("This doesn't look like a set of genomes. Try again?", type = "error")
        contigs <- "No genomes selected. Try again?"
        updateTextInput(session, "contig_file", value = contigs)
      }

    }

  })

  observeEvent(input$add_contigs, {
    if(mod_handle == ""){
      shinyalert("", "I can't add files unless you select or create a database first.", type = "error")
    }else{
    if(file.exists(input$contig_file)){
      mod_handle$set_genomes(input$contig_file, input$one_mag)
      initial_message <<- paste0(initial_message, "\nGenomes added!")

      output$message <- renderText(initial_message)

    }
    }
  })

  observeEvent(input$what_are_mags, {

    mag_info = "RecruitPlotEasy supports the visualization of Metagenome Assembled Genomes (MAGs) created through binning multiple assembled contigs.\n\nTo display multiple contigs as a single MAG, RecruitPlotEasy has to be aware of what contigs belong to what MAGs, which is done by adding MAGs to your database To add MAGs to a database, each MAG needs to have its contigs organized into separate FASTA format files. Add these files to a database using the 'add genomes' button with the MAG checkbox filled. You can add as many MAGs as you want, but each MAG has to be added one at a time."

    shinyalert("", mag_info, type = "info", imageWidth = 200)


  })

  #Add reads
  observeEvent(input$reads, {

    tryCatch({
      reads <- choose_file(caption = "Select Mapped Reads")
    },
    error = function(cond){
      reads <- "No reads selected. Try again?"
      updateTextInput(session, "read_file", value = reads)
    })

    if(length(reads) == 0){
      reads <- "No reads selected. Try again?"
      updateTextInput(session, "read_file", value = reads)
    }else{

      if(file.exists(reads)){
        fmt = recplot_py$detect_file_format(reads)

        if(fmt == "sam" || fmt == "bam" || fmt == "blast"){
          updateTextInput(session, "read_file", value = reads)
        }else{
          shinyalert("This doesn't look like mapped reads. Try again?", type = "error")
          reads <- "No reads selected. Try again?"
          updateTextInput(session, "read_file", value = reads)
        }
      }
    }



  })

  observeEvent(input$add_reads,{
    if(mod_handle == ""){
      shinyalert("", "I can't add files unless you select or create a database first.", type = "error")
    }else{

    if(file.exists(input$read_file)){
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Adding reads...", value = 0.33, detail = "Please be patient.")

      mod_handle$set_reads(input$read_file)

      initial_message <<- paste0(initial_message, "\nReads added!")

      output$message <- renderText(initial_message)

      progress$set(message = "Done!", value = 1)

      plot_handle$get_samples()

      samps = unlist(plot_handle$samples)

      samples_in_db <- unlist(samps)

      names(samples_in_db) = samps

      if(length(samps)==0){
        shinyalert("","I didn't find any samples in this database, but it looks like a RecruitPlotEasy database. It's possible that reads haven't yet been added to it. Reads must be added to the database before anything can be plotted.", type = "info")
      }else{
        updateSelectInput(session, "samples", choices = samples_in_db)
        updateSelectInput(session, "samples_interact", choices = samples_in_db)
      }

    }
    }
  })

  observeEvent(input$read_info, {
    shinyalert("", "RecruitPlotEasy allows you to add reads formatted in SAM, BAM, or tabular BLAST (generated by mapping reads using BLASTn with the '--outfmt 6' option).\n\nYou can plot after just adding reads, but you can't visualize genes and the plot may not cover the entire genome if you do.", type = "info")
  })

  #Add genes
  observeEvent(input$genes,{


    tryCatch({
      genes <- choose_file(caption = "Select Prodigal GFF genes")
    },
    error = function(cond){
      genes <- "No genes selected. Try again?"
      updateTextInput(session, "gene_file", value = genes)
    })

    if(length(genes) == 0){
      genes <- "No genes selected. Try again?"
      updateTextInput(session, "gene_file", value = genes)
    }else{

      if(file.exists(genes)){

        fmt = recplot_py$detect_file_format(genes)

        if(fmt == "genes"){
          updateTextInput(session, "gene_file", value = genes)
        }else{
          shinyalert("This doesn't look Prodigal GFF format genes. Try again?", type = "error")
          genes <- "No genes selected. Try again?"
          updateTextInput(session, "gene_file", value = genes)
        }
      }
    }


  })

  observeEvent(input$add_genes, {

    if(mod_handle == ""){
      shinyalert("", "I can't add files unless you select or create a database first.", type = "error")
    }else{
      if(file.exists(input$gene_file)){
        mod_handle$add_genes(input$gene_file)

        initial_message <<- paste0(initial_message, "\nGenes added!")
        output$message <- renderText(initial_message)

      }

    }


  })

  observeEvent(input$gene_info, {
    shinyalert("", "RecruitPlotEasy normally divides a genome into fixed-size windows spaced evenly along your genome. If you add genes for your genomes and select genes as the plotting task, then RecruitPlotEasy will instead use the start and stop loci of those genes as its windows, filling in the intergenic gaps as needed. The interactive plot will also list gene information with its annotations.\n\nGenes must be created by Prodigal and must be in the GFF format ('-f gff' in Prodigal's options).", type = "info")
  })

  #db management

  observeEvent(input$show_samples, {

  if(plot_handle != ""){
      plot_handle$get_samples()
      samples = unlist(plot_handle$samples)

      if(length(samples) > 0){
        text_out = paste(unlist(samples), collapse="\n")
      }else{
        text_out = "No samples found."
      }
      output$message2 <- renderText(text_out)

    }
  })

  observeEvent(input$current_args, {
    if(plot_handle != ""){
      #bin_width = 1000, bin_height = 0.5, in_group_minimum = 95, min_align_len = 50,
      #min_pct_align = 90, local_pct_id = True, only_one_alignment_per_read = True,
      #best_hit_criteria = "local_percent_ID", tad_percent = 80

      local_id = "No"
      if(plot_handle$pct_id_is_local){
        local_id = "Yes"
      }
      ma = "Yes"
      if(plot_handle$only_one_alignment_per_read){
        ma = "No"
      }

      writeout = paste0("Bin width: ", plot_handle$width,
                        "\nBin height: ", plot_handle$height,
                        "\n%ID minimum: ", plot_handle$in_grp,
                        "\nMinimum alignment length: ", plot_handle$min_align_length,
                        "\nMinimum percent alignment: ", plot_handle$min_pct_align,
                        "\n%ID is determined locally: ", local_id,
                        "\nMultiple alignments allowed: ", ma,
                        "\nBest hit criteria: ", plot_handle$best_hit_criteria,
                        "\nTAD-N: ", plot_handle$tad_middle)

      output$message2 <- renderText(writeout)

    }
  })


  observeEvent(input$local_id, {
    if(plot_handle != ""){
      plot_handle$pct_id_is_local <<- input$local_id
      }
  })

  observeEvent(input$multiple_aligns, {
    if(plot_handle != ""){

      plot_handle$only_one_alignment_per_read <<- input$multiple_aligns
    }
  })

  observeEvent(input$bh_crit, {
    if(plot_handle != ""){

      plot_handle$best_hit_criteria <<- input$bh_crit

    }
  })

  observeEvent(input$min_align_len, {
    if(plot_handle != ""){
      plot_handle$min_align_length <<- input$min_align_len
    }
  })

  observeEvent(input$min_pct_align, {
    if(plot_handle != ""){
      plot_handle$min_pct_align <<- input$min_pct_align
    }
  })

  observeEvent(input$tadn, {
    if(plot_handle != ""){
      plot_handle$tad_middle <<- input$tadn
    }
  })

  observeEvent(input$plot_genes, {
    if(plot_handle != ""){
      plot_handle$plotting_genes <<- input$plot_genes
    }
  })

  observeEvent(input$export_reads,{
    if(plot_handle != ""){
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Exporting your reads...", value = 0.5, detail = "Please be patient.")
      plot_handle$export_reads_to_file()
      progress$set(message = "Exporting your reads...", value = 1, detail = "Please be patient.")
    }
  })

  #Plotting args

  #Data output and plot saving.
  observeEvent(input$print_stat, {


    if(is.na(stat_plot[1])){

      shinyalert("Cannot print plot.", "There is no genome currently loaded. You need to hit 'View Selected Genome' first.", type = "error")

      return(NA)
    }else{

      if(input$pdf_name == ""){

        shinyalert(
          "The output plot needs a name.", "You can enter a name here, then hit Save PDF button again to print it.", type = "input",
          callbackR = function(x) {
            updateTextInput(session, "pdf_name", value = x)
          }
        )

        return(NA)
      }else{

        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Saving Recruitment Plot", value = 0.33, detail = "Please be patient.")

        print_name <- gsub(" ", "_", input$pdf_name)

        pdf(paste0(print_name, ".pdf"), height = 11, width = 17)

        suppressWarnings(print(stat_plot))

        dev.off()

        fwrite(main_dat, paste0(print_name, "_recruitment_data.tsv"), sep = "\t")
        fwrite(depth_dat, paste0(print_name, "_sequencing_depth.tsv"), sep = "\t")
        fwrite(hist_dat, paste0(print_name, "_bp_histogram.tsv"), sep = "\t")
        fwrite(data.table(sql_query = query), paste0(print_name, "_query.tsv"), sep = "\t")

        progress$set(message = "Saving Recruitment Plot", value = 1, detail = "Please be patient.")
      }

    }



  })

  observeEvent(input$print_interact, {

    if(is.na(interact_plot[1])){

      shinyalert("Cannot print plot.", "There is no genome currently loaded. You need to hit 'View Selected Genome' first.", type = "error")

      return(NA)

    }else{

      if(input$pdf_name_interact == ""){

        shinyalert(
          "The output interactive plot needs a name.", "You can enter a name here, then hit the Save Plot button again to save it.", type = "input",
          callbackR = function(x) {
            updateTextInput(session, "pdf_name_interact", value = x)
          }
        )

        return(NA)

      }else{


        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Saving Recruitment Plot", value = 0.33, detail = "Please be patient. Interactive plots take a bit.")

        print_name <- gsub(" ", "_", input$pdf_name_interact)

        saveWidget(interact_plot, paste0(print_name, ".html"))

        fwrite(main_dat, paste0(print_name, "_recruitment_data.tsv"), sep = "\t")
        fwrite(depth_dat, paste0(print_name, "_sequencing_depth.tsv"), sep = "\t")
        fwrite(hist_dat, paste0(print_name, "_bp_histogram.tsv"), sep = "\t")
        fwrite(data.table(sql_query = query), paste0(print_name, "_query.tsv"), sep = "\t")

        progress$set(message = "Saving Recruitment Plot", value = 1, detail = "Please be patient.")

      }

    }
    return(NA)
  })

  #Static plots

  #Select genomes
  observeEvent(input$samples, {
    if(plot_handle != ""){

    plot_handle$reset()

    plot_handle$set_sample(input$samples)
    plot_handle$link_to_file()

    if(!is.na(stat_plot) || !is.na(interact_plot)){
      shinyalert("", "You've switched samples. RecruitPlotEasy is no longer keeping track of your current plot and you cannot save it now.", type = "info")
    }

    stat_plot <<- NA
    interact_plot <<- NA
    main_dat <<- NA
    depth_dat <<- NA
    hist_dat <<- NA
    query <<- NA

    mags = names(plot_handle$mags_in_sample)
    mags_in_samp = mags

    names(mags_in_samp) = mags

    updateSelectInput(session, "mags_in_db", choices = mags_in_samp)
    updateSelectInput(session, "mags_in_db_interact", choices = mags_in_samp)
    }
  })

  observeEvent(input$samples_interact, {
    if(plot_handle != ""){
    plot_handle$reset()

    plot_handle$set_sample(input$samples_interact)
    plot_handle$link_to_file()

    if(!is.na(stat_plot) || !is.na(interact_plot)){
      shinyalert("", "You've switched samples. RecruitPlotEasy is no longer keeping track of your current plot and you cannot save it now.", type = "info")
    }

    stat_plot <<- NA
    interact_plot <<- NA
    main_dat <<- NA
    depth_dat <<- NA
    hist_dat <<- NA
    query <<- NA

    mags = names(plot_handle$mags_in_sample)
    mags_in_samp = mags

    names(mags_in_samp) = mags

    updateSelectInput(session, "mags_in_db", choices = mags_in_samp)
    updateSelectInput(session, "mags_in_db_interact", choices = mags_in_samp)
    }
  })

  observeEvent(input$mags_in_db, {
    if(plot_handle != ""){
    plot_handle$reset()
    plot_handle$set_mag(input$mags_in_db)
    }
  })

  observeEvent(input$mags_in_db_interact, {
    if(plot_handle != ""){
    plot_handle$reset()
    plot_handle$set_mag(input$mags_in_db_interact)
    }
  })
  observeEvent(input$width, {
    if(plot_handle != ""){
      plot_handle$width = input$width
    }
  })

  observeEvent(input$width_interact, {
    if(plot_handle != ""){
      plot_handle$width <<- input$width_interact
    }
  })

  #Height
  observeEvent(input$height, {
    if(plot_handle != ""){
      plot_handle$height <<- input$height
    }
  })

  observeEvent(input$height_interact, {
    if(plot_handle != ""){
      plot_handle$height <<- input$height_interact
    }
  })

  #blue box
  observeEvent(input$in_group_min_stat, {
    if(plot_handle != ""){
      plot_handle$in_grp <<- input$in_group_min_stat
    }
  })

  observeEvent(input$in_group_min_interact, {
    if(plot_handle != ""){
      plot_handle$in_grp <<- input$in_group_min_interact
    }
  })

  #static extended opts

  observeEvent(input$linear_stat, {
    if(plot_handle != ""){
      plot_handle$plot_linear <<- input$linear_stat
    }
  })

  # observeEvent(input$show_peaks, {
  #   if(plot_handle != ""){
  #     plot_handle$show_peaks <<- input$show_peaks
  #   }
  # })

  #data load

  observeEvent(input$get_a_mag, {
    if(plot_handle != ""){

      if(input$plot_genes){
        plot_handle$load_genes()
      }

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Creating Recruitment Plot", value = 0.33, detail = "Loading data. Please be patient.")

      plot_handle$create_frame()
      plot_handle$create_minimal_matrix()
      plot_handle$calculate_marginals()

      progress$set(message = "Creating Recruitment Plot", value = 0.66, detail = "Creating plot. Please be patient.")

      output$read_recruitment_plot <- renderPlot(suppressWarnings({

        if(plot_handle == ""){
          wp = warning_plot()
          return(wp)
        }

        if(length(unlist(plot_handle$samples)) == 0){
          wp = warning_plot()
          return(wp)
        }

        if(length(unlist(plot_handle$current_mag)) == 0){
          wp = warning_plot()
          return(wp)
        }

        if(length(unlist(plot_handle$contigs_in_mag)) == 0){
          wp = warning_plot()
          return(wp)
        }

        if(length(unlist(plot_handle$describe_x)) == 0){
          wp = warning_plot()
          return(wp)
        }

        static_plot_and_data = static_plot(plot_handle)

        stat_plot <<- static_plot_and_data[[1]]
        #interact_plot = NA
        main_dat <<- static_plot_and_data[[2]]
        depth_dat <<- static_plot_and_data[[3]]
        hist_dat <<- static_plot_and_data[[4]]
        query <<- static_plot_and_data[[5]]

        progress$set(message = "Creating Recruitment Plot", value = 1, detail = "Please be patient")

        return(stat_plot)

      }))
    }
  })

  #Create a reactive function to update the reactive list every time the plotly_relayout changes

  observeEvent(input$get_a_mag_interact, {
    if(plot_handle != ""){

      if(input$plot_genes){
        plot_handle$load_genes()
      }

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Creating Recruitment Plot", value = 0.33, detail = "Loading data. Please be patient.")

      plot_handle$create_frame()
      plot_handle$create_minimal_matrix()
      plot_handle$calculate_marginals()

      progress$set(message = "Creating interactive Recruitment Plot", value = 0.66, detail = "Please be patient. The interactive plots take longer.")

      output$Plotly_interactive <- renderPlotly(suppressWarnings({

        #Reset this each time to make checking for it more consistent.
        if(plot_handle == ""){
          wp = warning_plot()
          wp = ggplotly(wp)
          return(wp)
        }

        if(length(unlist(plot_handle$samples)) == 0){
          wp = warning_plot()
          wp = ggplotly(wp)
          return(wp)
        }

        if(length(unlist(plot_handle$current_mag)) == 0){
          wp = warning_plot()
          wp = ggplotly(wp)
          return(wp)
        }

        if(length(unlist(plot_handle$contigs_in_mag)) == 0){
          wp = warning_plot()
          wp = ggplotly(wp)
          return(wp)
        }

        if(length(unlist(plot_handle$describe_x)) == 0){
          wp = warning_plot()
          wp = ggplotly(wp)
          return(wp)
        }

        interactive_plot_and_data = interactive_plot(plot_handle)

        #stat_plot = NA
        interact_plot <<- interactive_plot_and_data[[1]]
        main_dat <<- interactive_plot_and_data[[2]]
        depth_dat <<- interactive_plot_and_data[[3]]
        hist_dat <<- interactive_plot_and_data[[4]]
        query <<- interactive_plot_and_data[[5]]

        return(interact_plot)

      }))

    }
  })

  observeEvent(input$addtl_stats_static, {
    if(plot_handle != ""){
      if(!is.null(plot_handle$describe_y)){
      stats = additional_stats(plot_handle)
      shinyalert("", stats, type = "info")
      }else{
        shinyalert("", "You need to load a plot first!", type = 'error')
      }
    }
  })

  observeEvent(input$addtl_stats_interact, {
    if(plot_handle != ""){
      if(!is.null(plot_handle$describe_y)){
      stats = additional_stats(plot_handle)
      shinyalert("", stats, type = "info")
    }else{
      shinyalert("", "You need to load a plot first!", type = 'error')
    }
    }
  })

  #Load reads

  observeEvent(input$add_cart_stat, {
    if(plot_handle != ""){
    if(!is.null(plot_handle$describe_y)){
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Adding reads to your cart", value = 0.33, detail = "Please be patient.")

      plot_handle$select_reads_for_export()

      progress$set(message = "Adding reads to your cart", value = 1 , detail = "Please be patient.")
  }else{
    shinyalert("", "You need to load a plot first!", type = 'error')
  }

    }

  })

  observeEvent(input$add_cart_interact, {
    if(plot_handle != ""){
      if(!is.null(plot_handle$describe_y)){
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Adding reads to your cart", value = 0.33, detail = "Please be patient.")

      plot_handle$select_reads_for_export()

      progress$set(message = "Adding reads to your cart", value = 1 , detail = "Please be patient.")
    }else{
      shinyalert("", "You need to load a plot first!", type = 'error')
    }
    }

  })

  #Consistency between plot panels

  observeEvent(input$tabs, {

    if(input$tabs == "Recruitment Plot"){
      if(!is.null(input$samples_interact)){
        updateSelectInput(session, "samples", selected = input$samples_interact)
      }
      if(!is.null(input$mags_in_db_interact)){
        updateSelectInput(session, "mags_in_db", selected = input$mags_in_db_interact)
      }
      updateNumericInput(session, "width", value = input$width_interact)
      updateNumericInput(session, "height", value = input$height_interact)
      #updateNumericInput(session, "low_bound", value = input$low_bound_interact)
      updateNumericInput(session, "in_group_min_stat", value = input$in_group_min_interact)
      updateSelectInput(session, "regions_stat", selected = input$regions_interact)
      updateTextInput(session, "pdf_name", value = input$pdf_name_interact)
    }
    if(input$tabs == "Interactive Plot"){
      if(!is.null(input$samples)){
        updateSelectInput(session, "samples_interact", selected = input$samples)
      }
      if(!is.null(input$mags_in_db)){
        updateSelectInput(session, "mags_in_db_interact", selected = input$mags_in_db)
      }
      updateNumericInput(session, "width_interact", value = input$width)
      updateNumericInput(session, "height_interact", value = input$height)
      #updateNumericInput(session, "low_bound_interact", value = input$low_bound)
      updateNumericInput(session, "in_group_min_interact", value = input$in_group_min_stat)
      updateSelectInput(session, "regions_interact", selected = input$regions_stat)
      updateTextInput(session, "pdf_name_interact", value = input$pdf_name)
    }

    if(input$tabs == "Database Management"){
      if(plot_handle == ""){
        shinyalert("", "Nothing went wrong, but nothing on this page works unless you've already selected a RecruitPlotEasy database. Go back to the database creation tab and make a database (it will automatically be selected if you do) or select an existing one, please.", type = 'info')
      }
    }



  })

  session$onSessionEnded(function() {

    cat("\nThank you for using Recruitment Plots!\n")
    cat("Any databases you created or plots you made are stored in: ")
    cat(getwd())
    cat("\n")

    stopApp()
  })

}


#Dev
#runApp(list(ui = recplot_UI(), server = recplot_server), launch.browser = T)

#' @export
RecruitPlotEasy <- function(){

  #Import python code.
  initiate()
  #Just make sure it's available.
  recplot_py <- get_python()

  cat("Recruitment plots ready.\n")

  if (interactive()) {
    runApp(list(ui = recplot_UI(), server = recplot_server), launch.browser = T)
  }

}





