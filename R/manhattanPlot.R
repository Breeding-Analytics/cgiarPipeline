### FUNCTION FOR INTERACTIVE MANHATTAN PLOT
manhattanPlot <- function(df = NULL, threshold = NULL) { # , field = NULL
  
  "%>%" <- dplyr::`%>%`
  y_min <- min(df$predictedValue)
  if(threshold < max(df$predictedValue)){y_max <- round(max(df$predictedValue)+1)}else{y_max <- round(threshold+1)} 
  # (max(df$predictedValue)+1)%/%1+1 # (threshold+1)%/%1+1
  axisdf <- df %>%
    dplyr::group_by(stdError) %>%
    dplyr::summarize(center=(max(BPcum) + min(BPcum))/2)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = BPcum, y = predictedValue)) + #, text = text
    ggplot2::geom_point(ggplot2::aes(color=as.factor(stdError)), alpha = 0.8, size = 1.5) +
    # ggplot2::scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    ggplot2::scale_color_manual(values = rep(c("grey", "skyblue"), .5*length(unique(df$stdError)))) +
    ggplot2::scale_x_continuous(label = axisdf$stdError, breaks= axisdf$center) +
    # ggplot2::ylim(y_min, y_max) +
    ggplot2::expand_limits(y=c(0,y_max)) +
    ggplot2::scale_y_continuous(breaks = seq(0,y_max)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dotted") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      title = ggplot2::element_text(colour = "black", size = 10),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    ) + ggpubr::border() # + ggplot2::labs(x = "Chromosome") # , y = expression(-log[10](italic("p"))) title = field, 
  p <- plotly::plotly_build(p)
  p <- p %>% plotly::layout(xaxis = list(showline = TRUE, mirror = TRUE, linecolor = 'black'),
                            yaxis = list(showline = TRUE, mirror = TRUE, linecolor = 'black'))
  # plotly::ggplotly(p)
}

