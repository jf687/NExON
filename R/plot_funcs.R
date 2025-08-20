# Functions to create network plots

#' @title create_layout()
#'
#' @description Creates a network layout using the Fruchterman–Reingold force-directed layout algorithm, resulting in a balanced, spaced layout.
#'
#' @param network The network w
#'
#' @return Matrix array that gives the coordinates of each node. This should then be parsed as the layout_coords argument in create_network_plots().
#' @export
#'
#' @examples
#' m <- matrix(sample(0:1, 25, replace=TRUE), nrow=5)
#' m[lower.tri(m)] <- t(m)[lower.tri(m)]
#' diag(m) <- 0
#' create_layout(m)
create_layout <- function(network){
  layout_coords <- sna::gplot.layout.fruchtermanreingold(
    network::as.matrix.network.adjacency(
      network::network(network > 0.5, directed = FALSE)
    ),
    layout.par = list()
  )
  return(layout_coords)
}



#' Title
#'
#' @param data_list
#' @param layout_coords
#' @param color
#' @param edge_color
#' @param title
#'
#' @return
#' @export
#'
#' @examples
create_network_plots <- function(data_list, layout_coords, color = "red",edge_color = "purple",
                                 title = "Network Plots") {
  N <- length(data_list)
  networks <- vector("list", N)
  plots <- vector("list", N)


  # Loop to generate networks and their corresponding plots
  for (i in seq_len(N)) {
    networks[[i]] <- network::network(data_list[[i]] > 0.5, directed = FALSE)
    plots[[i]] <- GGally::ggnet2(
      networks[[i]], edge.size = 0.3, alpha = 0.8, mode = layout_coords,
      color = color, label = F, size = 1, size.min = 0, edge.color = edge_color, title = "test"
    )
    #+
    #ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 6))

  }

  # Print all plots using gridExtra
  print(gridExtra::grid.arrange(grobs = plots, ncol = N,     top = grid::textGrob(
    label = title,
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )))
}
