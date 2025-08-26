# Functions to create network plots

#' @title create_layout()
#'
#' @description Creates a network layout using the Fruchterman–Reingold force-directed layout algorithm, resulting in a balanced, spaced layout.
#'
#' @param network Single matrix input of the network for which the layout will be optimised. Function is designed for a network, NOT a precision matrix.
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



#' @title create_network_plots
#'
#' @param net_list list of matrices of networks to plot
#' @param layout_coords Matrix array of coordinates which should have been generated using the function `create_layout`.
#' @param color Set colour of the nodes.
#' @param edge_color Set colour of the edges
#' @param title Set the title of the plot.
#'
#' @return Plots a 1xN array of networks, each with the same layout.
#' @export
#'
#' @examples
#' m <- matrix(sample(0:1, 25, replace=TRUE), nrow=5)
#' m2 <- matrix(sample(0:1, 25, replace=TRUE), nrow=5)
#' m3 <- matrix(sample(0:1, 25, replace=TRUE), nrow=5)
#' m[lower.tri(m)] <- t(m)[lower.tri(m)]
#' diag(m) <- 0
#' layout <- create_layout(m)
#' create_network_plots(list(m,m2,m3), layout_coords = layout)
create_network_plots <- function(net_list, layout_coords, color = "red",edge_color = "purple",
                                 title = "Network Plots") {
  N <- length(net_list)
  networks <- vector("list", N)
  plots <- vector("list", N)

  for (i in seq_len(N)) {
    networks[[i]] <- network::network(net_list[[i]] > 0.5, directed = FALSE)
    plots[[i]] <- GGally::ggnet2(
      networks[[i]], edge.size = 0.3, alpha = 0.8, mode = layout_coords,
      color = color, label = F, size = 1, size.min = 0, edge.color = edge_color, title = "test"
    )
  }

  # Print all plots using gridExtra
  print(gridExtra::grid.arrange(grobs = plots, ncol = N,     top = grid::textGrob(
    label = title,
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )))
}
