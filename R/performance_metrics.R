#' Internal Helper Function
#'
#' @noRd
#' @keywords internal
evaluate_network <- function(true_network, estimated_network,t=1e-9) {
  if (!all(dim(true_network) == dim(estimated_network))) {
    stop("The dimensions of the true network and the estimated network must match.")
  }



  LT <- lower.tri(true_network)
  TP <- sum(true_network[LT] >t & estimated_network[LT] > t)
  TN <- sum(true_network[LT] <t & estimated_network[LT] < t)
  FP <- sum(true_network[LT] <t & estimated_network[LT] > t)
  FN <- sum(true_network[LT] >t & estimated_network[LT] < t)

  return(list(
    true_positives = TP,
    true_negatives = TN,
    false_positives = FP,
    false_negatives = FN,
    conf.mat = matrix(c(TP, FN, FP, TN), nrow = 2, byrow = T)
  ))
}

## 'evaluate_network_list' takes multiple true and estimated networks and returns the confusion matrix components
## It also returns the confusion matrix in itself, and 'P' which is the dimension of the networks that it is calculating.

#' @title evaluate_network_list
#'
#' @description
#' Calculates the overall number of true positives, true negatives, false positives and false negatives attained
#'  by the network estimations and compiles into confusion matrix form
#'
#'
#' @param true_networks List of matrices containing the true (simulated) networks.
#' @param estimated_networks List of matrices containing the estimated networks.
#'
#' @return     $true_positives, $true_negatives, $false_positives, $false_negatives, $conf.mat
#' @export
#'
#' @examples .
evaluate_network_list <- function(true_networks, estimated_networks) {
  # Ensure both lists have the same length
  if (length(true_networks) != length(estimated_networks)) {
    stop("The lists of true networks and estimated networks must have the same length.")
  }

  total_TP <- 0
  total_TN <- 0
  total_FP <- 0
  total_FN <- 0

  conf.mat <- matrix(0, nrow = 2, ncol = 2)

  P <- ncol(true_networks[[1]])


  for (i in seq_along(true_networks)) {
    # Call the original evaluate_network function
    result <- evaluate_network(true_networks[[i]], estimated_networks[[i]] )

    total_TP <- total_TP + result$true_positives
    total_TN <- total_TN + result$true_negatives
    total_FP <- total_FP + result$false_positives
    total_FN <- total_FN + result$false_negatives
    conf.mat <- conf.mat + result$conf.mat
  }


  return(list(
    P = P,
    true_positives = total_TP,
    true_negatives = total_TN,
    false_positives = total_FP,
    false_negatives = total_FN,
    conf.mat = conf.mat)
  )
}

#' @title recall()
#'
#' @description Makes a calculation for precision based on the confusion matrix entries (TP/(TP+FN))
#'
#' @param conf.mat A 2x2 confusion matrix
#'
#' @return The recall corresponding to the confusion matrix.
#' @export
#'
#' @examples .
recall = function (conf.mat) {

  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else if (conf.mat[1, 1] == 0 & conf.mat[1, 2] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[1, 2]))
  }
}

#' @title precision()
#'
#' @description Makes a calculation for precision based on the confusion matrix entries (TP/(TP+FP))
#'
#' @param conf.mat A 2x2 confusion matrix
#'
#' @return The precision corresponding to the confusion matrix.
#' @export
#'
#' @examples .
precision = function (conf.mat) {

  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[2, 1]))
  }
}

#' @title sparsity()
#'
#' @description
#' Calculates the sparsity of a matrix
#'
#'
#' @param network A precision/adjacency matrix
#'
#' @return Sparsity calculation of the inputted matrix
#' @export
#'
#' @examples
#' m <- matrix(sample(0:1, 25, replace=TRUE), nrow=5)
#' sparsity(m)
sparsity <- function(network){
  return(sum(network[lower.tri(network)] != 0)/sum(lower.tri(network)))
}
