
#' @title simulate_networks()
#'
#' @description
#'  Implements an algorithm that generates symmetric positive definite precision matrices, where a fraction of entries are (linearly)
#'  dependent on the covariate ordinal value that corresponds to each matrix. The function also generates normally distributed data
#'  using the inverses of the precision matrices it generates, which can be used in simulations.
#'
#' @param P The number of variables in the networks. i.e. PxP precision matrices are generated.
#' @param A The number of precision matrices (/adjacency matrices) that you want to generate.
#' @param network_seed Sets the seed for the random generation of the initial, scale-free precision matrix.
#'  This seed will recursively increase by 1 until positive definite solutions are found. The final seed is stored
#'  in output and printed to help speed up future function usages.
#' @param data_seed Sets the seed for data generation. To generate several sets of data from the same precision matrices, keep `network_seed`
#'  constant and change `data_seed`.
#' @param frac_change The fraction of non-zero entries of the first precision matrix that will have negative correlation with the covariate
#'  the number of which is equal to the number that will appear. I.e. if the network has 100 edges and `frac_change = 0.2`, 20 edges will have
#'  negative correlation and 20 will have positive correlation.
#' @param Ns_sample For data generation, the set from which the number of samples taken using each precision matrix will be randomly drawn from.
#'
#'
#' @return
#' - $Omegas : The set of precision matrices.
#' - $Ys : List of simulated data matrices.
#' - $Ns : The number of drawn samples for each network.
#' - $neg_id : The indices of disappearing edges.
#' - $pos_id : The indices of appearing edges.
#' - $As : The adjacency matrices corresponding to the precision matrices.
#' @export
#'
#' @examples simulate_networks()
#' @examples simulate_networks(P = 100, A = 4, frac_change = 0.4)

simulate_networks <- function(P = 50, A = 3, network_seed = 123, data_seed = 123,
                               frac_change = 0.2, Ns_sample = c(150,150,150)){



  set.seed(network_seed)

  # Sample sizes for each dataset
  Ns <- sample(Ns_sample, A, replace = T)

  #generate an initial network
  net0 <- huge::huge.generator(n = Ns[1], d = P, graph = 'scale-free', v = 0.4, u = 0.05, verbose = F)
  diag.omega <- diag(net0$omega)
  net0$omega[as.numeric(net0$theta) != 1] <- 0
  diag(net0$omega) <- diag.omega

  pol <- matrix(data = sample(c(-1,1), P^2, replace = T), nrow = P)
  pol <- upper.tri(pol)*pol + t(upper.tri(pol))*t(pol)
  diag(pol) <- 1
  net0$omega <- net0$omega * pol

  # Replicating the initial precision matrix.
  Omegas <- replicate(A, net0$omega, simplify=FALSE)

  # Assert which vertices have an edge (edge_id) between them and which don't (non_edge_id)
  edge_id <- which(as.logical(net0$theta)  & upper.tri(net0$omega), arr.ind = T)
  non_edge_id <- which(!as.logical(net0$theta) & upper.tri(net0$omega), arr.ind = T)


  npos <- nneg <- ceiling(nrow(edge_id)*frac_change)


  # Of the absent edges, randomly select 'npos' of them to have a positive correlation with the covariate
  # Of the absent edges, randomly select 'nneg' of them to have a negative correlation with the covariate
  neg_id <- matrix(edge_id[sample(1:nrow(edge_id),nneg),], ncol=2)
  pos_id <- matrix(non_edge_id[sample(1:nrow(non_edge_id),npos),], ncol = 2)




  # Iterate over each disappearing edge and have them 'fade' out linearly
  for (i in 1:nneg) {
    val <- net0$omega[neg_id[i,1],neg_id[i,2]]
    #k <- sample(2:A, 1, replace = T)
    for (t in 2:A) {
      Omegas[[t]][neg_id[i,1],neg_id[i,2]] <- Omegas[[t]][neg_id[i,2],neg_id[i,1]] <- ((A-t)/(A-1))*val
    }
  }

  # Iterate over each appearing edge and have them 'fade' in linearly
  hi_val <- max(abs(net0$omega[lower.tri(net0$omega)]))
  for (i in 1:npos) {
    #k <- sample(2:A, 1, replace = T)
    val <- hi_val - sample(c(0,0.02,0.04,0.06,0.08,0.1), 1, replace = T)
    plus_min <- sample(c(-1,1),1,replace = T)
    for (t in 2:A) {
      Omegas[[t]][pos_id[i,1],pos_id[i,2]] <- Omegas[[t]][pos_id[i,2],pos_id[i,1]] <- plus_min*val*(t-1)/(A-1)
    }
  }


  Ys = list()

  for(i in seq_along(Omegas)){
    Omegas[[i]] <- (Omegas[[i]] + t(Omegas[[i]]))/2
  }
  for(i in seq_along(Omegas)){
    while(!matrixcalc::is.positive.definite(Omegas[[i]])){
      nets <- simulate_networks(P = P, A = A, network_seed = network_seed + 1, data_seed = data_seed,
                                 frac_change = frac_change)
      return(nets)
    }
  }
  set.seed(data_seed)

  for(y in 1:A){

    Ys[[y]] <- MASS::mvrnorm(n=Ns[y], mu = rep(0,P), Sigma = solve(Omegas[[y]]))
    Ys[[y]] <- scale(Ys[[y]])

  }
  #extract adjacencies

  As <- lapply(Omegas, function(omega) (abs(omega) > 0) * 1)

  cat("Solution found with network_seed = ",network_seed,"\n")

  return(list(Omegas = Omegas, Ys = Ys, Ns = Ns, neg_id = neg_id, As = As,
              net0_sig = net0$sigma, pos_id = pos_id, network_seed = network_seed, edge_id = edge_id))

}


