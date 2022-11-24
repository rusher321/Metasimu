#' Sim_diri
#'
#' @param pi, vector: the parameter of dirichlet distirbution
#' @param n, numeric , sample number
#' @param disp , numeric , shape of the simulation
#' @param beta , numeric , generate the causal taxa
#' @param sigN , numeric, number of causal taxa
#' @param have.bias , logistic, if add a bias factor to simualted data
#' @param disp2 , numeric, the other parameter of beta distribution
#'
#' @return matrix or list
#' @import dirmult
#' @export
#' @examples
Sim_diri <- function(pi = NULL, n, disp = 0.02, beta = 3, sigN = NULL, have.bias = FALSE, disp2 = 0.2){

  outlist <- list()

  alpha_sim <- (1-disp)/disp
  taxN <- length(pi)

  if(is.null(pi)){
      pi <- rbeta(n, shape1 = disp, shape2 = disp2)
      pi <- pi/sum(pi)
      otu.table.sim <- dirmult::rdirichlet( n = n, alpha = alpha_sim*pi)
  }else{
      otu.table.sim <- dirmult::rdirichlet( n = n, alpha = alpha_sim*pi)
  }

  if(is.null(sigN)){
    beta.otu.log <- log(rep(beta, n))
    Y <- c(rep(0, round(n/2), rep(1, n-round(n/2))))
    sigtax <- sample(1:length(pi), sigN)
    otu.table.sim[, sigtax] <- otu.table.sim[, sigtax] * exp(Y %*% t(beta.otu.log))
    otulist[["Y"]] <- Y
    otulist[["sigtax"]] <- sigtax
  }

  if (have.bias == TRUE){
    set.seed(0)
    bias.factor.log <- rep(NA, length(pi))
    bias.factor.log[sigtax] <- rnorm(length(sigtax), 1, 0.8)
    bias.factor.log[-sigtax] <- rnorm(length(pi)-length(sigtax), 0, 0.8)
    bias.factor <- exp(bias.factor.log)
    otu.table.sim <- t(t(otu.table.sim) * bias.factor)
  }

  otu.table.sim <- otu.table.sim/rowSums(otu.table.sim)
  otulist[["simdata"]] <- otu.table.sim

  return(otulist)

}
