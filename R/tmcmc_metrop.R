#' @title Simulate a Transformation based Markov Chain Monte Carlo (TMCMC) algorithm (Rcpp sped up version)
#'
#' @param target_pdf The log target density function from which the user wants to generate samples.
#' @param scale The proposal density scaling parameter. An approximation of the optimal scaling given the target_pdf is performed by OptimalScaling().
#'              The default scale is this estimated optimal scaling
#' @param base The starting value of the chain
#' @param nsamples The number of samples to be drawn using the TMCMC algorithm.
#' @param burn_in The number of samples assigned as burn-in period. The default burn-in is taken to be one-third of nsamples.
#'
#'
#' @description The function simulates a TMCMC chain of length nsamples using the scale, base and burn in taken optimally as default or specified by user. The function
#' outputs the full chain as well as the estimated posterior mean estimated from the samples drawn post burn-in.
#'
#
#'  @author  Kushal K Dey
#'
#'  @export


tmcmc_metrop <- function(target_pdf, scale, base, nsamples, burn_in=NULL)
{
  #Rcpp::sourceCpp('src/RcppExports.cpp')
  Rcpp::sourceCpp('src/utils.cpp')
  if(is.null(burn_in)) burn_in <- nsamples/3;
  chain <- matrix(0, nsamples, length(base))
  num=1;
  chain[1,] <- base;
  while(num <= nsamples) {
    eps <- rnorm(1,0,scale);
    b <- sample(c(-1,+1),length(base),replace=TRUE)
    chain[(num+1),] <- tmcmcUpdate(chain[num,],b,eps,target_pdf)$chain
    if(num %% 500 ==0)
      paste("The chain is at iteration:",num);
  }
  posterior_mean <- apply(chain[round(burn_in):nsamples,],2,mean);
  ll <- list("chain"=chain,"post.mean"=posterior_mean);
  return(ll)
}

