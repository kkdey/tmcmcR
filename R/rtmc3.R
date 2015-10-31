#' @title Simulate a Randomized Transformation based MC3 algorithm (Rcpp sped up version)
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
#'  @author  Kushal K Dey, Sourabh Bhattacharya
#'
#'  @export
#'



rtmc3 <- function(target_pdf, beta_set, scale, base, nsamples, burn_in=NULL)
{
  if(is.null(burn_in)) burn_in <- nsamples/3;
  rtmc3_chains <- vector("list", length(beta_set));
  num = 1
  while(num <= nsamples){
    rtmc3_chains <- mclapply(1:length(beta_set),
                  function(k){
                                if(num==1)
                                  chain <- base;
                                  out <- chain;
                                if(num > 1)
                                {
                                  chain <- rtmc3_chains[[k]][num,];
                                  eps <- rnorm(1,0,scale);
                                  b <- sample(c(-1,+1),length(base),replace=TRUE)
                                  chain <- tmcmcUpdate(chain,b,eps,target_pdf);
                                  out <- rbind(rtmc3_chains[[k]],chain);
                                 }
                                 return(out)
                              }
                          )

  for(k in 2:length(beta_set))
  {
    chain1 <- rtmc3_chains[[k]][num,];
    chain2 <- rtmc3_chains[[k-1]][num,];
    swap_rate <- min(1, exp((beta_set[k] - beta_set[(k-1)])*(target_pdf(chain2)-target_pdf(chain1))));
    w <- runif(0,1)
    if(w < swap_rate){
      rtmc3_chains[[k]][num,] <- chain2;
      rtmc3_chains[[k-1]][num,] <- chain1;
    }
  }
  }

  posterior_mean <- apply(rtmc3_chains[[1]], 2, mean);
  ll <- list("chain_set"=rtmc3_chains,"post.mean"=posterior_mean);
  return(ll)
}
