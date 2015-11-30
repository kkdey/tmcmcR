#' @title Simulate a Random Walk Metropolis Hastings algorithm (Rcpp sped up version)
#' @description The function simulates a RWMH chain of length nsamples using the scale, base and burn in taken optimally as default or specified by user. The function
#' outputs the full chain as well as the estimated posterior mean estimated from the samples drawn post burn-in.
#'
#' @param target_pdf The target density function from which the user wants to generate samples.
#' @param scale The proposal density scaling parameter. An approximation of the optimal scaling given the target_pdf is performed by OptimalScaling().
#'              The default scale is this estimated optimal scaling
#' @param base The starting value of the chain
#' @param nsamples The number of samples to be drawn using the RWMH algorithm.
#' @param burn_in The number of samples assigned as burn-in period. The default burn-in is taken to be one-third of nsamples.
#' @param verb logical parameter, if TRUE the function prints the progress of simulation.
#'
#' @return Returns a list containing the following items
#' \item{chain}{The full chain produced by the RWMH method.}
#' \item{post.mean}{The estimated posterior mean of the RWMH chain adjusting for burn-in.}
#
#'
#
#'  @author  Kushal K Dey
#'
#'  @useDynLib tmcmcR
#'  @export


rwmh_metrop <- function(target_pdf, scale, base, nsamples, burn_in=NULL, verb=TRUE)
{
  if(is.null(burn_in)) burn_in <- nsamples/3;
  chain <- matrix(0, nsamples, length(base))
  chain[1,] <- base;
  num=2;
  while(num <= nsamples) {
    eps <- rnorm(length(base),0,scale);
    chain[num,] <- rwmhUpdate(chain[(num-1),],eps,target_pdf)$chain
    if(verb){
    if(num %% 500 ==0)
      cat("The chain is at iteration:",num,"\n");
    }
    num <- num +1;
  }
  if(length(base) > 1)
    posterior_mean <- apply(chain[round(burn_in):nsamples,],2,mean);
  if(length(base)==1)
    posterior_mean <- mean(chain[round(burn_in):nsamples,]);
  ll <- list("chain"=chain,"post.mean"=posterior_mean);
  return(ll)
}

