#' @title Simulate adaptive TMCMC algorithm (Haario, Atchade and RAMA method)
#'
#' @param target_pdf The target density function from which the user wants to generate samples.
#' @param scale The proposal density scaling parameter. An approximation of the optimal scaling given the target_pdf is performed by OptimalScaling().
#'              The default scale is this estimated optimal scaling
#' @param base The starting value of the chain
#' @param nsamples The number of samples to be drawn using the TMCMC algorithm.
#' @param burn_in The number of samples assigned as burn-in period. The default burn-in is taken to be one-third of nsamples.
#' @param a_rama The scaling coefficients (a[n]) used in RAMA
#' @param b_rama The scaling coefficients (b[n]) used in RAMA
#' @param method The method of adaptation used for generating the chain. May be one of 3 types- Haario,
#'               Atchade and the RAMA methods.
#'
#' @description The function simulates adaptive TMCMC chain of length nsamples using one of the three methods of adaptation - Haario, Atchade and RAMA.
#'              The other inputs and the output same as tmcmc_metrop function.
#'
#
#'  @author  Kushal K Dey
#'
#'  @useDynLib tmcmcR
#'  @export


adapt_tmcmc_metrop <- function(target_pdf, base, nsamples, burn_in=NULL, a_rama=NULL,
                               b_rama=NULL, method=c("Atchade","Haario","Rama"))
{
  Rcpp::sourceCpp('src/RcppExports.cpp')
  Rcpp::sourceCpp('src/utils.cpp')
  if(is.null(burn_in)) burn_in <- nsamples/3;
  scale <- 1;
  chain <- matrix(0, nsamples, length(base))
  chain[1,] <- base;
  num=2;

  if(method=="Atchade"){
                          while(num <= nsamples) {
                          eps <- rnorm(1,0,scale);
                          b <- sample(c(-1,+1),length(base),replace=TRUE)
                          out <- tmcmcUpdate(chain[(num-1),],b,eps,target_pdf)
                          chain[num,] <- out$chain;
                          scale <- scale + (1/num) *(out$acc_rate - 0.439);
                          if(num %% 500 ==0)
                            paste("The chain is at iteration:",num);
                          num <- num + 1;
                          }
  }
  if(method=="Haario"){
                          store_eps <- array(0, nsamples)
                          store_eps[1] <- def.scale;
                          while(num <= nsamples) {
                            eps <- rnorm(1,0,scale);
                            b <- sample(c(-1,+1),length(base),replace=TRUE)
                            out <- tmcmcUpdate(chain[(num-1),],b,eps,target_pdf)
                            chain[num,] <- out$chain;
                            if(min(abs(out$chain-chain[num,]))==0){
                              store_eps[num] <- store_eps[(num-1)];
                            } else{
                              store_eps[num] <- eps;
                            }
                            scale <- def.scale * var(unique(store[1:num])) + 0.001*def.scale;
                            if(num %% 500 ==0)
                              paste("The chain is at iteration:",num);
                            num <- num + 1;
                            }
  }
  if(method=="Rama"){
                        while(num <= nsamples) {
                          if (norm(chain[num,],type='2') < length(base)){
                            scale <- exp(2*a_rama[num])
                          }else{
                            scale <- exp(2*b_rama[num])
                          }
                          eps <- rnorm(1,0,scale);
                          b <- sample(c(-1,+1),length(base),replace=TRUE)
                          out <- tmcmcUpdate(chain[(num-1),],b,eps,target_pdf)
                          chain[num,] <- out$chain;
                          if(num %% 500 ==0)
                            paste("The chain is at iteration:",num);
                          num <- num + 1;
                        }

  }
  posterior_mean <- apply(chain[round(burn_in):nsamples,],2,mean);
  ll <- list("chain"=chain,"post.mean"=posterior_mean);
  return(ll)
}
