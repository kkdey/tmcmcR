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
                               b_rama=NULL, M_rama =NULL, def.scale =1, verb=TRUE,
                               method=c("Atchade","SCAM","Rama"))
{
  if(is.null(burn_in)) burn_in <- nsamples/3;
  scale <- def.scale;
  chain <- matrix(0, nsamples, length(base))
  chain[1,] <- base;
  num=2;
  if(is.null(a_rama)) a_rama <- 0
  if(is.null(b_rama)) b_rama <- 0
  if(is.null(M_rama)) M_rama <- 100

  if(method=="Atchade"){
                          while(num <= nsamples) {
                          eps <- abs(rnorm(1,0,scale));
                          b <- sample(c(-1,+1),length(base),replace=TRUE)
                          out <- tmcmcUpdate(chain[(num-1),],b,eps,target_pdf)
                          chain[num,] <- out$chain;
                          scale <- scale + (1/num) *(out$acc_rate - 0.439);
                          if(num %% 500 ==0){
                            if(verb){
                              cat("The chain is at iteration:",num);
                            }
                          }
                          num <- num + 1;
                          }
  }
  if(method=="SCAM"){
                        scale <- 5;
                        store_eps <- array(0, nsamples);
                        while(num <= nsamples) {
                          eps <- abs(rnorm(1,0,scale));
                          b <- sample(c(-1,+1),length(base),replace=TRUE)
                          out <- tmcmcUpdate(chain[(num-1),],b,eps,target_pdf)
                          chain[num,] <- out$chain;
                          if(min(abs(out$chain-chain[num,]))==0){
                             store_eps[num] <- store_eps[(num-1)];
                            } else{
                              store_eps[num] <- eps;
                            }
                          if(num > 10){
                            scale <- sqrt((2.4)^2 * var(store_eps[1:num]) + 0.005);
                          }
                          if(num <= 10){
                            scale <- 5;
                          }
                          if(num %% 500 ==0){
                            if(verb){
                              cat("The chain is at iteration:",num);
                            }
                          }
                            num <- num + 1;
                          }
  }
  if(method=="Rama"){
                        while(num <= nsamples) {
                          if (norm(chain[(num-1),],type='2') < sqrt(length(base))){
                            flag <- 0;
                            scale <- exp(2*a_rama)
                          }else{
                            flag <- 1;
                            scale <- exp(2*b_rama)
                          }
                          eps <- abs(rnorm(1,0,scale));
                          b <- sample(c(-1,+1),length(base),replace=TRUE)
                          out <- tmcmcUpdate(chain[(num-1),],b,eps,target_pdf)
                          if(flag ==1)
                          {
                            b_rama <- b_rama + as.numeric(out$acc_rate > 0.44)*min(0.01, 1/sqrt(num));
                            if(b_rama > M_rama) b_rama <- M_rama
                            if(b_rama <= -M_rama) b_rama <- -M_rama
                          }
                          if(flag ==0)
                          {
                            a_rama <- a_rama + as.numeric(out$acc_rate > 0.44)*min(0.01, 1/sqrt(num));
                            if(a_rama > M_rama) a_rama <- M_rama
                            if(a_rama <= -M_rama) a_rama <- -M_rama
                          }
                          chain[num,] <- out$chain;
                          if(num %% 500 ==0){
                            if(verb){
                            cat("The chain is at iteration:",num);
                            }
                          }
                          num <- num + 1;
                        }

  }
  posterior_mean <- apply(chain[round(burn_in):nsamples,],2,mean);
  ll <- list("chain"=chain,"post.mean"=posterior_mean);
  return(ll)
}
