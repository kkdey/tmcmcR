#' @title Simulate adaptive RWMH algorithm (SCAM, Atchade and RAMA method)
#' @description The function simulates adaptive RWMH chain of length nsamples using one of the three
#'              methods of adaptation - SCAM, Atchade and RAMA.
#'
#' @param target_pdf The target density function from which the user wants to generate samples.
#' @param base The starting value of the chain
#' @param nsamples The number of samples to be drawn using the RWMH algorithm.
#' @param def.scale the initial scale of the proposal chosen, on which we update adaptively. Default is 1.
#' @param burn_in The number of samples assigned as burn-in period. The default burn-in is taken to be one-third of nsamples.
#' @param a_rama The scaling used in RAMA if norm of the chain at current iterate is less than dimension.
#' @param b_rama The scaling used in RAMA if norm of the chain at current iterate is greater than dimension.
#' @param M_rama The bounds on a_rama and b_rama, lower bound is -M_rama, upper bound is M_rama.
#' @param atchade_low The lower bound of scale for Atchade scheme
#' @param atchade_high The upper bound of scale for Atchade scheme
#' @param method The method of adaptation used for generating the chain. May be one of 3 types- SCAM,
#'               Atchade and the RAMA methods.
#' @param verb logical parameter, if TRUE the function prints the progress of simulation.
#'
#' @return Returns a list containing the following items
#' \item{chain}{The full chain produced by the adaptive RWMH method with user-chosen method of adaptation.}
#' \item{post.mean}{The estimated posterior mean of the adaptive RWMH chain adjusting for burn-in.}
#'  @author  Kushal K Dey
#'
#'  @useDynLib tmcmcR
#'  @export

adapt_rwmh_metrop <- function(target_pdf, base, nsamples, def.scale =1, burn_in=NULL, a_rama=NULL,
                               b_rama=NULL, M_rama =NULL,
                               atchade_low = NULL, atchade_high=NULL,
                               method=c("Atchade","SCAM","Rama"),
                               verb=TRUE)
{
  if(is.null(burn_in)) burn_in <- nsamples/3;
  scale <- def.scale;
  chain <- matrix(0, nsamples, length(base))
  chain[1,] <- base;
  num=2;

  if(method=="Atchade"){
                          if(is.null(atchade_low)) atchade_low <- 0.005;
                          if(is.null(atchade_high)) atchade_high <- 50;
                          while(num <= nsamples) {
                            eps <- rnorm(length(base),0,scale);
                            out <- rwmhUpdate(chain[(num-1),],eps,target_pdf)
                            chain[num,] <- out$chain;
                            scale <- scale + (1/num) *(out$acc_rate - 0.234);
                            if(scale < atchade_low) scale <- atchade_low;
                            if(scale > atchade_high) scale <- atchade_high;
                            if(num %% 500 ==0){
                              if(verb){
                                        cat("The chain is at iteration:",num,"\n");
                                      }
                            }
                            num <- num + 1;
                            }
                        }
  if(method=="SCAM"){
                        scale_vec <- rep(5, length(base));
                        while(num <= nsamples) {
                          eps <- rnorm(length(base),0,scale_vec);
                          out <- rwmhUpdate(chain[(num-1),],eps,target_pdf)
                          chain[num,] <- out$chain;
                          if(num > 10){
                              scale_vec <- sqrt((2.4)^2 * apply(chain[1:num,], 2, var) + 0.005);
                              }
                          if(num <= 10){
                              scale_vec <- rep(5, length(base));
                              }
                          if(num %% 500 ==0){
                              if(verb){
                                        cat("The chain is at iteration:",num,"\n");
                                      }
                              }
                          num <- num + 1;
                          }
                    }
  if(method=="Rama"){
                        if(is.null(a_rama)) a_rama <- 0
                        if(is.null(b_rama)) b_rama <- 0
                        if(is.null(M_rama)) M_rama <- 100
                        while(num <= nsamples) {
                              if (norm(chain[(num-1),],type='2') < sqrt(length(base))){
                                flag <- 0;
                                scale <- exp(2*a_rama)
                                }else{
                                flag <- 1;
                                scale <- exp(2*b_rama)
                                }
                                eps <- rnorm(length(base),0,scale);
                                out <- rwmhUpdate(chain[(num-1),],eps,target_pdf)
                              if(flag ==1)
                              {
                                  b_rama <- b_rama + as.numeric(out$acc_rate > 0.234)*min(0.01, 1/sqrt(num));
                                  if(b_rama > M_rama) b_rama <- M_rama
                                  if(b_rama <= -M_rama) b_rama <- -M_rama
                              }
                              if(flag ==0)
                              {
                                  a_rama <- a_rama + as.numeric(out$acc_rate > 0.234)*min(0.01, 1/sqrt(num));
                                  if(a_rama > M_rama) a_rama <- M_rama
                                  if(a_rama <= -M_rama) a_rama <- -M_rama
                              }
                              chain[num,] <- out$chain;
                              if(num %% 500 ==0){
                                  if(verb){
                                          cat("The chain is at iteration:",num,"\n");
                                          }
                                  }
                              num <- num + 1;
                        }

                  }
  if(length(base) > 1)
    posterior_mean <- apply(chain[round(burn_in):nsamples,],2,mean);
  if(length(base)==1)
    posterior_mean <- mean(chain[round(burn_in):nsamples,]);
  ll <- list("chain"=chain,"post.mean"=posterior_mean);
  return(ll)
}
