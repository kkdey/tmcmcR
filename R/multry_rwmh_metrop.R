#' @title Simulate a Multiple Try RandomWalk Metropolis Hastings (MT-RWMH) algorithm (Liu, Liang and Wong 2000)
#'
#' @description The function simulates a Multiple try RWMH (MT-RWMH) chain of length nsamples using the scale, base and
#' burn in taken optimally as default or specified by user. The function outputs the full chain as well as the estimated
#' posterior mean estimated from the samples drawn post burn-in.
#'
#' @param target_pdf The log target density function from which the user wants to generate samples.
#' @param scale The proposal density scaling parameter. An approximation of the optimal scaling given the target_pdf is performed by OptimalScaling().
#'              The default scale is this estimated optimal scaling
#' @param base The starting value of the chain
#' @param nmove The total number of moves proposed.
#' @param nsamples The number of samples to be drawn using the TMCMC algorithm.
#' @param burn_in The number of samples assigned as burn-in period. The default burn-in is taken to be one-third of nsamples.
#' @param verb logical parameter, if TRUE the function prints the progress of simulation.
#'
#' @return Returns a list containing the following items
#' \item{chain}{The full chain produced by the MT-TMCMC method.}
#' \item{post.mean}{The estimated posterior mean of the TMCMC chain adjusting for burn-in.}
#
#
#'  @author  Kushal K Dey
#'
#'  @useDynLib tmcmcR
#'  @import tmcmcR
#'  @export
#'
#'



mt_rwmh_metrop <- function(target_pdf, scale, base, nmove, nsamples, burn_in=NULL, verb=FALSE)
{
  if(is.null(burn_in)) burn_in <- nsamples/3;
  chain <- matrix(0, nsamples, length(base))
  chain[1,] <- base;
  num <- 2
  while(num <= nsamples) {
    eps <- matrix(rnorm(nmove*length(base),0,scale), nrow=nmove);
    trial_chains <- chain[(num-1),]+ eps;
    pi_trials <- unlist(lapply(1:nmove, function(x) target_pdf(trial_chains[x,])))
    pi_trials_norm <- exp(pi_trials - max(pi_trials));
    index <- sample(1:nmove, size=1,  prob=pi_trials_norm);
    move_candidate <- trial_chains[index,];
    eps_rev <- rbind(matrix(rnorm((nmove-1)*length(base),0,scale), nrow=(nmove-1)),eps[index,]);
    trial_chains_rev <- move_candidate + eps_rev;
    trial_chains_rev <- rbind(trial_chains_rev, chain[(num-1),]);
    pi_trials_rev <- unlist(lapply(1:nmove, function(x) target_pdf(trial_chains_rev[x,])))
    pi_trials_norm_rev <- exp(pi_trials_rev - max(pi_trials));

    acc_rate <- min(1, sum(pi_trials_norm)/sum(pi_trials_norm_rev));
    if(runif(1,0,1) < acc_rate) {chain[num,] <- move_candidate;}
    else{ chain[num,] <- chain[(num-1),];}

    if(verb){
      if(num %% 500 == 0)
        cat("The chain is at iteration:",num,"\n");
    }
    num <- num + 1;
  }
  if(length(base) > 1)
    posterior_mean <- apply(chain[round(burn_in):nsamples,],2,mean);
  if(length(base)==1)
    posterior_mean <- mean(chain[round(burn_in):nsamples,]);
  ll <- list("chain"=chain,"post.mean"=posterior_mean);
  return(ll)
}
