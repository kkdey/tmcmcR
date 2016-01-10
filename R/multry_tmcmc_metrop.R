#' @title Simulate a Multiple Try Transformation based Markov Chain Monte Carlo (MT-TMCMC) algorithm (Rcpp sped up version)
#'
#' @description The function simulates a Multiple try TMCMC (MT-TMCMC) chain of length nsamples using the scale, base and
#' burn in taken optimally as default or specified by user. The function outputs the full chain as well as the estimated
#' posterior mean estimated from the samples drawn post burn-in.
#'
#' @param target_pdf The log target density function from which the user wants to generate samples.
#' @param scale The proposal density scaling parameter. An approximation of the optimal scaling given the target_pdf is performed by OptimalScaling().
#'              The default scale is this estimated optimal scaling
#' @param base The starting value of the chain
#' @param nmove_size The number of distinct move sizes used.
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


mt_tmcmc_metrop <- function(target_pdf, scale, base, nmove_size, nmove, nsamples, burn_in=NULL, verb=TRUE)
{
  if(is.null(burn_in)) burn_in <- nsamples/3;
  chain <- matrix(0, nsamples, length(base))
  chain[1,] <- base;
  num <- 2
  while(num <= nsamples) {
    eps <- rep(rnorm(nmove_size,0,scale), each=floor(nmove/nmove_size));
    b <- t(sapply(1:nmove, function(l) sample(c(-1,+1),length(base),replace=TRUE)));
    trial_chains <- chain[(num-1),]+ b*eps;
    pi_trials <- unlist(lapply(1:nmove, function(x) target_pdf(trial_chains[x,])))
    pi_trials_norm <- exp(pi_trials - max(pi_trials));
    move_candidate <- trial_chains[sample(1:nmove, size=1,  prob=pi_trials_norm),];

    b_rev <- t(sapply(1:(nmove-1), function(l) sample(c(-1,+1),length(base),replace=TRUE)));
    trial_chains_rev <- move_candidate + b_rev*eps[1:(nmove-1)];
    trial_chains_rev <- rbind(trial_chains_rev, chain[(num-1),]);
    pi_trials_rev <- unlist(lapply(1:nmove, function(x) target_pdf(trial_chains_rev[x,])))
    pi_trials_norm_rev <- exp(pi_trials_rev - max(pi_trials_rev));

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

