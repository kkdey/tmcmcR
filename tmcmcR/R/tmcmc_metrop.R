
tmcmc_metrop <- function(target_pdf, scale, base, nsamples, burn_in=NULL)
{
  if(is.null(burn_in)) burn_in <- nsamples/3;
  chain <- matrix(0, nsamples, length(base))
  num=1;
  chain[1,] <- base;
  while(num <= nsamples) {
    eps <- rnorm(1,0,scale);
    b <- sample(c(-1,+1),length(base),replace=TRUE)
    chain[(num+1),] <- tmcmcUpdate(chain[num,],b,eps,target_pdf)
  }
  posterior_mean <- apply(chain[round(burn_in):nsamples,],2,mean);
  ll <- list("chain"=chain,"post.mean"=posterior_mean);
  return(ll)
}

