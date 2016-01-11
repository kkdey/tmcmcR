## Running tmcmcR on some simple illustrative examples

library(devtools)
#install_github('kkdey/tmcmcR')
library(tmcmcR)
library(mcmc)
d=30;  ##  dimension of the simulated variable
L=30; ###   the number of replications we use for finding KS statistic
nsamples <- 5000;
Mult_Mattingly=array(0,c(2,L,nsamples,d));


mu_target=rep(0,d);
Sigma_target = 0.01*diag(1/(1:(d))*d);

L=30; ###   the number of replications we use for finding KS statistic

Mult_Mattingly=array(0,c(2,L,nsamples,d));


Mattingly_matrix <- 100*(diag(1-0.7,d)+0.7*rep(1,d)%*%t(rep(1,d)));

library(mvtnorm)

pdf = function(x)
{
  return (dmvnorm(x,mu_target,Sigma_target,log=TRUE)-t(x)%*%Mattingly_matrix%*%x)
}

base=rnorm(d,0,1);

for ( l in 1:L)
{
  Mult_Mattingly[1,l,,] <- tmcmcR:::mt_tmcmc_metrop(pdf, base=base, scale=1,nmove_size=5, nmove=50, nsamples=1000,burn_in = NULL)$chain;
  Mult_Mattingly[2,l,,] <- tmcmcR:::tmcmc_metrop(pdf, base=base, scale=1, nsamples=1000,burn_in = NULL)$chain;
  cat("We are at iter:",l, "\n")
}


KSval_TMCMC=array(0,nsamples);
KSval_mtry_TMCMC=array(0,nsamples);

for(d in 1:30)
{
  for(n in 1:nsamples)
  {
    simulate.vec <- rnorm(L,mu_target[d],Sigma_target[d,d]/(sqrt(2*Sigma_target[d,d]^2+1)));
    KSval_mtry_TMCMC[n]=ks.test(Mult_Mattingly[1,,n,d],simulate.vec)$statistic;
    KSval_TMCMC[n]=ks.test(Mult_Mattingly[2,,n,d],simulate.vec)$statistic;

  }

  plot(1:nsamples,KSval_mtry_TMCMC,col="red",type="l",lwd=1,pch=2,xlab="",ylab="")
  lines(1:nsamples,KSval_TMCMC,col="blue",lwd=1,pch=3)
  title(xlab="Time step of run");
  title(ylab="KS test distance");
  title(main="KS plot comparison");
  legend("topright",c("mtry_TMCMC","TMCMC"),fill=c("red","blue"),border="black");
}



