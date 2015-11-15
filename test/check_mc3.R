## Running tmcmcR on some simple illustrative examples

## RTMC-3 (Randomized Transformation based Maetropolis Couple Markov Chain Monte Carlo) application

library(devtools)
#install_github('kkdey/tmcmcR')
library(tmcmcR)

library(mcmc)
d=50;  ##  dimension of the simulated variable
L=20; ###   the number of replications we use for finding KS statistic
nsamples <- 5000;
Mult_Mattingly=array(0,c(3,L,nsamples,d));


mu_target=rep(0,d);
Sigma_target = 0.01*diag(1/(1:(d))*d);

Mult_Mattingly=array(0,c(2,L,nsamples,d));

Mattingly_matrix <- 100*(diag(1-0.7,d)+0.7*rep(1,d)%*%t(rep(1,d)));

library(mnormt)
library(fMultivar)
library(mvtnorm)
library(parallel)

pdf = function(x)
{
  return (dmvnorm(x,mu_target,Sigma_target,log=TRUE)-t(x)%*%Mattingly_matrix%*%x)
}

base=rnorm(d,0,1);
beta_set <- seq(1,0.05,length.out=10);

for ( l in 1:L)
{
  Mult_Mattingly[1,l,,] <- rtmc3(pdf,beta_set=beta_set,base=base, scale=1, cycle=100, nsamples=2000, verb=FALSE)$chain_set[[1]];
  Mult_Mattingly[2,l,,] <- rmc3(pdf,beta_set=beta_set, base=base, scale=1, cycle=100, nsamples=1000, verb=FALSE)$chain_set[[1]];
  cat("We are at iter:",l, "\n")
}

KSval1=array(0,nsamples);
KSval2=array(0,nsamples);

for(d in 1:40)
{
  for(n in 1:nsamples)
  {
    simulate.vec <- rnorm(L,mu_target[d],Sigma_target[d,d]/(sqrt(2*Sigma_target[d,d]^2+1)));
    KSval1[n]=ks.test(Mult_Mattingly[1,,n,d],simulate.vec)$statistic;
    KSval2[n]=ks.test(Mult_Mattingly[2,,n,d],simulate.vec)$statistic;
  }

  plot(1:nsamples,KSval1,col="red",type="l",lwd=1,pch=2,xlab="",ylab="")
  lines(1:nsamples,KSval2,col="blue",lwd=1,pch=3)
  title(xlab="Time step of run");
  title(ylab="KS test distance");
  title(main="KS plot comparison");
  legend("topright",c("RTMC3","RMC3"),fill=c("red","blue"),border="black");
}


