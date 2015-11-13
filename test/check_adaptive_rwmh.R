## Running tmcmcR on some simple illustrative examples

## Adaptive set up

library(devtools)
#install_github('kkdey/tmcmcR')
library(tmcmcR)

library(mcmc)
d=50;  ##  dimension of the simulated variable
L=30; ###   the number of replications we use for finding KS statistic
nsamples <- 5000;
Mult_Mattingly=array(0,c(3,L,nsamples,d));


mu_target=rep(0,d);
Sigma_target = 0.01*diag(1/(1:(d))*d);

Mult_Mattingly=array(0,c(3,L,nsamples,d));

Mattingly_matrix <- 100*(diag(1-0.7,d)+0.7*rep(1,d)%*%t(rep(1,d)));

library(mnormt)
library(fMultivar)
library(mvtnorm)

pdf = function(x)
{
  return (dmvnorm(x,mu_target,Sigma_target,log=TRUE)-t(x)%*%Mattingly_matrix%*%x)
}

base=rnorm(d,0,1);

for ( l in 1:L)
{
  Mult_Mattingly[1,l,,] <- adapt_rwmh_metrop(pdf,base=base, nsamples=5000, method="Atchade", verb=FALSE)$chain;
  Mult_Mattingly[2,l,,] <- adapt_rwmh_metrop(pdf,base=base, nsamples=5000, method="SCAM", verb=FALSE)$chain;
  Mult_Mattingly[3,l,,] <- adapt_rwmh_metrop(pdf,base=base, nsamples=5000, method="Rama", verb=FALSE)$chain
  cat("We are at iter:",l, "\n")
}


KSval1=array(0,nsamples);
KSval2=array(0,nsamples);
KSval3=array(0,nsamples);

for(d in 1:40)
{
  for(n in 1:nsamples)
  {
    simulate.vec <- rnorm(L,mu_target[d],Sigma_target[d,d]/(sqrt(2*Sigma_target[d,d]^2+1)));
    KSval1[n]=ks.test(Mult_Mattingly[1,,n,d],simulate.vec)$statistic;
    KSval2[n]=ks.test(Mult_Mattingly[2,,n,d],simulate.vec)$statistic;
    KSval3[n]=ks.test(Mult_Mattingly[3,,n,d],simulate.vec)$statistic;
  }

  plot(1:nsamples,KSval1,col="red",type="l",lwd=1,pch=2,xlab="",ylab="")
  lines(1:nsamples,KSval2,col="blue",lwd=1,pch=3)
  lines(1:nsamples,KSval3,col="green",lwd=1,pch=3)
  title(xlab="Time step of run");
  title(ylab="KS test distance");
  title(main="KS plot comparison");
  legend("topright",c("Atchade","SCAM", "Rama"),fill=c("red","blue","green"),border="black");
}

