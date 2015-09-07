
mean <- c(0,0,0);
Sd <- diag(3);
normpdf <- function(x)
{
  return (dmvnorm(x,mean,Sd,log=TRUE))
}

tmcmcUpdate(c(0,0,0),c(1,1,1),0.5,normpdf)

normpdf(c(-400,-400,-400))
normpdf(c(100,100,100))
