.rand_generate <- function(unipdf, method=c("TMCMC", "RWMH"))
{
  base <- 0
  scale <- 1
  nsamples <- 1000;
  if(method=="TMCMC")
    out <- tmcmc_metrop(unipdf, scale, base, nsamples)$chain[1000,];
  if(method=="RMWH")
    out <- rwmh_metrop(unipdf, scale, base, nsamples)$chain[1000,];
  return(out)
}
