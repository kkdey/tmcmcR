.rand_generate <- function(unipdf, method=c("TMCMC", "RWMH"), nsamples=1000, base=0, scale=1)
{
  if(method=="TMCMC")
    out <- tmcmc_metrop(unipdf, scale, base, nsamples)$chain[1000,];
  if(method=="RMWH")
    out <- rwmh_metrop(unipdf, scale, base, nsamples)$chain[1000,];
  return(out)
}

