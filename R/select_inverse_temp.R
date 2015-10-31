#' @title Select inverse temperatures using Stochatic approximation for swapping in product of iid components distribution (RMC3 and RTMC3)
#'
#'
#' @param pdf_component The marginal distribution of the iid product component target distribution
#' @param minbeta the minimum cut off of inverse temperature we allow in RMC3 or RTMC3 (default is 0.05)
#' @param L_iter The number of sub-iterations required to fix each iterate of the inverse temperatures
#' @param method The method used to simulate from the marginal pdf component. Choices include TMCMC and RWMH.
#'
#' @description The function simulates a TMCMC chain of length nsamples using the scale, base and burn in taken optimally as default or specified by user. The function
#' outputs the full chain as well as the estimated posterior mean estimated from the samples drawn post burn-in.
#'
#'  @author  Kushal K Dey
#'
#'  @export
#'

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

select_inverse_temp <- function(pdf_component, minbeta=0.05, L_iter =50,
                                method=c("RWMH","TMCMC"))
{
  beta_array <- c();
  counter <- 1
  current_beta = 1;
  while(current_beta > minbeta)
  {
    rho <- 0;
    for (n in 1:L_iter)
    {
      temp_beta <- current_beta * (1/(1 + exp(rho)));
      pdf_1 <- function(x, current_beta = current_beta) { return(pdf_component(x)*current_beta)};
      pdf_2 <- function(x, temp_beta=temp_beta) { return(pdf_component(x)*temp_beta)};

      x_curr <- .rand_generate(pdf_1, method=method);
      x_temp <- .rand_generate(pdf_2, method=method);

      B <- -(temp_beta - current_beta)* (pdf_2(x_temp,temp_beta) - pdf1(x_curr, current_beta));
      alpha <- min(1, exp(B));
      rho <- rho + (1/counter) * (alpha - 0.44);
    }

    current_beta <- current_beta * (1/(1 + exp(rho)));
    beta_array <- c(beta_array, current_beta);
  }
}
