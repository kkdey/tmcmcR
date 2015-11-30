#' @title Select inverse temperatures in RMC3 and RTMC3
#' @description The function selects the inverse temperatures using a Stochastic approximation
#' algorithm for RMC3 or RTMC3 chains.
#'
#' @param pdf_component The univariate marginal distribution of the iid product component target distribution
#' @param minbeta the cut off of inverse temperature which when crossed, we stop at that value of inverse
#'        temperature. Helps to check the lower bound of the inverse temperatures. Default is 0.01.
#' @param L_iter The number of sub-iterations required to fix each iterate of the inverse temperatures
#' @param sim_method The method used to simulate from the marginal pdf component. Choices include TMCMC and RWMH.
#' @param inv_temp_scheme The selection of inverse temperatures scheme followed, affects the stochastic approximation
#'        driving the inverse temperature updates. Choices include "fixed" and "randomized".
#' @param rho_start The scale update for the inverse temperature selection. Default is 0.
#' @param scale The stochastic approximation scaling parameter. Default is 0.1. The smaller this value,
#'        more inverse temperatures will be selected and algorithm will be more fine-grained.
#'
#' @return Returns a vector of inverse temperatures starting from 1 and ending at the first value of inverse
#'        temperature smaller than minbeta.
#'
#' @author  Kushal K Dey
#'
#' @export
#'


select_inverse_temp <- function(pdf_component, minbeta=0.01, L_iter =50,
                                sim_method=c("RWMH","TMCMC"),
                                inv_temp_scheme = c("randomized","fixed"),
                                rho_start=0,
                                scale=0.1)
{
  beta_array <- 1;
  counter <- 1
  current_beta = 1;
  while(current_beta > minbeta)
  {
    rho <- rho_start;
    for (l in 1:L_iter)
    {
      temp_beta <- current_beta * (1/(1 + exp(rho)));
      pdf_1 <- function(x) { return(pdf_component(x)*current_beta)};
      pdf_2 <- function(x) { return(pdf_component(x)*temp_beta)};

      x_curr <- .rand_generate(pdf_1, method=sim_method);
      x_temp <- .rand_generate(pdf_2, method=sim_method);

      B <- -(temp_beta - current_beta)* (pdf_2(x_temp) - pdf_1(x_curr));
      alpha <- min(1, exp(B));
      if(inv_temp_scheme=="randomized")
        rho <- rho + (scale/l) * (alpha - 0.44);
      if(inv_temp_scheme=="fixed")
        rho <- rho + (scale/l) * (alpha - 0.234);
    }

    current_beta <- current_beta * (1/(1 + exp(rho)));
    paste("The inverse temperature selected:", counter);
    beta_array <- c(beta_array, current_beta);
    counter <- counter + 1;
  }
  return(beta_array)
}
