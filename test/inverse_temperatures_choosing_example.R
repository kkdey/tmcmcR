
## Select inverse temperatures

unipdf <- function(x) return (dnorm(x,0,1,log=TRUE))
d <- 1
nsamples <- 5000;

base <- 0

temp <- tmcmc_metrop(unipdf,base=base, scale=1,nsamples=5000,burn_in = NULL)$chain;

temp1 <- .rand_generate(unipdf, method="TMCMC")

inv_temp_fixed_tmcmc <- select_inverse_temp(unipdf, minbeta=0.05, L_iter =50, sim_method="TMCMC", inv_temp_scheme="fixed")
inv_temp_rand_tmcmc <- select_inverse_temp(unipdf, minbeta=0.05, L_iter =50, sim_method="TMCMC", inv_temp_scheme="randomized")
inv_temp_fixed_rwmh <- select_inverse_temp(unipdf, minbeta=0.05, L_iter =50, sim_method="RWMH", inv_temp_scheme="fixed")
inv_temp_rand_rwmh <- select_inverse_temp(unipdf, minbeta=0.05, L_iter =50, sim_method="RWMH", inv_temp_scheme="randomized")

print(inv_temp_fixed_tmcmc)
print(inv_temp_rand_tmcmc)
print(inv_temp_fixed_rwmh)
print(inv_temp_rand_rwmh)
