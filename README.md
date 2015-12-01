# tmcmcR

A R/C++ package for fast implementataion of Transformation based Markov Chain Monte Carlo (TMCMC) Methods, which has higher acceptance rate and faster
convergence in many scenarios, compared to the standard Random Walk Metropolis Hastings (RWMH) approach. The package implements R/C++ implementtaions
for both the standard TMCMC and the RWMH algorithms, but can also be used for adaptive versions of both methods ([SCAM](http://link.springer.com/article/10.1007%2FBF02789703), 
[RAMA](http://probability.ca/jeff/ftpdir/adaptex.pdf) and 
[Atchade](https://projecteuclid.org/euclid.bj/1130077595) methods). The package also implements Metropolis coupling variations 
of the TMCMC and the RWMH approaches (applicable for simulating from multimodal densities), together with a inverse temperature selection scheme with both fixed or randomized 
step designs.

### Description

- Version: 0.1.2
- Authors: [Kushal K Dey](http://kkdey.github.io/), [Sourabh Bhattacharya](http://www.isical.ac.in/~biru/sb.html).
- Maintainer: Kushal K. Dey
- License: GPL (>=2)

### Installation

To install the **tmcmcR** package, 

```
library(devtools)
install_github('kkdey/tmcmcR')
```

To load the package,

```
require(tmcmcR)
```
### Functions

The chains we implement in this package, together with the exact functionality, include

- **rwmh_metrop( )** : Standard Random Walk Metropolis Hastings (RWMH)
- **tmcmc_metrop( )**: Standard Transformation based Markov Chain Monte Carlo (TMCMC)
- **adapt_rwmh_metrop( )** : Adaptive Random Walk Metropolis Hastings, method: *SCAM, RAMA, Atchade*
- **adapt_tmcmc_metrop( )**: Adaptive Transformation based Markov Chain Monte Carlo, method: *SCAM, RAMA, Atchade*
- **rmc3( )**: Metropolis Coupling Markov Chain Monte Carlo with RWMH update 
- **rtmc3( )**: Metropolis Coupling Markov Chain Monte Carlo with TMCMC update (TMC3)

The last two methods need selection of inverse temperatures to run multiple chains and use these parallel chains for swapping
at regular intervals. For target densities with product of known iid components, we have a selection scheme for inverse temperatures
which may be of two types, depending on if the step size is fixed or randomized.

- **select_inverse_temp( )**: Inverse temperature selection scheme, method: *TMCMC*, *RWMH*, scheme: *fixed*, *randomized*.

### Vignettes

For checking example usage of functions under the **tmcmcR** package, check out the [ vignette](https://rpubs.com/kkdey/132076).
Also check the Github folder [test](https://github.com/kkdey/tmcmcR/tree/master/test) for simulation examples for each function
in the package (we use these codes to generate the figures in the vignette).

### Citation

If you are using the **tmcmcR** R package, please cite 

KK Dey, S Bhattacharya. **tmcmcR: R package for MCMC with improved acceptance rate and coverage**. [![RPubs doc](icons/R-icon.png)](http://rpubs.com/kkdey/132076).

For the TMCMC methods and materials, please check

S Dutta, S Bhattacharya. **Markov Chain Monte Carlo Based on Deterministic Transformations**. [![pdf](icons/pdf-icon.png)](http://www.sciencedirect.com/science/article/pii/S1572312713000683)

KK Dey, S Bhattacharya. **On Optimal Scaling of Additive Transformation Based Markov Chain Monte Carlo**. Submitted *Brazilian Journal of Probability and Statistics*. [![pdf](icons/pdf-icon.png)](http://arxiv.org/pdf/1307.1446v4.pdf)

KK Dey, S Bhattacharya. **On Single Variable Transformation Approach to Markov Chain Monte Carlo**. *JSM Proceedings*. 2014.  [![pdf](icons/pdf-icon.png)](http://arxiv.org/abs/1408.6667)

## Contact

For any queries, contact [kkdey@uchicago.edu](kkdey@uchicago.edu)

## Acknowledgements

- Gao Wang, University of Chicago
- Nan Xiao, University of Chicago





