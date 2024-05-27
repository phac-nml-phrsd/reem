
The name of this R package, `reem`, stands for: **R**enewal **E**quation
based **E**pidemic **M**odel

This package simulates and fits infectious disease epidemic to clinical
and wastewater data.

To install: `devtools::install_github("phac-nml-phrsd/reem")`

**Note:** This package implements the model as a __Reference Class__. This is a way to access the Object Oriented Programming functionalities in R. 
Developers, if you are not familiar with Reference Classes in R, please see this [very short introduction](https://www.datamentor.io/r-programming/reference-class) and [Hadley Wickham's](http://adv-r.had.co.nz/R5.html) course for more details. 

## Model description

The epidemic model is a traditional SIR model but implemented as a
renewal equation instead of the more populate ordinary differential
equations (ODE). It has been shown that the renewal equation
implementation is equivalent to the ODE one. See for example: D. Fargue,
Reducibilite des systemes hereditaires, Int. J. Nonlinear Mech., 9
(1974), pp. 331–338, D. Breda et al On the formulation of epidemic
models (an appraisal of Kermack and McKendrick), J. Biol. Dyn.,6 (2012),
pp. 103–117, and Champredon et al. Equivalence of the Erlang-Distributed
SEIR Epidemic Model and the Renewal Equation, SIAM J. Appl. Math., 78
(2018).

The renewal equation of the pathogen transmission process is as follows:

![i(t) = \left(\frac{S(t)}{N}\right)^{\alpha}\\ \mathcal{R}\_0 \\ B(t) \sum\_{k=1}^\ell g(k)i(t-k)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%28t%29%20%3D%20%5Cleft%28%5Cfrac%7BS%28t%29%7D%7BN%7D%5Cright%29%5E%7B%5Calpha%7D%5C%2C%20%5Cmathcal%7BR%7D_0%20%5C%2C%20B%28t%29%20%5Csum_%7Bk%3D1%7D%5E%5Cell%20g%28k%29i%28t-k%29 "i(t) = \left(\frac{S(t)}{N}\right)^{\alpha}\, \mathcal{R}_0 \, B(t) \sum_{k=1}^\ell g(k)i(t-k)")

![S(t) = \max(0, S(t-1) - i(t))](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%28t%29%20%3D%20%5Cmax%280%2C%20S%28t-1%29%20-%20i%28t%29%29 "S(t) = \max(0, S(t-1) - i(t))")

where
![i(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%28t%29 "i(t)")
is the incidence at time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t"),
![\mathcal{R}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BR%7D_0 "\mathcal{R}_0")
is the basic reproduction number,
![S_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_t "S_t")
are the number of susceptible individuals at time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t"),
![N](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N "N")
is the total population size,
![g](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g "g")
is the intrinsic generation interval distribution and
![B](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B "B")
a function to change the transmission rate with respect to time (e.g.,
to implement public health measures, vaccination campaigns, etc.).

## Simulation example

The code below simulates an epidemic, tracking the spread in the
population and the pathogen concentration in the wastewater.

``` r
library(reem)
library(snowfall)
```

``` r
# Define model parameters

prms = list(
  horizon  = 90,  # horizon of the simulation
  last.obs = 88,  # last observation time (must be < horizon)
  B        = rep(1,90), # time dependent multiplicative factor for transmission
  i0prop  = 1e-3,  # initial proportion of the population infected
  date.start = lubridate::ymd('2022-01-01'), # start date of the epidemic
  start.delta = 0, 
  R0      = 2.0, # Basic reproduction number
  N       = 1e4, # population size
  alpha   = 0, # transmission heterogeneity (alpha=0: homogeneous)
  I.init  = c(1,1,3,5), # initial incidence (overwritten in fit ABC)
  lag     = 7,   # Aggregation window for clinical reports
  rho     = 0.1, # mean reporting ratio
  g       = get_gi(), # Generation interval distribution
  fec     = get_fecalshed(), # fecal shedding kinetics
  kappa   = 0.18, # pathogen decay in ww
  psi     = get_psi(), # plug flow simulation,
  shed.mult = 0.2 # deposited fecal shedding multiplier  
)

# Create the model object instance
obj = new('reem', 
          name = 'foo', 
          prms = prms, 
          is.fitted = FALSE)

# Print (formatted) the model parameters
obj$print_prms()
```

    ## 
    ## --- Parameters for REEM ` foo `
    ## horizon  =  90 
    ## last.obs  =  88 
    ## B  =  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
    ## i0prop  =  0.001 
    ## date.start  =  18993 
    ## start.delta  =  0 
    ## R0  =  2 
    ## N  =  10000 
    ## alpha  =  0 
    ## I.init  =  1 1 3 5 
    ## lag  =  7 
    ## rho  =  0.1 
    ## g  =  0 0.02432345 0.3100143 0.3955746 0.1967229 0.05837831 0.01249735 0.002135546 0.000309429 3.952736e-05 4.571663e-06 
    ## fec  =  0.001341379 0.04085834 0.1509359 0.2350638 0.2280526 0.1640111 0.09593837 0.04824391 0.02161408 0.008840312 0.003358844 0.001200831 0.0004079245 0.0001326712 
    ## kappa  =  0.18 
    ## psi  =  0.85 0.1 0.05 
    ## shed.mult  =  0.2 
    ##  --------------------------

``` r
# Simulate an epidemic based on the parameters
simepi  = obj$simulate_epi(deterministic = FALSE)
```

We can use the built-in function `plot_epi()` to plot the epidemic

``` r
g = plot_epi(simepi)
plot(patchwork::wrap_plots(g, ncol = 1))
```

![](README_files/figure-gfm/plot_epi-1.png)<!-- -->

## Fit example

If we attach observation data (from clinical and/or wastewater
surveillance), to a `reem` object, then it is possible to fit its model
parameters to these data. The example below shows how to do this.

First, let’s create synthetic observation data from the simulation
example above. We use synthetic data (as opposed to real data) because
we to assess the accuracy of the fit. The synthetic data are generated
from the vey same model `reem` where we know all the parameters.

``` r
# we use the example above and make a copy of its
# simulated epidemics. The time series will be the "observations".

sim.data = simepi

# Set the date for "today"
asof = lubridate::ymd('2023-03-01') 

# Retrieve the simulated observations
obs.cl = dplyr::filter(sim.data$obs.cl, date <= asof)
obs.ww = dplyr::filter(simepi$obs.ww, date <= asof)

# shift the dates such that they 
# do not perfectly align with the simulation
obs.cl$t <- obs.cl$t + 1
obs.cl$date <- obs.cl$date + 1

# Attached (simulated) data to new `reem` object:
prms$t.obs.cl <- NULL
prms$t.obs.ww <- NULL
obj  = new('reem', 
           name = 'foo2', 
           prms = prms, 
           obs.cl = obs.cl,
           obs.ww = obs.ww,
           is.fitted = FALSE)
```

Now, we set up the fitting algorithm. This algorithm is an Approximate
Bayesian Computation (ABC), which is a very straightforward and robust
Bayesian method but not very efficient.

``` r
# Define the parameters for the ABC fitting algorithm
prm.abc = list(
  n.abc = 500,   # Total number of ABC iterations (the larger the better)
  n.sim = 0,     # Number of simulation for a given set of prior parameters. `0` for deterministic, else`8` should be enough
  p.abc = 0.02,  # Acceptance probability (the lower the better)
  n.cores = 4,   # number of cores used for parallel computing
  use.cl = TRUE, # use clinical observations in the fit?
  use.ww = TRUE, # use wastewater observation in the fit?
  err.type = 'L2'# Type of error calculated during the fit.
)

# Define the priors of the parameters to be fitted:
prms.to.fit = list(
  # Gamma distrib. with mean 2.5 and variance 0.25:
  R0          = list('gamma', 2.5, 0.251), 
  # Other prior distributions definition self-explanatory:
  alpha       = list('norm', 2, 3),
  i0prop      = list('unif', -5, -2),
  start.delta = list('unif_int', -7, 7) # `unif_int` = uniform with integers
)
```

Now that all the parameters for the ABC fit have been defined, we can
launch the actual fit. Note that we use the object `obj` already built
and we call its *member function* named `fit_abc()`:

``` r
#  Launch the fit 
thefit = obj$fit_abc(prm.abc, prms.to.fit)  
```

    ## 
    ## ----- ABC FIT -----
    ## 
    ## Target data sources :
    ##   clinical   = yes
    ##   wastewater = yes
    ## 
    ## Number of priors     : 500
    ## Number of posteriors : 10 (accept ratio = 0.02)
    ## Number of cores      : 4 (125 iters per core)
    ## 
    ## Data horizon : 92 (days) 
    ## 
    ## ---------------------

    ## Warning in searchCommandline(parallel, cpus = cpus, type = type, socketHosts =
    ## socketHosts, : Unknown option on commandline:
    ## rmarkdown::render('/Users/davidchampredon/GitHub/reem/README.Rmd',~+~~+~encoding~+~

    ## R Version:  R version 4.2.0 (2022-04-22)

    ## snowfall 1.84-6.3 initialized (using snow 0.4-4): parallel execution on 4 CPUs.

    ## Library dplyr loaded.

    ## Library dplyr loaded in cluster.

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Library purrr loaded.

    ## Library purrr loaded in cluster.

    ## 
    ## Stopping cluster

Once the ABC fit is completed, we can call built-in functions to assess
the goodness of fit.

First, lets visualize the posterior time series trajectories for the the
clinical cases and wastewater concentrations.

``` r
gg = obj$plot_fit()

plot(gg$traj.cl)
```

![](README_files/figure-gfm/plot_fit_results-1.png)<!-- -->

``` r
plot(gg$traj.ww)
```

![](README_files/figure-gfm/plot_fit_results-2.png)<!-- -->

We can also plot the posterior distributions of the fitted parameters
(with their prior distribution overlayed in grey):

``` r
plot(gg$post.prms)
```

![](README_files/figure-gfm/plot%20posteriors-1.png)<!-- -->

This is plotting the pairwise 2D posterior densities, to check for any
strong correlations between fitted parameters

``` r
plot(gg$post.prms.2d)
```

![](README_files/figure-gfm/plot%202D%20posteriors-1.png)<!-- -->

And finally, we can assess the efficiency of the ABC algorithm by
plotting the posterior distance values as a function of the sorted
iterations

``` r
plot(gg$dist)
```

![](README_files/figure-gfm/plot%20fitted%20distances-1.png)<!-- -->
