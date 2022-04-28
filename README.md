# Financial modeling

A package for financial modeling, and simulation of stochastic processes such as encountered in finance, e.g. Geometric Brownian Motions (Black-Scholes asset modeling) and Jump Decision Processes (Merton asset modeling).

This R-cpp package mostly allows to reproduce the graphics of the blogpost [Financial modeling: a simple statistical point of view](https://medium.com/towards-data-science/financial-modeling-a-clean-short-and-simple-statistical-point-of-view-853dc29efb00). This blogpost explains how to derive the main results of financial modeling such as option pricing under the Black-Scholes model (European and American exercise styles), hedging, portfolio optimization, jump assets models, exotic options pricing using simulation ...


## Install instructions

To install the R dependencies, run (if needed) from the R command line:
```{r}
install.packages(c("Rcpp", "RcppArmadillo", "cpp11"))
```

You need cuda installed (and a GPU of course). 
Check the CUDA_PATH on your computer in the makefile "./src/Makevars", as well as the R_HOME (this one is left unspecified).

Once you have correct paths specified, install the package from the console:
```console
cd PACKAGE_SOURCE_FOLDER
R CMD INSTALL -l /INSTALL_PATH ./ --preclean
```


## Demos

Demos scripts are in the folder ./demos.
There are 9 demo scripts which allow to reproduce the graphics of the blogpost. Plus two additional demos showing how to use cuda kernels to perform simulation of asset evolution under Black-Scholes and Merton models.

I refer you to the blog for the figures.

1. ./demos/gen_fig1.R: Simulations of Brownian Motion over a definite time horizon.
![Alt text](./results/fig1.png?raw=true "Simulations of Brownian Motion")

2. ./demos/gen_fig2.R: Black-Scholes surface obtained by making vary both the exercise price and the spot price.
![Alt text](./results/fig0.png?raw=true "Black-Scholes surface")

3. ./demos/gen_fig3.R: Determining numerically the optimal exercise price for American options.
![Alt text](./results/fig5.png?raw=true "Opimal exercise price for american options")

4. ./demos/gen_fig4.R: Realization of a jump decision process over a definitie time horizon.
![Alt text](./results/fig2.png?raw=true "Realization of a jump decision process over a definitie time horizon")

5. ./demos/gen_fig5.R: Log-normal law from trajectory simulations.
![Alt text](./results/fig3.png?raw=true "Log-normal law from trajectory simulations")

7. ./demos/gen_fig6.R: Sampling several Geometric Brownian Motion (GBM) trajectories starting from the same spot price.
![Alt text](./results/fig6.png?raw=true "Sampling GBM trajectories")

6. ./demos/gen_fig7.R: Estimating the value of a Barrier option under a jumping asset price model using simulation.
![Alt text](./results/fig4.png?raw=true "Barrier option value for three different barrier values")

8. ./demos/sim_cuda_BS.R: 
Comparison between the European call option value as computed using the BS analytic formula, and the estimated value obtained by GPU simulations.

9. ./demos/sim_cuda_JDP.R: 
GPU sampling of trajectories of a Jump decision process (Merton models), and plotting the distribution of realizations at exercise time.


