# JAGS-examples
Some example Bayesian models and simulations in the JAGS language

Step 1: First download latest version of JAGS from http://mcmc-jags.sourceforge.net/

Step 2: Open R or RStudio and install the relevant MCMC and JAGS packages: coda, stats4, runjags, and segemented (for simulation 2).

Other useful packages for this include ggplot2, MASS, and rjags.

Step 3: read the introductory powerpoint file that loosely goes over Bayesian, MCMC, and applications

Step 4: Run the various .R simulation files

  - Simulation 1.R goes over simple linear regression

  - Simulation 2.R goes over a breakpoint analysis

  - Simulation 3.R goes over a hierarchical model

Step 5: have fun. 

*Note* - Many people will code the JAGS model template in R or other text editor and then save that template as a .JAGS file. The user then uploads that file and calls that model in another R script function that uses the 'rjags' or 'runjags' functions to run the model.

*Note* - The JAGS user manual has a list of other syntax options, the various distributions included and their parameterizations.
