This repository contains the code for the LSHTM contribution to the [Ebola forecasting challenge](http://www.ebola-challenge.org).

It contains the following files:

* ```disease_params_challenge.r``` reads the line lists and transmission trees to estimate parameters and save them in a file.
* ```ebola_fit_challenge.bi``` is the transmission model, written in [libbi](http://libbi.org).
* ```fit_challenge.r``` fits the transmission data to the data.
* ```analyse_challenge.r``` analyses the results of the fits and produces CSV files with the projections and parameter estimates.
* ```script/*``` are scripts for running things in parallel on a HPC cluster.
