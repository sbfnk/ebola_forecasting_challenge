#!/bin/bash
#$ -N ebofit_nb
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -m n
#$ -l mem_free=1000M,h_vmem=2000M
#$ -pe smp 4
#$ -q parallel.q
#$ -R y

scenario=$1 ## to be passed to script: any of 2, 3, 4

source ~/.bashrc
source ~/perl5/perlbrew/etc/bashrc

Rscript ~/code/ebola_forecasting_challenge/fit_challenge.r -n 10000 -p 256 -t 4 -s $scenario -i 2

