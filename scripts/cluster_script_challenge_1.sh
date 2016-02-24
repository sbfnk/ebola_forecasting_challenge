#!/bin/bash
#$ -N ebofit_nb
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -m n
#$ -l mem_free=1000M,h_vmem=2000M
#$ -pe smp 4
#$ -q parallel.q
#$ -R y
#$ -t 1-6

region=${SGE_TASK_ID}
scenario=1

source ~/.bashrc
source ~/perl5/perlbrew/etc/bashrc

Rscript ~/code/ebola/libbi/fit_challenge.r -r $region -n 10000 -p 256 -t 4 -s $scenario

