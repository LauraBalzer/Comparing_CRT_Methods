# Comparing_CRT_Methods

This GitHub project contains the R code to implement the simulation studies and real data analysis in "Defining and Estimating Effects in Cluster Randomized Trials: A Methods Comparison" by Alejandra Benitez, Maya L. Petersen, Mark J. van der Laan, Nicole Santos, Elizabeth Butrick, Dilys Walker, Rakesh Ghosh, Phelgona Otieno, Peter Waiswa, Laura B. Balzer.

The corresponding technical report is available at https://arxiv.org/abs/2110.09633.

Questions and comments can be addressed to Laura at lbalzer@umass.edu.

#-----#-----#-----

sim1.R: code for Simulation Study I to compare common CRT estimators for their natural target of inference with an effect and under the null

sim1_supp.R: code for Supplementary Approaches in Simulation Study I to compare TMLEs for estimation of the cluster- and indv-level when breaking and keeping the matches

sim2.R: code for Simulation Study 2 to compare TMLEs for estimation of the cluster- and indv-level effects with highly informative cluster size

gen_data.R: code to generate data for the simulation studies

ptbi_analysis.R: run the analysis for the real-data example (PTBi study) 

MainFunctions.R: wrapper functions to do estimation and inference 

Stage2_Functions.R: code to use TMLE with and without Adaptive Prespecification to estimate intervention effects; Modified from https://github.com/LauraBalzer/SEARCH_Analysis_Adults

Adapt_Functions.R: code to implement Adaptive Prespecification; Modified from https://github.com/LauraBalzer/SEARCH_Analysis_Adults

Estimators.R: code to implement a t-test, CARE, GEE, and Aug-GEE. 

