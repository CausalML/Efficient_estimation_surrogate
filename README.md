This repo contains code for the paper [On the Role of Surrogates in the Efficient Estimation of
Treatment Effects with Limited Outcome Data](https://arxiv.org/abs/2003.12408).

There are two folders in this repo: california-gain and offset-regression, corresponding to the experiments in Sections 5.1 and 5.2 of the paper, respectively. In each folder, the experiment.R file implements the main function used in the experiments, the utility.R file contains all auxiliary functions, and the analysis. R file analyzes the results and generates the figures and tables in the paper. The .rds files contain all experiment results (e.g., estimation errors, coverage) used to generate the figures and tables.  
