# A framework for simulating genotype by environment interaction using multiplicative models

## Introduction

This repository contains R scripts that demonstrate simulation of genotye by environemnt interaction. This collection of R scripts serves as supplementary material for the manuscript:

Bancic et al. (2024) bioRxiv TODO

    @article{Bancic2024,
      title = {A framework for simulating genotype by environment interaction using multiplicative models},
      author = {Bančič, Jon and Gorjanc, Gregor and Daniel, Tolhurst},
      journal = {bioRxiv},
      year = {2024},
      doi = {TODO}
    }

The manuscript is currently under revision in Applied and Theoretical Genetics, so changes are expected. 

## Repository contents

  * `README.md` is this file.

  * `jbancic_GEIsim.Rproj` is the RStudio/Posit project file (to set the working directory etc.).

  * `00SimulateFunctions` script contains functions used by scripts `01-03`.

  * `01SimulateGxEinteraction` script demonstrates simulation of between-environment genetic correlation matrices.

  * `02SimulateMET_ToyExampleAppendix` script demonstrates simulation of phenotypic data to construct a MET dataset as shown in Appendix 1 of the manuscript.

  * `03SimulateMET` script demonstrates simulation of phenotypic data to construct a MET dataset using independent genotype slopes or correlated genotype slopes with AlphaSimR. This script uses FieldSimR to simulate plot errors and agricolae to simulate a completely randomized experimental design.
  
  * `04ModelComparison` folder contains scripts to run a single replicate of statistical model comparison. The scripts are:
  
    * `00SimulateTPEsamples` script creates `samplesModerateGEI.RData` file which is used by `04ModelComparison.R`.
  
    * `samplesModerateGEI.RData` RData file contains pre-sampled true between-environment genetic correlation and variance matrices of dimensions 1000x1000, true genetic effects and variance parameters, a simulated MET dataset with 1000 environments, and 1000 samples for constructing MET datasets of sizes 5, 10, 20, and 50.
    
    * `04ModelComparison.R` script runs a single replicate of statistical model comparison with a MET dataset of size 10 environments. Models include phenotypic averages (Pheno), main effect only (Main), compound symmetric (Comp), main plus diagonal (MDiag), diagonal only (Diag) and factor analytic 1-4 (FA1-4) models.
  
    * `Results_ModerateGEInoSpatial10` summary file of results from `04ModelComparison.R`.
  
  * `05BreedingSimulation` folder contains scripts to run a single replicate of the phenotypic and genomic selection line breeding program simulation without and with moderate level of genotype by environment interaction. The files are:
  
    * `00OPrepareTPE.R` script prepares TPE genetic variance matrix and samples 600 environments to be used during a 30-year simulation period.
  
    * `01CreateHaplotypes.R` script creates founder haplotypes to be used by all `RUNME` scripts.
    
    * `02-05RUNME` scripts run a single replicate of phenotypic (`_pheno`) and genomic selection (`_GS`) line breeding programs without (`_noGEI`) or with moderate GEI (`_withGEI`) for 30 years.
    
    * `PlotResults.R` takes outputs from simulations and plots them.
    
    * Remaining scripts are necessary to run breeding programs and are called by `RUNME` scripts. 
    

## How to work with the provided R scripts

  1) Download this whole repository
  
  2) In RStudio/Posit, open project file `jbancic_GEIsim.Rproj`
  
  3) Actions:
      i) Run scripts `01-03` independently to learn about simulating genotype by environment interaction and a MET dataset. 
      ii) To run a model comparison example, enter `04ModelComparison` folder. First run script `00SimulateTPEsamples` followed by `04ModelComparison.R`.
      iii) To run a breeding simulation example, enter `05BreedingSimulation` folder. First run script `00OPrepareTPE` followed by `01CreateHaplotypes` and then different            scenarios using `RUNME` scripts.
