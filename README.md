# A framework for simulating genotype by environment interaction using multiplicative models

## Introduction

This repository contains supplementary R scripts that demonstrate the simulation of genotype by environment interaction using the R package FieldSimR as a part of a manuscript "A framework for simulating genotype by environment interaction using multiplicative models". This manuscript is currently under revision in Applied and Theoretical Genetics, so changes are expected. 

Citation:

@article{Bancic2024,
    title = {A framework for simulating genotype by environment interaction using multiplicative models},
    author = {Bančič, Jon and Gorjanc, Gregor and Tolhurst, Daniel},
    journal = {Research Square},
    year = {2024},
    note = {Preprint version 1, 28 January 2024},
    doi = {10.21203/rs.3.rs-3855188/v1},
    url = {https://doi.org/10.21203/rs.3.rs-3855188/v1}
}


## Repository contents

  * `README.md` is this file.

  * `jbancic_GEIsim.Rproj` is the RStudio/Posit project file (to set the working directory etc.).
  
  * `S0_WorkingExample` demonstrates how to simulate a small example TPE with 100 environments from which a MET dataset is sampled with 10 environments using base R functions. This script is equivalent to Appendix C in the manuscript.

  * `S1_SimulateGe` demonstrates how to simulate a reduced rank between-environment genetic variance matrix, Ge, that captures genotype by environment interaction.

  * `S2_SampleMET` demonstrates how to simulate a MET dataset that captures genotype by environment interaction and spatial variation within environments.

  * `S3_ModelComparison` demonstrates model comparison using ASReml-R with a simulated MET dataset generated in `S2_SampleMET`.

  * `S4_BreedingProgrammeComparison` demonstrates how to simulate a phenotypic line breeding program in the presence of genotype by environment interaction.
