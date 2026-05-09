# benchmark-integrative-methods

**Introduction**

This repository contains scripts and methods for the comparative analysis of integrative classification methods for multi-omics data.

**Usage**

Scripts make_example.R & run_simulation.R are meant to be used at the beginning of the analysis workflow to simulate data.

**Installation of Methods**

After running scripts make_example.R & run_simulation.R, all the methods should be installed from Cran/Bioconductor/github as suggested by the authors of each method. 

**Main Scripts**

Each method has a corresponding main script, with the method's name specified in the filename. These scripts are used to run the method and compute the performances.

**Additional Files**

In some cases, there are additional files accompanying the main scripts to aid in the analysis process. There are prior R code and posterior R code available for the PIMKL method, which should be executed as part of the analysis workflow.
