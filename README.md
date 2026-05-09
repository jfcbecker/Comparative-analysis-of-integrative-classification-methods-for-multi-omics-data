# benchmark-integrative-methods

**Introduction**

This repository contains scripts and methods for the comparative analysis of integrative classification methods for multi-omics data. It provides a thorough comparison of six methods, representative of the main families of intermediate integrative approaches for supervised multi-omics classification, evaluated across 15 simulation scenarios and three real-world datasets spanning infectious diseases, oncology, and vaccine research.

**Reference**

Novoloaca A, Broc C, Beloeil L, Yu W-H, Becker J. *Comparative analysis of integrative classification methods for multi-omics data*. Briefings in Bioinformatics, 2024, 25(4):bbae331. [https://doi.org/10.1093/bib/bbae331](https://doi.org/10.1093/bib/bbae331)

A `CITATION.cff` file is provided at the root of the repository for citation managers.

**Keywords**

benchmark, data integration, multi-omics data, prediction models, supervised analysis

**Methods Compared**

The benchmark covers integrative classification methods including DIABLO, SIDA, PIMKL, netDx, block forest, and stacked generalization, with random forest used as a non-integrative baseline.

**Usage**

Scripts `make_example.R` and `run_simulation.R` are meant to be used at the beginning of the analysis workflow to simulate data.

**Installation of Methods**

After running scripts `make_example.R` and `run_simulation.R`, all the methods should be installed from CRAN/Bioconductor/GitHub as suggested by the authors of each method.

**Main Scripts**

Each method has a corresponding main script, with the method's name specified in the filename. These scripts are used to run the method and compute the performances.

**Additional Files**

In some cases, there are additional files accompanying the main scripts to aid in the analysis process. There are prior R code and posterior R code available for the PIMKL method, which should be executed as part of the analysis workflow.
