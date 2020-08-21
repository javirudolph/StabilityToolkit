
<!-- README.md is generated from README.Rmd. Please edit that file -->

# StabilityToolkit

<!-- badges: start -->

<!-- badges: end -->

The goal of StabilityToolkit is to …

This package contains functions to estimate plasmid cost, segregation
rate and horizontal transfer rates from stability assays data. The
stability asays follow periodically (usually every 24 hours) the changes
in the fraction of plasmid free cells in serial passage cultures. At the
end of each cycle, the experimenter takes a sample of the liquid culture
and plates it after dilution. Using an antibiotic marker present in the
plasmid, the experimenter can count, out of a total number of colonies
plated, how many of these are plasmid-carrying. Here, we provide an
array of statistical and mathematical functions to estimate the
parameters of a series of dynamic models consisting of a pair of
difference equations that follow the abundance of plasmid-free and
plasmid-carrying cells under various biological hypotheses. The package
also includes functions to simulate data under each model. Estimation is
done via Maximum Likelihood. Hierarchical stochastic models are also
considered. Further information can be found in the author papers made
in collaboration with Eva Top, see:
<http://people.clas.ufl.edu/josemi/papers/>

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("javirudolph/StabilityToolkit")
```
