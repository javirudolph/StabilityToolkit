
<!-- README.md is generated from README.Rmd. Please edit that file -->

# StabilityToolkit

<!-- badges: start -->

<!-- badges: end -->

The goal of StabilityToolkit is to provide an array of statistical and
mathematical functions to estimate parameters of a series of dynamic
models consisting of a pair of difference equations that follow the
abundance of plasmid-free and plasmid-carrying cells under various
biological hypotheses. The package also includes functions to simulate
data under each model. Estimation is done via Maximum Likelihood.
Hierarchical stochastic models are also considered. Further information
can be found in the author
[papers](http://people.clas.ufl.edu/josemi/papers/) made in
collaboration with [Eva Top](http://www.thetoplab.org/).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("javirudolph/StabilityToolkit")
```

Once the package is installed, you can start using it by loading it
with:

``` r
library(StabilityToolkit)
```

## Preamble

If you have questions and run into problems, please e-mail me at
<josemi@ufl.edu>. Mathematics and statistics are fantastic tools to
translate fundamental questions in biology into testable predictions.
It’s the best mechanism that we have in science to formulate relevant
biological questions in a way that these can be directly confronted with
the data sets at hand. If you have multiple tentative explanations for
your results and you speak even just a little of the language of math,
you have what you need to translate those tentative explanations into
statistical hypotheses that can all be simultaneously tested against
your data. So if you want to make the most of this package, you should
familiarize yourself with the models in the [Top
lab](http://www.thetoplab.org/). A lot of time and effort went into
writing the models in a way that they are easily understood. A great
place to start is the paper by [De Gelder et
al 2004](https://www.genetics.org/content/168/3/1131) and the model
section of the paper by [Ponciano et
al 2007](https://doi.org/10.1534/genetics.106.061937). There you will
find all the models used in this package. Reading those papers is a
pre-requisite to use this package smartly, and not just as a
plug-and-play tool. Even if you don’t fully understand all of the
statistical details, you won’t regret spending time on these papers. It
is my hope that after reading these you will have what you need to
de-bug your analyses.
