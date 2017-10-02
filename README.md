
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`medltmle`
============

[![Travis-CI Build Status](https://travis-ci.org/podTockom/medltmle.svg?branch=master)](https://travis-ci.org/podTockom/medltmle) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/podTockom/medltmle?branch=master&svg=true)](https://ci.appveyor.com/project/podTockom/medltmle) [![Coverage Status](https://img.shields.io/codecov/c/github/podTockom/medltmle/master.svg)](https://codecov.io/github/podTockom/medltmle?branch=master) [![CRAN](http://www.r-pkg.org/badges/version/medltmle)](http://www.r-pkg.org/pkg/medltmle) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/medltmle)](https://CRAN.R-project.org/package=medltmle) [![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Machine learning-based summary of association with multivariate outcomes

**Authors:** [Ivana Malenica](https://github.com/podTockom) and Wenjing Zheng

Description
-----------

`medltmle` estimates the natural mediation effect for a longitudinal setting with time-varying mediators. This R package implements several estimators of the data dependent parameter and non-data dependent parameter (fully conditional on the past) for direct and indirect mediation effects over multiple time points, adjusting for measured time-varying confounding and informative right-censoring. The theoretical justifications for using either of the aforementioned parameters are outlined in the vignette.

Currently available estimators include

1.  TMLE for longitudinal data
2.  The longitudinal G-computation
3.  Inverse Probability of Treatment Weighted (IPTW)

Future releases will support longitudinal data in long format and will integrate with the [`stremr` package](https://github.com/osofr/stremr) in order to handle more elaborate longitudinal data structures.

------------------------------------------------------------------------

Installation
------------

You can install the most recent *stable release* from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("podTockom/medltmle")
```

------------------------------------------------------------------------

Issues
------

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/podTockom/medltmle/issues).

------------------------------------------------------------------------

Example
-------

To get an idea of how `medltmle` works, let's try using it with a simple simulated data set:

``` r
# setup
suppressMessages(library(medltmle))
set.seed(67394)

# simulation parameters
end.time = 2
n <- 400

# simulate data
data <- GenerateData(n = n, end.time = end.time)
```

Next, we can generate simple models for conditional densities and iterative expectations, and define counterfactual exposures:

``` r
# define models
spec <- make.sim.spec(2)

# define counterfactuals
abar <- 1
abar.prime <- 0
```

Having gone through the above steps, we can now obtain IPTW and TMLE estimates of the *natural mediation effect*:

``` r
# let's fit the longitudinal TMLE for the fully conditional on the past mediation parameter:
result_10 <- suppressMessages(
              medltmle(data = data,
                       Anodes = names(data)[grep("^A", names(data))],
                       Cnodes = names(data)[grep("^C", names(data))],
                       Znodes = names(data)[grep("^Z", names(data))],
                       Lnodes = names(data)[grep("^L", names(data))],
                       Ynodes = names(data)[grep("^Y", names(data))],
                       survivalOutcome = TRUE,
                       QLform = spec$QL.c,
                       QZform = spec$QZ.c,
                       gform = spec$g.c,
                       qzform = spec$qz.c,
                       qLform = spec$qL.c,
                       abar = rep(abar, end.time),
                       abar.prime = rep(abar.prime, end.time),
                       CSE=TRUE,
                       time.end = end.time
                      )
              )
              
#>"<0x10f267b40>"
#>tracemem[0x10f267b40 -> 0x10f27c490]: MainCalcsMediation LtmleMediationMSMFromInputs ltmleMediation medltmle #>withCallingHandlers suppressMessages 
#>tracemem[0x10f27c490 -> 0x10c900350]: EstimateG MainCalcsMediation LtmleMediationMSMFromInputs ltmleMediation #>medltmle withCallingHandlers suppressMessages 
#>tracemem[0x10f27c490 -> 0x10c948570]: EstimateMultiDens MainCalcsMediation LtmleMediationMSMFromInputs #>ltmleMediation medltmle withCallingHandlers suppressMessages 
#>tracemem[0x10f27c490 -> 0x1088f6ce0]: EstimateMultiDens MainCalcsMediation LtmleMediationMSMFromInputs #>ltmleMediation medltmle withCallingHandlers suppressMessages 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  1 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  2 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  3 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  4 
#>tracemem[0x10f27c490 -> 0x101b3d780]: EstimateG MainCalcsMediation LtmleMediationMSMFromInputs ltmleMediation #>medltmle withCallingHandlers suppressMessages 
#>tracemem[0x10f27c490 -> 0x10c99f6a0]: EstimateMultiDens MainCalcsMediation LtmleMediationMSMFromInputs #>ltmleMediation medltmle withCallingHandlers suppressMessages 
#>tracemem[0x10f27c490 -> 0x1088b31b0]: EstimateMultiDens MainCalcsMediation LtmleMediationMSMFromInputs #>ltmleMediation medltmle withCallingHandlers suppressMessages 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  1 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  2 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  3 
#>Wed Jul 26 23:35:43 2017 EstimateLYnodes node  4 
#>[1] "2017-07-26 23:35:43 PDT"

# let's examine the estimates:
result_10$estimates
#>     tmle      iptw 
#>0.8903333 0.9439669
```

------------------------------------------------------------------------

License
-------

© 2017 Ivana Malenica

The contents of this repository are distributed under the MIT license. See below for details:

    The MIT License (MIT)

    Copyright (c) 2017 Ivana Malenica and Wenjing Zheng

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
