`medltmle`
================

## Description

`medltmle` estimates natural and stochastic mediation effect for a longitudinal setting with time-varying mediators. The package implements several estimators of the direct and indirect mediation effect over multiple time-points, adjusting for measured time-varying confounding and informative right-censoring. Currently available estimators include TMLE for longitudinal data, IPW, and longitudinal G-formula.

---

## Installation

Install the most recent _stable release_ from GitHub:
  ```
  devtools::install_github("podTockom/medltmle")
  ```
  
---

## Example with Simulated Data
  ```
  require(medltmle)
  set.seed(2)
  
  end.time=2
  ```

Generate some simulated data. 
  ```
   data<-GenerateData(n=400, end.time=end.time, abar=NULL,abar.prime=NULL)
  ```

Generate simple models for conditional densities and iterative expectations:
  ```
   spec<-make.sim.spec(2)
  ```

Define counterfactual exposures.
  ```
  abar <- 1
  abar.prime <- 0
  ```

IPTW and TMLE estimate of the natural mediation effect.
  ```
  result.c <- medltmle(data=data,
                           Anodes=names(data)[grep('^A',names(data))],
                           Cnodes=names(data)[grep('^C',names(data))],
                           Znodes=names(data)[grep('^Z',names(data))],
                           Lnodes=names(data)[grep('^L',names(data))],
                           Ynodes=names(data)[grep('^Y',names(data))],
                           survivalOutcome = T,
                           QLform=spec$QL.c,
                           QZform=spec$QZ.c,
                           gform=spec$g.c,
                           qzform=spec$qz.c,
                           qLform=spec$qL.c,
                           abar=rep(abar,end.time),
                           abar.prime=rep(abar.prime,end.time),
                           estimand="NE",
                           time.end=end.time
                           )

  ```
---

## License
&copy; 2017-2018 Ivana Malenica & Wenjing Zheng

The contents of this repository are distributed under the MIT license. See file LICENSE for details.
