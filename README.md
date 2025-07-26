# bayesLifeHIV

[![R Build Status](https://github.com/PPgp/bayesLifeHIV/workflows/R-CMD-check/badge.svg)](https://github.com/PPgp/bayesLifeHIV/actions?workflow=R-CMD-check)

Extension of the **bayesLife** R package that takes into account HIV/AIDS, as described in Godwin and Raftery (2017). It requires **bayesLife** version 4.0-2 or higher.

### Installation

You can either install it the traditional way of cloning GitHub repositories and installing the packages from local directories. Or you can use the **devtools** package as shown below.

* Install **bayesLifeHIV** from GitHub:

	```
	library(devtools)
	install_github("PPgP/bayesLifeHIV")
	```
	

### Usage

```
library(bayesLifeHIV)
```
Since **bayesLifeHIV** shares (and overwrites) some of the global objects of **bayesLife**, it is recommended not to mix simulations from **bayesLifeHIV** and **bayesLife**  in the same session. Or, if you work on objects created with **bayesLifeHIV**, start with 

```
using.bayesLifeHIV()
```

When you switch to a **bayesLife** simulation, do 

```
using.bayesLife()
```

The two functions above reset the global settings to the defaults of the respective package.
 
#### MCMC estimation

Here is an example of a toy simulation:

```
sim.dir <- "e0simHIV"
m <- run.e0hiv.mcmc(nr.chains = 2, iter = 50, thin = 1, 
			output.dir = sim.dir, replace.output = TRUE,
			verbose = TRUE)

```
The function above estimates MCMCs for ALL countries, non-epidemic as  well as epidemic. Argument ``annual`` can be set to ``TRUE`` if this should be an estimation on annual data (available for ``wpp.year`` >= 2022).

Countries considered as epidemic in the estimation can be viewed via

```
hiv.countries.est(m$meta)
```

Summaries, DL curves etc. can be viewed as with usual **bayesLife** objects, see ``?bayesLife``. E.g. 

```
e0.partraces.plot(m)
```

#### Projections

To generate predictions for all countries for the toy simulation above, run the following command:
 
```
pred <- e0hiv.predict(sim.dir = sim.dir, burnin = 10, 
			replace.output=TRUE)
```
Again, this generates predictions for non-epidemic as well as epidemic countries. 
Countries considered as epidemic in the prediction can be viewed via

```
hiv.countries.pred(m$meta)
```

Summaries, trajectory plots, maps etc. can be viewed as with ususal **bayesLife** objects, see ``?bayesLife``. E.g.


```
e0.map.gvis(pred)
```

#### MCMC settings

In this version of **bayesLife** (version 4.0 and higher) and thus **bayesLifeHIV**, the ``run.e0.mcmc()`` and ``run.e0hiv.mcmc()`` functions were made cleaner in a way that many of the arguments were moved into global options. They can be viewed using 

```
e0mcmc.options()
```

These options are reset by ``using.bayesLifeHIV()`` (HIV setup) and ``using.bayesLife()`` (normal setup). 

To change some of these values, one can use the  ``e0mcmc.options()`` command prior to ``run.e0.mcmc()``. E.g. to change the upper limit of the ``Triangle``'s prior, do 

```
e0mcmc.options(Triangle = list(
					prior.up = c(110, 110, 110, 110))
				)
```

The run function stores the current options into the mcmc object (``meta$mcmcm.options``). Thus, if we estimate with the above settings, we can see that the resulting object contains the changes:

```
sim.dir <- "e0simHIVopttest1"
m1 <- run.e0hiv.mcmc(nr.chains = 1, iter = 10, thin = 1, 
              output.dir = sim.dir, replace.output = TRUE,
              verbose = TRUE)
m1$meta$mcmc.options$Triangle$prior.up
# compare to the original mcmc object 
m$meta$mcmc.options$Triangle$prior.up
```

Another way of changing the setting is to pass it directly into the run function, namely in the ``mcmc.options`` argument. Let's now change the lower bound of ``Triangle``:

```
sim.dir <- "e0simHIVopttest2"
m2 <- run.e0hiv.mcmc(nr.chains = 1, iter = 10, thin = 1, 
              output.dir = sim.dir, replace.output = TRUE,
              mcmc.options = list(Triangle = list(
                         prior.low = c(5, 5, 0, 5))),
              verbose = TRUE)
m2$meta$mcmc.options$Triangle$prior.low
``` 

The difference is that this time the global options were not changed, see ``e0mcmc.options("Triangle")``.

### Required Datasets

The package requires three additional datasets that are beyond of what is needed in **bayesLife**. Currently, there are two versions of these three datasets included in the package, one for a 5-year simulation (``annual = FALSE``, default) and one for an annual simulation (``annual = TRUE``). The latter (have "1y" at the end of their respective names) are newer - used for WPP 2024. Which version is used by default
depends on the value of the ``annual`` argument. 

When creating those datasets, they should all be tab delimited. 
 
#### Estimation

  * **HIVprevalence**, **HIVprevalence1y**: contain historical and projected HIV prevalence for countries with past or/and future epidemics. Rows correspond to countries; columns correspond to time periods, spanning from 1970 to 2100
  (note that columns for future years are only used for scaling, see below). It can also have a column "include\_code". Only countries are considered as epidemic in the estimation, if their "include\_code" is 1. View the datasets via
  
    ```
    data(HIVprevalence, HIVprevalence1y)
    head(HIVprevalence)   # 5-year data
    head(HIVprevalence1y) # annual data
    ```
    The default dataset can be overwritten by passing the name of the new file as the argument ``my.hiv.file`` in the ``run.e0hiv.mcmc()`` function.
    
  * **ARTcoverage**, **ARTcoverage1y**: contain historical and projected values of ART coverage for countries with past or/and future epidemics. It has the same structure as HIVprevalence(1y), including the optional "include\_code" column. If the column is present, countries considered as epidemic must have their "include\_code" set to 1 in both datasets, HIVprevalence(1y) and ARTcoverage(1y). If the "include\_code" column is missing, all countries included in these files are considered as epidemic. Note that in the estimation, only columns for historical time periods are used, while future time periods of ARTcoverage(1y) are used in the projections. View the dataset via
  
    ```
    data(ARTcoverage, ARTcoverage1y)
    head(ARTcoverage)    # 5-year data
    head(ARTcoverage1y)  # annual data
    ```
    The default dataset can be overwritten by passing the name of the new file as the argument ``my.art.file`` in the ``run.e0hiv.mcmc()`` function.

During the estimation, the HIV prevalence and ART coverage datasets are merged and converted into a delta(nonART) dataset, which is delta(hiv * (100 - art)/100), containing all countries. It can be accessed from an mcmc object, e.g. for an object ``m``:

```
m$meta$dlt.nart[, 1:10]
```

#### Projections

  * **ARTcoverage**, **ARTcoverage1y**: as mentioned above, columns that correspond to future years are used in the projections. The dataset is again filtered using the "include\_code" column, if available. The default ARTcoverage(1y) can be overwritten by passing the name of the new file as the argument ``my.art.file`` in the ``e0hiv.predict()`` function.

  * **HIVprevTrajectories**, **HIVprevTrajectories1y**: these are probabilistic trajectories of HIV prevalence for countries considered as epidemic in the future. It should have columns "country\_code", "Trajectory", and time periods "2010-2015", ..., "2095-2100" (for the 5-year version) or "2010", "2011" ,..., "2100" (for the annual version). The default dataset contains 1000 trajectories for the countries that are epidemic (scaled to the HIVprevalence(1y) dataset, see below). View the dataset by 
  
    ```
    data(HIVprevTrajectories, HIVprevTrajectories1y)
    head(HIVprevTrajectories)
    head(HIVprevTrajectories1y)
    ``` 
    The default dataset can be overwritten by passing the name of the new file as the argument ``my.hivtraj.file`` in the ``e0hiv.predict()`` function. It has to have the same format as the default dataset. If scaling to the HIVprevalence(1y) dataset is desired on the user-defined data, set the argument ``scale.hivtraj`` to ``TRUE``.
    
    
Similarly to the estimation, the HIV and ART datasets are merged and converted into a dataset of delta(nonART) trajectories, which are then used in the prediction step.

By default, countries considered as epidemic in the projection are those that have "include\_code" set to 3 in the ``include_{wpp}`` dataset in the **bayesLife** package. This can be overwritten by the argument ``hiv.countries`` in ``e0hiv.predict()``, which should be a vector of country codes. However, in both cases, such countries are required to have records in both, the ARTcoverage(1y) and HIVprevTrajectories(1y) datasets.

To view countries considered as epidemic in the projection in WPP 2024, one can do

```
data(include_2024, package = "bayesLife")
data(ARTcoverage1y, package = "bayesLifeHIV")
merge(subset(include_2024, include_code == 3), 
      ARTcoverage1y[, c("country_code", "name")])
```

##### Trajectories Scaling 

We have scaled our original HIV trajectories so that the median for each country aligns with the values in the **HIVprevalence** dataset (columns corresponding to future time periods). The **HIVprevTrajectories** is a result of such scaling. 

The scaling (which uses adjusted logit) is implemented in the function ``e0hiv.scale.trajectories()``. One can pass a dataset of trajectories (in the same format as HIVprevTrajectories) and a dataset of one time series per country (in the same format as HIVprevalence). Alternatively, it can be performed on the fly within the prediction function ``e0hiv.predict`` by setting the argument ``scale.hivtraj`` to ``TRUE``. Argument ``scale.hivtraj.tofile`` in ``e0hiv.predict`` can be used to pass the file name containing the dataset to be scaled to. The HIVprevalence dataset is used as the default.


### References

[Godwin, J. and Raftery, A.E. (2017): Bayesian projection of life expectancy accounting for the HIV/AIDS epidemic. Demographic Research 37(48):1549–1610.](http://www.demographic-research.org/Volumes/Vol37/48)
