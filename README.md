# bayesLifeHIV

[![Travis-CI Build Status](https://travis-ci.org/PPgp/bayesLifeHIV.svg?branch=master)](https://travis-ci.org/PPgp/bayesLifeHIV)

Extension of the **bayesLife** R package that takes into account HIV/AIDS. It requires a version of **bayesLife** from the [bL4 branch](https://github.com/PPgp/bayesLife/tree/bL4).

### Installation

You can either install it the traditional way of cloning GitHub repositories and installing the packages from local directories. Or you can use the **devtools** package as shown below.

1. Since **bayesLifeHIV** is currently a private repository, create an access token from [here](https://github.com/settings/tokens):
	1. Click on "Generate New Token".
	2. Enter a description (e.g. "installing bayesLifeHIV") and check the "repo" button.
	3. Click "Generate token".
	4. Save the token somewhere.

2. Install **bayesLifeHIV** from GitHub:

	```
	install_github("PPgP/bayesLifeHIV", auth_token = "my_token")
	```
	Replace ``my_token`` above with the token you saved in the previous step. This command should install **bayesLife** from the bL4 branch.

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

Those two functions reset the global settings to the defaults of the respective package.
 
#### MCMC estimation

Here is an example of a toy simulation:

```
sim.dir <- "e0simHIV"
m <- run.e0hiv.mcmc(nr.chains = 2, iter = 50, thin = 1, 
			output.dir = sim.dir, replace.output = TRUE,
			verbose = TRUE)

```

Countries considered as epidemic in the estimation can be viewed via

```
hiv.countries.est(m$meta)
```

Summaries, DL curves etc. can be viewed as with ususal **bayesLife** objects, see ``?bayesLife``. E.g. 

```
e0.partraces.plot(m)
```

#### Projections

To generate predictions for the toy simulation above, run the following command:
 
```
pred <- e0hiv.predict(sim.dir = sim.dir, burnin = 10, 
			replace.output=TRUE)
```

Countries considered as epidemic in the prediction can be viewed via

```
hiv.countries.pred(m$meta)
```
Summaries, trajectory plots, maps etc. can be viewed as with ususal **bayesLife** objects, see ``?bayesLife``. E.g.


```
e0.map.gvis(pred)
```

#### MCMC settings

In this version of **bayesLife** (branch bL4) and thus **bayesLifeHIV**, the ``run.e0.mcmc()`` and ``run.e0hiv.mcmc()`` functions were made cleaner in a way that many of the arguments were moved into global options. They can be viewed using 

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

