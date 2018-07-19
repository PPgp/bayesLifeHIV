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
#### MCMC estimation

