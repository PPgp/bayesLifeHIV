library(bayesLifeHIV)
source('test_functions.R')

cran <- TRUE
wpp <- 2019

test.options()
if(cran)
    test.simulate5y(wpp.year = wpp)

## Time-expensive tests
if(!cran) {
	for (wpp in c(2019, 2024)) test.simulate5y(wpp.year = wpp)
    for (wpp in c(2022, 2024)) test.simulate1y(wpp.year = wpp)
}
