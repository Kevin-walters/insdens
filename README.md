# insdens
This package uses MCMC to calculate the posterior probability, for each bacterial gene in the data supplied, that it is essential using the insertion density derived from TraDIS data.
For full details on the format of input files, how to run the MCMC and what to do with the output, please read the vignette contained within the package.

To install the package first install the devtools package:
`install.packages("devtools")`
Then install the insdens package
`devtools::install_github("Kevin-walters/insdens", build_vignettes = T)`
Make sure you include the `build_vignettes = T` argument otherwise there will be no vignette in your package.
