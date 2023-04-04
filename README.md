# A composite zero-dose vulnerability index for equity assessment and spatial prioritization in low and middle-income countries
<sub>Last edited: 04/04/2023</sub>

## Data availability

## R Scripts

`covrank.R`

This R script contains the code to the `covrank()` function. This function is necessary in order to run the code within the `index.R` script. The `covrank()` function ranks the risk factors using their predictive $R^2$ statistics by fitting simple binomial regression modelled to cluster level data using the coverage indicators.

`index.R`

This R script contains the code used to calculate the vulnerability indices proposed in this research. The code in this R script is dependent on the high-resolution outputs produced from the `indicators.R` script. The code in this R script is also dependent on the `covrank()` function within the `covrank.R` script.

`indicators.R`

This R script contains the code used to construct point-referenced spatial binomial models for the zero-dose indicators extracted from the Demographic Health Survey (DHS) datasets. The constructed models are then fitted in the Bayesian context with INLA. Finally, the model is used to produce high-resolution 1x1 km prediction and uncertainty surfaces.
