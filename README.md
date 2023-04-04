# A composite zero-dose vulnerability index for equity assessment and spatial prioritization in low and middle-income countries
<sub>Last edited: 04/04/2023</sub>

## Data availability

## R Scripts

`covrank.R`

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc sit amet ligula euismod, sollicitudin mauris vitae, auctor lectus. Aliquam gravida justo vitae erat ullamcorper, quis consectetur erat tempor. Cras sit amet ligula ac massa dictum interdum ut ac mi. Phasellus eu velit maximus, pretium mi at, auctor felis. Donec rutrum volutpat lobortis. Cras ultricies, tellus vitae cursus elementum, elit metus faucibus nisl, eu congue felis turpis sed metus. Etiam tortor libero, viverra at leo eget, condimentum suscipit ligula. Integer tincidunt mollis hendrerit. Curabitur dignissim dictum pulvinar. Fusce eget ligula id ante auctor luctus. Morbi ac orci in elit tincidunt auctor et in orci. Morbi ut leo at mauris dignissim sodales. Maecenas eget quam eu lectus accumsan malesuada vel ac nisl. Etiam eu posuere urna, vitae mollis arcu. Fusce sit amet libero at nisi porttitor faucibus.


`index.R`

This R script contains the code used to calculate the vulnerability indices proposed in this research. The code in this R script is dependent on the high-resolution outputs produced from the `indicators.R` script. The code in this R script is also dependent on the function within the `covrank.R` script.

`indicators.R`

This R script contains the code used to construct point-referenced spatial binomial models for the zero-dose indicators extracted from the Demographic Health Survey (DHS) datasets. The constructed models are then fitted in the Bayesian context with INLA. Finally, the model is used to produce high-resolution 1x1 km prediction and uncertainty surfaces.
