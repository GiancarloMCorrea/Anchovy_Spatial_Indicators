# Anchovy_Spatial_Indicators
Indicators found in:

> Moron, G., Galloso, P., Gutierrez, D., Torrejon-Magallanes, J., 2019. Temporal changes in mesoscale aggregations and spatial distribution scenarios of the Peruvian anchovy (Engraulis ringens). Deep Sea Research Part II: Topical Studies in Oceanography 159, 75–83. [https://doi.org/10.1016/j.dsr2.2018.11.009](https://doi.org/10.1016/j.dsr2.2018.11.009)

Create a folder called 'data' to save all acoustic files (e.g. 'Acoustic_Data_170304.csv'). Colnames are lon_m, lat_m, ancp.

The main code is 'Nongaussian_lognormal.R'. Run the INLA model using this code. A folder per survey will be created. Then, save those folders in a folder called 'surveys'.

Use the function 'getIndices.R' to calculate spatial indicators found in Moron et al (2019).
