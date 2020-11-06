# Anchovy_Spatial_Indicators
Indicators found in Moron et al (2019)

Create a folder called 'data' to save all acoustic files (e.g. 'Acoustic_Data_170304.csv'). Colnames are lon_m, lat_m, ancp.

The main code is 'Nongaussian_lognormal.R'. Run the INLA model using this code. A folder per survey will be created. Then, save those folders in a folder called 'surveys'.

Use the function 'getIndices.R' to calculate spatial indicators found in Moron et al (2019).
