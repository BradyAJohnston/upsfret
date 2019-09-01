# Analysis of upsFRET images
Repository for attempting to analyse upsFRET images. Images are 3-channel images of a a 96 well plate. Image channels are Transfer, Cy3 and Cy5. 
### Steps in analysis.
1. Import channeled .tif images
2. Combine to single, 3-channel image.
3. Do spot analysis on the 96 wells.
4. Background subtract and remove noisey pixels that are contaminants.
5. Compute mean, median etc of the plate. 

Export data into table format for plotting and analysis in R or other package.

Trying to get a jupyter notebook running from github to do analysis using `Juno` on the iPad.
