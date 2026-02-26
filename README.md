# MebukiITO_Riskmap
Sample R code for Article "Spatial distribution of ticks and tick-borne pathogens in central Hokkaido, Japan and associated ecological factors revealed by intensive short-term survey in 2024"

This repository contains the representative R script for the ecological and spatial analyses conducted in the study: "Spatial distribution of ticks and tick-borne pathogens in central Hokkaido, Japan and associated ecological factors revealed by intensive short-term survey in 2024"


Overview

Due to ethical policies regarding the exact GPS coordinates of the sampling sites, the original spatial dataset cannot be publicly shared. To ensure the transparency and reproducibility of our methodology, we provide this representative R script. This script demonstrates the analytical code used in the study, allowing researchers to apply the same methods to their own ecological data by replacing the generic data object with their own data frame.


Contents

GitHub_MebukiITO_Riskmap.R: The R script containing the analytical workflow (Spatial statistics, GLM/GAM modeling, and Niche overlap calculations).


Requirements

The script was developed using R version 4.4.1.
The script was run in R. The following packages are required:

pROC, MuMIn, dismo, dplyr, ggplot2, sf, spdep, mgcv, caret


Usage

Clone or download this repository.

Open the R script sample_analysis.R.

Load your own spatial dataset and name it data. Ensure your dataset contains the necessary target binary response columns (e.g., species), coordinate columns (Longitude, Latitude), and environmental variables.

