#### Running order of code for generating figures and tables in the paper: "Frequency Domain Resampling for Gridded Spatial Data" by Bera, Nordman and Bandyopadhyay. 


### Basic setup is for HPC

## Covariance parameter and variogram model
# Use the following execution order for all code except the isotropy code:
# precompute.sh >> subsample.sh
#
# Process parameters, lag settings, and bandwidths are specified in the .R file.
# The block length can be adjusted in the "subsample.sh" file.

## Isotropy code
# Run the .slurm file directly.
# Parameters can be modified in the corresponding .R file.
