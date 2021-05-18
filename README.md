# ALT-FISH
Scripts for automated analysis of ALT-FISH microscopy data with or without additional co-staining (e.g, cytoplasm marker, immunolabeling).
The scripts require 16-bit .TIF image z-stacks with at least 2 (and maximum 4) channels as input.

## Telosegment toolkit
### What does this script do?
* automatically reads in and processes .TIF image stacks for extraction of image features
* optimal focus volume determination to remove empty frames and identify out-of-focus positions 
* data reduction by z-projection for subsequent 2D segmentation procedures
* segmentation of cell nuclei (comprising adjustable size filter and removal of border nuclei)
* segmentation of spots within the segmented nuclear areas
* quantification of:
  * whole image parameters (optimal z-volume, background estimates) 
  * size, xy-positions and shape parameters for the generated nucleus and spot masks
  * mask fluorescence intensity parameters in up to 4 channels (mean, median, maximum, standard deviation, etc.)
  * number of spots and aggregated spot area per nucleus
* customizable batch processing of image stacks from multiple folders

### Script outputs:
* each run  generates a results folder with date-time stamp containing all the output
* .csv results files contain the above-mentioned quantification results for each cell and/or image position plus run-specific parameters
  
    "_positions" : list with a numerical position id (unique identifier) assigned to each processed image stack 
  
    "_nucleus_shape_features" : nucleus mask xy-positions and shape features sorted by position and nucleus id 
    
    "_nucleus_spotmask_npmask_area_intensity_results" : spot numbers, aggregated spot areas and intensity features sorted by position and nucleus id 
    
    "_parameters_used_in_run" : summary of user-adjustable parameters used in this run
    
    "_position_results" : per-position image features (identified z-focus range, type of projection, background estimates)
    
    "_single_spot_shape_features_and_intensities" : shape and intensity features for each segmented spot sorted by position, nucleus and spot id
  
* images with the created nuclear and spot masks for visual inspection, publication or optional 
* fit parameters (only in R command line)
* optional: an image with the created nuclear and aggregate masks for visual inspection (for each cell)




## Telosegment-cytosegment toolkit: What does the script do?
The script is limited to 3 channel images, where one channel constitutes a nuclear staining (for nucleus segmentation), another one provides a cytoplasm marker (for cell segmentation) and a third channel corresponds to the channel where spot segmentation is carried out (e.g., ALT-FISH signal).

This script extracts the same image features as the Telosegment toolkit script, but is   


# Script 1: Telosegment toolkit



# Script 2: Telosegment-Cytosegment toolkit


# Filtering and import of segmentation data
