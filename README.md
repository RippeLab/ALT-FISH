# ALT-FISH
Collection of scripts for automated analysis of ALT-FISH stained microscopy data without or with cytoplasm staining.
The scripts require 16-bit .TIF image z-stacks as input, with at least 2 channels (one for nucleus segmentation, one for spot segmentation), but the analysis is extensible for up to 4-channel images.

# Script 1: Telosegment toolkit
* simple segmentation of nuclei
* segmentation of spots within the segmented nuclei
* quantification of:
  * whole image  
  * size, xy-position and shape parameters for the generated nuclei and spot masks
  * mask fluorescence intensity parameters in up to 4 channels (mean, median, maximum, standard deviation, etc.)





# Script 2: Telosegment-Cytosegment toolkit


# Filtering and import of segmentation data
