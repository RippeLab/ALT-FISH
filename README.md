# ALT-FISH
Scripts for automated analysis of ALT-FISH microscopy data with or without additional co-staining (e.g, cytoplasm marker, immunolabeling).
The scripts require 16-bit .TIF image z-stacks with at least 2 (and maximum 4) channels as input.

# Script 1: Telosegment toolkit
### What does this script do?
* importing and processing of .TIF image stacks for the quantitative extraction of image features
* optimal focus volume determination ("BestSlices") to remove empty frames and identify out-of-focus positions 
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
Each run generates a folder with unique date-time stamp containing all the results
* results files (.csv) containing quantification results for each cell and/or image position plus run-specific parameters:
    - "_positions" : list with a numerical position id (unique identifier) assigned to each processed image stack 
    - "_nucleus_shape_features" : nucleus mask xy-positions and shape features sorted by position and nucleus id 
    - "_nucleus_spotmask_npmask_area_intensity_results" : spot numbers, aggregated spot areas and intensity features sorted by position and nucleus id 
    - "_parameters_used_in_run" : summary of user-adjustable parameters used in this run
    - "_position_results" : per-position image features (identified z-focus range, type of projection, background estimates)
    - "_single_spot_shape_features_and_intensities" : shape and intensity features for each segmented spot sorted by position, nucleus and spot id
* images with the created masks for visual inspection, publication or optional downstream analyses :
    - "_BestSlices_projections_with_masks_binary" : BestSlices z-projections used for segmentation with binarized masks
    - "_BestSlices_projections_with_masks_outlines" : as above, but with mask outlines and including the nucleoplasm mask
    - "_BestSlices_projections_with_masks_scaledRGB" : scaled RGB image with colored mask representation for visual inspection
* carbon copy of the script used in this run ("_Telosegment_toolkit_script_used_in_this_run.R") – for traceability purposes
* optional: R workspace after run completion ("_R_workspace_after_run_completion.RData") – useful for trouble-shooting

### Requirements
* R version 3.6.0 or higher (tested with 4.0.2)
* R packages:
  * [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html) version 4.32.0
* associated custom [functions](https://github.com/lfra/ALT-FISH/tree/main/functions) for image processing and segmentation



# Script 2: Telosegment-cytosegment toolkit
The working principle and output format is similar to script 1 (see above). However, the functionality is limited to 3 channel images, where one channel constitutes a nuclear staining (for nucleus segmentation), another one provides a cytoplasm marker (for cell segmentation) and a third channel corresponds to the channel where spot segmentation is carried out (e.g., ALT-FISH signal). The script segments the nucleus and cytoplasm area of each cell in 2D, followed by the quantification of spot numbers and features in these areas. Spot detection and segmentation stringenies for nuclei and cytoplasmic areas can be adjusted separately (see comments within the script for details). Moreover, the script provides means for automatically removing border cells (recommended) and segmentation artifacts (e.g.,faulty matching of cytoplasmic and nuclear masks).

### Additional Requirements


# Script 1: Telosegment toolkit



# Script 2: Telosegment-Cytosegment toolkit


# Filtering and import of segmentation data
