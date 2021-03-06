# ALT-FISH script collection
Scripts for automated analysis of ALT-FISH microscopy data with or without additional co-staining (e.g, cytoplasm marker, immunolabeling).
The scripts have been written to process 1024x1024px 16-bit .TIF image z-stacks with at least 2 (and maximum 4) channels as input.

# Script 1: [Telosegment toolkit](Telosegment_toolkit_multirun.R)
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
* customizable batch processing of image stacks from multiple folders (multirun mode)

### Outputs:
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
* carbon copy of the script used in this run ("_Telosegment_toolkit_script_used_in_this_run.R") ??? for traceability purposes
* optional: R workspace after run completion ("_R_workspace_after_run_completion.RData")????? useful for trouble-shooting

### Requirements
* R version 3.6.0 or higher (tested with 4.0.2)
* R packages:
  * [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html) version 4.32.0
* custom functions :
  - [findBestSlices()](functions/findBestSlices_v2.R): returns the optimal z-range based on the nuclear staining (e.g., DAPI) intensity 
  - [makeNucMask()](functions/makeNucMask_v1.R): creates nuclear masks 
  - [makeSpotMask()](functions/makeSpotMask_v1.R): creates spot masks 


# Script 2: [Telosegment-cytosegment toolkit](Telosegment_cytosegment_multirun.R)
### What does this script do?
The working principle and output format is similar to script 1 (see above). However, the functionality is limited to 3 channel images, where one channel constitutes a nuclear staining (for nucleus segmentation), another one provides a cytoplasm marker (for cell segmentation) and a third channel corresponds to the channel where spot segmentation is carried out (e.g., ALT-FISH signal). The script segments the nucleus and cytoplasm area of each cell in 2D, followed by the quantification of spot numbers and features in these areas separately. Spot detection and segmentation stringenies for nuclei and cytoplasmic areas can be adjusted individually (see comments within the script for details). Moreover, the script provides means for the removal of border cells (recommended) and segmentation artifacts (e.g.,faulty matching of cytoplasmic and nuclear masks).

### Outputs:
Additional files generated by this script:
* results files (.csv) containing quantification results for the cytoplasm area
    - "_cytoplasm_spotmask_spotfreemask_area_intensity_results" : cytoplasmic spot numbers, aggregated spot areas and intensity features sorted by position and    cell id
    - "_single_spot_shape_features_and_intensities_cytoplasm" : shape and intensity features for each segmented cytoplasmic spot sorted by position, cell and spot id

### Additional Requirements
* custom functions:
  - [makeCytoMask()](functions/makeCytoMask_v1.R): returns masks corresponding to the entire cell area using both nucleus and cytoplasm channel images
  - [matchCellMasks()](functions/matchCellMasks_v5.R): required for matching of nuclei and cytoplasm masks


# Script 3: [Telosegment toolkit suspension](Telosegment_toolkit_suspension_multirun.R)
### What does this script do?
The workflow and output format is identical to script 1. However, this script additionally performs a focus check on each single segmented nucleus using the custom function [evaluateZProfile()](functions/evaluateZProfile_v1.R). It thereby enables the later removal of cells that are only partially within the recorded z-volume. This is particularly important when analyzing imaging data recorded from, for example, stained cells in suspension. The focus check result (in focus = "y" or "n") is appended as an additional column to the "_nucleus_spotmask_npmask_area_intensity_results.csv" results table. In addition, nucleus marker intensity profiles along all z-slices are exported for each nucleus and position and can be found in the folder "_nuc_focuscheck".

# Script 4: [CollectAndFilter Telosegment results](CollectAndFilter_Telosegment_results.R)
### What does this script do?
The script imports data (.csv tables) obtained as output from Telosegment toolkit scripts (scripts 1-3) into the current R Studio environment for further processing or plotting. The imported data is formatted as lists of data frames, thereby preserving the original folder and table name structure. In addition to importing the raw data, the script offers two filtering steps (both recommended):

1) Removal of data from out-of-focus image positions: Both Telosegment and Telosegment-cytosegment toolkit scripts label image positions that do not meet the specified focus-criteria. During this step, these positions are identified and associated data is removed from the corresponding data frames. 

2) Removal of data from cells whose nuclear masks do not meet certain shape criteria. Segmentation artifacts arising from autofluorescent objects (dirt particles, locally altered fluorescence intensity background, etc.), abberantly shaped or poorly-separated nuclei (doublets, triplets) are identified based on a nuclear shape quality score (NSQS), see script comments for details. NSQS is a size-independent measure of circularity/irregularity and a default NSQS cutoff of > 3 was found to efficiently remove most of the above-mentioned artifacts (score become higher the more irregular the objects are). Although this empirically determined cutoff value is relatively stringent, we recommend to inspect the results after filtering and adjust the NSQS cutoff until the desired filtering efficiency is reached. It is helpful to retrieve NSQS values for individual nuclei/objects from the "_nucleus_shape_features" table/data frame after the first filtering step and compare them to the segmented images (mask overlays) in the Telosegment result folder in order to fine-tune the filtering stringency.

### Output
The script's output consists of R objects containing the raw data, out-of-focus filtered data and out-of-focus+nsqs filtered data. Each object is a list in which each element corresponds to one complete Telosegment toolkit run result. Folder paths are used as element names. Each list element is itself a list of data frames containing the data from all the .csv tables found in this specific results folder (see comments in script for details). 
