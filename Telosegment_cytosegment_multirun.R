############################################ Segmentation of spots in nucleus and cytoplasm ############################################################
# Thresholding-based segmentation of cell nuclei and spots within nuclei 
# Thresholding-based segmentation of cell cytoplasm and spots within cytoplasm 
# Quantification of number, intensity, area and shape features of cell masks and single spots
# This script requires at least 3-channel image stacks as input (nucleus marker, cytoplasm marker, spots)
# MULTIRUN MODE: This allows to loop over several specified .tif data folders sequentially using the same script version

# load libraries and functions
library(EBImage)
library(abind)
source("makeNucMask_v1.R")
source("makeSpotMask_v1.R")
source("findBestSlices_v2.R")
source("makeCytoMask_v1.R")
source("matchCellMasks_v5.R")


####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################
# specify a folder list as character vector (vector of file paths)
# target folders should contain the .tif files to be analyzed

folderList <- c("/home/user/folder/tif_stack_dataset_1/",
                "/home/user/folder/tif_stack_dataset_2/")

###### loop through folder list ######
for(f in 1:length(folderList)){
  
# specify location of the currently opened R script version
# this is required for saving a .txt copy of this script in the  results folder for documentation
scriptversion <- list.files(getwd(),pattern = "Telosegment_cytosegment_multirun.R")

# specify location of the .tif input files 
# .tif files are required to be multi-channel z-stacks (2-4 channels)
folder <- folderList[f]

# specify channel names, order & assign to desired operation
# set non-existing channels to 0 if less than 4 channels are needed
# non-existing channels will be replaced by placeholders (pixel value 0) at the level of projections to simplify downstream processing
# fill in the channel definition in the order of the image stack, avoid putting empty channels "in between"
# e.g if image stack is DAPI/telo/marker1 use channels <- c("dapi","telo","marker1",0) and not c("dapi","telo",0,"marker1") 
channels <- c("telo","phalloidin","dapi",0)

# define channel used for nucleus segmentation
nuc_ch <- "dapi"
# define channel used for spot segmentation within the nuclei
spot_seg_ch <- "telo"
# define marker for cytoplasm segmentation -only use ch1 for that! -, if not used set to 0
ch1 <- "phalloidin"
# define additional marker for quantification, if not used set to 0
ch2 <- 0

#number of z-slices per channel
zperch <- 51

# remove nucleus or cytoplasm masks touching the image borders? "y" = yes, "n" = no
# default setting should always be "y"
remove_bordercells <- "y"

# adjustable adaptive thresholding offset for nucleus segmentation, default = 0.0008, see EBImage function thresh() for details
nuc_threshold_offset <- 0.0008

# parameter for fine-tuning the two-pass cytoplasm segmentation done by makeCytoMask (details in the function)
# default value is 2
cyto_threshold <- 2

# thresholding parameter for spot segmentation relative to the nuclear median in the spot channel
# default value is 2.5 meaning that pixel groups with intensities above 2.5*median(nuclear intensity in the spot channel) are handed over to the makeSpotMask() function
spot_cutoff <- 2.5 
# same parameter for spot segmentation in the cytoplasm area
# default value is 2
spot_cutoff_cytoplasm <- 2

# size-filter for nucleus segmentation (pixels)
# only nucleus masks with areas larger than this value will be returned by makeNucMask()
# default is 4000
nucSize_cutoff <- 4000 

# site-filter for cytoplasm segmentation (pixels)
# only cytoplasm masks with areas larger than this value will be returned by makeCytoMask()
# default is 4000
cytoSize_cutoff <- 35000

# maximum number of nuclei expected per image
# required for variable initialization
# if mn < actual number of segmented nuclei per image problems can arise when saving results in the data tables
# default is 30
mn <- 30 

# passfilter for detecting positions without nuclei
# is handed over to findBestSlices (v2)
# position is considered empty if the maximum of the mean distribution is not at least x-fold the minimum of the mean distribution
# default is 1.2 which also works for samples where the DAPI staining was weak
passfilter <- 1.2

####################################################################################################################
############################################### FIXED ##############################################################
####################################################################################################################

# create timestamp variable for this run
# will be attached to results folder and all result files for the current run
date_time <- gsub(" ","_",gsub(":","-",substr(Sys.time(),1,19)))

# create results folders
dir.create(paste0(folder,date_time,"_Telosegment_cytosegment_toolkit_results"))
results_folder <- paste0(folder,date_time,"_Telosegment_cytosegment_toolkit_results")
dir.create(paste0(results_folder,"/",date_time,"_nuc_shape_features"))
dir.create(paste0(results_folder,"/",date_time,"_spot_shape_features"))
dir.create(paste0(results_folder,"/",date_time,"_spot_shape_features_cytoplasm"))
  
# document adjustable segmentation parameters used in this run
parameters_df <- data.frame(channels[1], channels[2], channels[3], channels[4], nuc_ch, spot_seg_ch, ch1, ch2, remove_bordercells, nuc_threshold_offset,cyto_threshold, spot_cutoff,spot_cutoff_cytoplasm, nucSize_cutoff,cytoSize_cutoff, mn, passfilter)
colnames(parameters_df) <- c("channel_1","channel_2","channel_3","channel_4","nuc_ch","spot_seg_ch", "ch1", "ch2", "remove_bordercells", "nuc_threshold_offset","cyto_threshold", "spot_cutoff","spot_cutoff_cytoplasm","nucSize_cutoff","cytoSize_cutoff","mn","findBestSlices_passfilter")
write.table(file=paste0(results_folder,"/",date_time,"_parameters_used_in_run.csv"), parameters_df, sep = "\t", append = F, row.names = F, quote = F)

# generate a position list to document the files that were processed in this run
positions <- list.files(path=folder, pattern="tif")

# write position list
write.table(file=paste0(results_folder,"/",date_time,"_positions.csv"), cbind(seq(1,length(positions),by=1) ,positions), sep = "\t", append = F, row.names = F, quote = F)

# export the currently run script as text file for logging purposes 
write(readLines(scriptversion),file=paste0(results_folder,"/",date_time,"_Telosegment_cytosegment_toolkit_script_used_in_this_run.R"))

##################################################################################################################################
##################################################################################################################################
#                                                     INITIALIZE VARIABLES
##################################################################################################################################
##################################################################################################################################

### nuclei masks and segmented spot area within ####
# background intensities
nuc_ch_median_intensity <- rep(0, length(positions))
nuc_ch_bgmask_mean_intensity <- rep(0, length(positions))
spot_seg_ch_median_intensity <- rep(0, length(positions))
spot_seg_ch_bgmask_mean_intensity <- rep(0, length(positions))
ch1_median_intensity <- rep(0, length(positions))
ch1_bgmask_mean_intensity <- rep(0, length(positions))
ch2_median_intensity <- rep(0, length(positions))
ch2_bgmask_mean_intensity <- rep(0, length(positions))

# intensities in nuclear masks (mean, median, sd)
mean_nuclear_nuc_ch_intensity <- rep(0, mn*length(positions))
median_nuclear_nuc_ch_intensity <- rep(0, mn*length(positions)) 
sd_nuclear_nuc_ch_intensity <- rep(0, mn*length(positions)) 
mean_nuclear_spot_ch_intensity <- rep(0, mn*length(positions))
median_nuclear_spot_ch_intensity <- rep(0, mn*length(positions))
sd_nuclear_spot_ch_intensity <- rep(0, mn*length(positions))
mean_nuclear_ch1_intensity <- rep(0, mn*length(positions))
median_nuclear_ch1_intensity <- rep(0, mn*length(positions))
sd_nuclear_ch1_intensity <- rep(0, mn*length(positions))
mean_nuclear_ch2_intensity <- rep(0, mn*length(positions))
median_nuclear_ch2_intensity <- rep(0, mn*length(positions))
sd_nuclear_ch2_intensity <- rep(0, mn*length(positions))

# mask areas
nuc_area <- rep(0, mn*length(positions))
area_all_spotmask <- rep(0, mn*length(positions))
area_npmask <- rep(0, mn*length(positions))

# number of spots
num_spots <- rep(0, mn*length(positions))

# intensitites in spot and np masks
mean_nuc_ch_all_spotmask_intensity <- rep(0, mn*length(positions))
sd_nuc_ch_all_spotmask_intensity <- rep(0, mn*length(positions))
mean_spot_ch_all_spotmask_intensity <- rep(0, mn*length(positions))
sd_spot_ch_all_spotmask_intensity <- rep(0, mn*length(positions))
mean_ch1_all_spotmask_intensity <- rep(0, mn*length(positions))
sd_ch1_all_spotmask_intensity <- rep(0, mn*length(positions))
mean_ch2_all_spotmask_intensity <- rep(0, mn*length(positions))
sd_ch2_all_spotmask_intensity <- rep(0, mn*length(positions))

mean_nuc_ch_npmask_intensity <- rep(0, mn*length(positions))
sd_nuc_ch_npmask_intensity <- rep(0, mn*length(positions))
mean_spot_ch_npmask_intensity <- rep(0, mn*length(positions))
sd_spot_ch_npmask_intensity <- rep(0, mn*length(positions))
mean_ch1_npmask_intensity <- rep(0, mn*length(positions))
sd_ch1_npmask_intensity <- rep(0, mn*length(positions))
mean_ch2_npmask_intensity <- rep(0, mn*length(positions))
sd_ch2_npmask_intensity <- rep(0, mn*length(positions))

# other variables
# total nucleus counter
cnt <- 1 
#total cpmask counter
cnt_cp <- 1
# total position counter vector
pos <- rep(0, mn*length(positions))
# vector for storing BestSlices ranges determined for each position
slices_range_used <- rep(0, length(positions))

### cytoplasm masks and and segmented spot area within ####
# intensities in cytoplasm masks (mean, median, sd)
mean_cp_nuc_ch_intensity <- rep(0, mn*length(positions))
median_cp_nuc_ch_intensity <- rep(0, mn*length(positions)) 
sd_cp_nuc_ch_intensity <- rep(0, mn*length(positions)) 
mean_cp_spot_ch_intensity <- rep(0, mn*length(positions))
median_cp_spot_ch_intensity <- rep(0, mn*length(positions))
sd_cp_spot_ch_intensity <- rep(0, mn*length(positions))
mean_cp_ch1_intensity <- rep(0, mn*length(positions))
median_cp_ch1_intensity <- rep(0, mn*length(positions))
sd_cp_ch1_intensity <- rep(0, mn*length(positions))
mean_cp_ch2_intensity <- rep(0, mn*length(positions))
median_cp_ch2_intensity <- rep(0, mn*length(positions))
sd_cp_ch2_intensity <- rep(0, mn*length(positions))

# cytoplasm mask areas
cp_area <- rep(0, mn*length(positions))
area_all_spotmask_cp <- rep(0, mn*length(positions))
area_spotfreemask_cp <- rep(0, mn*length(positions)) #this corresponds to the cp mask area minus the sum of all spot areas in the cytoplasm

# number of spots in cytoplasm area
num_spots_cp <- rep(0, mn*length(positions))

# intensities in the summed area of cytoplasmic spot masks and spot-free cytoplasm area
mean_nuc_ch_all_cpspotmask_intensity <- rep(0, mn*length(positions))
sd_nuc_ch_all_cpspotmask_intensity <- rep(0, mn*length(positions))
mean_spot_ch_all_cpspotmask_intensity <- rep(0, mn*length(positions))
sd_spot_ch_all_cpspotmask_intensity <- rep(0, mn*length(positions))
mean_ch1_all_cpspotmask_intensity <- rep(0, mn*length(positions))
sd_ch1_all_cpspotmask_intensity <- rep(0, mn*length(positions))
mean_ch2_all_cpspotmask_intensity <- rep(0, mn*length(positions))
sd_ch2_all_cpspotmask_intensity <- rep(0, mn*length(positions))

mean_nuc_ch_spotfree_cp_intensity <- rep(0, mn*length(positions))
sd_nuc_ch_spotfree_cp_intensity <- rep(0, mn*length(positions))
mean_spot_ch_spotfree_cp_intensity <- rep(0, mn*length(positions))
sd_spot_ch_spotfree_cp_intensity <- rep(0, mn*length(positions))
mean_ch1_spotfree_cp_intensity <- rep(0, mn*length(positions))
sd_ch1_spotfree_cp_intensity <- rep(0, mn*length(positions))
mean_ch2_spotfree_cp_intensity <- rep(0, mn*length(positions))
sd_ch2_spotfree_cp_intensity <- rep(0, mn*length(positions))

#######################################################################################################################################
##############################################    POSITION LOOP i 1: length(positions)  ###############################################
#######################################################################################################################################

for(i in 1:length(positions)){
  # read-in multi-channel image stack i
  data <- readImage(paste0(folder,"/",positions[i]))
  
  # determine the number of channels defined for this run
  number_of_ch <- length(which(c(channels[]!=0)==TRUE))
  
  # split image data into the respective channels
  # empty channels are filled with 0
  if(number_of_ch<2){print("WARNING: You need at least 2 channel images for the analysis. Please check if the channel definition or the format of your input images is correct.")}
  if(number_of_ch==2){
    nuc_data <- data[,,seq(which(channels[]==nuc_ch),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    spot_data <- data[,,seq(which(channels[]==spot_seg_ch),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    ch1_data <- 0
    ch2_data <- 0
  }
  if(number_of_ch==3){
    nuc_data <- data[,,seq(which(channels[]==nuc_ch),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    spot_data <- data[,,seq(which(channels[]==spot_seg_ch),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    ch1_data <- data[,,seq(which(channels[]==ch1),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    ch2_data <- 0
  }
  if(number_of_ch==4){
    nuc_data <- data[,,seq(which(channels[]==nuc_ch),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    spot_data <- data[,,seq(which(channels[]==spot_seg_ch),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    ch1_data <- data[,,seq(which(channels[]==ch1),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
    ch2_data <- data[,,seq(which(channels[]==ch2),zperch*length(which(c(channels[]!=0)==TRUE)),by=length(which(c(channels[]!=0)==TRUE)))]
  }
  
  ################### Data reduction 1: select z-slices of focus volume (BestSlices) ###############################
  # This step confines the z-range used for creating z-projections to the nuclear z-volume
  # findBestSlices() uses the nuc_channel intensity distribution along z to select the z-range containing the full nuclear volume
  # e.g. black z-slices with no information above or below the imaged nuclei are removed
  # findBestSlices() returns 0 if a position is out of focus or contains no cells (details in the function)
  
  # determine best_slices_range used for z-reduction
  best_slices_range <- findBestSlices(nuc_data,zperch,passfilter)
  
  # additional filtering step 
  # only z-ranges greater than 10 will be considered
  # if below or 0 all slices will be used for the projection and slices_range_used will be marked with "1-zperch_NA"
  if (length(best_slices_range)<10){
    best_slices_range <- c(1:zperch)
    slices_range_used[i] <- paste0(min(best_slices_range),"-", max(best_slices_range),"_NA")
  }else{
  # initialize variable to document slices range for each position
  slices_range_used[i] <- paste0(min(best_slices_range),"-", max(best_slices_range))
  }
  
  ################### Data reduction 2: dimensionality reduction by z-projection ##################################
  # generate maximum z-projection for best z-slices of the nuc channel
  # maximum projection is recommended for nucleus segmentation
  # define which type of projections were used for each channel by changing projection_type in the order nuc_ch, spot_ch, ch1, ch2
  # NOTE that projection_type only documents the projection type used
  # for changing, the function apply()'d to the images needs to be manually changed in the subsequent subsections 

  projection_type <- c("max","max", "max", "max")
  nuc_data_best <- nuc_data[,,best_slices_range]
  nuc_data_projection <- apply(nuc_data_best, c(1:2),max)

  # generate maximum projections for spot and other channels using the best slices z-range
  if(number_of_ch<2){print("WARNING: You need at least 2 channel images for the analysis. Please check if the channel definition or the format of your input images is correct.")}
  if(number_of_ch==2){
    spot_data_best <- spot_data[,,best_slices_range]
    spot_data_best_projection <- apply(spot_data_best, c(1,2),max)
    # ch1 placeholder projection (all pixels 0)
    ch1_data_best_projection <- Image(0,dim=c(dim(nuc_data)[1],dim(nuc_data)[2]))
    # ch2 placeholder projection (all pixels 0)
    ch2_data_best_projection <- Image(0,dim=c(dim(nuc_data)[1],dim(nuc_data)[2]))
    }
  if(number_of_ch==3){
    spot_data_best <- spot_data[,,best_slices_range]
    spot_data_best_projection <- apply(spot_data_best, c(1,2),max)
    ch1_data_best <- ch1_data[,,best_slices_range]
    ch1_data_best_projection <- apply(ch1_data_best, c(1,2),max)
    # ch2 placeholder projection (all pixels zero)
    ch2_data_best_projection <- Image(0,dim=c(dim(nuc_data)[1],dim(nuc_data)[2]))
  }
  if(number_of_ch==4){
    spot_data_best <- spot_data[,,best_slices_range]
    spot_data_best_projection <- apply(spot_data_best, c(1,2),max)
    ch1_data_best <- ch1_data[,,best_slices_range]
    ch1_data_best_projection <- apply(ch1_data_best, c(1,2),max)
    ch2_data_best <- ch2_data[,,best_slices_range]
    ch2_data_best_projection <- apply(ch2_data_best, c(1,2),max)
  }
  
  ################################### Segmentation step 1A: Nuclei segmentation ##################################################
  
  # generate mask for all nuclei in the current position
  # makeNucMask() uses a gauss-blurred best slices z-projection of the nuc_channel and the above-defined adjustable parameters
  nucmask <- makeNucMask(gblur(nuc_data_projection, sigma=1), nuc_threshold_offset, nucSize_cutoff, remove_bordercells)
  
  # if no best_slices_range could be determined for this position and was thus set to 1:zperch or no nucleus was found:
  # create a dummy image with one squared "nucleus" of value 1 in the top left corner (40000 px area) and background intensity 0
  # the same is applied to the cytoplasm mask, which will receive the same placeholder image
  # all parameters that would be needed for downstream processing will be set for this placeholder dataset
  if(length(best_slices_range)==zperch | max(nucmask)==0){
    dummy_image <- Image(0,dim(nuc_data)[1:2])
    dummy_image[50:249,50:249] <- 1
    nucmask <- makeNucMask(dummy_image, nuc_threshold_offset, nucSize_cutoff, remove_bordercells)
    cytomask <- Image(0,dim(nuc_data)[1:2]) 
    cytomask[30:269,30:269] <- 1
    cpmask <- matchCellMasks(nucmask,cytomask, returnMaskType = "cytoplasm")
    matchingtable <- data.frame(final_nuc_ids=1,final_cpmask_ids=1)
    rois <- as.numeric(matchingtable$final_nuc_ids)
    cp_rois <- as.numeric(matchingtable$final_cpmask_ids)
    allspotmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
    allnpmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
    allspotmasks_cp <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
    allspotfreemasks_cp <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  }else{
  
  # initialize variables for storing whole-image spot and nucleoplasm masks (npmasks)
  # whole-image spot and npmasks are later only used for visualization of segmentation results in RGB images (see below)
  allspotmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  allnpmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  
  # initialize variables for storing whole-image spot and spotfree masks of cytoplasm
  # whole-image spot and npmasks are later only used for visualization of segmentation results in RGB images (see below)
  allspotmasks_cp <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  allspotfreemasks_cp <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])

  ################################### Segmentation step 1B: Cell body and cytoplasm area segmentation ##################################################

  # generate mask for all cell bodies (cytoplasm+nucleus are combined) in the current position
  # makeCytoMask() requires both the nucleus and the cytoplasm marker (ch1) channel as well as three adjustable parameters as input
  cytomask <- makeCytoMask(nuc_data_projection,ch1_data_best_projection,cytoSize_cutoff,cyto_threshold,remove_bordercells)
  # display masks (optional)
  # display(cytomask)
  # display(colorLabels(cytomask), all=TRUE)
  
  #if no suitable cytomask was found, repeat the placeholder procedure from above, else proceed
  if(max(cytomask)==0){
    dummy_image <- Image(0,dim(nuc_data)[1:2])
    dummy_image[50:249,50:249] <- 1
    nucmask <- makeNucMask(dummy_image, nuc_threshold_offset, nucSize_cutoff, remove_bordercells)
    cytomask <- Image(0,dim(nuc_data)[1:2]) 
    cytomask[30:269,30:269] <- 1
    cpmask <- matchCellMasks(nucmask,cytomask, returnMaskType = "cytoplasm")
    matchingtable <- data.frame(final_nuc_ids=1,final_cpmask_ids=1)
    rois <- as.numeric(matchingtable$final_nuc_ids)
    cp_rois <- as.numeric(matchingtable$final_cpmask_ids)
    allspotmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
    allnpmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
    allspotmasks_cp <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
    allspotfreemasks_cp <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  }else{
  ################################### "Mask matching" and filtering step ##################################################
  
  # In this part the segmented nucleus and cell body masks are compared and matched
  # to allow subsequent removal of nuclei masks or cell body masks without matching counterpart (e.g. nuclei masks without matching cytoplasm mask).
  # Moreover, cell body masks are simultaneously transformed to cytoplasm masks, which do not comprise the nucleus area.
  # The function matchCellMasks() requires nucmask and cytomask objects from the same image in the order indicated below (details in the function).
  # Setting the parameter returnMaskType = "nuclei" returns only nuclei masks having a matching cytoplasm mask.
  # Setting the parameter returnMaskType = "cytoplasm" returns cytoplasm masks having a matching nuclear mask.
  # Setting the parameter returnMaskType = "matchingmask" returns a data frame containing nucleus mask id matched to the corresponding cpmask id (after having filtered out unmatched masks or duplicates)
  # matchCellMasks() also filters out masks with unambiguous matching, for example cytoplasm masks with two matching nuclei, which occasionally arise during segmentation.
  # In brief, the function generates a binary image in which each nucleus is represented by a circular mask of value 1 on a dark background of value 0.
  # The circles are drawn with the nucleus mask object coordinates (center of mass) as center and radius = min.radius + 15 px.
  # This circle mask image is then used to retrieve the IDs (bwlabel indices) of the cytoplasm mask that belongs to these nuclei for matching and subsequent filtering.
  
  # return matching nucleus masks and set this as new nucmask for further segmentation
  nucmask_match <- matchCellMasks(nucmask,cytomask, returnMaskType = "nuclei")
  nucmask <- nucmask_match
  
  # return matching cytoplasm masks
  cpmask <-matchCellMasks(nucmask,cytomask, returnMaskType = "cytoplasm")
  # return data frame with matching information (nucleus id --> cpmask id)
  matchingtable <-matchCellMasks(nucmask,cytomask, returnMaskType = "matchingtable")
  # retrieve indices of single nucleus masks for looping over nuclei
  rois <- as.numeric(matchingtable$final_nuc_ids)
  
  # Optional: display the masks for visual inspection
  #display(colorLabels(nucmask_match), all=TRUE)
  #display(colorLabels(cytomask_match), all=TRUE)
  #display(nuc_data_projection^0.5)
  }
  }
  ######################################  Quantify whole-image background estimates ################################################
  
  # This step measures intensities of whole projections in all channels that are potentially informative for image background estimation
  # 1. computes the median intensity of the whole image projection
  # 2. computes the mean intensity in the negative of the nuclear masks
  # the first value is a robust estimate for the background of the image, however will be biased for positions with high cell density 
  # the second value is more precisely the extra-nuclear background, however will be biased by border nuclei contribution when choosing remove_bordercells <- "y" or in cases where many nuclei are removed due to failed cpmask matches
  
  nuc_ch_median_intensity[i] <- median(nuc_data_projection)
  nuc_ch_bgmask_mean_intensity[i] <- mean(nuc_data_projection[nucmask==0])
  spot_seg_ch_median_intensity[i] <- median(spot_data_best_projection)
  spot_seg_ch_bgmask_mean_intensity[i] <- mean(spot_data_best_projection[nucmask==0])
  ch1_median_intensity[i] <- median(ch1_data_best_projection)
  ch1_bgmask_mean_intensity[i] <- mean(ch1_data_best_projection[nucmask==0])
  ch2_median_intensity[i] <- median(ch2_data_best_projection)
  ch2_bgmask_mean_intensity[i] <- mean(ch2_data_best_projection[nucmask==0])

#######################################################################################################################################
##############################################    NUCLEUS LOOP j in 1: length(rois)  ##################################################
#######################################################################################################################################  

for(j in 1:length(rois)){
    # whole nucleus quantification for current nucleus j (mean, median, sd in all channels)
    mean_nuclear_nuc_ch_intensity[cnt] <- mean(nuc_data_projection[nucmask==rois[j]])
    median_nuclear_nuc_ch_intensity[cnt] <- median(nuc_data_projection[nucmask==rois[j]])
    sd_nuclear_nuc_ch_intensity[cnt] <- sd(nuc_data_projection[nucmask==rois[j]])
    
    mean_nuclear_spot_ch_intensity[cnt] <- mean(spot_data_best_projection[nucmask==rois[j]])
    median_nuclear_spot_ch_intensity[cnt] <- median(spot_data_best_projection[nucmask==rois[j]])
    sd_nuclear_spot_ch_intensity[cnt] <- sd(spot_data_best_projection[nucmask==rois[j]])
    
    mean_nuclear_ch1_intensity[cnt] <- mean(ch1_data_best_projection[nucmask==rois[j]])
    median_nuclear_ch1_intensity[cnt] <- median(ch1_data_best_projection[nucmask==rois[j]])
    sd_nuclear_ch1_intensity[cnt] <- sd(ch1_data_best_projection[nucmask==rois[j]])
    
    mean_nuclear_ch2_intensity[cnt] <- mean(ch2_data_best_projection[nucmask==rois[j]])
    median_nuclear_ch2_intensity[cnt] <- median(ch2_data_best_projection[nucmask==rois[j]])
    sd_nuclear_ch2_intensity[cnt] <- sd(ch2_data_best_projection[nucmask==rois[j]])

    # document nucleus area
    # cnt is the continuous nucleus counter over all positions analyzed
    nuc_area[cnt] <- sum(nucmask==rois[j])
    
    # postion number assignment
    pos[cnt] <- i

    # compute nucleus shape features and xy-positions and export them to a file
    # the temporary results files created during the run for each nucleus (nuc_shape_features_pos_1_cell_1 etc.) will later be deleted and combined into a unified data table
    nuc_shape_features <- cbind(i,j,computeFeatures.moment(nucmask==rois[j]),computeFeatures.shape(nucmask==rois[j])) # cbind(i,j,...) and nucmask==rois[j] for both
    colnames(nuc_shape_features)[1:2] <- c("position", "cell")
    write.table(file=paste0(paste0(results_folder,"/",date_time,"_nuc_shape_features"),"/nuc_shape_features_pos_",i,"_cell_",j), nuc_shape_features, sep = "\t", append = F, row.names = F, quote = F)

    ################################### Segmentation step 2A: Spot segmentation in nuclei ##################################################

    #generate spot masks in the current nucleus j
    spot_shape_features <- NULL
    mask <- nucmask
    mask[mask!=rois[j]] <- 0 #set all pixels of nucmask to 0 that do not have the mask intensity value of the current nucleus j
    nuc <- gblur(spot_data_best_projection, sigma=1) #store a gauss-blurred image of the spot projection reduced to the current nucleus mask
    nuc[mask!=rois[j]] <- NA #set all regions in nuc to NA that are outside of the current nuclear mask (mask)
    
    # hand over image to the spot segmentation function with a threshold that depends on the nuclear median and user-defined spot_cutoff
    # alternatively, an absolute thresholding can also be handed over spotmask <- makeSpotMask(nuc, 0.05, mask) 
    # makeSpotMask() generates a binary mask where all pixels above the threshold are set to value 1 and all others to 0
    # the mask is further processed to smoothen connected pixel groups
    # connected pixel groups (spots) are identified as individual objects and assigned an id to enable downsteam single spot analyses 
    spotmask <- makeSpotMask(nuc, median(nuc, na.rm=T)*spot_cutoff, mask)

    #convert spotmask into an Image object to calculate spot position and shape features for each individual spot in the current nucleus j
    spotmask_image <-Image(spotmask, dim=dim(spotmask))
    spot_shape_features <- NULL
    spot_shape_features <- cbind(computeFeatures.moment(spotmask_image),computeFeatures.shape(spotmask_image))
    
    #if spot_shape_features cannot be computed because no spots were segmented in the current nucleus j, set all values to 0 but keep position and cell id
    if(length(spot_shape_features)<11){
      spot_shape_features <- array(dim=c(1,13))
      colnames(spot_shape_features) <- c("position","cell","m.cx","m.cy","m.majoraxis","m.eccentricity","m.theta","s.area","s.perimeter",	"s.radius.mean", "s.radius.sd",	"s.radius.min",	"s.radius.max")
      spot_shape_features[,"position"] <- i
      spot_shape_features[,"cell"] <- j
      spot_shape_features[,3:13] <- 0
    }else{
    # format spot shape features data table to contain position and cell id
    spot_shape_features <- cbind(i,j,spot_shape_features) #cbind(i,j,spot_shape_features)
    colnames(spot_shape_features) <- c("position","cell","m.cx","m.cy","m.majoraxis","m.eccentricity","m.theta","s.area","s.perimeter",	"s.radius.mean", "s.radius.sd",	"s.radius.min",	"s.radius.max")
    }
    # create a combined spotmask where all values that have been identified as spots (value larger than 0) are combined into one unified region with the value 1
    all_spotmask <- spotmask
    all_spotmask[all_spotmask!=0] <- 1
    # add to whole-image spotmask
    allspotmasks[all_spotmask>0] <- 1
    
    # create a nucleoplasm mask with pixel value 1 for the current nucleus j - uses the current nucleus mask "mask"
    mask[mask!=0] <- 1
    npmask <- mask-all_spotmask
    # add to whole-image npmask
    allnpmasks[npmask>0] <- 1

    ############################### Quantify all_spotmask, npmask features and intensities for the current nucleus j #########################################

    # compute the sum of all spot areas in the current nucleus
    area_all_spotmask[cnt] <- sum(all_spotmask>0)
    area_npmask[cnt] <- sum(npmask>0)

    # determine number of segmented spots (spotmask objects) per nucleus
    # -1 is required because image regions in spotmask always comprise the background region (value 0) in addition to the spots
    num_spots[cnt] <- length(table(spotmask))-1

    # compute mean intensity and standard deviation for all_spotmask region in all channels
    mean_nuc_ch_all_spotmask_intensity[cnt] <- mean(nuc_data_projection[all_spotmask>0])
    sd_nuc_ch_all_spotmask_intensity[cnt] <- sd(nuc_data_projection[all_spotmask>0])
    mean_spot_ch_all_spotmask_intensity[cnt] <- mean(spot_data_best_projection[all_spotmask>0])
    sd_spot_ch_all_spotmask_intensity[cnt] <- sd(spot_data_best_projection[all_spotmask>0])
    mean_ch1_all_spotmask_intensity[cnt] <- mean(ch1_data_best_projection[all_spotmask>0])
    sd_ch1_all_spotmask_intensity[cnt] <- sd(ch1_data_best_projection[all_spotmask>0])
    mean_ch2_all_spotmask_intensity[cnt] <- mean(ch2_data_best_projection[all_spotmask>0])
    sd_ch2_all_spotmask_intensity[cnt] <- sd(ch2_data_best_projection[all_spotmask>0])

    # compute mean intensity and standard deviation for npmask region in all channels
    mean_nuc_ch_npmask_intensity[cnt] <- mean(nuc_data_projection[npmask>0])
    sd_nuc_ch_npmask_intensity[cnt] <- sd(nuc_data_projection[npmask>0])
    mean_spot_ch_npmask_intensity[cnt] <- mean(spot_data_best_projection[npmask>0])
    sd_spot_ch_npmask_intensity[cnt] <- sd(spot_data_best_projection[npmask>0])
    mean_ch1_npmask_intensity[cnt] <- mean(ch1_data_best_projection[npmask>0])
    sd_ch1_npmask_intensity[cnt] <- sd(ch1_data_best_projection[npmask>0])
    mean_ch2_npmask_intensity[cnt] <- mean(ch2_data_best_projection[npmask>0])
    sd_ch2_npmask_intensity[cnt] <- sd(ch2_data_best_projection[npmask>0])

    # get numbering of spots for looping over single spots in current nucleus j
    rois_spots <- as.numeric(names(table(spotmask)[2:length(table(spotmask))]))
    

    
    ###########################################################################################################################################
    ##############################################    SINGLE SPOT LOOP s in 1: length(rois_spots) (Nucleus) ###################################
    ###########################################################################################################################################  
    
    # if no spots were segmented in the current nucleus, fill single spot intensities with zeros
    # also compute various other measures in the spot channel: minimum, maximum, median intensity & 20% / 95% percentiles of the spot mask pixel group intensities 
    if (mean(spot_shape_features[,3:13])==0){
      mean_nuc_ch_single_spot_intensity <- 0
      mean_spot_ch_single_spot_intensity <- 0
      min_spot_ch_single_spot_intensity <- 0
      max_spot_ch_single_spot_intensity <- 0
      median_spot_ch_single_spot_intensity <- 0
      twentieth_quant_spot_ch_single_spot_intensity <- 0
      ninetyfifth_quant_spot_ch_single_spot_intensity <- 0
      
      mean_ch1_single_spot_intensity <- 0
      mean_ch2_single_spot_intensity <- 0
    }else{
    # initialize single spot variables
    mean_nuc_ch_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    mean_spot_ch_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    min_spot_ch_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    max_spot_ch_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    median_spot_ch_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    twentieth_quant_spot_ch_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    ninetyfifth_quant_spot_ch_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    
    mean_ch1_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))
    mean_ch2_single_spot_intensity <- rep(0, length(rownames(spot_shape_features)))

    ############################### Quantify mean intensities in all channels for the current spot s ###########################################
    for(s in 1:length(rois_spots)){
      mean_nuc_ch_single_spot_intensity[s] <- mean(nuc_data_projection[spotmask==rois_spots[s]])
      mean_spot_ch_single_spot_intensity[s] <- mean(spot_data_best_projection[spotmask==rois_spots[s]])
      min_spot_ch_single_spot_intensity[s] <- min(spot_data_best_projection[spotmask==rois_spots[s]])
      max_spot_ch_single_spot_intensity[s] <- max(spot_data_best_projection[spotmask==rois_spots[s]])
      median_spot_ch_single_spot_intensity[s] <- median(spot_data_best_projection[spotmask==rois_spots[s]])
      twentieth_quant_spot_ch_single_spot_intensity[s] <- quantile(spot_data_best_projection[spotmask==rois_spots[s]],0.2)
      ninetyfifth_quant_spot_ch_single_spot_intensity[s] <- quantile(spot_data_best_projection[spotmask==rois_spots[s]],0.95)
      
      mean_ch1_single_spot_intensity[s] <- mean(ch1_data_best_projection[spotmask==rois_spots[s]])
      mean_ch2_single_spot_intensity[s] <- mean(ch2_data_best_projection[spotmask==rois_spots[s]])
    }
    #############################################################################################################################################
    ######################################################## END SINGLE SPOT LOOP s #############################################################
    #############################################################################################################################################
    }
    # compute spot shape features and xy-positions for all spots in the current nucleus and export them to a file
    # the temporary results files created during the run for each nucleus (spot_shape_features_pos_1_cell_1 etc.) will later be deleted and combined into a unified data table
    spot_shape_features <- cbind(spot_shape_features,mean_nuc_ch_single_spot_intensity,mean_spot_ch_single_spot_intensity,min_spot_ch_single_spot_intensity,max_spot_ch_single_spot_intensity,median_spot_ch_single_spot_intensity,twentieth_quant_spot_ch_single_spot_intensity,ninetyfifth_quant_spot_ch_single_spot_intensity,mean_ch1_single_spot_intensity,mean_ch2_single_spot_intensity)
    write.table(file=paste0(results_folder,"/",date_time,"_spot_shape_features","/spot_shape_features_pos_",i,"_cell_",j), spot_shape_features, sep = "\t", append = F, row.names = F, quote = F)

 
    #increase total nucleus counter
    cnt <- cnt +1
  }
#######################################################################################################################################
######################################################## END NUCLEUS LOOP j ###########################################################
#######################################################################################################################################

#######################################################################################################################################
##############################################    CPMASK LOOP jj in 1: length(cp_rois)  ################################################
#######################################################################################################################################   
  # retrieve indices of single cpmasks for looping
  cp_rois <- as.numeric(matchingtable$final_cpmask_ids)
  
  for(jj in 1:length(cp_rois)){
    # whole cytoplasm quantification for current cpmask jj (mean, median, sd in all channels)
    # cnt_cp is the continuous nucleus/cpmask counter over all positions analyzed
    mean_cp_nuc_ch_intensity[cnt_cp] <- mean(nuc_data_projection[cpmask==cp_rois[jj]])
    median_cp_nuc_ch_intensity[cnt_cp] <- median(nuc_data_projection[cpmask==cp_rois[jj]])
    sd_cp_nuc_ch_intensity[cnt_cp] <- sd(nuc_data_projection[cpmask==cp_rois[jj]])
    
    mean_cp_spot_ch_intensity[cnt_cp] <- mean(spot_data_best_projection[cpmask==cp_rois[jj]])
    median_cp_spot_ch_intensity[cnt_cp] <- median(spot_data_best_projection[cpmask==cp_rois[jj]])
    sd_cp_spot_ch_intensity[cnt_cp] <- sd(spot_data_best_projection[cpmask==cp_rois[jj]])
    
    mean_cp_ch1_intensity[cnt_cp] <- mean(ch1_data_best_projection[cpmask==cp_rois[jj]])
    median_cp_ch1_intensity[cnt_cp] <- median(ch1_data_best_projection[cpmask==cp_rois[jj]])
    sd_cp_ch1_intensity[cnt_cp] <- sd(ch1_data_best_projection[cpmask==cp_rois[jj]])
    
    mean_cp_ch2_intensity[cnt_cp] <- mean(ch2_data_best_projection[cpmask==cp_rois[jj]])
    median_cp_ch2_intensity[cnt_cp] <- median(ch2_data_best_projection[cpmask==cp_rois[jj]])
    sd_cp_ch2_intensity[cnt_cp] <- sd(ch2_data_best_projection[cpmask==cp_rois[jj]])
    
    # document cpmask area
    cp_area[cnt_cp] <- sum(cpmask==cp_rois[jj])
    
    # no shape features or coordinates are computed for the cytoplasm masks
    
    ################################### Segmentation step 2B: Spot segmentation in cytoplasm masks ######################################
    # generate spot masks in the current cytoplasm jj
    # empty spot_shape_features_cp variable
    spot_shape_features_cp <- NULL
    mask_cp <- cpmask
    mask_cp[mask_cp!=cp_rois[jj]] <- 0 #set all pixels of cpmask to 0 that do not have the mask id value of the current cpmask jj
    cp <- gblur(spot_data_best_projection, sigma=1) #store a gauss-blurred image of the spot projection reduced to the current cpmask jj
    cp[mask_cp!=cp_rois[jj]] <- NA #set all regions in cp to NA that are outside of the current cpmask jj
    
    # hand over image to the spot segmentation function with a threshold that depends on the nuclear median and user-defined spot_cutoff
    # alternatively, an absolute thresholding can also be handed over spotmask <- makeSpotMask(nuc, 0.05, mask) 
    # makeSpotMask() generates a binary mask where all pixels above the threshold are set to value 1 and all others to 0
    # the mask is further processed to smoothen connected pixel groups
    # connected pixel groups (spots) are identified as individual objects and assigned an id using bwlabel() to enable downsteam single spot analyses 
    spotmask_cp <- makeSpotMask(cp, median(cp, na.rm=T)*spot_cutoff_cytoplasm, mask_cp)
    # renumerate spotmask_cp 
    spotmask_cp <- bwlabel(spotmask_cp)
    
    #convert spotmask_cp into an Image object to calculate spot position and shape features for each individual spot in the current cpmask jj
    spotmask_image_cp <-Image(spotmask_cp, dim=dim(spotmask_cp))
    spot_shape_features_cp <- NULL
    spot_shape_features_cp <- cbind(computeFeatures.moment(spotmask_image_cp),computeFeatures.shape(spotmask_image_cp))
    
    #if spot_shape_features_cp cannot be computed because no spots were segmented in the current cpmask jj, set all values to 0 but keep position and cell id
    if(length(spot_shape_features_cp)<11){
      spot_shape_features_cp <- array(dim=c(1,13))
      colnames(spot_shape_features_cp) <- c("position","cell","m.cx","m.cy","m.majoraxis","m.eccentricity","m.theta","s.area","s.perimeter",	"s.radius.mean", "s.radius.sd",	"s.radius.min",	"s.radius.max")
      spot_shape_features_cp[,"position"] <- i
      spot_shape_features_cp[,"cell"] <- jj
      spot_shape_features_cp[,3:13] <- 0
    }else{
      # format spot shape features data table to contain position and cell id
      spot_shape_features_cp <- cbind(i,jj,spot_shape_features_cp) #cbind(i,j,spot_shape_features)
      colnames(spot_shape_features_cp) <- c("position","cell","m.cx","m.cy","m.majoraxis","m.eccentricity","m.theta","s.area","s.perimeter",	"s.radius.mean", "s.radius.sd",	"s.radius.min",	"s.radius.max")
    }
    # create a combined spotmask where all values that have been identified as spots (value larger than 0) are combined into one unified region with the value 1
    all_spotmask_cp <- spotmask_cp
    all_spotmask_cp[all_spotmask_cp!=0] <- 1
    # add to whole-image spotmask
    allspotmasks_cp[all_spotmask_cp>0] <- 1
    
    # create a spotfree cytoplasm mask (inverse of the spot mask for the current cpmask) with pixel value 1 for the current cpmask j - uses the current cpmask mask_cp
    mask_cp[mask_cp!=0] <- 1
    spotfreemask <- mask_cp-all_spotmask_cp
    # add to whole-image npmask
    allspotfreemasks_cp[spotfreemask>0] <- 1
    
    ############################### Quantify all_spotmask_cp, spotfreemask features and intensities for the current cpmask jj #########################################
    
    # compute the sum of all spot areas in the current cpmask
    area_all_spotmask_cp[cnt_cp] <- sum(all_spotmask_cp>0)
    area_spotfreemask_cp[cnt_cp] <- sum(spotfreemask>0)
    
    # determine number of segmented spots (spotmask objects) for the current cpmask
    # -1 is required because image regions in spotmask always comprise the background region (value 0) in addition to the spots
    num_spots_cp[cnt_cp] <- length(table(spotmask_cp))-1
    
    # compute mean intensity and standard deviation for all_spotmask_cp region in all channels
    mean_nuc_ch_all_cpspotmask_intensity[cnt_cp] <- mean(nuc_data_projection[all_spotmask_cp>0])
    sd_nuc_ch_all_cpspotmask_intensity[cnt_cp] <- sd(nuc_data_projection[all_spotmask_cp>0])
    mean_spot_ch_all_cpspotmask_intensity[cnt_cp] <- mean(spot_data_best_projection[all_spotmask_cp>0])
    sd_spot_ch_all_cpspotmask_intensity[cnt_cp] <- sd(spot_data_best_projection[all_spotmask_cp>0])
    mean_ch1_all_cpspotmask_intensity[cnt_cp] <- mean(ch1_data_best_projection[all_spotmask_cp>0])
    sd_ch1_all_cpspotmask_intensity[cnt_cp] <- sd(ch1_data_best_projection[all_spotmask_cp>0])
    mean_ch2_all_cpspotmask_intensity[cnt_cp] <- mean(ch2_data_best_projection[all_spotmask_cp>0])
    sd_ch2_all_cpspotmask_intensity[cnt_cp] <- sd(ch2_data_best_projection[all_spotmask_cp>0])
    
    # compute mean intensity and standard deviation for spot-free cytoplasm mask region in all channels
    mean_nuc_ch_spotfree_cp_intensity[cnt_cp] <- mean(nuc_data_projection[spotfreemask>0])
    sd_nuc_ch_spotfree_cp_intensity[cnt_cp] <- sd(nuc_data_projection[spotfreemask>0])
    mean_spot_ch_spotfree_cp_intensity[cnt_cp] <- mean(spot_data_best_projection[spotfreemask>0])
    sd_spot_ch_spotfree_cp_intensity[cnt_cp] <- sd(spot_data_best_projection[spotfreemask>0])
    mean_ch1_spotfree_cp_intensity[cnt_cp] <- mean(ch1_data_best_projection[spotfreemask>0])
    sd_ch1_spotfree_cp_intensity[cnt_cp] <- sd(ch1_data_best_projection[spotfreemask>0])
    mean_ch2_spotfree_cp_intensity[cnt_cp] <- mean(ch2_data_best_projection[spotfreemask>0])
    sd_ch2_spotfree_cp_intensity[cnt_cp] <- sd(ch2_data_best_projection[spotfreemask>0])
    
    # get numbering of spots for looping over single spots in current nucleus j
    rois_spots_cp <- as.numeric(names(table(spotmask_cp)[2:length(table(spotmask_cp))]))
    
    ###########################################################################################################################################
    ####################################  SINGLE SPOT LOOP ss in 1: length(rois_spots_cp) (Cytoplasm) #########################################
    ###########################################################################################################################################  
    # if no spots were segmented in the current cpmask, fill single spot intensities with zeros
    # also compute various other measures in the spot channel: minimum, maximum, median intensity & 20% / 95% percentiles of the spot mask pixel group intensities 
    if (mean(spot_shape_features_cp[,3:13])==0){
      mean_nuc_ch_single_spot_cp_intensity <- 0
      mean_spot_ch_single_spot_cp_intensity <- 0
      min_spot_ch_single_spot_cp_intensity <- 0
      max_spot_ch_single_spot_cp_intensity <- 0
      median_spot_ch_single_spot_cp_intensity <- 0
      twentieth_quant_spot_ch_single_spot_cp_intensity <- 0
      ninetyfifth_quant_spot_ch_single_spot_cp_intensity <- 0
      
      mean_ch1_single_spot_cp_intensity <- 0
      mean_ch2_single_spot_cp_intensity <- 0
    }else{
      # initialize single spot variables
      mean_nuc_ch_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      mean_spot_ch_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      min_spot_ch_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      max_spot_ch_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      median_spot_ch_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      twentieth_quant_spot_ch_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      ninetyfifth_quant_spot_ch_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      
      mean_ch1_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      mean_ch2_single_spot_cp_intensity <- rep(0, length(rownames(spot_shape_features_cp)))
      
      ############################### Quantify intensities in all channels for the current spot ss ###########################################
      for(ss in 1:length(rois_spots_cp)){
        mean_nuc_ch_single_spot_cp_intensity[ss] <- mean(nuc_data_projection[spotmask_cp==rois_spots_cp[ss]])
        mean_spot_ch_single_spot_cp_intensity[ss] <- mean(spot_data_best_projection[spotmask_cp==rois_spots_cp[ss]])
        min_spot_ch_single_spot_cp_intensity[ss] <- min(spot_data_best_projection[spotmask_cp==rois_spots_cp[ss]])
        max_spot_ch_single_spot_cp_intensity[ss] <- max(spot_data_best_projection[spotmask_cp==rois_spots_cp[ss]])
        median_spot_ch_single_spot_cp_intensity[ss] <- median(spot_data_best_projection[spotmask_cp==rois_spots_cp[ss]])
        twentieth_quant_spot_ch_single_spot_cp_intensity[ss] <- quantile(spot_data_best_projection[spotmask_cp==rois_spots_cp[ss]],0.2)
        ninetyfifth_quant_spot_ch_single_spot_cp_intensity[ss] <- quantile(spot_data_best_projection[spotmask_cp==rois_spots_cp[ss]],0.95)
        
        mean_ch1_single_spot_cp_intensity[ss] <- mean(ch1_data_best_projection[spotmask_cp==rois_spots_cp[ss]])
        mean_ch2_single_spot_cp_intensity[ss] <- mean(ch2_data_best_projection[spotmask_cp==rois_spots_cp[ss]])
      }
      ###########################################################################################################################################
      ######################################################## END SINGLE SPOT LOOP ss ##########################################################
      ###########################################################################################################################################
      }
    # compute spot shape features and xy-positions for all spots in the current cpmask and export them to a file
    # the temporary results files created during the run for each nucleus (spot_shape_features_pos_1_cell_1 etc.) will later be deleted and combined into a unified data table
    spot_shape_features_cp <- cbind(spot_shape_features_cp,mean_nuc_ch_single_spot_cp_intensity,mean_spot_ch_single_spot_cp_intensity,min_spot_ch_single_spot_cp_intensity,max_spot_ch_single_spot_cp_intensity,median_spot_ch_single_spot_cp_intensity,twentieth_quant_spot_ch_single_spot_cp_intensity,ninetyfifth_quant_spot_ch_single_spot_cp_intensity,mean_ch1_single_spot_cp_intensity,mean_ch2_single_spot_cp_intensity)
    write.table(file=paste0(results_folder,"/",date_time,"_spot_shape_features_cytoplasm","/spot_shape_features_cp_pos_",i,"_cell_",jj), spot_shape_features_cp, sep = "\t", append = F, row.names = F, quote = F)
    
    #increase total cpmask counter
    cnt_cp <- cnt_cp +1
  }
#######################################################################################################################################
######################################################## END CPMASK LOOP jj ###########################################################
#######################################################################################################################################
  
  # export projections with generated nucmask and all spotmasks as .tif for potential re-analysis of data
  # empty channels will contain all pixels with value 0
  writeImage(combine(nuc_data_projection,nucmask, spot_data_best_projection, allspotmasks, ch1_data_best_projection, cpmask, allspotmasks_cp, ch2_data_best_projection), bits.per.sample=16, files=paste0(results_folder,"/",gsub(".tif","_BestSlices_projections_with_masks_binary.tif",positions[i])))
  
  ############################### Create images for visualization of segmentation results #########################################
  
  # RGB images scaled to the maximum image intensity for quick visual inspection of results
  # write RGB image overlays of masks onto projections of current position
  # projections are divided by their maximum pixel value for scaling
  best_nuc_projection_with_masks <- toRGB(nuc_data_projection/max(nuc_data_projection))
  best_spot_projection_with_masks_nonpmasks <- toRGB(spot_data_best_projection/max(spot_data_best_projection))
  best_spot_projection_with_masks <- toRGB(spot_data_best_projection/max(spot_data_best_projection))
  best_ch1_projection_with_masks <- toRGB(ch1_data_best_projection/max(ch1_data_best_projection))
  best_ch2_projection_with_masks <- toRGB(ch2_data_best_projection/max(ch2_data_best_projection))
  best_ch1_projection_with_cpmasks <- toRGB(ch1_data_best_projection/max(ch1_data_best_projection))
  best_ch1_projection_with_cp_spot_masks <- toRGB(spot_data_best_projection/max(spot_data_best_projection))
  best_ch1_projection_with_cp_spot_and_spotfree_masks<- toRGB(spot_data_best_projection/max(spot_data_best_projection))
  
  # nucleus mask outlines = yellow
  # spot mask outlines = blue
  # nucleoplasm mask = transparent red shading
  best_nuc_projection_with_masks <- paintObjects(nucmask, best_nuc_projection_with_masks, opac=c(1,0), col=c("yellow","yellow"))
  best_spot_projection_with_masks <- paintObjects(allspotmasks, best_spot_projection_with_masks, opac=c(1,0), col=c("blue","blue"))
  best_spot_projection_with_masks <- paintObjects(allnpmasks, best_spot_projection_with_masks, opac=c(0,0.1), col=c("red","red"))
  best_spot_projection_with_masks_nonpmasks <- paintObjects(allspotmasks, best_spot_projection_with_masks_nonpmasks, opac=c(1,0), col=c("blue","blue"))
  best_ch1_projection_with_masks <- paintObjects(allspotmasks, best_ch1_projection_with_masks, opac=c(1,0), col=c("blue","blue"))
  best_ch1_projection_with_masks <- paintObjects(allnpmasks, best_ch1_projection_with_masks, opac=c(0,0.1), col=c("red","red"))
  best_ch2_projection_with_masks <- paintObjects(allspotmasks, best_ch2_projection_with_masks, opac=c(1,0), col=c("blue","blue"))
  best_ch2_projection_with_masks <- paintObjects(allnpmasks, best_ch2_projection_with_masks, opac=c(0,0.1), col=c("red","red"))
  
  #cytoplasm masks and cp spot masks
  # cytoplasm mask outlines = yellow
  # cytoplasmic spot mask outlines = blue
  # cytoplasm spot-free region = transparent red shading
  best_ch1_projection_with_cpmasks <- paintObjects(cpmask, best_ch1_projection_with_cpmasks, opac=c(1,0), col=c("yellow","yellow"))
  best_ch1_projection_with_cp_spot_masks <- paintObjects(allspotmasks_cp, best_ch1_projection_with_cp_spot_masks, opac=c(1,0), col=c("blue","blue"))
  best_ch1_projection_with_cp_spot_and_spotfree_masks <- paintObjects(allspotmasks_cp, best_ch1_projection_with_cp_spot_and_spotfree_masks, opac=c(1,0), col=c("blue","blue"))
  best_ch1_projection_with_cp_spot_and_spotfree_masks <- paintObjects(allspotfreemasks_cp, best_ch1_projection_with_cp_spot_and_spotfree_masks, opac=c(0,0.1), col=c("red","red"))
  writeImage(combine(best_nuc_projection_with_masks,best_spot_projection_with_masks,best_spot_projection_with_masks_nonpmasks,best_ch1_projection_with_masks,best_ch2_projection_with_masks,best_ch1_projection_with_cpmasks,best_ch1_projection_with_cp_spot_masks,best_ch1_projection_with_cp_spot_and_spotfree_masks), bits.per.sample=8, files=paste0(results_folder,"/",gsub(".tif","_BestSlices_projections_with_masks_scaledRGB.tif",positions[i])))
  
  ##########
  # write out non-scaled 16-bit grayscale images with mask outlines as separate channels for creating customizable composites etc. in FIJI
  # the mask outlines are set to pixel value 65636 (nucleoplasm mask is not drawn as outline)
  blackimage <- Image(0,dim=c(dim(nuc_data)[1],dim(nuc_data)[2]))
  blackimage_nucmask <- channel(paintObjects(nucmask, blackimage, opac=c(1,0), col=c("white","white")), mode = "gray")
  blackimage_allspotmasks <- channel(paintObjects(allspotmasks, blackimage, opac=c(1,0), col=c("white","white")), mode ="gray")
  blackimage_allnpmasks <- channel(paintObjects(allnpmasks, blackimage, opac=c(0,1), col=c("white","white")), mode="gray")
  blackimage_cpmask <- channel(paintObjects(cpmask, blackimage, opac=c(1,0), col=c("white","white")), mode = "gray")
  blackimage_allspotmasks_cp <- channel(paintObjects(allspotmasks_cp, blackimage, opac=c(1,0), col=c("white","white")), mode ="gray")
  blackimage_allspotfreemasks_cp <- channel(paintObjects(allspotfreemasks_cp, blackimage, opac=c(0,1), col=c("white","white")), mode="gray")
  writeImage(combine(nuc_data_projection,blackimage_nucmask,spot_data_best_projection,blackimage_allspotmasks,blackimage_allnpmasks,ch1_data_best_projection,blackimage_cpmask,blackimage_allspotmasks_cp,blackimage_allspotfreemasks_cp,ch2_data_best_projection), bits.per.sample=16, files=paste0(results_folder,"/",gsub(".tif","_BestSlices_projections_with_masks_outlines.tif",positions[i])))
}
#######################################################################################################################################
######################################################## END POSITION LOOP i ##########################################################
#######################################################################################################################################


##################################################################################################################################################
################################################_____ Collect NUCLEUS results____#################################################################
# collect, order and summarize nucleus shape features and xy positions
nuc_shape_folder <- paste0(results_folder,"/",list.files(results_folder)[grep("nuc_shape_features",list.files(results_folder))])
nuc_shape_files <- list.files(nuc_shape_folder)
nucleus_features <- data.frame()
for(b in 1:length(nuc_shape_files)){
  nucleus_features[b,1:length(colnames(read.table(file=paste0(nuc_shape_folder,"/",nuc_shape_files[b]),header=T)))] <- read.table(file=paste0(nuc_shape_folder,"/",nuc_shape_files[b]),header=T,)
}
#re-order nucleus_features according to columns position and cell (i and j) and export to results file, delete temporary nuc_shape_features folder
nucleus_features <- as.data.frame(nucleus_features)
nucleus_features <- nucleus_features[order(nucleus_features$position,nucleus_features$cell),]
write.table(file=paste0(results_folder,"/",date_time, "_nucleus_shape_features.csv"), nucleus_features, sep = "\t", append = F, row.names = F, quote = F)
if (file.exists(nuc_shape_folder)){unlink(nuc_shape_folder, recursive = T)}

#collect, order and summarize single spot shape, xy-position and intensity features
spot_shape_folder <- paste0(results_folder,"/",list.files(results_folder)[grep("spot_shape_features$",list.files(results_folder))])
spot_shape_files <- list.files(spot_shape_folder)
spot_features <- data.frame()
for(r in 1:length(spot_shape_files)){
  spot_features <- rbind(spot_features,read.table(file=paste0(spot_shape_folder,"/",spot_shape_files[r]),header=T,))
}
spot_features <- spot_features[order(spot_features$position,spot_features$cell),]

# introduce spot_id
# the spot_id is useful to uniquely track back single spots to an image by position id + cell id + spot_id
spot_ids_per_cell <- vector()
for(h in 1:length(positions)){
  for(e in 1:length(as.numeric(table(spot_features[which(spot_features$position==h),]$cell)))){
    spot_ids_per_cell <- c(spot_ids_per_cell,1:as.numeric(table(spot_features[which(spot_features$position==h),]$cell))[e])
  }
}

# add spot_id column with values to spot_features data frame
# cells with no segmented spot temporarily get the spot_id 1
spot_features_with_ids <- cbind(spot_features,spot_ids_per_cell)
colnames(spot_features_with_ids) <- c(colnames(spot_features_with_ids)[1:length(colnames(spot_features))],"spot_id")

# re-order table to have spot_id as third column
spot_features_with_ids <- spot_features_with_ids[c("position","cell","spot_id",colnames(spot_features)[3:length(colnames(spot_features))])]

# check if there are cells with no segmented spots
# this is the case if the spot x-position has the value 0 (instead cells with only one segmented spot will have a value > 0)
# if yes, set the spot_id of these cells with no segmented spots to 0
if(0 %in% spot_features_with_ids$m.cx==TRUE){spot_features_with_ids[which(spot_features_with_ids$spot_id==1 & spot_features_with_ids$m.cx==0),]$spot_id <- 0}

# export spot_features_with_ids to results file, delete temporary spot_shape_features folder
spot_features_with_ids <- as.data.frame(spot_features_with_ids)
write.table(file=paste0(results_folder,"/",date_time, "_single_spot_shape_features_and_intensities.csv"), spot_features_with_ids, sep = "\t", append = F, row.names = F, quote = F)
if (file.exists(spot_shape_folder)){unlink(spot_shape_folder, recursive = T)}

##################################################################################################################################################
################################################_____Summarize and write NUCLEUS results____######################################################

# results/parameters per analyzed position e.g. type of z-projection, estimates for background intensity, best slices range etc.
# summarize types of projection used for each channel
nuc_ch_proj_type <- rep(projection_type[1],length(positions))
spot_ch_proj_type <- rep(projection_type[2],length(positions))
ch1_proj_type <- rep(projection_type[3],length(positions))
ch2_proj_type <- rep(projection_type[4],length(positions))
pos_results <- cbind(c(as.data.frame(table(pos[which(pos!=0)]))$Var1),positions,slices_range_used,nuc_ch_proj_type,spot_ch_proj_type,ch1_proj_type,ch2_proj_type,nuc_ch_bgmask_mean_intensity,nuc_ch_median_intensity,spot_seg_ch_bgmask_mean_intensity, spot_seg_ch_median_intensity,ch1_bgmask_mean_intensity,ch1_median_intensity,ch2_bgmask_mean_intensity,ch2_median_intensity)
pos_results <- as.data.frame(pos_results)
colnames(pos_results)[1] <- "position"
colnames(pos_results)[2] <- "position_name"
write.table(file=paste0(results_folder,"/",date_time, "_position_results.csv"), pos_results, sep = "\t", append = F, row.names = F, quote = F)

# results/parameters per analyzed cell e.g. mean, median and sd of nuclear intensities, mean intensities, areas and sd of spot and npmasks, etc.
int_area_res <- cbind(pos[1:(cnt-1)],nucleus_features$cell,nucleus_features$s.area, area_all_spotmask[1:(cnt-1)],area_npmask[1:(cnt-1)], num_spots[1:(cnt-1)],
                      mean_nuclear_nuc_ch_intensity[1:(cnt-1)], median_nuclear_nuc_ch_intensity[1:(cnt-1)],sd_nuclear_nuc_ch_intensity[1:(cnt-1)],mean_nuc_ch_all_spotmask_intensity[1:(cnt-1)],sd_nuc_ch_all_spotmask_intensity[1:(cnt-1)],mean_nuc_ch_npmask_intensity[1:(cnt-1)],sd_nuc_ch_npmask_intensity[1:(cnt-1)],
                      mean_nuclear_spot_ch_intensity[1:(cnt-1)], median_nuclear_spot_ch_intensity[1:(cnt-1)], sd_nuclear_spot_ch_intensity[1:(cnt-1)],mean_spot_ch_all_spotmask_intensity[1:(cnt-1)],sd_spot_ch_all_spotmask_intensity[1:(cnt-1)],mean_spot_ch_npmask_intensity[1:(cnt-1)],sd_spot_ch_npmask_intensity[1:(cnt-1)],
                      mean_nuclear_ch1_intensity[1:(cnt-1)], median_nuclear_ch1_intensity[1:(cnt-1)], sd_nuclear_ch1_intensity[1:(cnt-1)],mean_ch1_all_spotmask_intensity[1:(cnt-1)],sd_ch1_all_spotmask_intensity[1:(cnt-1)],mean_ch1_npmask_intensity[1:(cnt-1)],sd_ch1_npmask_intensity[1:(cnt-1)],  
                      mean_nuclear_ch2_intensity[1:(cnt-1)], median_nuclear_ch2_intensity[1:(cnt-1)], sd_nuclear_ch2_intensity[1:(cnt-1)],mean_ch2_all_spotmask_intensity[1:(cnt-1)],sd_ch2_all_spotmask_intensity[1:(cnt-1)],mean_ch2_npmask_intensity[1:(cnt-1)],sd_ch2_npmask_intensity[1:(cnt-1)])
int_area_res <- data.frame(int_area_res)
colnames(int_area_res) <- c("position","cell", "area_nucmask", "area_all_spotmask", "area_npmask","num_spots",
                            "mean_nuclear_nuc_ch_intensity", "median_nuclear_nuc_ch_intensity", "sd_nuclear_nuc_ch_intensity", "mean_nuc_ch_all_spotmask_intensity","sd_nuc_ch_all_spotmask_intensity","mean_nuc_ch_npmask_intensity","sd_nuc_ch_npmask_intensity",
                            "mean_nuclear_spot_ch_intensity", "median_nuclear_spot_ch_intensity", "sd_nuclear_spot_ch_intensity","mean_spot_ch_all_spotmask_intensity","sd_spot_ch_all_spotmask_intensity","mean_spot_ch_npmask_intensity","sd_spot_ch_npmask_intensity",
                            "mean_nuclear_ch1_intensity", "median_nuclear_ch1_intensity", "sd_nuclear_ch1_intensity","mean_ch1_all_spotmask_intensity","sd_ch1_all_spotmask_intensity","mean_ch1_npmask_intensity","sd_ch1_npmask_intensity",  
                            "mean_nuclear_ch2_intensity", "median_nuclear_ch2_intensity", "sd_nuclear_ch2_intensity","mean_ch2_all_spotmask_intensity","sd_ch2_all_spotmask_intensity","mean_ch2_npmask_intensity","sd_ch2_npmask_intensity")
                    
write.table(file=paste0(results_folder,"/",date_time, "_nucleus_spotmask_npmask_area_intensity_results.csv"), int_area_res, sep = "\t", append = F, row.names = F, quote = F)


##################################################################################################################################################
################################################_____ Collect CYTOPLASM results____###############################################################

#collect, order and summarize single spot shape, xy-position and intensity features for cytoplasmic spots
spot_shape_folder_cp <- paste0(results_folder,"/",list.files(results_folder)[grep("spot_shape_features_cytoplasm",list.files(results_folder))])
spot_shape_files_cp <- list.files(spot_shape_folder_cp)
spot_features_cp <- data.frame()
for(rr in 1:length(spot_shape_files_cp)){
  spot_features_cp <- rbind(spot_features_cp,read.table(file=paste0(spot_shape_folder_cp,"/",spot_shape_files_cp[rr]),header=T,))
}
spot_features_cp <- spot_features_cp[order(spot_features_cp$position,spot_features_cp$cell),]

# introduce spot_id
# the spot_id is useful to uniquely track back single spots to an image by position id + cell id + spot_id
spot_ids_per_cell_cp <- vector()
for(hh in 1:length(positions)){
  for(ee in 1:length(as.numeric(table(spot_features_cp[which(spot_features_cp$position==hh),]$cell)))){
    spot_ids_per_cell_cp <- c(spot_ids_per_cell_cp,1:as.numeric(table(spot_features_cp[which(spot_features_cp$position==hh),]$cell))[ee])
  }
}

# add spot_id column with values to spot_features_cp data frame
# cells with no segmented spot temporarily get the spot_id 1
spot_features_with_ids_cp <- cbind(spot_features_cp,spot_id=spot_ids_per_cell_cp)

# re-order table to have spot_id as third column
spot_features_with_ids_cp <- spot_features_with_ids_cp[c("position","cell","spot_id",colnames(spot_features_cp)[3:length(colnames(spot_features_cp))])]

# check if there are cells with no segmented spots
# this is the case if the spot x-position has the value 0 (instead cells with only one segmented spot will have a value > 0)
# if yes, set the spot_id of these cells with no segmented spots to 0
if(0 %in% spot_features_with_ids_cp$m.cx==TRUE){spot_features_with_ids_cp[which(spot_features_with_ids_cp$spot_id==1 & spot_features_with_ids_cp$m.cx==0),]$spot_id <- 0}

# export spot_features_with_ids_cp to results file, delete temporary spot_shape_features_cytoplasm folder
spot_features_with_ids_cp <- as.data.frame(spot_features_with_ids_cp)
write.table(file=paste0(results_folder,"/",date_time, "_single_spot_shape_features_and_intensities_cytoplasm.csv"), spot_features_with_ids_cp, sep = "\t", append = F, row.names = F, quote = F)
if (file.exists(spot_shape_folder_cp)){unlink(spot_shape_folder_cp, recursive = T)}

##################################################################################################################################################
################################################_____Summarize and write CYTOPLASM results____####################################################

# results/parameters per analyzed cell e.g. mean, median and sd of nuclear intensities, mean intensities, areas and sd of spot and npmasks, etc.
int_area_res_cp <- cbind(pos[1:(cnt_cp-1)],nucleus_features$cell,cp_area[1:(cnt_cp-1)], area_all_spotmask_cp[1:(cnt_cp-1)],area_spotfreemask_cp[1:(cnt_cp-1)], num_spots_cp[1:(cnt_cp-1)],
                         mean_cp_nuc_ch_intensity[1:(cnt_cp-1)], median_cp_nuc_ch_intensity[1:(cnt_cp-1)],sd_cp_nuc_ch_intensity[1:(cnt_cp-1)],mean_nuc_ch_all_cpspotmask_intensity[1:(cnt_cp-1)],sd_nuc_ch_all_cpspotmask_intensity[1:(cnt_cp-1)],mean_nuc_ch_spotfree_cp_intensity[1:(cnt_cp-1)],sd_nuc_ch_spotfree_cp_intensity[1:(cnt_cp-1)],
                         mean_cp_spot_ch_intensity[1:(cnt_cp-1)], median_cp_spot_ch_intensity[1:(cnt_cp-1)], sd_cp_spot_ch_intensity[1:(cnt_cp-1)],mean_spot_ch_all_cpspotmask_intensity[1:(cnt_cp-1)],sd_spot_ch_all_cpspotmask_intensity[1:(cnt_cp-1)],mean_spot_ch_spotfree_cp_intensity[1:(cnt_cp-1)],sd_spot_ch_spotfree_cp_intensity[1:(cnt_cp-1)],
                         mean_cp_ch1_intensity[1:(cnt_cp-1)], median_cp_ch1_intensity[1:(cnt_cp-1)], sd_cp_ch1_intensity[1:(cnt_cp-1)],mean_ch1_all_cpspotmask_intensity[1:(cnt_cp-1)],sd_ch1_all_cpspotmask_intensity[1:(cnt_cp-1)],mean_ch1_spotfree_cp_intensity[1:(cnt_cp-1)],sd_ch1_spotfree_cp_intensity[1:(cnt_cp-1)],  
                         mean_cp_ch2_intensity[1:(cnt_cp-1)], median_cp_ch2_intensity[1:(cnt_cp-1)], sd_cp_ch2_intensity[1:(cnt_cp-1)],mean_ch2_all_cpspotmask_intensity[1:(cnt_cp-1)],sd_ch2_all_cpspotmask_intensity[1:(cnt_cp-1)],mean_ch2_spotfree_cp_intensity[1:(cnt_cp-1)],sd_ch2_spotfree_cp_intensity[1:(cnt_cp-1)])
int_area_res_cp <- data.frame(int_area_res_cp)
colnames(int_area_res_cp) <- c("position","cell", "cp_area", "area_all_spotmask_cp", "area_spotfreemask_cp","num_spots_cp",
                            "mean_cp_nuc_ch_intensity", "median_cp_nuc_ch_intensity", "sd_cp_nuc_ch_intensity", "mean_nuc_ch_all_cpspotmask_intensity","sd_nuc_ch_all_cpspotmask_intensity","mean_nuc_ch_spotfree_cp_intensity","sd_nuc_ch_spotfree_cp_intensity",
                            "mean_cp_spot_ch_intensity", "median_cp_spot_ch_intensity", "sd_cp_spot_ch_intensity","mean_spot_ch_all_cpspotmask_intensity","sd_spot_ch_all_cpspotmask_intensity","mean_spot_ch_spotfree_cp_intensity","sd_spot_ch_spotfree_cp_intensity",
                            "mean_cp_ch1_intensity", "median_cp_ch1_intensity", "sd_cp_ch1_intensity","mean_ch1_all_cpspotmask_intensity","sd_ch1_all_cpspotmask_intensity","mean_ch1_spotfree_cp_intensity","sd_ch1_spotfree_cp_intensity",  
                            "mean_cp_ch2_intensity", "median_cp_ch2_intensity", "sd_cp_ch2_intensity","mean_ch2_all_cpspotmask_intensity","sd_ch2_all_cpspotmask_intensity","mean_ch2_spotfree_cp_intensity","sd_ch2_spotfree_cp_intensity")

write.table(file=paste0(results_folder,"/",date_time, "_cytoplasm_spotmask_spotfreemask_area_intensity_results.csv"), int_area_res_cp, sep = "\t", append = F, row.names = F, quote = F)


# save R workspace after the run for easier access to data & potential troubleshooting
save.image(file=paste0(results_folder,"/",date_time, "_R_workspace_after_run_completion.RData"))

# clear workspace, keeping only the variable folderList which is required for looping
rm(list = setdiff(ls(),"folderList"))

# (re-)load packages and functions
library(EBImage)
library(abind)
source("makeNucMask_v1.R")
source("makeSpotMask_v1.R")
source("findBestSlices_v2.R")
source("makeCytoMask_v1.R")
source("matchCellMasks_v5.R")


# end of folderList loop
}
