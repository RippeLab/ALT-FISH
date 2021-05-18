############################################ Telosegment toolkit MULTIRUN MODE ############################################################
# Thresholding-based segmentation of cell nuclei and spots within nuclei 
# Quantification of number, intensity, area and shape features of nuclei and single spots
# MULTIRUN MODE: This allows to loop over several specified .tif data folders sequentially using the same script version

# load libraries and functions
library(EBImage)
library(abind)
source("makeNucMask_v1.R")
source("makeSpotMask_v1.R")
source("findBestSlices_v2.R")

####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################

# specify a folder list as character vector (vector of file paths)
# target folders should contain the .tif files to be analyzed

folderList <- c("/home/user/folder/tif_stack_dataset_1/",
                "/home/user/folder/tif_stack_dataset_2/")

###### loop through folder list ######
for(f in 1:length(folderList)){

# specify location of the currently opened Telosegment_toolkit R script version
# this is required for saving a .txt copy of this script in the Telosegment_toolkit results folder for documentation
scriptversion <- list.files(getwd(),pattern = "Telosegment_toolkit_multirun.R")

# specify location of the .tif input files 
# .tif files are required to be multi-channel z-stacks (2-4 channels)
folder <- folderList[f]

# specify channel names, order & assign to desired operation
# set non-existing channels to 0 if less than 4 channels are needed
# non-existing channels will be replaced by placeholders (pixel value 0) at the level of projections to simplify downstream processing
# fill in the channel definition in the order of the image stack, avoid putting empty channels "in between"
# e.g if image stack is DAPI/telo/marker1 use channels <- c("dapi","telo","marker1",0) and not c("dapi","telo",0,"marker1") 
channels <- c("telo","dapi",0,0)

# define channel used for nucleus segmentation
nuc_ch <- "dapi"
# define channel used for spot segmentation within the nuclei
spot_seg_ch <- "telo"
# define additional marker for quantification 1, if not used set to 0
ch1 <- 0
# define additional marker for quantification 2, if not used set to 0
ch2 <- 0

#number of z-slices per channel
zperch <- 51

# remove incomplete nuclei at the image borders? "y" = yes, "n" = no
# default setting is "y"
remove_bordercells <- "y"

# adjustable adaptive thresholding offset for nucleus segmentation, default = 0.0008, see EBImage function thresh() for details
nuc_threshold_offset <- 0.0008

# thresholding parameter for spot segmentation relative to the nuclear median in the spot channel
# default value is 2.5 meaning that pixel groups with intensities above 2.5*median(nuclear intensity in the spot channel) are handed over to the makeSpotMask() function
spot_cutoff <- 1.8

# size-filter for nucleus segmentation (pixels)
# only nucleus masks with areas larger than this value will be returned by makeNucMask()
# default is 4000
nucSize_cutoff <- 4000 

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
dir.create(paste0(folder,date_time,"_Telosegment_toolkit_results"))
results_folder <- paste0(folder,date_time,"_Telosegment_toolkit_results")
dir.create(paste0(results_folder,"/",date_time,"_nuc_shape_features"))
dir.create(paste0(results_folder,"/",date_time,"_spot_shape_features"))
  
# document adjustable segmentation parameters used in this run
parameters_df <- data.frame(channels[1], channels[2], channels[3], channels[4], nuc_ch, spot_seg_ch, ch1, ch2, remove_bordercells, spot_cutoff,nuc_threshold_offset, nucSize_cutoff, mn, passfilter)
colnames(parameters_df) <- c("channel_1","channel_2","channel_3","channel_4","nuc_ch","spot_seg_ch", "ch1", "ch2", "remove_bordercells", "spot_cutoff","nuc_threshold_offset", "nucSize_cutoff", "mn","findBestSlices_passfilter")
write.table(file=paste0(results_folder,"/",date_time,"_parameters_used_in_run.csv"), parameters_df, sep = "\t", append = F, row.names = F, quote = F)

# generate a position list to document the files that were processed in this run
#positions <- list.files(path=folder, pattern="tif")
positions <- list.files(path=folder, pattern="tif")

# write position list
write.table(file=paste0(results_folder,"/",date_time,"_positions.csv"), cbind(seq(1,length(positions),by=1) ,positions), sep = "\t", append = F, row.names = F, quote = F)

# export the currently run script as text file for logging purposes 
write(readLines(scriptversion),file=paste0(results_folder,"/",date_time,"_Telosegment_toolkit_script_used_in_this_run.R"))

##################################################################################################################################
##################################################################################################################################
#                                                     INITIALIZE VARIABLES
##################################################################################################################################
##################################################################################################################################
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
# total position counter vector
pos <- rep(0, mn*length(positions))
# vector for storing BestSlices ranges determined for each position
slices_range_used <- rep(0, length(positions))

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
  
  ################################### Segmentation step 1: Nuclei segmentation ##################################################
  
  # generate mask for all nuclei in the current position
  # makeNucMask() uses a gauss-blurred best slices z-projection of the nuc_channel and the above-defined adjustable parameters
  nucmask <- makeNucMask(gblur(nuc_data_projection, sigma=1), nuc_threshold_offset, nucSize_cutoff, remove_bordercells)
  
  # if no best_slices_range could be determined for this position and was thus set to 1:zperch or no nucleus was found:
  # create a dummy image with one squared "nucleus" of value 1 in the top left corner (40000 px area) and background intensity 0
  # this nuc_data_projection dummy will be processed as all other images, thereby maintaining results for all positions and marking out-of-focus positions for later filtering
  if(length(best_slices_range)==zperch | max(nucmask)==0){
    dummy_image <- Image(0,dim(nuc_data)[1:2])
    dummy_image[50:249,50:249] <- 1
    nucmask <- makeNucMask(dummy_image, nuc_threshold_offset, nucSize_cutoff, remove_bordercells)
  }
  
  # initialize variables for storing whole-image spot and nucleoplasm masks (npmasks)
  # whole-image spot and npmasks are later only used for visualization of segmentation results in RGB images (see below)
  allspotmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  allnpmasks <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  
  # retrieve indices of single nucleus masks for looping
  rois <- as.numeric(names(table(nucmask)[2:length(table(nucmask))]))

  ######################################  Quantify whole-image background estimates ################################################
  
  # This step measures intensities of whole projections in all channels that are potentially informative for image background estimation
  # 1. computes the median intensity of the whole image projection
  # 2. computes the mean intensity in the negative of the nuclear masks
  # the first value is a robust estimate for the background of the image, however will be biased for positions with high cell density 
  # the second value is more precisely the extra-nuclear background, however will be biased by border nuclei contribution when choosing remove_bordercells <- "y"
  
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

    ################################### Segmentation step 2: Spot segmentation ##################################################

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
    ##############################################    SINGLE SPOT LOOP s in 1: length(rois_spots)  ############################################
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
  
  # export projections with generated nucmask and all spotmasks as .tif for potential re-analysis of data
  # empty channels will contain all pixels with value 0
  writeImage(combine(nuc_data_projection,nucmask, spot_data_best_projection, allspotmasks, ch1_data_best_projection, ch2_data_best_projection), bits.per.sample=16, files=paste0(results_folder,"/",gsub(".tif","_BestSlices_projections_with_masks_binary.tif",positions[i])))
  
  ############################### Create images for visualization of segmentation results #########################################
  
  # RGB images scaled to the maximum image intensity for quick visual inspection of results
  # write RGB image overlays of masks onto projections of current position
  # projections are divided by their maximum pixel value for scaling
  best_nuc_projection_with_masks <- toRGB(nuc_data_projection/max(nuc_data_projection))
  best_spot_projection_with_masks_nonpmasks <- toRGB(spot_data_best_projection/max(spot_data_best_projection))
  best_spot_projection_with_masks <- toRGB(spot_data_best_projection/max(spot_data_best_projection))
  best_ch1_projection_with_masks <- toRGB(ch1_data_best_projection/max(ch1_data_best_projection))
  best_ch2_projection_with_masks <- toRGB(ch2_data_best_projection/max(ch2_data_best_projection))
  
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
  writeImage(combine(best_nuc_projection_with_masks,best_spot_projection_with_masks,best_spot_projection_with_masks_nonpmasks,best_ch1_projection_with_masks,best_ch2_projection_with_masks), bits.per.sample=8, files=paste0(results_folder,"/",gsub(".tif","_BestSlices_projections_with_masks_scaledRGB.tif",positions[i])))
  
  ##########
  # write out non-scaled 16-bit grayscale images with mask outlines as separate channels for creating customizable composites etc. in FIJI
  # the mask outlines are set to pixel value 65636 (nucleoplasm mask is not drawn as outline)
  blackimage <- Image(0,dim=c(dim(nuc_data)[1],dim(nuc_data)[2]))
  blackimage_nucmask <- channel(paintObjects(nucmask, blackimage, opac=c(1,0), col=c("white","white")), mode = "gray")
  blackimage_allspotmasks <- channel(paintObjects(allspotmasks, blackimage, opac=c(1,0), col=c("white","white")), mode ="gray")
  blackimage_allnpmasks <- channel(paintObjects(allnpmasks, blackimage, opac=c(0,1), col=c("white","white")), mode="gray")
  writeImage(combine(nuc_data_projection,blackimage_nucmask,spot_data_best_projection,blackimage_allspotmasks,blackimage_allnpmasks,ch1_data_best_projection,ch2_data_best_projection), bits.per.sample=16, files=paste0(results_folder,"/",gsub(".tif","_BestSlices_projections_with_masks_outlines.tif",positions[i])))
}
#######################################################################################################################################
######################################################## END POSITION LOOP i ##########################################################
#######################################################################################################################################

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
spot_shape_folder <- paste0(results_folder,"/",list.files(results_folder)[grep("spot_shape_features",list.files(results_folder))])
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

#######################################################################################################################################
################################################_____Summarize and write main results____##############################################

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

# save R workspace after the run for easier access to data & potential troubleshooting
#save.image(file=paste0(results_folder,"/",date_time, "_R_workspace_after_run_completion.RData"))

# clear workspace, keeping only the variable folderList which is required for looping
rm(list = setdiff(ls(),"folderList"))

# (re-)load packages and functions
library(EBImage)
library(abind)
source("makeNucMask_v1.R")
source("makeSpotMask_v1.R")
source("findBestSlices_v2.R")

# end of folderList loop
}