####################################################################################################################################################################
################################ Read and filter data from Telosegment toolkit and Telosegment cytosegment toolkit results folders #################################
####################################################################################################################################################################
# This script retrieves all .csv files in folders named "Telosegment_toolkit_results" or "Telosegment_cytosegment_toolkit_results" within a specified main folder.
# The files are read in and organized into a list containing the path to the result folders as element name (for easier access).
# Each list element itself is a list of data frames with the content of the csv files.
# These sublists contain the file name as element names.
# Telosegment_toolkit_results list elements contain 6 csv tables/elements
# Telosegment_cytosegment toolkit_results list elements contain 8 csv tables/elements (additional data for cytoplasm)

# Once the data is retrieved, it is filtered in two steps
# 1. Removal of out of focus positions and calculation of the nucleus shape quality scores (nsqs) for each nucleus.
# 2. Filter out data corresponding to nuclei (spot, nucleus and corresp. cytoplasm data) with a nsqs-score greater than a defined nsqs_cutoff (default = 3).
# Data is exported as Rdata object in three versions:
# * raw data
# * data filtered for out-of-focus positions with computed nsqs score added
# * data filtered for out-of-focus positions with computed nsqs score added, and filtered for a user-defined nsqs_cutoff

# This data formats can afterwards be imported for plotting.

############################################################################################################################################################
############################ Generate list containing paths to and additional information about the results folders (completeness, run number) #############
############################################################################################################################################################
# browse all subfolders, identify folders named "_toolkit_results" and check folder completeness
# if multiple results folders are found within the same folder (e.g. multiple runs conducted on the same dataset)
# duplicate folders receive a run_id with the newest (datetime stamp) results folder having the highest run_id number
##########################################################################################################################################################
# specify location of the main folder 
# Ideally the results folders are located in the corresponding condition-sorted imaging data folders e.g. probe/cell_line/replicate1/...
all_results_folder <- "/home/user/folder_containing_Telosegment_results_folders"

# specify a pre-fix used for naming during data export/storage
data_prefix <- "Experiment_Set1_"

##########################################################################################################################################################
##############################################   Identify results folder paths and completeness of runs   ################################################
##########################################################################################################################################################
# fetch paths of all Telosegment_toolkit_results or Telosegment_cytosegment_toolkit_results folders contained within the specified main folder 
ttkit_results_folderpaths <- list.dirs(all_results_folder)[grep("_toolkit_results",list.dirs(all_results_folder))]
# the list is further filtered removing Telosegment_toolkit_results subfolder paths that are fetched for aborted Telosegment toolkit runs (e.g. /nuc_shape_features)
if(length(grep("results/",ttkit_results_folderpaths))==0){
  ttkit_results_folderpaths_filtered <- ttkit_results_folderpaths
}else{ttkit_results_folderpaths_filtered <- ttkit_results_folderpaths[-grep("results/",ttkit_results_folderpaths)]}

# check if the identified folders are complete (contain all expected results files) and annotate accordingly by creating a new data.frame
ttkit_results_folderpaths_filtered <- cbind(ttkit_results_folderpaths_filtered,c(1:length(ttkit_results_folderpaths_filtered)))
colnames(ttkit_results_folderpaths_filtered) <- c("folderpath","status")

for(v in 1:length(ttkit_results_folderpaths_filtered[,"folderpath"])){
  if (length(list.files(ttkit_results_folderpaths_filtered[v,"folderpath"], pattern=".csv"))>=6 & length(list.files(ttkit_results_folderpaths_filtered[v,"folderpath"], pattern="used_in_this_run.R"))==1){
    ttkit_results_folderpaths_filtered[v,"status"] <- "complete"
  }else{ttkit_results_folderpaths_filtered[v,"status"]<- "incomplete"}
}
ttkit_results_folderpaths_annotated <- as.data.frame(ttkit_results_folderpaths_filtered)

# only keep complete results folders
ttkit_results_folderpaths_complete <- ttkit_results_folderpaths_annotated[which(ttkit_results_folderpaths_annotated$status=="complete"),]

# Assign numeric datetime stamp variable based on folder datetime stamp
ttkit_results_folderpaths_complete <- cbind(ttkit_results_folderpaths_complete,rep(0,length(ttkit_results_folderpaths_complete$folderpath)))
colnames(ttkit_results_folderpaths_complete) <- c(colnames(ttkit_results_folderpaths_complete)[1:2],"numeric_datetime")
folder_str <- vector()
datetime_str <- vector()
pathtofolder <- vector()
for(q in 1:length(ttkit_results_folderpaths_complete$folderpath)){
  folder_str[q] <- substr(as.character(ttkit_results_folderpaths_complete[q,"folderpath"]),nchar(as.character(ttkit_results_folderpaths_complete[q,"folderpath"]))-46,nchar(as.character(ttkit_results_folderpaths_complete[q,"folderpath"])))
  pathtofolder[q] <- substr(as.character(ttkit_results_folderpaths_complete[q,"folderpath"]),1,nchar(as.character(ttkit_results_folderpaths_complete[q,"folderpath"]))-47)
  datetime_str[q] <- substr(folder_str[q],1,19)
  ttkit_results_folderpaths_complete[q,"numeric_datetime"] <- gsub("_","",gsub("-","",datetime_str[q]))
}

# Check if there is more than one Telosegment_toolkit_results folder per condition subfolder (for examples when multiple runs were conducted for one sample set with different parameters)
# If this is the case, different runs on the same dataset (condition subfolder) will be assigned a run_id
# If there is only one run per dataset, it will receive the run_id 1
# If there are multiple runs, the newest run (datetime stamp) will receive the highes id
if (length(duplicated(pathtofolder))!=0){
  print("More than one Telosegment toolkit results folder was found per condition subfolder. Duplicate or more Telosegment toolkit runs per condition subfolder will be assigned a run_id for subsequent filtering")
  duplicated_folder_positions <- vector()
  run_id_per_condition_folder <- as.vector(rep(0,length(ttkit_results_folderpaths_complete$folderpath)))
  ttkit_results_folderpaths_final <- ttkit_results_folderpaths_complete[order(ttkit_results_folderpaths_complete$numeric_datetime),]
  for(y in 1:length(which(duplicated(pathtofolder)))){
    duplicated_folder_positions <-  grep(pathtofolder[which(duplicated(pathtofolder))[y]],ttkit_results_folderpaths_final$folderpath)
    for(e in 1:length(duplicated_folder_positions)){
      run_id_per_condition_folder[duplicated_folder_positions[e]] <- e
    }
  }
  run_id_per_condition_folder[which(run_id_per_condition_folder==0)] <- 1
  ttkit_results_folderpaths_final <- cbind(ttkit_results_folderpaths_final,run_id_per_condition_folder)
  colnames(ttkit_results_folderpaths_final) <- c(colnames(ttkit_results_folderpaths_final)[1:3],"run_id")
}else{
  run_id_per_condition_folder <- rep(1,length(ttkit_results_folderpaths_complete$folderpath))
  ttkit_results_folderpaths_final <- ttkit_results_folderpaths_complete[order(ttkit_results_folderpaths_complete$numeric_datetime),]
  ttkit_results_folderpaths_final <- cbind(ttkit_results_folderpaths_final,run_id_per_condition_folder)
  colnames(ttkit_results_folderpaths_final) <- c(colnames(ttkit_results_folderpaths_final)[1:3],"run_id")
}
##########################################################################################################################################################
##########################################################################################################################################################

# NOTE: The ttkit_results_folderpaths_final object contains all the information about the retrieved results folder (paths, completness, run number)
# and can be examined for overview or manual pre-filtering, if desired.

##########################################################################################################################################################
############################      Generate main data list containing all csv-file results as lists of data frames     ####################################
##########################################################################################################################################################
# generate vector list containing all .csv filenames associated to a given results folder path
csvlist <- list()
for(aa in 1:length(ttkit_results_folderpaths_final$folderpath)){
  csvlist[[aa]]<-list.files(as.character(ttkit_results_folderpaths_final$folderpath[aa]), pattern=".csv")
}
# name list elements with associated Telosegment_toolkit_results folder path
names(csvlist) <- as.character(ttkit_results_folderpaths_final$folderpath)

# generate a main list containing all csv data as data frames
# the elements of the main list correspond to one results folder dataset
# one element is itself a list of data frames with different dimensions (6 data frames: nucleus_shape_features.csv, positions.csv, etc...)

# initialize main list variable
list_of_csvdatalists <- list()

# loop over all paths in of csvlist and read in the corresponding dataset into a temporary csvdatalist (one for each Telosegment_toolkit_results folder)
# csvdatalist will be added to the main list list_of_csvdatalists
for(o in 1:length(ttkit_results_folderpaths_final$folderpath)){
  nuc_shape_features <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("nucleus_shape_features.csv",csvlist[[o]])]), header=T)
  nuc_area_int_results <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("nucleus_spotmask_npmask_area_intensity_results.csv",csvlist[[o]])]), header=T)
  run_parameters <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("parameters_used_in_run.csv",csvlist[[o]])]), header=T)
  position_results <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("position_results.csv",csvlist[[o]])]), header=T)
  positions <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("positions.csv",csvlist[[o]])]), header=T)
  single_spot_shape_int_results <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("single_spot_shape_features_and_intensities.csv",csvlist[[o]])]), header=T)
  
  # in case of Telosegment_cytosegment_toolkit_results, also retrieve the additional csv files containing the cytoplasm quantifications
  if(length(csvlist[[o]][grep("cytoplasm_spotmask_spotfreemask_area_intensity_results.csv",csvlist[[o]])])==1 & length(csvlist[[o]][grep("single_spot_shape_features_and_intensities_cytoplasm.csv",csvlist[[o]])])==1){
    cp_area_int_results <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("cytoplasm_spotmask_spotfreemask_area_intensity_results.csv",csvlist[[o]])]), header=T)
    single_spot_shape_int_results_cp <- read.table(file=paste0(names(csvlist[o]),"/",csvlist[[o]][grep("single_spot_shape_features_and_intensities_cytoplasm.csv",csvlist[[o]])]), header=T)
    
    csvdatalist <- list(nuc_shape_features,nuc_area_int_results,run_parameters,position_results,positions,single_spot_shape_int_results,cp_area_int_results,single_spot_shape_int_results_cp)
    names(csvdatalist) <- c(csvlist[[o]][grep("nucleus_shape_features.csv",csvlist[[o]])],csvlist[[o]][grep("nucleus_spotmask_npmask_area_intensity_results.csv",csvlist[[o]])],csvlist[[o]][grep("parameters_used_in_run.csv",csvlist[[o]])],csvlist[[o]][grep("position_results.csv",csvlist[[o]])],csvlist[[o]][grep("positions.csv",csvlist[[o]])],csvlist[[o]][grep("single_spot_shape_features_and_intensities.csv",csvlist[[o]])],csvlist[[o]][grep("cytoplasm_spotmask_spotfreemask_area_intensity_results.csv",csvlist[[o]])],csvlist[[o]][grep("single_spot_shape_features_and_intensities_cytoplasm.csv",csvlist[[o]])])
    list_of_csvdatalists[[o]] <- csvdatalist
    csvdatalist <- NULL
  }else{
    csvdatalist <- list(nuc_shape_features,nuc_area_int_results,run_parameters,position_results,positions,single_spot_shape_int_results)
    names(csvdatalist) <- c(csvlist[[o]][grep("nucleus_shape_features.csv",csvlist[[o]])],csvlist[[o]][grep("nucleus_spotmask_npmask_area_intensity_results.csv",csvlist[[o]])],csvlist[[o]][grep("parameters_used_in_run.csv",csvlist[[o]])],csvlist[[o]][grep("position_results.csv",csvlist[[o]])],csvlist[[o]][grep("positions.csv",csvlist[[o]])],csvlist[[o]][grep("single_spot_shape_features_and_intensities.csv",csvlist[[o]])])
    
    list_of_csvdatalists[[o]] <- csvdatalist
    csvdatalist <- NULL
  }
}
# naming of elements in list_of_csvdatalists for easier access (since the path contains information about the experimental condition)
# the unnamed list_of_csvdatalists is preserved 
list_of_csvdatalists_named <- list_of_csvdatalists
names(list_of_csvdatalists_named) <- as.character(ttkit_results_folderpaths_final$folderpath)

##########################################################################################################################################################
#########################################################     SAVE RAW DATA (1/3)      ###################################################################

# NOTE: The list_of_csvdatalists_named object contains all the raw data 
# save to file
saveRDS(list_of_csvdatalists_named,paste0(data_prefix,"raw_data.rds"))

############################################################################################################################################################
#######################################               DATA FILTERING STEP 1 (out-of-focus positions)               #########################################
############################################################################################################################################################
# filter out data from out-of-fous (oof)/empty image positions 
# such positions are marked by a squared dummy nucleus in the upper left corner with the area of 40000px
# this is taken as a criterion to remove data rows in the relevant data.frames of list_of_csvdatalists_named
# for Telosegment_cytosegment_toolkit results, the same dummy nuclei can be used for filtering
# for filtering the cytoplasm_spotmask_spotfreemask_area_intensity_results.csv data, positions with dummy nuclei receive a corresponding dummy cytoplasm with the area value of 17600
# The % of oof positions filtered out per data set is printed for overview.

list_of_csvdatalists_named_filtered <- list_of_csvdatalists_named
for(c in 1:length(list_of_csvdatalists_named)){
  # if a dataset does not contain any dummy nuclei, proceed to the next
  if (length(which(list_of_csvdatalists_named_filtered[[c]][[1]]$s.area==40000))!=0){
    
    # remove affected data.frame rows in single_spot_shape_features_and_intensities.csv
    # this is done by finding the position and cell id of 40000px-area nuclei and removing spot data rows with matching position and cell id
    for(cc in 1:length(list_of_csvdatalists_named_filtered[[c]][[2]][which(list_of_csvdatalists_named_filtered[[c]][[2]]$area_nucmask==40000),][,1])){
      list_of_csvdatalists_named_filtered[[c]][[6]] <- list_of_csvdatalists_named_filtered[[c]][[6]][-which(list_of_csvdatalists_named_filtered[[c]][[6]]$position==list_of_csvdatalists_named_filtered[[c]][[2]][which(list_of_csvdatalists_named_filtered[[c]][[2]]$area_nucmask==40000),][,"position"][cc] & list_of_csvdatalists_named_filtered[[c]][[6]]$cell==list_of_csvdatalists_named_filtered[[c]][[2]][which(list_of_csvdatalists_named_filtered[[c]][[2]]$area_nucmask==40000),][,"cell"][cc]),]
      ### for Telosegment_cytosegment_toolkit_results, also remove affected data.frame rows in single_spot_shape_features_and_intensities_cytoplasm.csv and cytoplasm_spotmask_spotfreemask_area_intensity_results.csv
      ### dummy cytoplasm masks have an area value of 17600
      if (length(list_of_csvdatalists_named[[c]])==8){
        list_of_csvdatalists_named_filtered[[c]][[8]] <- list_of_csvdatalists_named_filtered[[c]][[8]][-which(list_of_csvdatalists_named_filtered[[c]][[8]]$position==list_of_csvdatalists_named_filtered[[c]][[2]][which(list_of_csvdatalists_named_filtered[[c]][[2]]$area_nucmask==40000),][,"position"][cc] & list_of_csvdatalists_named_filtered[[c]][[8]]$cell==list_of_csvdatalists_named_filtered[[c]][[2]][which(list_of_csvdatalists_named_filtered[[c]][[2]]$area_nucmask==40000),][,"cell"][cc]),]
        list_of_csvdatalists_named_filtered[[c]][[7]] <- list_of_csvdatalists_named_filtered[[c]][[7]][-which(list_of_csvdatalists_named_filtered[[c]][[7]]$cp_area==17600),]
      }
    }
    ### for both Telosegment and Telosegment_cytosegment type of results
    # remove affected data.frame rows in _nucleus_shape_features.csv and _nucleus_spotmask_npmask_area_intensity_results.csv 
    list_of_csvdatalists_named_filtered[[c]][[1]] <- list_of_csvdatalists_named_filtered[[c]][[1]][-which(list_of_csvdatalists_named_filtered[[c]][[1]]$s.area==40000),]
    list_of_csvdatalists_named_filtered[[c]][[2]] <- list_of_csvdatalists_named_filtered[[c]][[2]][-which(list_of_csvdatalists_named_filtered[[c]][[2]]$area_nucmask==40000),]
    print(paste0(round((length(list_of_csvdatalists_named[[c]][[1]]$position)-length(list_of_csvdatalists_named_filtered[[c]][[1]]$position))/length(list_of_csvdatalists_named[[c]][[1]]$position)*100, digits = 1)," % dummy nuclei for dataset: ",names(list_of_csvdatalists_named_filtered[c])))
    
    # introduce nucleus shape quality score (nsqs) as column to nucleus_shape_features data frame (see below for details about nsqs)
    nsqs_numerator <- list_of_csvdatalists_named_filtered[[c]][[1]]$s.radius.max*list_of_csvdatalists_named_filtered[[c]][[1]]$s.perimeter*list_of_csvdatalists_named_filtered[[c]][[1]]$s.radius.sd
    nsqs_denominator <- list_of_csvdatalists_named_filtered[[c]][[1]]$s.radius.min*list_of_csvdatalists_named_filtered[[c]][[1]]$s.area
    nsqs <- nsqs_numerator/nsqs_denominator 
    list_of_csvdatalists_named_filtered[[c]][[1]] <- cbind(list_of_csvdatalists_named_filtered[[c]][[1]], nsqs)
    
  }else {
    print(paste0("0 % of dummy nuclei for dataset: ",names(list_of_csvdatalists_named_filtered[c])))
    # introduce nucleus shape quality score (nsqs) as column to nucleus_shape_features data frame (see below for details about nsqs)
    nsqs_numerator <- list_of_csvdatalists_named_filtered[[c]][[1]]$s.radius.max*list_of_csvdatalists_named_filtered[[c]][[1]]$s.perimeter*list_of_csvdatalists_named_filtered[[c]][[1]]$s.radius.sd
    nsqs_denominator <- list_of_csvdatalists_named_filtered[[c]][[1]]$s.radius.min*list_of_csvdatalists_named_filtered[[c]][[1]]$s.area
    nsqs <- nsqs_numerator/nsqs_denominator 
    list_of_csvdatalists_named_filtered[[c]][[1]] <- cbind(list_of_csvdatalists_named_filtered[[c]][[1]], nsqs)
  }
}
##########################################################################################################################################################
####################################################     SAVE OOF-FILTERED DATA (2/3)      ###############################################################

# NOTE: the list_of_csvdatalists_named_filtered object contains all the data after filtering out out-of-focus (oof) positions
# Moreover, the _nucleus_shape_features.csv data frame of each result data set receives an additional column containing the computed nucleus shape quality score (nsqs)
# save to file
saveRDS(list_of_csvdatalists_named_filtered,paste0(data_prefix,"oof_filtered_data.rds"))


############################################################################################################################################################
#######################################               DATA FILTERING STEP 2 (nsqs score)               #####################################################
############################################################################################################################################################
# filter out data coming from abberantly shaped nuclei, DAPI artifacts, attached nuclei (doublets, triplets) 
# this is done by using the "nucleus shape quality score" (nsqs) already calculated above and added to the _nucleus_shape_features table ($nsqs)
# nsqs is given by the term (s.radius.max*s.perimeter*s.radius.sd)/(s.radius.min*s.area)
# and becomes high for irregularly shaped nuclei, independent of their absolute area/size
# default nsqs_cutoff value is 3 (>3: bad quality, < 3: good quality)
# it may however be adjusted, depending on the cell type
# IMPORTANT: in some very rare cases, Telosegment_cytosegment_toolkit can produce very small (few pixels) nuclear masks as a result of the cytoplasm/nucleus
# segmentation matching process. Usually these nuclei have nsqs > 3, but occasionally this results in "NaN" for nsqs (because shape features cannot 
# be computed for e.g. only 1px). This is also accounted for in the filtering procedure below.

nsqs_cutoff <- 3
list_of_csvdatalists_named_filtered_nsqs <- list_of_csvdatalists_named_filtered

# initialize row index vectors for filtering
goodspotrows <- vector()
goodspotrows_cp <- vector()

# loop over each dataset
for(d in 1:length(list_of_csvdatalists_named_filtered_nsqs)){
  # remove all nuclei from _nucleus_spotmask_npmask_area_intensity_results that have nsqs > nsqs_cutoff
  list_of_csvdatalists_named_filtered_nsqs[[d]][[2]] <- list_of_csvdatalists_named_filtered_nsqs[[d]][[2]][which(list_of_csvdatalists_named_filtered_nsqs[[d]][[1]]$nsqs<nsqs_cutoff & list_of_csvdatalists_named_filtered_nsqs[[d]][[1]]$nsqs!="NaN"),]
  # remove all nuclei from _nuc_shape_features table that have nsqs > nsqs_cutoff
  list_of_csvdatalists_named_filtered_nsqs[[d]][[1]] <- list_of_csvdatalists_named_filtered_nsqs[[d]][[1]][which(list_of_csvdatalists_named_filtered_nsqs[[d]][[1]]$nsqs<nsqs_cutoff & list_of_csvdatalists_named_filtered_nsqs[[d]][[1]]$nsqs!="NaN"),]
  # remove all cytoplasm rows from cytoplasm_spotmask_spotfreemask_area_intensity_results.csv whose corresponding nucleus has nsqs > nsqs_cutoff
  if (length(list_of_csvdatalists_named_filtered_nsqs[[d]])==8){
    list_of_csvdatalists_named_filtered_nsqs[[d]][[7]] <- list_of_csvdatalists_named_filtered_nsqs[[d]][[7]][which(list_of_csvdatalists_named_filtered_nsqs[[d]][[1]]$nsqs<nsqs_cutoff & list_of_csvdatalists_named_filtered_nsqs[[d]][[1]]$nsqs!="NaN"),]
  }
  
  # only retain nuclei single_spot_shape_features_and_intensity rows with position and cell matching to the nsqs-filtered data
  for(dd in 1:length(list_of_csvdatalists_named_filtered_nsqs[[d]][[1]][,1])){
    # the vector goodspotrows documents the row indices of the current single_spot_shape_features_and_intensity data frame where position and cell id match the nsqs-filtered data (=rows to keep)
    goodspotrows <- c(goodspotrows, which(list_of_csvdatalists_named_filtered_nsqs[[d]][[6]][,"position"]==list_of_csvdatalists_named_filtered_nsqs[[d]][[1]][,"position"][dd] & list_of_csvdatalists_named_filtered_nsqs[[d]][[6]][,"cell"]==list_of_csvdatalists_named_filtered_nsqs[[d]][[1]][,"cell"][dd]))
    # do the same for cytoplasmic spots (Telosegment_cytosegment_toolkit_results)
    if (length(list_of_csvdatalists_named_filtered_nsqs[[d]])==8){
      goodspotrows_cp <- c(goodspotrows_cp, which(list_of_csvdatalists_named_filtered_nsqs[[d]][[8]][,"position"]==list_of_csvdatalists_named_filtered_nsqs[[d]][[1]][,"position"][dd] & list_of_csvdatalists_named_filtered_nsqs[[d]][[8]][,"cell"]==list_of_csvdatalists_named_filtered_nsqs[[d]][[1]][,"cell"][dd]))
    }
  }
  # subset the single_spot_shape_features_and_intensity data frame, keeping only "good" rows
  list_of_csvdatalists_named_filtered_nsqs[[d]][[6]] <- list_of_csvdatalists_named_filtered_nsqs[[d]][[6]][goodspotrows,]
  # empty goodspotrows vector
  goodspotrows <- NULL
  if (length(list_of_csvdatalists_named_filtered_nsqs[[d]])==8){
    # do the same for the single_spot_shape_features_and_intensities_cytoplasm data frames
    list_of_csvdatalists_named_filtered_nsqs[[d]][[8]] <- list_of_csvdatalists_named_filtered_nsqs[[d]][[8]][goodspotrows_cp,]
  }
  # empty goodspotrows_cp vector
  goodspotrows_cp <- NULL
}

##########################################################################################################################################################
####################################################     SAVE OOF-NSQS-FILTERED DATA (3/3)      ##########################################################

# NOTE: the list_of_csvdatalists_named_filtered_nsqs object contains all the data after further filtering out data from nuclei with nsqs > nsqs_cutoff
# save to file
saveRDS(list_of_csvdatalists_named_filtered_nsqs,paste0(data_prefix,"oof_nsqs_filtered_data.rds"))

##########################################################################################################################################################
#################################     SAVE R workspace in case something needs to be adjusted/pre-filtered  ##############################################
save.image(file=paste0(data_prefix,"R_workspace.RData"))
