evaluateZProfile <- function(im_nuc,nucmask,z_num) {
  
  #collect nucleus ids
  rois_nuc <- as.numeric(names(table(nucmask)))[2:length(as.numeric(names(table(nucmask))))]
  
  #create empty data frame with each column representing one nucleus and each row the z coordinate/z slice
  mean_z_nuc_ch_intensity <- data.frame(matrix(ncol = length(rois_nuc), nrow = z_num))
  colnames(mean_z_nuc_ch_intensity) <- paste0("nuc_",1:length(rois_nuc))
  
  #create empty vector for storing the focus evaluation results (one value per nucleus)
  nuc_infocus <- vector(length=length(rois_nuc))
  
  # compute mean dapi signal for each nucleus mask region in each z-slice and summarize in data frame
  for(j in 1:length(rois_nuc)){
    for(zcoord in 1:zperch){
      nuc <- nuc_data[,,zcoord]
      #nuc[mask!=rois_nuc[j]] <- NA
      mean_z_nuc_ch_intensity[zcoord,j] <- mean(nuc[nucmask==rois_nuc[j]])
    }
    # evaluate nucleus as in-focus or not-in-focus based on the DAPI signal distribution along z
    # find indices of the slices where the mean dapi signal is higher than the median of the dapi signal distribution along z
    best_slices <- c(which(mean_z_nuc_ch_intensity[,j]>median(mean_z_nuc_ch_intensity[,j])))
    if(length(best_slices)==0){best_slices_range <-0}
    best_slices_range <- c(best_slices[1]:best_slices[length(best_slices)])
    # consider position out-of-focus and set best_slices_range to zero when the max of the mean distribution is too close to the z-stack ends
    if(which.max(mean_z_nuc_ch_intensity[,j])>36 | which.max(mean_z_nuc_ch_intensity[,j])<10){
      best_slices_range <- 0
    }
    # if a suitable best_slices_range could be found, annotate nucleus focus as "y", else as "n"
    if(length(best_slices_range)==1){
      nuc_infocus[j] <- "n"
    }else{nuc_infocus[j] <- "y"}
  }
  
  #combine results: dapi intensity profiles + annotation
  #annotation will be stored in the last row (row number: z_num + 1)
  mean_z_nuc_ch_intensity_annotated <- rbind(mean_z_nuc_ch_intensity,nuc_infocus)
  # return index range of best slices
  return(mean_z_nuc_ch_intensity_annotated)
}
