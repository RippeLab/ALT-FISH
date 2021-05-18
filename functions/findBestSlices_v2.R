findBestSlices <- function(dapi_im, z_num, passfilter) {
  
  # compute mean dapi signal for each z-slice
  dapi_mean_distribution <- vector()
  for(u in 1:z_num){
      dapi_mean_distribution <- c(dapi_mean_distribution, mean(dapi_im[,,u]))
  }
  
  # find indices of the slices where the mean dapi signal is higher than the median of the mean dapi signal distribution (=focus volume)
  best_slices <- c(which(dapi_mean_distribution>median(dapi_mean_distribution)))
  if(length(best_slices)==0){best_slices_range <-0}
  best_slices_range <- c(best_slices[1]:best_slices[length(best_slices)])
  
  # consider position out-of-focus and set best_slices_range to zero when the max of the mean distribution is too close to the z-stack ends
  if(which.max(dapi_mean_distribution)>46 | which.max(dapi_mean_distribution)<5){
    best_slices_range <- 0
  }
  
  # consider frame empty if the maximum of the mean distribution is not at least x-fold the minimum of the mean distribution 
  # (even frames with only a single cell usually do not meet this criterion)
  if(max(dapi_mean_distribution)/min(dapi_mean_distribution) < passfilter){
    best_slices_range <- 0
  }
  
  # return index range of best slices
  return(best_slices_range)
}
