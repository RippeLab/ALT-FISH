#' This function segments the area occupied by the cell nucleus
#'
#' @param im Image series
#' @return Segmentation result

makeNucMask <- function(im, threshold, nucSize_cutoff, remove_bordercells, dilate_disc=3, erode_disc=3) {
  mask <- thresh(im,100,100,threshold) #uses adaptive thresholding function thresh()
  mask <- dilate(mask, makeBrush(dilate_disc, shape='disc'))
  mask <- erode(mask, makeBrush(erode_disc, shape='disc'))
  mask <- fillHull(mask) # fill holes in the mask
  mask <- bwlabel(mask) # find connected sets of pixels
  
  if(remove_bordercells=="y"){
  ###removal of nuclei at image borders###
  # subset for boundary pixels
  dims <- dim(mask)
  border <- c(mask[1:dims[1],1], mask[1:dims[1],dims[2]],mask[1,1:dims[2]], mask[dims[1],1:dims[2]])
  # extract object identifiers at the boundary
  ids <- unique(border[which(border != 0)])
  mask <- rmObjects(mask, ids) #remove objects touching the image border
  ###----------------------------------###
  # select the largest regions of interest
  regions <- table(mask)
  largest_regions <- which(regions[2:length(regions)]>nucSize_cutoff)
  mask[!(mask %in% largest_regions)] <- 0
  
  # return segmentation result
  return(mask)
  
  }else
  {
    # select the largest regions of interest
    regions <- table(mask)
    largest_regions <- which(regions[2:length(regions)]>nucSize_cutoff)
    mask[!(mask %in% largest_regions)] <- 0
    
    # return segmentation result
    return(mask)}
}