makeSpotMask <- function(im, threshold, nuc_mask, dilate_disc=3, erode_disc=3, ignore_boundary=3, minsize=5) {
  # find regions of interest
  mask <- im>threshold
  mask <- dilate(mask, makeBrush(dilate_disc, shape='disc'))
  mask <- erode(mask, makeBrush(erode_disc, shape='disc'))
  mask[erode(nuc_mask, makeBrush(ignore_boundary, shape='disc'))==0] <- 0 #removes spots that touch the nucleus boundary
  mask <- bwlabel(mask) # find connected sets of pixels
  mask <- fillHull(mask) # fill holes in the mask
  
  # filter out too small regions and remove background region (largest region)
  regions <- table(mask)
  largest_region <- which(regions==max(regions))
  small_regions <- which(regions<minsize)-1 #-1 is required for correctly identifying the region indices
  mask[mask %in% small_regions] <- 0 #set all too small spots to 0
  mask[mask==as.numeric(names(largest_region)[1])] <- 0 #set background(largest region) to 0
  
  # return mask
  return(mask)
}