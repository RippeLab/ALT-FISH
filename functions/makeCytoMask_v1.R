#' This function segments the area occupied by the cell cytoplasm
#'
#' @param im Image series
#' @return Segmentation result

makeCytoMask <- function(im_nuc, im_cyto, cytoSize_cutoff, cyto_threshold, remove_bordercells) {
  #make combined and intensity-transformed image of nucleus and cytoplasm channel
  im <- gblur(im_nuc+im_cyto^0.2,2)
  #rough 1st pass thresholding on combined image to generate preliminary mask
  mask <- thresh(im,w=500,h=500)
  #use inverse of preliminary mask to estimate background outside the mask area
  #2nd pass thresholding using the background estimate * cyto_threshold
  mask <- im_cyto > (median(im_cyto[which(mask==0)])*cyto_threshold)
  #fill holes and erode
  mask <- erode(mask, makeBrush(5, shape='diamond'))
  mask <- fillHull(mask)
  #use watershed function on distancemap of the mask to separate object
  #distmap() creates an image object where each foreground pixel (value 1) contains the distance to the closest background pixel(value 0) as value
  #distance maps are required by the watershed function
  mask <- watershed(distmap(mask),10)
  
  if(remove_bordercells=="y"){
  ###removal of nuclei at image borders###
  # subset for boundary pixels
  dims <- dim(mask)
  border <- c(mask[1:dims[1],1], mask[1:dims[1],dims[2]],mask[1,1:dims[2]], mask[dims[1],1:dims[2]])
  # extract object identifiers at the boundary
  ids <- unique(border[which(border != 0)])
  # account for objects that only slightly touch the borders (100 px "border touch cutoff")
  ids_filtered <- c(as.integer(names(which(table(border)[2:length(table(border))]>100))))
  
  mask <- rmObjects(mask, ids_filtered) #remove objects touching the image border
  ###----------------------------------###
  # select the largest regions of interest
  regions <- table(mask)
  largest_regions <- which(regions[2:length(regions)]>cytoSize_cutoff)
  mask[!(mask %in% largest_regions)] <- 0
  
  # return segmentation result
  return(mask)
  
  }else
  {
    # select the largest regions of interest
    regions <- table(mask)
    largest_regions <- which(regions[2:length(regions)]>cytoSize_cutoff)
    mask[!(mask %in% largest_regions)] <- 0
    
    # return segmentation result
    return(mask)}
}