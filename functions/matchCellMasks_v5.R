#' This function filters out cytoplasm masks without matching nucleus mask
#'
#' @param (bw)labeled nucleus and cytoplasm masks from the same image
#' @return "cytoplasm-only" masks that have a matching nucleus in the image or nucleus masks which have matching cytoplasm
#' this is defined by returnMaskType parameter. Importantly, also masks with unambiguous assignments (e.g. cytoplasmic mask with two nuclei) are filtered out.
matchCellMasks <- function(nucmask, cytomask, returnMaskType) {
  #compute nucleus features to extract coordinates and minimum radii
  nucfeatures <- data.frame(cbind(computeFeatures.moment(nucmask),computeFeatures.shape(nucmask)))
  nuccoord_minrad <- data.frame("m.cx"=c(round(nucfeatures$m.cx)),"m.cy"=c(round(nucfeatures$m.cy)), "s.radius.min"=c(nucfeatures$s.radius.min), "matchingmask"=as.logical(c(rep("FALSE",length(nucfeatures$m.cx)))))
  #generate binary cytomask
  #note that cytomask is a mask comprising the cytoplasm + nucleus area
  cytomask_binary <- cytomask
  cytomask_binary[cytomask_binary>1] <- 1
  #and check if cytomask_binary value is 1 (mask) at extracted nucleus coordinates which is documented in the matchingmask column as TRUE/FALSE
  #create image with values 0 of the same dimensions as cpmask
  #this image is used to draw circular objects of pixel value 1 with radius = min.radius + 15 px at the center of mass of each nucleus 
  blackframe <- Image(0, dim=dim(cytomask))
  for(w in 1:length(rownames(nuccoord_minrad))){
    nuccoord_minrad$matchingmask[w] <- isTRUE(cytomask_binary[nuccoord_minrad$m.cx[w],nuccoord_minrad$m.cy[w]]==1)
    blackframe <- drawCircle(blackframe, nuccoord_minrad$m.cx[w], nuccoord_minrad$m.cy[w], nuccoord_minrad$s.radius.min[w]+15, col=1, fill=TRUE)
  }
  #generate "cytoplasm only" mask by setting all pixels 0 which correspond to nuclei identified in the nucmask
  cpmask <- cytomask
  cpmask[nucmask>=1] <- 0
  #make vector of object IDs of the single cpmasks, counting from 2:end to exclude the first object (background) 
  cpmask_ids <-as.numeric(names(table(cpmask)[2:length(table(cpmask))]))
  
  ## procedure to identify cells with more than one nucleus per cytoplasm mask
  #label nucleus representation objects in blackframe
  blackframe_labelled <- bwlabel(blackframe)
  #retrieve object ids from blackframe representation
  repres_ids <- as.numeric(names(table(blackframe_labelled)[2:length(table(blackframe_labelled))]))
  #loop over each nucleus representation object in blackframe and retrieve it's corresponding cpmask id
  #initialize vector to store corresponding cpmask ids
  corresp_cpmask_ids <- rep(0,length(repres_ids))
  for(v in 1:length(repres_ids)){
    # extract correspondng cpmask id for current nucleus representation v by taking the maximum pixel value in the cpmask image in the area confined by the current nucleus representation
    corresp_cpmask_ids[v] <- max(cpmask[blackframe_labelled==repres_ids[v]])
  }
  
  #identify object IDs of cpmasks by extracting pixel values >0 of cpmask in the area confined by circles in blackframe
  #pixel values corresponds to the object ID
  #zero pixels are excluded because they are background
  matchedcpmasks_ids <- as.numeric(names(table(cpmask[which(blackframe==1)])))[which(as.numeric(names(table(cpmask[which(blackframe==1)])))>0)]
  #identify differences between cpmask_ids and matchedcpmask_ids and remove objects with unmatched ids from cpmask
  unmatched_ids <- setdiff(cpmask_ids,matchedcpmasks_ids)
  #change storage mode of cpmask to "integer" (is required by rmObjects()) and is later changed back to "double"
  storage.mode(cpmask) <- "integer"
  if(length(unmatched_ids)==0){
    #remove duplicate ids (e.g. two nuclei which have been assigned the same cpmask)
    cpmask_filtered <- rmObjects(cpmask,corresp_cpmask_ids[which(duplicated(corresp_cpmask_ids))], reenumerate = FALSE)
  }else{
    #remove duplicate ids (e.g. two nuclei which have been assigned the same cpmask)
    cpmask_filtered <- rmObjects(cpmask,corresp_cpmask_ids[which(duplicated(corresp_cpmask_ids))],reenumerate = FALSE)
    #remove unmatched ids
    cpmask_filtered <- rmObjects(cpmask, unmatched_ids,reenumerate = FALSE)
  }
  #return filtered result
  if(returnMaskType=="nuclei"){
    #return filtered nucmask
    #the filtering removes nuclei without matching cpmask
    #this is done by generating a "filled" and dilated cpmask_filtered and simply setting every pixel value to 0 in the nucmask which is outside of this filled mask
    cpmask_filtered_filled <- dilate(cpmask_filtered, makeBrush(61, shape='diamond'))
    cpmask_filtered_filled <- fillHull(cpmask_filtered_filled)
    nucmask_filtered <- nucmask
    nucmask_filtered[cpmask_filtered_filled==0] <- 0
    #re-label/renumerate the filtered nucmask by bwlabel
    nucmask_filtered <- bwlabel(nucmask_filtered)
    return(nucmask_filtered)
  }else{
    if(returnMaskType=="cytoplasm"){
      #return filtered cytoplasm mask
      storage.mode(cpmask_filtered) <- "double"
      return(cpmask_filtered)
    }
    if(returnMaskType=="matchingtable"){
      #same procedure as above (removal of lone nuclei)
      cpmask_filtered_filled <- dilate(cpmask_filtered, makeBrush(61, shape='diamond'))
      cpmask_filtered_filled <- fillHull(cpmask_filtered_filled)
      nucmask_filtered <- nucmask
      nucmask_filtered[cpmask_filtered_filled==0] <- 0
      #re-label/renumerate the filtered nucmask by bwlabel
      nucmask_filtered <- bwlabel(nucmask_filtered)
      
      # collect the final ids (bwlabel pixel values) of cpmasks with corresponding nucleus in the order of the filtered nucmask ids (bwlabel pixel values)
      # remove potential duplicate ids which arise during unambiguous matching but preserving the original order
      final_cpmask_ids <- corresp_cpmask_ids
      final_cpmask_ids <- final_cpmask_ids[!final_cpmask_ids %in% unique(final_cpmask_ids[which(duplicated(final_cpmask_ids))])]
      # collect the final ids (bwlabel pixel values) of nucmasks remaining after filtering
      final_nuc_ids <- as.numeric(names(table(nucmask_filtered)[2:length(table(nucmask_filtered))]))
      # combine both into a data frame and return
      return(as.data.frame(cbind(final_nuc_ids,final_cpmask_ids)))
    }else{
      print("Please specify the returnMaskType parameter as one of the following: nuclei, cytoplasm, matchingtable")
    }
  }
  
}