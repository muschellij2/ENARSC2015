#' @title N3 or N4 Correct
#' @description This function wraps ANTsR bias field corrections and returns
#' nifti objects
#' @param file (character) image to be manipulated
#' @param correction (character) N3 or N4 correction?
#' @param outfile (character) resultant image name (optional)
#' @param retimg (logical) return image of class nifti
#' @param reorient (logical) If retimg, should file be reoriented when read in?
#' Passed to \code{\link{readNIfTI}}.
#' @param shrinkfactor Shrink factor passed to 
#' \code{N3BiasFieldCorrection}
#' @param ... additional arguments passed to 
#' \code{N3BiasFieldCorrection}
#' @return If \code{retimg} then object of class nifti.  Otherwise,
#' Result from system command, depends if intern is TRUE or FALSE.
#' @export
bias_correct = function(
  file,
  correction = c("N3", "N4"),
  outfile=NULL, 
  retimg = FALSE,
  reorient = FALSE,
  shrinkfactor = "4",
  ...){
  require(ANTsR)
  
  correction = match.arg(correction, c("N3", "N4"))
  func = paste0(correction, "BiasFieldCorrection")
  
  if (retimg){
    if (is.null(outfile)) {
      outfile = tempfile()
      outfile = paste0(nii.stub(outfile), '.nii.gz')      
    }
  } else {
    stopifnot(!is.null(outfile))
  }

  if (inherits(file, "antsImage")){
    img = file
  } else {
    img <- antsImageRead(file, 3)
  }
  imgn3 <- antsImageClone(img)
  funclist = list(d=img@dimension, i=img, o=imgn3, shrinkfactor, ...)
  res = do.call(func, funclist)

  antsImageWrite( imgn3, filename = outfile)

  if (retimg){
    x = readNIfTI(outfile, reorient = reorient)
  } else {
    x = outfile
  }
  return(x)
}